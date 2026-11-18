"""Write variants present in one VCF tree and absent from another.

Two ways to match the trees:

* by file name (default): ``chrN/variants_<start>_<end>.vcf`` in the first tree is compared with the
  file of the same name in the second. Right for trees split by genomic position.
* per chromosome (``concatenate=True``): every variant of a chromosome in the second tree removes
  that variant from the first, whatever files the two trees split it into. Needed for older trees
  split by variant count (``variants_00000.vcf``, ``variants_00001.vcf``, ...), whose segment
  boundaries differ between runs. chrX/chrY files are matched within their group (``females``,
  ``males_par1``, ``males_par2``, ``males_nonPar``), not across groups.

Output files always keep the first tree's names.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from inherited.filter import PREFIXES
from inherited.params import write_params
from inherited.xchrom import X_BUCKETS, Y_BUCKETS

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
_BUCKETS = sorted(set(X_BUCKETS) | set(Y_BUCKETS), key=len, reverse=True)
_VCF_NAME = re.compile(
    r"^variants(?:_(?P<bucket>" + "|".join(_BUCKETS) + r"))?"
    r"(?:_(?P<first>\d+)(?:_(?P<second>\d+))?)?\.vcf$"
)
TSV_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"

POSITION, SEGMENT, WHOLE = "position", "segment", "whole"


@dataclass(frozen=True)
class VcfFile:
    """One ``variants*.vcf`` file: its group (chrX/chrY bucket or "") and naming style."""

    path: Path
    bucket: str
    naming: str  # POSITION (_<start>_<end>), SEGMENT (_<index>) or WHOLE (no number)
    order: int

    @property
    def name(self) -> str:
        return self.path.name


def parse_vcf_name(path: Path) -> VcfFile | None:
    match = _VCF_NAME.fullmatch(path.name)
    if match is None:
        return None
    bucket = match["bucket"] or ""
    if match["second"] is not None:
        return VcfFile(path, bucket, POSITION, int(match["first"]))
    if match["first"] is not None:
        return VcfFile(path, bucket, SEGMENT, int(match["first"]))
    return VcfFile(path, bucket, WHOLE, 0)


def subtract_variants(
    input_dir1: Path,
    input_dir2: Path,
    output_dir: Path,
    *,
    sites: bool = True,
    concatenate: bool = False,
) -> dict:
    """Subtract the variants of ``input_dir2`` from ``input_dir1``.

    A row of the first tree is written when its ``CHROM``, ``POS``, ``REF`` and ``ALT`` are absent
    from the second tree: from the file of the same name by default, or from the whole chromosome
    (and chrX/chrY group) with ``concatenate``. Output files keep the first tree's names; a file
    with nothing to compare against is copied. When ``sites`` is true, the same variants are
    written as short-format TSVs under ``{output_dir}-sites``, with ``PATIENTS`` copied from
    ``{input_dir1}-sites``, plus a ``patients.txt``. When ``sites`` is false, only the VCF tree is
    written and no sites directory is read or created. Parameters are recorded in
    ``{output_dir}/params.json``. Returns the per-group file counts written there.
    """
    if not input_dir1.is_dir():
        raise ValueError(f"input directory not found: {input_dir1}")
    if not input_dir2.is_dir():
        raise ValueError(f"input directory not found: {input_dir2}")
    sites_dir1 = input_dir1.with_name(input_dir1.name + "-sites")
    if sites and not sites_dir1.is_dir():
        raise ValueError(f"sites directory not found: {sites_dir1}")
    sites_out = output_dir.with_name(output_dir.name + "-sites")

    plan = _plan(input_dir1, input_dir2, concatenate)
    patients: set[str] = set()
    groups: dict[str, dict[str, int]] = {}
    for chrom, bucket, files1, blocked_sources in plan:
        label = f"{chrom}/{bucket}" if bucket else chrom
        groups[label] = {"files1": len(files1), "files2": len(blocked_sources)}
        if concatenate:
            group_blocked = _keys_of(blocked_sources)
        for vcf in files1:
            if concatenate:
                blocked = group_blocked
            else:
                blocked = _keys_of(blocked_sources[vcf.name] if vcf.name in blocked_sources else [])
            destination = output_dir / chrom / vcf.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not sites:
                _write_vcf_difference(vcf.path, destination, blocked)
                continue
            sites_path = _sites_file(sites_dir1, chrom, vcf.name)
            site_patients = _read_site_patients(sites_path)
            tsv_destination = sites_out / chrom / sites_path.name
            tsv_destination.parent.mkdir(parents=True, exist_ok=True)
            patients.update(
                _write_difference(
                    vcf.path, destination, tsv_destination, blocked, site_patients, sites_path
                )
            )
    if sites:
        _write_patient_list(sites_out / "patients.txt", patients)
    _write_params(input_dir1, input_dir2, output_dir, sites=sites, concatenate=concatenate,
                  groups=groups)
    return groups


def _plan(input_dir1: Path, input_dir2: Path, concatenate: bool):
    """List (chrom, bucket, files of tree 1, what to block) for every group of tree 1.

    By file name, "what to block" maps a tree-1 file name to the tree-2 file of that name;
    with ``concatenate`` it is the list of all tree-2 files of the group.
    """
    plan = []
    for chrom_dir in _chrom_dirs(input_dir1):
        chrom = chrom_dir.name
        files1 = _variant_files(chrom_dir)
        dir2 = input_dir2 / chrom
        files2 = _variant_files(dir2) if dir2.is_dir() else []
        if not concatenate:
            for tree, files in ((input_dir1, files1), (input_dir2, files2)):
                segmented = [f.name for f in files if f.naming == SEGMENT]
                if segmented:
                    raise ValueError(
                        f"{tree / chrom} has segment-numbered files ({segmented[0]}, ...), whose "
                        "boundaries differ between runs and cannot be paired by name; "
                        "use --concatenate"
                    )
            files1 = [f for f in files1 if f.naming == POSITION]
            by_name = {f.name: [f] for f in files2 if f.naming == POSITION}
            for bucket, group in _by_bucket(files1).items():
                plan.append((chrom, bucket, group,
                             {f.name: by_name[f.name] for f in group if f.name in by_name}))
            continue
        groups2 = _by_bucket(files2)
        for tree, groups in ((input_dir1, _by_bucket(files1)), (input_dir2, groups2)):
            for bucket, group in groups.items():
                styles = sorted({f.naming for f in group})
                if len(styles) > 1:
                    where = f"{tree / chrom}" + (f" ({bucket})" if bucket else "")
                    raise ValueError(
                        f"{where} mixes file naming styles ({', '.join(styles)}): "
                        + ", ".join(f.name for f in group[:4])
                    )
        for bucket, group in _by_bucket(files1).items():
            plan.append((chrom, bucket, group, groups2.get(bucket, [])))
    if not plan:
        expected = "variants*.vcf" if concatenate else "variants_<start>_<end>.vcf"
        raise ValueError(f"no {expected} files found under {input_dir1}/chr*")
    return plan


def _by_bucket(files: list[VcfFile]) -> dict[str, list[VcfFile]]:
    groups: dict[str, list[VcfFile]] = defaultdict(list)
    for vcf in files:
        groups[vcf.bucket].append(vcf)
    return dict(groups)


def _keys_of(files: list[VcfFile]) -> set[str]:
    """Compact keys (CHROM, POS, REF, ALT joined by tabs) of every row in ``files``."""
    keys: set[str] = set()
    for vcf in files:
        for key, _line in _iter_vcf_rows(vcf.path):
            keys.add(_compact(key))
    return keys


def _compact(key: tuple[str, str, str, str]) -> str:
    return "\t".join(key)


def _write_params(
    input_dir1: Path,
    input_dir2: Path,
    output_dir: Path,
    *,
    sites: bool,
    concatenate: bool,
    groups: dict[str, dict[str, int]],
) -> Path:
    payload = {
        "command": "subtract",
        "concatenate": concatenate,
        "files": groups,
        "input_dir1": str(input_dir1.resolve()),
        "input_dir2": str(input_dir2.resolve()),
        "matching": (
            "per chromosome (chrX/chrY per group): all input_dir2 files of a chromosome together"
            if concatenate
            else "per file: input_dir1 file against the input_dir2 file of the same name"
        ),
        "output_dir": str(output_dir.resolve()),
        "sites": sites,
    }
    return write_params(output_dir / "params.json", payload)


def _chrom_dirs(input_dir: Path) -> list[Path]:
    dirs = [
        path
        for path in input_dir.iterdir()
        if path.is_dir() and _CHROM_DIR.fullmatch(path.name)
    ]
    return sorted(dirs, key=lambda path: _chrom_sort_key(path.name))


def _chrom_sort_key(name: str) -> tuple[int, int | str]:
    suffix = name[3:]
    if suffix == "X":
        return (1, "X")
    if suffix == "Y":
        return (2, "Y")
    return (0, int(suffix))


def _variant_files(chrom_dir: Path) -> list[VcfFile]:
    files = [
        parsed
        for path in chrom_dir.iterdir()
        if path.is_file() and (parsed := parse_vcf_name(path)) is not None
    ]
    return sorted(files, key=lambda f: (f.bucket, f.order, f.name))


def _iter_vcf_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                raise ValueError(f"{path} expected at least 5 tab-separated columns")
            yield (fields[0], fields[1], fields[3], fields[4]), line


def _sites_file(sites_dir: Path, chrom: str, vcf_name: str) -> Path:
    rest = Path(vcf_name).stem.removeprefix("variants")
    chrom_dir = sites_dir / chrom
    matches = [
        chrom_dir / f"{prefix}{rest}.tsv"
        for prefix in PREFIXES
        if (chrom_dir / f"{prefix}{rest}.tsv").is_file()
    ]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError(f"sites file not found for {chrom}/{vcf_name} in {sites_dir}")
    names = ", ".join(path.name for path in matches)
    raise ValueError(f"ambiguous sites files for {chrom}/{vcf_name}: {names}")


def _read_site_patients(path: Path) -> dict[tuple[str, str, str, str], str]:
    patients: dict[tuple[str, str, str, str], str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                if line.startswith("#") and line.split("\t")[-1] == "TRIO_CALLS":
                    raise ValueError(f"{path} requires input in short format")
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path} expected 6 tab-separated columns")
            if "=" in fields[5]:
                raise ValueError(f"{path} requires input in short format")
            patients[(fields[0], fields[1], fields[3], fields[4])] = fields[5]
    return patients


def _write_difference(
    source: Path,
    destination: Path,
    tsv_destination: Path,
    blocked: set[str],
    site_patients: dict[tuple[str, str, str, str], str],
    sites_path: Path,
) -> set[str]:
    found: set[str] = set()
    with destination.open("w", encoding="utf-8") as vcf, tsv_destination.open(
        "w", encoding="utf-8"
    ) as tsv:
        tsv.write(TSV_HEADER)
        for key, line in _iter_vcf_rows(source):
            if _compact(key) in blocked:
                continue
            if key not in site_patients:
                chrom, pos, ref, alt = key
                raise ValueError(
                    f"{chrom}:{pos}:{ref}:{alt} from {source.name} is missing from {sites_path}"
                )
            patients = site_patients[key]
            vcf.write(f"{line}\n")
            tsv.write(f"{line}\t{patients}\n")
            found.update(patient for patient in patients.split(";") if patient)
    return found


def _write_vcf_difference(
    source: Path, destination: Path, blocked: set[str]
) -> None:
    with destination.open("w", encoding="utf-8") as vcf:
        for key, line in _iter_vcf_rows(source):
            if _compact(key) not in blocked:
                vcf.write(f"{line}\n")


def _write_patient_list(path: Path, patients: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for patient in sorted(patients):
            handle.write(f"{patient}\n")
