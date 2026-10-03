"""Write variants present in one VCF tree and absent from the paired file in another."""

from __future__ import annotations

import re
from pathlib import Path

from inherited.filter import PREFIXES

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
_VCF_NAME = re.compile(r"^variants_\d+_\d+\.vcf$")
TSV_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"


def subtract_variants(input_dir1: Path, input_dir2: Path, output_dir: Path) -> None:
    """Subtract paired ``variants_${start}_${end}.vcf`` files.

    A row from the first directory is written when its ``CHROM``, ``POS``,
    ``REF``, and ``ALT`` are absent from the same filename in the second
    directory. Files that exist only in the first directory are copied.
    The same variants are written as short-format TSVs under
    ``{output_dir}-sites``, with ``PATIENTS`` copied from ``{input_dir1}-sites``.
    """
    if not input_dir1.is_dir():
        raise ValueError(f"input directory not found: {input_dir1}")
    if not input_dir2.is_dir():
        raise ValueError(f"input directory not found: {input_dir2}")
    sites_dir1 = input_dir1.with_name(input_dir1.name + "-sites")
    if not sites_dir1.is_dir():
        raise ValueError(f"sites directory not found: {sites_dir1}")
    sites_out = output_dir.with_name(output_dir.name + "-sites")

    patients: set[str] = set()
    for chrom_dir in _chrom_dirs(input_dir1):
        for vcf_path in _variant_files(chrom_dir):
            paired = input_dir2 / chrom_dir.name / vcf_path.name
            blocked = _variant_keys(paired) if paired.is_file() else set()
            sites_path = _sites_file(sites_dir1, chrom_dir.name, vcf_path.name)
            site_patients = _read_site_patients(sites_path)
            destination = output_dir / chrom_dir.name / vcf_path.name
            tsv_destination = sites_out / chrom_dir.name / sites_path.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            tsv_destination.parent.mkdir(parents=True, exist_ok=True)
            patients.update(
                _write_difference(
                    vcf_path,
                    destination,
                    tsv_destination,
                    blocked,
                    site_patients,
                    sites_path,
                )
            )
    _write_patient_list(sites_out / "patients.txt", patients)


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


def _variant_files(chrom_dir: Path) -> list[Path]:
    files = [
        path
        for path in chrom_dir.iterdir()
        if path.is_file() and _VCF_NAME.fullmatch(path.name)
    ]
    return sorted(files, key=lambda path: path.name)


def _variant_keys(path: Path) -> set[tuple[str, str, str, str]]:
    return {key for key, _line in _iter_vcf_rows(path)}


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
    blocked: set[tuple[str, str, str, str]],
    site_patients: dict[tuple[str, str, str, str], str],
    sites_path: Path,
) -> set[str]:
    found: set[str] = set()
    with destination.open("w", encoding="utf-8") as vcf, tsv_destination.open(
        "w", encoding="utf-8"
    ) as tsv:
        tsv.write(TSV_HEADER)
        for key, line in _iter_vcf_rows(source):
            if key in blocked:
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


def _write_patient_list(path: Path, patients: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for patient in sorted(patients):
            handle.write(f"{patient}\n")
