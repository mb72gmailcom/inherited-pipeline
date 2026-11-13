"""Invert short-format site TSVs (variant -> patients) into patient -> variants TSVs."""

from __future__ import annotations

import re
from array import array
from collections.abc import Iterator
from pathlib import Path

from inherited.filter import PREFIXES
from inherited.params import write_params

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
PATIENTS_HEADER = "#PATIENT\tN_VARIANTS\tVARIANTS\n"
PATIENTS_FILE = "patients.tsv"
COUNTS_FILE = "counts.tsv"


def invert_sites(
    input_dir: Path,
    output_dir: Path,
    *,
    prefix: str = "inherited",
    chroms: list[str] | None = None,
) -> dict:
    """Write patient -> variants TSVs from ``{input_dir}/chrN/{prefix}_*.tsv``.

    Input rows are ``#CHROM POS ID REF ALT PATIENTS`` with ``;``-separated patients.
    Every shard file of a chromosome is read, whatever its naming
    (``inherited_00000.tsv`` or ``inherited_<start>_<end>.tsv``).

    Outputs in ``output_dir``:
        chrN/patients.tsv  one line per patient: PATIENT, N_VARIANTS, VARIANTS
        counts.tsv         variants per patient: PATIENT, TOTAL, one column per chromosome
        params.json
    VARIANTS are ``chrom_pos_ref_alt`` keys joined by ``;``, ordered by position,
    REF and ALT. No genome-wide patient file is written. Patients are sorted by ID. A variant listed in more
    than one input row is written once per patient.
    """
    if prefix not in PREFIXES:
        raise ValueError(f"prefix must be one of {', '.join(PREFIXES)}, got {prefix!r}")
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")

    chrom_dirs = _chrom_dirs(input_dir)
    if chroms:
        wanted = [_normalize_chrom(chrom) for chrom in chroms]
        available = {path.name: path for path in chrom_dirs}
        missing = [chrom for chrom in wanted if chrom not in available]
        if missing:
            raise ValueError(f"chromosome folders not found under {input_dir}: {', '.join(missing)}")
        chrom_dirs = sorted({available[c] for c in wanted}, key=lambda p: _chrom_sort_key(p.name))
    if not chrom_dirs:
        raise ValueError(f"no chrN folders under {input_dir}")
    _check_output_dir(output_dir)

    summary: dict[str, dict] = {}
    counts: dict[str, dict[str, int]] = {}
    for chrom_dir in chrom_dirs:
        files = _site_files(chrom_dir, prefix)
        if not files:
            raise ValueError(f"no {prefix}_*.tsv files in {chrom_dir}")
        destination = output_dir / chrom_dir.name / PATIENTS_FILE
        stats = _invert_chromosome(chrom_dir.name, files, destination)
        for patient, n in stats.pop("per_patient").items():
            counts.setdefault(patient, {})[chrom_dir.name] = n
        summary[chrom_dir.name] = stats
        print(
            f"{chrom_dir.name}: {stats['variants']} variants, {stats['pairs']} patient-variant "
            f"pairs, {stats['patients']} patients -> {destination}",
            flush=True,
        )

    names = [path.name for path in chrom_dirs]
    _write_counts(output_dir / COUNTS_FILE, counts, names)

    totals = {
        "variants": sum(s["variants"] for s in summary.values()),
        "pairs": sum(s["pairs"] for s in summary.values()),
        "patients": len(counts),
    }
    write_params(
        output_dir / "params.json",
        {
            "command": "invert",
            "input_dir": str(input_dir.resolve()),
            "output_dir": str(output_dir.resolve()),
            "prefix": prefix,
            "chromosomes": names,
            "variant_key": "chrom_pos_ref_alt",
            "layout": "per chromosome: chrN/patients.tsv",
            "per_chromosome": summary,
            "totals": totals,
        },
    )
    return totals


def _invert_chromosome(chrom: str, files: list[Path], destination: Path) -> dict:
    rows: list[tuple[int, str, str, str]] = []
    for path in files:
        rows.extend(_read_rows(path, chrom))
    rows.sort(key=lambda row: (row[0], row[1], row[2]))

    keys: list[str] = []
    per_patient: dict[str, array] = {}
    duplicate_rows = 0
    previous = None
    for pos, ref, alt, patients in rows:
        if (pos, ref, alt) == previous:
            duplicate_rows += 1
        else:
            keys.append(f"{chrom}_{pos}_{ref}_{alt}")
            previous = (pos, ref, alt)
        index = len(keys) - 1
        for patient in patients.split(";"):
            if not patient:
                continue
            indices = per_patient.get(patient)
            if indices is None:
                per_patient[patient] = array("I", [index])
            elif indices[-1] != index:
                indices.append(index)
    del rows

    destination.parent.mkdir(parents=True, exist_ok=True)
    pairs = 0
    with destination.open("w", encoding="utf-8") as handle:
        handle.write(PATIENTS_HEADER)
        for patient in sorted(per_patient):
            indices = per_patient[patient]
            pairs += len(indices)
            variants = ";".join(keys[i] for i in indices)
            handle.write(f"{patient}\t{len(indices)}\t{variants}\n")
    return {
        "files": len(files),
        "variants": len(keys),
        "duplicate_rows": duplicate_rows,
        "pairs": pairs,
        "patients": len(per_patient),
        "per_patient": {patient: len(indices) for patient, indices in per_patient.items()},
    }


def _read_rows(path: Path, chrom: str) -> Iterator[tuple[int, str, str, str]]:
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.split("\t")[-1] == "TRIO_CALLS":
                    raise ValueError(f"{path} requires input in short format (PATIENTS column)")
                continue
            fields = line.split("\t")
            if len(fields) != 6:
                raise ValueError(f"{path}:{number} expected 6 tab-separated columns, got {len(fields)}")
            if fields[0] != chrom:
                raise ValueError(f"{path}:{number} has chromosome {fields[0]}, expected {chrom}")
            if "=" in fields[5]:
                raise ValueError(f"{path} requires input in short format (PATIENTS column)")
            try:
                pos = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{number} has a non-integer POS {fields[1]!r}") from exc
            yield pos, fields[3], fields[4], fields[5]


def _write_counts(path: Path, counts: dict[str, dict[str, int]], chroms: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(["#PATIENT", "TOTAL", *chroms]) + "\n")
        for patient in sorted(counts):
            row = [counts[patient].get(chrom, 0) for chrom in chroms]
            handle.write("\t".join([patient, str(sum(row)), *map(str, row)]) + "\n")


def _check_output_dir(output_dir: Path) -> None:
    """Refuse an output tree whose chrN folders hold anything but our patients.tsv.

    This protects VCF trees (chrN/variants_*.vcf), site trees (chrN/inherited_*.tsv)
    and other results from a mistyped --output-dir; re-running into a previous
    invert output is allowed.
    """
    if not output_dir.is_dir():
        return
    for chrom_dir in _chrom_dirs(output_dir):
        foreign = sorted(p.name for p in chrom_dir.iterdir() if p.name != PATIENTS_FILE)
        if foreign:
            shown = ", ".join(foreign[:3]) + (" ..." if len(foreign) > 3 else "")
            raise ValueError(
                f"output directory {output_dir} is not an invert output: "
                f"{chrom_dir.name}/ already contains {shown}"
            )


def _site_files(chrom_dir: Path, prefix: str) -> list[Path]:
    return sorted(path for path in chrom_dir.glob(f"{prefix}_*.tsv") if path.is_file())


def _chrom_dirs(input_dir: Path) -> list[Path]:
    dirs = [p for p in input_dir.iterdir() if p.is_dir() and _CHROM_DIR.fullmatch(p.name)]
    return sorted(dirs, key=lambda path: _chrom_sort_key(path.name))


def _chrom_sort_key(name: str) -> tuple[int, int | str]:
    suffix = name[3:]
    if suffix == "X":
        return (1, "X")
    if suffix == "Y":
        return (2, "Y")
    return (0, int(suffix))


def _normalize_chrom(text: str) -> str:
    return text if text.startswith("chr") else f"chr{text}"
