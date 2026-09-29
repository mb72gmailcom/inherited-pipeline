"""Filter result TSVs by how common a variant is and how many variants a person carries."""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

PREFIXES = ("inherited", "denovo", "mendelian_bad")

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")

VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n"
TSV_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
_HIST_HEADER = "count\tn\n"


def filter_results(
    input_dir: Path,
    output_dir: Path,
    prefix: str,
    *,
    variant_cap: int | None = None,
    patient_cap: int | None = None,
) -> Path | None:
    """Count variants and, when both caps are set, write filtered VCF and TSV trees.

    Returns the ``-sites`` directory when caps are applied, otherwise ``None``.
    """
    if prefix not in PREFIXES:
        raise ValueError(
            f"prefix must be one of {', '.join(PREFIXES)}, got {prefix!r}"
        )
    _validate_caps(variant_cap, patient_cap)
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")

    tsv_files = _discover_tsv_files(input_dir, prefix)
    person_counts, carrier_counts = _count_variants(tsv_files)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_histogram(
        output_dir / f"{prefix}_variants_per_person_hist.tsv",
        Counter(person_counts.values()),
    )
    _write_histogram(
        output_dir / f"{prefix}_patients_per_variant_hist.tsv",
        carrier_counts,
    )
    if variant_cap is None or patient_cap is None:
        return None

    sites_dir = output_dir.with_name(output_dir.name + "-sites")
    sites_dir.mkdir(parents=True, exist_ok=True)
    _write_filtered(
        tsv_files,
        output_dir,
        sites_dir,
        person_counts,
        variant_cap,
        patient_cap,
    )
    return sites_dir


def _validate_caps(variant_cap: int | None, patient_cap: int | None) -> None:
    if (variant_cap is None) != (patient_cap is None):
        raise ValueError("--variant-cap and --patient-cap must be given together")
    for name, value in (("--variant-cap", variant_cap), ("--patient-cap", patient_cap)):
        if value is not None and value < 0:
            raise ValueError(f"{name} must be an integer >= 0")


def _discover_tsv_files(input_dir: Path, prefix: str) -> list[tuple[str, Path]]:
    found: list[tuple[str, Path]] = []
    for chrom_dir in _chrom_dirs(input_dir):
        matches = [
            path
            for path in chrom_dir.glob(f"{prefix}*.tsv")
            if path.is_file()
        ]
        for path in sorted(matches, key=lambda item: item.name):
            found.append((chrom_dir.name, path))
    return found


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


def _count_variants(tsv_files: list[tuple[str, Path]]) -> tuple[Counter[str], Counter[int]]:
    person_counts: Counter[str] = Counter()
    carrier_counts: Counter[int] = Counter()
    for _chrom, path in tsv_files:
        for _columns, patients in _iter_short_rows(path):
            carrier_counts[len(patients)] += 1
            person_counts.update(patients)
    return person_counts, carrier_counts


def _iter_short_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.split("\t")[-1] == "TRIO_CALLS":
                    raise ValueError(f"{path} requires input in short format")
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path} expected 6 tab-separated columns")
            patients_field = fields[5]
            if "=" in patients_field:
                raise ValueError(f"{path} requires input in short format")
            patients = [patient for patient in patients_field.split(";") if patient]
            yield fields[:5], patients


def _write_histogram(path: Path, counts: Counter[int]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(_HIST_HEADER)
        for count in sorted(counts):
            handle.write(f"{count}\t{counts[count]}\n")


def _write_filtered(
    tsv_files: list[tuple[str, Path]],
    output_dir: Path,
    sites_dir: Path,
    person_counts: Counter[str],
    variant_cap: int,
    patient_cap: int,
) -> None:
    for chrom, path in tsv_files:
        vcf_path = output_dir / chrom / f"{path.stem}.vcf"
        tsv_path = sites_dir / chrom / path.name
        vcf_path.parent.mkdir(parents=True, exist_ok=True)
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        with vcf_path.open("w", encoding="utf-8") as vcf, tsv_path.open(
            "w", encoding="utf-8"
        ) as tsv:
            vcf.write(VCF_HEADER)
            tsv.write(TSV_HEADER)
            for columns, patients in _iter_short_rows(path):
                if len(patients) > patient_cap:
                    continue
                kept = [
                    patient
                    for patient in patients
                    if person_counts[patient] <= variant_cap
                ]
                if not kept:
                    continue
                site = "\t".join(columns)
                vcf.write(f"{site}\n")
                tsv.write(f"{site}\t{';'.join(kept)}\n")
