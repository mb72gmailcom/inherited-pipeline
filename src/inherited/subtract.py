"""Write variants present in one VCF tree and absent from the paired file in another."""

from __future__ import annotations

import re
from pathlib import Path

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
_VCF_NAME = re.compile(r"^variants_\d+_\d+\.vcf$")


def subtract_variants(input_dir1: Path, input_dir2: Path, output_dir: Path) -> None:
    """Subtract paired ``variants_${start}_${end}.vcf`` files.

    A row from the first directory is written when its ``CHROM``, ``POS``,
    ``REF``, and ``ALT`` are absent from the same filename in the second
    directory. Files that exist only in the first directory are copied.
    """
    if not input_dir1.is_dir():
        raise ValueError(f"input directory not found: {input_dir1}")
    if not input_dir2.is_dir():
        raise ValueError(f"input directory not found: {input_dir2}")

    for chrom_dir in _chrom_dirs(input_dir1):
        for vcf_path in _variant_files(chrom_dir):
            paired = input_dir2 / chrom_dir.name / vcf_path.name
            blocked = _variant_keys(paired) if paired.is_file() else set()
            destination = output_dir / chrom_dir.name / vcf_path.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            _write_difference(vcf_path, destination, blocked)


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


def _write_difference(
    source: Path,
    destination: Path,
    blocked: set[tuple[str, str, str, str]],
) -> None:
    with destination.open("w", encoding="utf-8") as handle:
        for key, line in _iter_vcf_rows(source):
            if key in blocked:
                continue
            handle.write(f"{line}\n")
