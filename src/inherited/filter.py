"""Filter result TSVs by how common a variant is and how many variants a person carries."""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from inherited.params import write_params

PREFIXES = ("inherited", "denovo", "mendelian_bad")

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")

TSV_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
_TRIO_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tTRIO_CALLS\n"
_PATIENTS = "PATIENTS"
_TRIO_CALLS = "TRIO_CALLS"
_COLUMNS = {_PATIENTS, _TRIO_CALLS}
_HIST_HEADER = "count\tn\n"
_CAP_DISABLED = "disabled"


def filter_results(
    input_dir: Path,
    output_dir: Path,
    prefix: str,
    *,
    variant_cap: int | None = None,
    patient_cap: int | None = None,
) -> Path | None:
    """Count variants and, when a cap is set, write filtered VCF and TSV trees.

    The last column may be ``PATIENTS`` or ``TRIO_CALLS``. A directory that
    mixes the two stops before any output is written. Returns the ``-sites``
    directory when at least one cap is applied, otherwise ``None``.
    """
    if prefix not in PREFIXES:
        raise ValueError(
            f"prefix must be one of {', '.join(PREFIXES)}, got {prefix!r}"
        )
    _validate_caps(variant_cap, patient_cap)
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")

    tsv_files = _discover_tsv_files(input_dir, prefix)
    person_counts, carrier_counts, column = _count_variants(input_dir, tsv_files)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_histogram(
        output_dir / f"{prefix}_variants_per_person_hist.tsv",
        Counter(person_counts.values()),
    )
    _write_histogram(
        output_dir / f"{prefix}_patients_per_variant_hist.tsv",
        carrier_counts,
    )
    sites_dir = None
    if variant_cap is not None or patient_cap is not None:
        sites_dir = output_dir.with_name(output_dir.name + "-sites")
        sites_dir.mkdir(parents=True, exist_ok=True)
        _write_filtered(
            tsv_files,
            output_dir,
            sites_dir,
            person_counts,
            prefix,
            column,
            variant_cap,
            patient_cap,
        )
    _write_params(
        input_dir,
        output_dir,
        prefix,
        variant_cap=variant_cap,
        patient_cap=patient_cap,
    )
    return sites_dir


def _validate_caps(variant_cap: int | None, patient_cap: int | None) -> None:
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


def _count_variants(
    input_dir: Path,
    tsv_files: list[tuple[str, Path]],
) -> tuple[Counter[str], Counter[int], str]:
    person_counts: Counter[str] = Counter()
    carrier_counts: Counter[int] = Counter()
    chosen: str | None = None
    chosen_path: Path | None = None
    for _chrom, path in tsv_files:
        column, rows = _iter_calls(path)
        for _columns, calls in rows:
            carrier_counts[len(calls)] += 1
            person_counts.update(person for person, _piece in calls)
        if column[0] is None:
            continue
        if chosen is None:
            chosen = column[0]
            chosen_path = path
            continue
        if column[0] != chosen:
            assert chosen_path is not None
            raise ValueError(
                _mixed_message(input_dir, chosen, chosen_path, column[0], path)
            )
    return person_counts, carrier_counts, chosen or _PATIENTS


def _iter_calls(path: Path):
    column: list[str | None] = [None]

    def rows():
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                line = line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("#"):
                    last = line.split("\t")[-1]
                    if last in _COLUMNS:
                        _set_column(column, path, last)
                    continue
                fields = line.split("\t")
                if len(fields) < 6:
                    raise ValueError(f"{path} expected 6 tab-separated columns")
                payload = fields[5]
                if payload:
                    row_column = _TRIO_CALLS if "=" in payload else _PATIENTS
                    _set_column(column, path, row_column)
                elif column[0] is None:
                    column[0] = _PATIENTS
                if column[0] == _TRIO_CALLS:
                    calls = _parse_trio_calls(path, payload)
                else:
                    calls = [(person, person) for person in payload.split(";") if person]
                yield fields[:5], calls

    return column, rows()


def _set_column(column: list[str | None], path: Path, found: str) -> None:
    current = column[0]
    if current is None:
        column[0] = found
        return
    if current != found:
        raise ValueError(f"{path} mixes PATIENTS and TRIO_CALLS")


def _parse_trio_calls(path: Path, payload: str) -> list[tuple[str, str]]:
    calls: list[tuple[str, str]] = []
    for part in payload.split(";"):
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"{path} TRIO_CALLS entry {part!r} is missing '='")
        person, _genotypes = part.split("=", 1)
        if not person:
            raise ValueError(f"{path} TRIO_CALLS entry is missing a child id")
        calls.append((person, part))
    return calls


def _mixed_message(
    input_dir: Path,
    left_column: str,
    left_path: Path,
    right_column: str,
    right_path: Path,
) -> str:
    if left_column != _PATIENTS:
        left_column, left_path, right_column, right_path = (
            right_column,
            right_path,
            left_column,
            left_path,
        )
    return (
        f"mixed last columns in {input_dir}: "
        f"{left_path.relative_to(input_dir).as_posix()} has {left_column} and "
        f"{right_path.relative_to(input_dir).as_posix()} has {right_column}"
    )


def _write_params(
    input_dir: Path,
    output_dir: Path,
    prefix: str,
    *,
    variant_cap: int | None,
    patient_cap: int | None,
) -> Path:
    payload = {
        "command": "filter",
        "input_dir": str(input_dir.resolve()),
        "output_dir": str(output_dir.resolve()),
        "patient_cap": _cap_param(patient_cap),
        "prefix": prefix,
        "variant_cap": _cap_param(variant_cap),
    }
    return write_params(output_dir / "params.json", payload)


def _cap_param(value: int | None) -> int | str:
    if value is None:
        return _CAP_DISABLED
    return value


def _write_histogram(path: Path, counts: Counter[int]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(_HIST_HEADER)
        for count in sorted(counts):
            handle.write(f"{count}\t{counts[count]}\n")


def _vcf_name(stem: str, prefix: str) -> str:
    return f"variants{stem.removeprefix(prefix)}.vcf"


def _write_filtered(
    tsv_files: list[tuple[str, Path]],
    output_dir: Path,
    sites_dir: Path,
    person_counts: Counter[str],
    prefix: str,
    column: str,
    variant_cap: int | None,
    patient_cap: int | None,
) -> None:
    header = TSV_HEADER if column == _PATIENTS else _TRIO_HEADER
    for chrom, path in tsv_files:
        vcf_path = output_dir / chrom / _vcf_name(path.stem, prefix)
        tsv_path = sites_dir / chrom / path.name
        vcf_path.parent.mkdir(parents=True, exist_ok=True)
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        wrote = False
        with vcf_path.open("w", encoding="utf-8") as vcf, tsv_path.open(
            "w", encoding="utf-8"
        ) as tsv:
            _column, rows = _iter_calls(path)
            for columns, calls in rows:
                if patient_cap is not None and len(calls) > patient_cap:
                    continue
                if variant_cap is None:
                    kept = calls
                else:
                    kept = [
                        call
                        for call in calls
                        if person_counts[call[0]] <= variant_cap
                    ]
                if not kept:
                    continue
                if not wrote:
                    tsv.write(header)
                    wrote = True
                site = "\t".join(columns)
                vcf.write(f"{site}\n")
                tsv.write(f"{site}\t{';'.join(piece for _person, piece in kept)}\n")
        if not wrote:
            vcf_path.unlink()
            tsv_path.unlink()
            _remove_empty_dir(vcf_path.parent)
            _remove_empty_dir(tsv_path.parent)
    _remove_empty_dir(sites_dir)


def _remove_empty_dir(path: Path) -> None:
    if path.is_dir() and not any(path.iterdir()):
        path.rmdir()
