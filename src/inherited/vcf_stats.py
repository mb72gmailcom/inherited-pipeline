from __future__ import annotations

import gzip
import json
import re
from pathlib import Path
from typing import TextIO

from inherited.alleles import trim_alleles
from inherited.shards import discover_vcf_shards

STATS_KIND = "vcf-stats"
STATS_VERSION = 1
SUMMARY_NAME = "summary.json"
INDEL_STATS_NAME = "indels-stats.json"
LENGTH_DIFF_STATS_NAME = "insertions-deletions-stats.json"

_CHROM_DIR_NAME = re.compile(r"^chr(?:\d+|X|Y|M)$")
_SHARD_VCF_NAME = re.compile(
    r"^(?P<stem>.+)\.(?P<chrom>chr(?:\d+|X|Y|M))_(?P<start>\d+)_(?P<end>\d+)\.vcf(?:\.gz)?$"
)
_WHOLE_VCF_NAME = re.compile(
    r"^(?P<stem>.+)\.(?P<chrom>chr(?:\d+|X|Y|M))\.vcf(?:\.gz)?$"
)
_SHARD_SUMMARY_NAME = re.compile(r"^summary_(\d+)_(\d+)\.json$")

_COUNT_FIELDS = (
    "bi_line_count",
    "multi_line_count",
    "bi_del_count",
    "bi_ins_count",
    "bi_mnv_count_raw",
    "bi_mnv_count",
    "bi_snv_count_raw",
    "bi_snv_count",
    "multi_del_count",
    "multi_ins_count",
    "multi_mnv_count_raw",
    "multi_mnv_count",
    "multi_snv_count_raw",
    "multi_snv_count",
    "symbolic_count",
    "ref_equals_alt",
    "no_alt",
    "skipped_other_chrom",
    "skipped_malformed",
)


def normalize_chrom_label(value: str) -> str:
    """Return ``chr1``, ``chrX``, ``chrY``, or ``chrM``."""
    text = value.strip()
    if not text:
        raise ValueError("--chr must be a chromosome number, X, Y, or M")
    rest = text[3:] if text.lower().startswith("chr") else text
    upper = rest.upper()
    if upper == "MT":
        upper = "M"
    if upper in {"X", "Y", "M"}:
        return f"chr{upper}"
    if rest.isdigit() and int(rest) > 0:
        return f"chr{int(rest)}"
    raise ValueError(f"--chr must be a chromosome number, X, Y, or M, not {value!r}")


def discover_chromosome_vcfs(input_dir: Path, prefix: str, chrom: str) -> list[Path]:
    """Find one whole-chromosome VCF or the position shards for that chromosome.

    Whole file: ``{input_dir}/{prefix}.{chrom}.vcf.gz``.
    Shards: ``{input_dir}/{chrom}/{prefix}.{chrom}_{start}_{end}.vcf.gz``.
    """
    if not prefix or prefix.strip() != prefix or "/" in prefix or "\\" in prefix:
        raise ValueError("--prefix must be a non-empty filename prefix")
    if not input_dir.is_dir():
        raise FileNotFoundError(f"VCF directory not found: {input_dir}")

    whole = _whole_chromosome_files(input_dir, prefix, chrom)
    shards = _chromosome_shard_files(input_dir / chrom, prefix, chrom)
    if whole and shards:
        raise ValueError(
            f"Both a whole-chromosome VCF and position shards match {prefix} for {chrom} "
            f"in {input_dir}"
        )
    if not whole and not shards:
        raise FileNotFoundError(
            f"No VCF for {chrom} matching {prefix}.{{chrom}}.vcf.gz in {input_dir} "
            f"or {prefix}.{{chrom}}_{{start}}_{{end}}.vcf.gz in {input_dir / chrom}"
        )
    return whole or shards


def parse_vcf_filename(path: Path) -> tuple[str, str, tuple[int, int] | None]:
    """Return ``(chrom, prefix, span)``. ``span`` is ``(start, end)`` for a shard."""
    shard = _SHARD_VCF_NAME.fullmatch(path.name)
    if shard is not None:
        start = int(shard.group("start"))
        end = int(shard.group("end"))
        if start > end:
            raise ValueError(f"Invalid shard coordinates in {path.name}: start {start} > end {end}")
        return shard.group("chrom"), shard.group("stem"), (start, end)
    whole = _WHOLE_VCF_NAME.fullmatch(path.name)
    if whole is not None:
        return whole.group("chrom"), whole.group("stem"), None
    raise ValueError(
        f"VCF name must be {{prefix}}.chrN.vcf.gz or "
        f"{{prefix}}.chrN_{{start}}_{{end}}.vcf.gz, not {path.name}"
    )


def collect_vcf_file(path: Path) -> tuple[dict[str, object], tuple[int, int] | None]:
    """Stream one VCF. Sample columns are not read."""
    if not path.is_file():
        raise FileNotFoundError(f"VCF not found: {path}")
    chrom, prefix, span = parse_vcf_filename(path)
    payload = _stats_payload(chrom, prefix, [path], [path])
    if span is not None:
        payload["start"] = span[0]
        payload["end"] = span[1]
    return payload, span


def collect_vcf_stats(input_dir: Path, chrom: str, prefix: str) -> dict[str, object]:
    """Stream one chromosome and return allele counts. Sample columns are not read."""
    label = normalize_chrom_label(chrom)
    paths = discover_chromosome_vcfs(input_dir, prefix, label)
    return _stats_payload(label, prefix, paths, paths)


def aggregate_chromosome_stats(input_dir: Path, chrom: str) -> dict[str, object]:
    """Add ``summary_<start>_<end>.json`` files for one chromosome. Does not read VCFs."""
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Stats directory not found: {input_dir}")
    label = normalize_chrom_label(chrom)
    chrom_dir = input_dir / label
    if not chrom_dir.is_dir():
        raise FileNotFoundError(f"No stats directory for {label}: {chrom_dir}")
    shards = _shard_summaries(chrom_dir)
    totals = {name: 0 for name in _COUNT_FIELDS}
    indel_lengths: dict[int, int] = {}
    length_diffs: dict[int, int] = {}
    files: list[str] = []
    for summary, indel_path, diff_path in shards:
        for name in _COUNT_FIELDS:
            totals[name] += int(summary.get(name, 0))
        _add_distribution(indel_lengths, _read_distribution(indel_path))
        _add_distribution(length_diffs, _read_distribution(diff_path))
        recorded = summary.get("files")
        if isinstance(recorded, list):
            files.extend(str(name) for name in recorded)
    return {
        "kind": STATS_KIND,
        "version": STATS_VERSION,
        "sum": False,
        "chrom": label,
        "files": files,
        **totals,
        "d_indel": _distribution(indel_lengths),
        "d_dl": _distribution(length_diffs),
    }


def sum_vcf_stats(input_dir: Path) -> dict[str, object]:
    """Add per-chromosome stats directories. Does not read VCFs."""
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Stats directory not found: {input_dir}")
    chrom_dirs = sorted(
        (
            path
            for path in input_dir.iterdir()
            if path.is_dir() and _CHROM_DIR_NAME.fullmatch(path.name)
        ),
        key=lambda path: _chrom_sort_key(path.name),
    )
    if not chrom_dirs:
        raise FileNotFoundError(f"No chrN stats directories in {input_dir}")

    totals = {name: 0 for name in _COUNT_FIELDS}
    indel_lengths: dict[int, int] = {}
    length_diffs: dict[int, int] = {}
    chroms: list[str] = []
    for chrom_dir in chrom_dirs:
        summary = _read_summary(chrom_dir)
        chroms.append(str(summary["chrom"]))
        for name in _COUNT_FIELDS:
            totals[name] += int(summary.get(name, 0))
        _add_distribution(indel_lengths, _read_distribution(chrom_dir / INDEL_STATS_NAME))
        _add_distribution(length_diffs, _read_distribution(chrom_dir / LENGTH_DIFF_STATS_NAME))
    return {
        "kind": STATS_KIND,
        "version": STATS_VERSION,
        "sum": True,
        "chroms": chroms,
        **totals,
        "d_indel": _distribution(indel_lengths),
        "d_dl": _distribution(length_diffs),
    }


def write_chromosome_stats(
    output_dir: Path,
    payload: dict[str, object],
    span: tuple[int, int] | None = None,
) -> Path:
    """Write ``output_dir/chrN/`` with the summary and the two indel distributions.

    A shard span adds ``_<start>_<end>`` to each filename. The unsuffixed names are the
    whole-chromosome totals.
    """
    chrom = str(payload["chrom"])
    destination = output_dir / chrom
    destination.mkdir(parents=True, exist_ok=True)
    suffix = "" if span is None else f"_{span[0]}_{span[1]}"
    summary = {key: value for key, value in payload.items() if key not in {"d_indel", "d_dl"}}
    _write_json(destination / f"summary{suffix}.json", summary)
    _write_json(destination / f"indels-stats{suffix}.json", payload["d_indel"])
    _write_json(destination / f"insertions-deletions-stats{suffix}.json", payload["d_dl"])
    return destination


def write_sum_stats(output_dir: Path, payload: dict[str, object]) -> Path:
    """Write the three stats files directly into ``output_dir``."""
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {key: value for key, value in payload.items() if key not in {"d_indel", "d_dl"}}
    _write_json(output_dir / SUMMARY_NAME, summary)
    _write_json(output_dir / INDEL_STATS_NAME, payload["d_indel"])
    _write_json(output_dir / LENGTH_DIFF_STATS_NAME, payload["d_dl"])
    return output_dir


def _whole_chromosome_files(input_dir: Path, prefix: str, chrom: str) -> list[Path]:
    found = [
        path
        for suffix in (".vcf.gz", ".vcf")
        if (path := input_dir / f"{prefix}.{chrom}{suffix}").is_file()
    ]
    if len(found) > 1:
        names = ", ".join(path.name for path in found)
        raise ValueError(f"Both compressed and plain VCFs are present for {chrom}: {names}")
    return found


def _chromosome_shard_files(shard_dir: Path, prefix: str, chrom: str) -> list[Path]:
    if not shard_dir.is_dir():
        return []
    try:
        shards = discover_vcf_shards(shard_dir, prefix)
    except FileNotFoundError:
        return []
    foreign = sorted({shard.chrom for shard in shards if shard.chrom != chrom})
    if foreign:
        raise ValueError(
            f"Shard directory {shard_dir} is for {chrom} but contains {foreign}"
        )
    return [shard.path for shard in shards]


def _accumulate_vcf(
    path: Path,
    chrom: str,
    counts: dict[str, int],
    indel_lengths: dict[int, int],
    length_diffs: dict[int, int],
) -> None:
    with _open_text(path) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            parsed = _parse_site(line)
            if parsed is None:
                counts["skipped_malformed"] += 1
                continue
            site_chrom, pos, ref, alt_field = parsed
            if site_chrom != chrom:
                counts["skipped_other_chrom"] += 1
                continue
            alleles = [allele.strip() for allele in alt_field.split(",") if allele.strip()]
            if not alleles or alleles == ["."]:
                counts["no_alt"] += 1
                continue
            prefix = "bi" if len(alleles) == 1 else "multi"
            counts[f"{prefix}_line_count"] += 1
            for alt in alleles:
                if _is_symbolic(alt):
                    counts["symbolic_count"] += 1
                    continue
                _count_allele(
                    pos, ref, alt, prefix, counts, indel_lengths, length_diffs
                )


def _parse_site(line: str) -> tuple[str, int, str, str] | None:
    parts = line.rstrip("\n").split("\t", 4)
    if len(parts) < 5 or not parts[3]:
        return None
    alt_field = parts[4].split("\t", 1)[0]
    if not alt_field:
        return None
    try:
        pos = int(parts[1])
    except ValueError:
        return None
    try:
        chrom = normalize_chrom_label(parts[0])
    except ValueError:
        return None
    return chrom, pos, parts[3], alt_field


def _count_allele(
    pos: int,
    ref: str,
    alt: str,
    group: str,
    counts: dict[str, int],
    indel_lengths: dict[int, int],
    length_diffs: dict[int, int],
) -> None:
    trimmed = trim_alleles(pos, ref, alt)
    if trimmed is None:
        counts["ref_equals_alt"] += 1
        return
    _add_class(counts, f"{group}_", _classify(ref, alt), raw=True)
    _, trimmed_ref, trimmed_alt = trimmed
    kind = _classify(trimmed_ref, trimmed_alt)
    _add_class(counts, f"{group}_", kind, raw=False)
    if kind in {"del", "ins"}:
        longer = max(len(trimmed_ref), len(trimmed_alt))
        diff = len(trimmed_ref) - len(trimmed_alt)
        indel_lengths[longer] = indel_lengths.get(longer, 0) + 1
        length_diffs[diff] = length_diffs.get(diff, 0) + 1


def _add_class(counts: dict[str, int], prefix: str, kind: str, *, raw: bool) -> None:
    if raw:
        if kind in {"snv", "mnv"}:
            counts[f"{prefix}{kind}_count_raw"] += 1
        return
    counts[f"{prefix}{kind}_count"] += 1


def _classify(ref: str, alt: str) -> str:
    ref_length = len(ref)
    alt_length = len(alt)
    if ref_length == alt_length == 1:
        return "snv"
    if ref_length == alt_length:
        return "mnv"
    if ref_length > alt_length:
        return "del"
    return "ins"


def _is_symbolic(alt: str) -> bool:
    return alt == "*" or alt.startswith("<") or "[" in alt or "]" in alt


def _chrom_sort_key(name: str) -> tuple[int, int | str]:
    rest = name[3:]
    if rest.isdigit():
        return (0, int(rest))
    return (1, rest)


def _stats_payload(
    chrom: str,
    prefix: str,
    record_paths: list[Path],
    vcf_paths: list[Path],
) -> dict[str, object]:
    counts = {name: 0 for name in _COUNT_FIELDS}
    indel_lengths: dict[int, int] = {}
    length_diffs: dict[int, int] = {}
    for path in vcf_paths:
        _accumulate_vcf(path, chrom, counts, indel_lengths, length_diffs)
    return {
        "kind": STATS_KIND,
        "version": STATS_VERSION,
        "sum": False,
        "chrom": chrom,
        "prefix": prefix,
        "files": [path.name for path in record_paths],
        **counts,
        "d_indel": _distribution(indel_lengths),
        "d_dl": _distribution(length_diffs),
    }


def _shard_summaries(chrom_dir: Path) -> list[tuple[dict[str, object], Path, Path]]:
    found: list[tuple[int, int, Path]] = []
    for path in chrom_dir.iterdir():
        if not path.is_file():
            continue
        match = _SHARD_SUMMARY_NAME.fullmatch(path.name)
        if match is None:
            continue
        found.append((int(match.group(1)), int(match.group(2)), path))
    found.sort()
    if not found:
        raise FileNotFoundError(f"No shard summary_<start>_<end>.json files in {chrom_dir}")
    for prev, current in zip(found, found[1:]):
        if current[0] <= prev[1]:
            raise ValueError(
                f"Overlapping shard summaries in {chrom_dir}: "
                f"{prev[0]}_{prev[1]} and {current[0]}_{current[1]}"
            )
    rows: list[tuple[dict[str, object], Path, Path]] = []
    for start, end, summary_path in found:
        summary = _read_summary_file(summary_path, expected_chrom=chrom_dir.name)
        indel_path = chrom_dir / f"indels-stats_{start}_{end}.json"
        diff_path = chrom_dir / f"insertions-deletions-stats_{start}_{end}.json"
        rows.append((summary, indel_path, diff_path))
    return rows


def _read_summary(chrom_dir: Path) -> dict[str, object]:
    path = chrom_dir / SUMMARY_NAME
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {SUMMARY_NAME} in {chrom_dir}. "
            "Aggregate shard files with --chr before --sum."
        )
    return _read_summary_file(path, expected_chrom=chrom_dir.name)


def _read_summary_file(path: Path, *, expected_chrom: str) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("kind") != STATS_KIND:
        raise ValueError(f"{path} is not a vcf-stats summary")
    if payload.get("version") != STATS_VERSION:
        raise ValueError(f"Unsupported vcf-stats version in {path}: {payload.get('version')}")
    if payload.get("chrom") != expected_chrom:
        raise ValueError(
            f"{path} records chromosome {payload.get('chrom')!r}, not {expected_chrom}"
        )
    return payload


def _read_distribution(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path.name} in {path.parent}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must be an object of length counts")
    return payload


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _distribution(counts: dict[int, int]) -> dict[str, int]:
    return {str(key): counts[key] for key in sorted(counts)}


def _add_distribution(total: dict[int, int], payload: object) -> None:
    if not isinstance(payload, dict):
        raise ValueError("Indel distribution must be an object")
    for key, value in payload.items():
        length = int(key)
        total[length] = total.get(length, 0) + int(value)


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")
