from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import TextIO

from inherited.shards import discover_vcf_shards

STATS_KIND = "vcf-stats"
STATS_VERSION = 1

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


def collect_vcf_stats(input_dir: Path, chrom: str, prefix: str) -> dict[str, object]:
    """Stream one chromosome and return allele counts. Sample columns are not read."""
    label = normalize_chrom_label(chrom)
    paths = discover_chromosome_vcfs(input_dir, prefix, label)
    counts = {name: 0 for name in _COUNT_FIELDS}
    indel_lengths: dict[int, int] = {}
    length_diffs: dict[int, int] = {}
    for path in paths:
        _accumulate_vcf(path, label, counts, indel_lengths, length_diffs)
    return {
        "kind": STATS_KIND,
        "version": STATS_VERSION,
        "sum": False,
        "chrom": label,
        "prefix": prefix,
        "files": [path.name for path in paths],
        **counts,
        "d_indel": _distribution(indel_lengths),
        "d_dl": _distribution(length_diffs),
    }


def sum_vcf_stats(input_dir: Path) -> dict[str, object]:
    """Add per-chromosome ``vcf-stats`` JSON files. Does not read VCFs."""
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Stats directory not found: {input_dir}")
    payloads = []
    for path in sorted(input_dir.iterdir()):
        if not path.is_file() or path.suffix != ".json":
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict) or payload.get("kind") != STATS_KIND:
            continue
        if payload.get("sum") is True or "chrom" not in payload:
            continue
        if payload.get("version") != STATS_VERSION:
            raise ValueError(
                f"Unsupported vcf-stats version in {path.name}: {payload.get('version')}"
            )
        payloads.append(payload)
    if not payloads:
        raise FileNotFoundError(f"No per-chromosome vcf-stats JSON files in {input_dir}")

    totals = {name: 0 for name in _COUNT_FIELDS}
    indel_lengths: dict[int, int] = {}
    length_diffs: dict[int, int] = {}
    chroms: list[str] = []
    for payload in payloads:
        chroms.append(str(payload["chrom"]))
        for name in _COUNT_FIELDS:
            totals[name] += int(payload.get(name, 0))
        _add_distribution(indel_lengths, payload.get("d_indel", {}))
        _add_distribution(length_diffs, payload.get("d_dl", {}))
    return {
        "kind": STATS_KIND,
        "version": STATS_VERSION,
        "sum": True,
        "chroms": chroms,
        **totals,
        "d_indel": _distribution(indel_lengths),
        "d_dl": _distribution(length_diffs),
    }


def write_vcf_stats(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def trim_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str] | None:
    """Trim a shared prefix and suffix. A prefix trim shifts pos.

    Returns None when the remaining REF and ALT are identical.
    """
    ref_bases = ref.upper()
    alt_bases = alt.upper()
    while len(ref_bases) > 1 and len(alt_bases) > 1 and ref_bases[-1] == alt_bases[-1]:
        ref_bases = ref_bases[:-1]
        alt_bases = alt_bases[:-1]
    while len(ref_bases) > 1 and len(alt_bases) > 1 and ref_bases[0] == alt_bases[0]:
        ref_bases = ref_bases[1:]
        alt_bases = alt_bases[1:]
        pos += 1
    if ref_bases == alt_bases:
        return None
    return pos, ref_bases, alt_bases


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
