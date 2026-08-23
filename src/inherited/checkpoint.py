from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


CHECKPOINT_FILENAME = "checkpoint.json"
STATS_CUMULATIVE_FILENAME = "stats_cumulative.json"


def validate_file_suffix(suffix: str) -> str:
    """Require a single path-safe filename token."""
    if not suffix or suffix.strip() != suffix:
        raise ValueError("--file-suffix must be a non-empty filename token")
    if suffix in {".", ".."} or Path(suffix).name != suffix:
        raise ValueError("--file-suffix must not contain a path separator")
    return suffix


def suffixed_filename(stem: str, ext: str, file_suffix: str | None = None) -> str:
    if file_suffix:
        return f"{stem}_{file_suffix}{ext}"
    return f"{stem}{ext}"


@dataclass
class CumulativeStats:
    variants_seen: int = 0
    alleles_tested: int = 0
    inherited_entries: int = 0
    inherited_variants: int = 0
    mendelian_bad_entries: int = 0
    mendelian_bad_variants: int = 0
    denovo_entries: int = 0
    denovo_variants: int = 0
    inherited_per_person: dict[str, int] = field(default_factory=dict)
    denovo_per_person: dict[str, int] = field(default_factory=dict)
    mendelian_bad_per_gt: dict[str, int] = field(default_factory=dict)

    def to_dict(self, *, include_details: bool = True) -> dict[str, Any]:
        data: dict[str, Any] = {
            "variants_seen": self.variants_seen,
            "alleles_tested": self.alleles_tested,
            "inherited_entries": self.inherited_entries,
            "inherited_variants": self.inherited_variants,
            "mendelian_bad_entries": self.mendelian_bad_entries,
            "mendelian_bad_variants": self.mendelian_bad_variants,
            "denovo_entries": self.denovo_entries,
            "denovo_variants": self.denovo_variants,
        }
        if include_details:
            data["inherited_per_person"] = self.inherited_per_person
            data["denovo_per_person"] = self.denovo_per_person
            data["mendelian_bad_per_gt"] = self.mendelian_bad_per_gt
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> CumulativeStats:
        return cls(
            variants_seen=int(data.get("variants_seen", 0)),
            alleles_tested=int(data.get("alleles_tested", 0)),
            inherited_entries=int(data.get("inherited_entries", 0)),
            inherited_variants=int(data.get("inherited_variants", 0)),
            mendelian_bad_entries=int(data.get("mendelian_bad_entries", 0)),
            mendelian_bad_variants=int(data.get("mendelian_bad_variants", 0)),
            denovo_entries=int(data.get("denovo_entries", 0)),
            denovo_variants=int(data.get("denovo_variants", 0)),
            inherited_per_person=dict(data.get("inherited_per_person", {})),
            denovo_per_person=dict(data.get("denovo_per_person", {})),
            mendelian_bad_per_gt=dict(data.get("mendelian_bad_per_gt", {})),
        )


@dataclass
class Checkpoint:
    chrom: str
    last_pos: int
    segment_index: int
    cumulative: CumulativeStats
    completed: bool = False
    details_external: bool = False
    shard_start: int | None = None
    shard_end: int | None = None

    def to_dict(self, *, include_details: bool = True) -> dict[str, Any]:
        data: dict[str, Any] = {
            "chrom": self.chrom,
            "last_pos": self.last_pos,
            "segment_index": self.segment_index,
            "completed": self.completed,
            "details_external": not include_details,
            "cumulative": self.cumulative.to_dict(include_details=include_details),
        }
        if self.shard_start is not None:
            data["shard_start"] = self.shard_start
        if self.shard_end is not None:
            data["shard_end"] = self.shard_end
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Checkpoint:
        shard_start = data.get("shard_start")
        shard_end = data.get("shard_end")
        return cls(
            chrom=str(data.get("chrom", "")),
            last_pos=int(data.get("last_pos", 0)),
            segment_index=int(data.get("segment_index", 0)),
            completed=bool(data.get("completed", False)),
            cumulative=CumulativeStats.from_dict(data.get("cumulative", {})),
            details_external=bool(data.get("details_external", False)),
            shard_start=int(shard_start) if shard_start is not None else None,
            shard_end=int(shard_end) if shard_end is not None else None,
        )


def checkpoint_path(output_dir: Path, file_suffix: str | None = None) -> Path:
    return output_dir / suffixed_filename("checkpoint", ".json", file_suffix)


def load_checkpoint(
    output_dir: Path, file_suffix: str | None = None
) -> Checkpoint | None:
    path = checkpoint_path(output_dir, file_suffix)
    if not path.is_file():
        return None
    with path.open(encoding="utf-8") as handle:
        return Checkpoint.from_dict(json.load(handle))


def save_checkpoint(
    output_dir: Path,
    checkpoint: Checkpoint,
    *,
    include_details: bool = True,
    file_suffix: str | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = checkpoint_path(output_dir, file_suffix)
    tmp = path.with_suffix(".json.tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(
            checkpoint.to_dict(include_details=include_details),
            handle,
            indent=2,
            sort_keys=True,
        )
    tmp.replace(path)


def write_json_atomic(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    tmp.replace(path)
