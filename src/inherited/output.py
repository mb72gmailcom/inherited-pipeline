from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from inherited.checkpoint import (
    Checkpoint,
    CumulativeStats,
    save_checkpoint,
    suffixed_filename,
    validate_file_suffix,
    write_json_atomic,
)
from inherited.constants import DEFAULT_BLOCK_SIZE, DEFAULT_SEGMENT_SIZE
from inherited.xchrom import (
    CHROM_MODE_AUTOSOMAL,
    CHROM_MODE_X,
    CHROM_MODE_Y,
    X_BUCKETS,
    Y_BUCKETS,
)

HitRecord = tuple[str, ...]


def serialize_patient_ids(hits: dict[str, HitRecord]) -> str:
    """Serialize affected patient IDs only."""
    return ";".join(sorted(hits))


def serialize_trio_calls(hits: dict[str, HitRecord]) -> str:
    """Serialize genotype payloads for one or more children.

    Trio records are ``(mGT, fGT, cGT, cGQ)``.
    Mother-child pair records (male nonPAR) are ``(mGT, cGT, cGQ)``.
    """
    parts = []
    for person_id in sorted(hits):
        record = hits[person_id]
        if len(record) == 3:
            mother_gt, child_gt, child_gq = record
            parts.append(f"{person_id}={mother_gt}|{child_gt}|{child_gq}")
        else:
            mother_gt, father_gt, child_gt, child_gq = record
            parts.append(
                f"{person_id}={mother_gt}|{father_gt}|{child_gt}|{child_gq}"
            )
    return ";".join(parts)


def parse_trio_calls(payload: str) -> dict[str, HitRecord]:
    hits: dict[str, HitRecord] = {}
    if not payload:
        return hits
    for part in payload.split(";"):
        person_id, genotypes = part.split("=", 1)
        fields = genotypes.split("|")
        hits[person_id] = tuple(fields)
    return hits


def serialize_payload(
    hits: dict[str, HitRecord],
    *,
    short_format: bool,
) -> str:
    if short_format:
        return serialize_patient_ids(hits)
    return serialize_trio_calls(hits)


def parse_patient_ids(payload: str) -> list[str]:
    if not payload:
        return []
    return payload.split(";")


def _gt_key(record: HitRecord) -> str:
    if len(record) == 3:
        mother_gt, child_gt, _child_gq = record
        return f"{mother_gt}:{child_gt}"
    mother_gt, father_gt, child_gt, _child_gq = record
    return f"{mother_gt}:{father_gt}:{child_gt}"


class JsonObjectStreamWriter:
    """Write a JSON object incrementally without retaining its entries."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.entries_written = 0
        self._handle = path.open("w", encoding="utf-8")
        self._handle.write("{")

    def append(self, key: str, value: int) -> None:
        if self.entries_written:
            self._handle.write("\n,")
        else:
            self._handle.write("\n")
        self._handle.write(f"{json.dumps(key)}: {json.dumps(value)}")
        self.entries_written += 1

    def append_encoded(self, entry: str) -> None:
        """Append one already JSON-encoded ``key: value`` entry."""
        if self.entries_written:
            self._handle.write("\n,")
        else:
            self._handle.write("\n")
        self._handle.write(entry)
        self.entries_written += 1

    def close(self) -> None:
        if self._handle.closed:
            return
        if self.entries_written:
            self._handle.write("\n")
        self._handle.write("}\n")
        self._handle.close()


class BlockWriter:
    """Buffer tab-separated result lines and flush to disk in blocks."""

    HEADER_SHORT = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
    HEADER_FULL = "#CHROM\tPOS\tID\tREF\tALT\tTRIO_CALLS\n"

    def __init__(
        self,
        path: Path,
        block_size: int = DEFAULT_BLOCK_SIZE,
        *,
        short_format: bool = True,
    ) -> None:
        self.path = path
        self.block_size = block_size
        self.block: list[str] = []
        self.lines_written = 0
        self._handle = path.open("w", encoding="utf-8")
        header = self.HEADER_SHORT if short_format else self.HEADER_FULL
        self._handle.write(header)

    def append(self, chrom: str, pos: str, ref: str, alt: str, payload: str) -> None:
        self.block.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{payload}")
        if len(self.block) >= self.block_size:
            self.flush()

    def flush(self) -> None:
        if not self.block:
            return
        self._handle.write("\n".join(self.block) + "\n")
        self.lines_written += len(self.block)
        self.block.clear()

    def close(self) -> None:
        self.flush()
        self._handle.close()


@dataclass
class BucketStats:
    inherited_entries: int = 0
    inherited_variants: int = 0
    mendelian_bad_entries: int = 0
    mendelian_bad_variants: int = 0
    denovo_entries: int = 0
    denovo_variants: int = 0
    inherited_per_person: dict[str, int] = field(default_factory=dict)
    denovo_per_person: dict[str, int] = field(default_factory=dict)
    mendelian_bad_per_gt: dict[str, int] = field(default_factory=dict)
    detail_inherited_entries: int = 0
    detail_inherited_variants: int = 0
    detail_mendelian_bad_entries: int = 0
    detail_mendelian_bad_variants: int = 0
    detail_denovo_entries: int = 0
    detail_denovo_variants: int = 0
    detail_inherited_per_person: dict[str, int] = field(default_factory=dict)
    detail_denovo_per_person: dict[str, int] = field(default_factory=dict)
    detail_mendelian_bad_per_gt: dict[str, int] = field(default_factory=dict)


@dataclass
class ResultWriter:
    output_dir: Path
    block_size: int = DEFAULT_BLOCK_SIZE
    segment_size: int = DEFAULT_SEGMENT_SIZE
    short_format: bool = True
    chrom_mode: str = CHROM_MODE_AUTOSOMAL
    cumulative: CumulativeStats = field(default_factory=CumulativeStats)
    segment_index: int = 0
    last_chrom: str = ""
    last_pos: int = 0
    inherited_segment_lines: int = 0
    mendelian_bad_segment_lines: int = 0
    denovo_segment_lines: int = 0
    detail_inherited_per_person: dict[str, int] = field(default_factory=dict)
    detail_denovo_per_person: dict[str, int] = field(default_factory=dict)
    detail_mendelian_bad_per_gt: dict[str, int] = field(default_factory=dict)
    bucket_stats: dict[str, BucketStats] = field(default_factory=dict)
    resume_mode: bool = False
    shard_mode: bool = False
    shard_start: int | None = None
    shard_end: int | None = None
    file_suffix: str | None = None
    _writers_ready: bool = False
    _inherited: dict[str, BlockWriter] = field(default_factory=dict, repr=False)
    _mendelian_bad: dict[str, BlockWriter] = field(default_factory=dict, repr=False)
    _denovo: dict[str, BlockWriter] = field(default_factory=dict, repr=False)
    _inherited_per_variant: dict[str, JsonObjectStreamWriter] = field(
        default_factory=dict, repr=False
    )
    _denovo_per_variant: dict[str, JsonObjectStreamWriter] = field(
        default_factory=dict, repr=False
    )

    def __post_init__(self) -> None:
        if self.file_suffix is not None:
            self.file_suffix = validate_file_suffix(self.file_suffix)
            self.segment_size = 0
            if self.shard_mode:
                raise ValueError("--file-suffix cannot be used with --vcf-dir")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        if not self.resume_mode:
            self._detail_deltas_path().write_text("", encoding="utf-8")
        self._init_bucket_stats()
        if self.resume_mode and not self.shard_mode:
            self._open_segment_writers()
            self._writers_ready = True

    @classmethod
    def from_checkpoint(
        cls,
        output_dir: Path,
        checkpoint: Checkpoint,
        *,
        block_size: int = DEFAULT_BLOCK_SIZE,
        segment_size: int = DEFAULT_SEGMENT_SIZE,
        short_format: bool = True,
        chrom_mode: str = CHROM_MODE_AUTOSOMAL,
        shard_mode: bool = False,
        file_suffix: str | None = None,
    ) -> ResultWriter:
        writer = cls(
            output_dir=output_dir,
            block_size=block_size,
            segment_size=segment_size,
            short_format=short_format,
            chrom_mode=chrom_mode,
            cumulative=checkpoint.cumulative,
            segment_index=checkpoint.segment_index + 1,
            last_chrom=checkpoint.chrom,
            last_pos=checkpoint.last_pos,
            resume_mode=True,
            shard_mode=shard_mode,
            file_suffix=file_suffix,
        )
        if checkpoint.details_external:
            writer._load_detail_deltas(checkpoint.segment_index)
        else:
            writer._write_detail_baseline(checkpoint.segment_index)
        return writer

    def _init_bucket_stats(self) -> None:
        for bucket in self._bucket_names():
            if bucket:
                self.bucket_stats.setdefault(bucket, BucketStats())

    def set_chrom_mode(self, mode: str) -> None:
        """Set autosomal / chrX / chrY output layout before the first write."""
        if mode not in (CHROM_MODE_AUTOSOMAL, CHROM_MODE_X, CHROM_MODE_Y):
            raise ValueError(f"Unknown chromosome mode: {mode!r}")
        if self._writers_ready:
            if mode != self.chrom_mode:
                raise ValueError("Cannot change chromosome mode after writers are open")
            return
        self.chrom_mode = mode
        self._init_bucket_stats()
        self._open_segment_writers()
        self._writers_ready = True

    def set_x_mode(self, enabled: bool) -> None:
        """Backward-compatible wrapper around :meth:`set_chrom_mode`."""
        self.set_chrom_mode(CHROM_MODE_X if enabled else CHROM_MODE_AUTOSOMAL)

    @property
    def uses_named_buckets(self) -> bool:
        return self.chrom_mode in (CHROM_MODE_X, CHROM_MODE_Y)

    @property
    def x_mode(self) -> bool:
        return self.chrom_mode == CHROM_MODE_X

    @property
    def y_mode(self) -> bool:
        return self.chrom_mode == CHROM_MODE_Y

    def _bucket_names(self) -> tuple[str, ...]:
        if self.chrom_mode == CHROM_MODE_X:
            return X_BUCKETS
        if self.chrom_mode == CHROM_MODE_Y:
            return Y_BUCKETS
        return ("",)

    def _kind_stem(self, kind: str, bucket: str) -> str:
        if self.uses_named_buckets:
            return f"{kind}_{bucket}"
        return kind

    def _output_path(self, stem: str, ext: str, *, hidden: bool = False) -> Path:
        name = suffixed_filename(stem, ext, self.file_suffix)
        if hidden:
            name = f".{name}"
        return self.output_dir / name

    def _result_path(self, kind: str, bucket: str, segment_index: int) -> Path:
        stem = self._kind_stem(kind, bucket)
        if self.shard_mode:
            return self.output_dir / f"{stem}_{self.shard_start}_{self.shard_end}.tsv"
        if self.segment_size <= 0:
            return self._output_path(stem, ".tsv")
        return self.output_dir / f"{stem}_{segment_index:05d}.tsv"

    def _per_variant_path(
        self, kind: str, bucket: str, segment_index: int
    ) -> Path:
        stem = self._kind_stem(f"{kind}_per_variant", bucket)
        if self.shard_mode:
            return self.output_dir / f"{stem}_{self.shard_start}_{self.shard_end}.json"
        if self.segment_size <= 0:
            return self._output_path(stem, ".json.part", hidden=True)
        return self.output_dir / f"{stem}_seg{segment_index:05d}.json"

    def _inherited_per_variant_path(self, bucket: str, segment_index: int) -> Path:
        return self._per_variant_path("inherited", bucket, segment_index)

    def _denovo_per_variant_path(self, bucket: str, segment_index: int) -> Path:
        return self._per_variant_path("denovo", bucket, segment_index)

    def _detail_deltas_path(self) -> Path:
        return self._output_path("cumulative_detail_deltas", ".jsonl")

    def _close_writers(self) -> None:
        for writer in self._inherited.values():
            writer.close()
        for writer in self._mendelian_bad.values():
            writer.close()
        for writer in self._denovo.values():
            writer.close()
        for writer in self._inherited_per_variant.values():
            writer.close()
        for writer in self._denovo_per_variant.values():
            writer.close()
        self._inherited.clear()
        self._mendelian_bad.clear()
        self._denovo.clear()
        self._inherited_per_variant.clear()
        self._denovo_per_variant.clear()

    def _open_segment_writers(self) -> None:
        if self.shard_mode and (self.shard_start is None or self.shard_end is None):
            raise ValueError("begin_shard() must be called before opening shard writers")
        self._close_writers()
        for bucket in self._bucket_names():
            self._inherited[bucket] = BlockWriter(
                self._result_path("inherited", bucket, self.segment_index),
                self.block_size,
                short_format=self.short_format,
            )
            self._mendelian_bad[bucket] = BlockWriter(
                self._result_path("mendelian_bad", bucket, self.segment_index),
                self.block_size,
                short_format=self.short_format,
            )
            self._denovo[bucket] = BlockWriter(
                self._result_path("denovo", bucket, self.segment_index),
                self.block_size,
                short_format=self.short_format,
            )
            self._inherited_per_variant[bucket] = JsonObjectStreamWriter(
                self._inherited_per_variant_path(bucket, self.segment_index)
            )
            self._denovo_per_variant[bucket] = JsonObjectStreamWriter(
                self._denovo_per_variant_path(bucket, self.segment_index)
            )
        self.inherited_segment_lines = 0
        self.mendelian_bad_segment_lines = 0
        self.denovo_segment_lines = 0

    def begin_shard(self, start: int, end: int) -> None:
        """Rotate output files to a genomic shard labeled ``{start}_{end}``."""
        if start > end:
            raise ValueError(f"Invalid shard coordinates: start {start} > end {end}")
        if self._writers_ready:
            self._finish_segment(open_next=False)
        self.shard_start = start
        self.shard_end = end
        if self._writers_ready:
            self._open_segment_writers()

    def write_inherited(
        self,
        chrom: str,
        pos: str,
        ref: str,
        alt: str,
        variant_key: str,
        hits: dict[str, HitRecord],
        *,
        bucket: str = "",
    ) -> None:
        if not hits:
            return
        if not self._writers_ready:
            self.set_chrom_mode(CHROM_MODE_AUTOSOMAL)
        if bucket not in self._inherited:
            raise KeyError(f"Unknown output bucket: {bucket!r}")
        self._inherited[bucket].append(
            chrom, pos, ref, alt, serialize_payload(hits, short_format=self.short_format)
        )
        self.inherited_segment_lines += 1
        self.cumulative.inherited_variants += 1
        self.cumulative.inherited_entries += len(hits)
        self._inherited_per_variant[bucket].append(variant_key, len(hits))

        bucket_state = self.bucket_stats.get(bucket)
        if bucket_state is not None:
            bucket_state.inherited_variants += 1
            bucket_state.inherited_entries += len(hits)
            bucket_state.detail_inherited_variants += 1
            bucket_state.detail_inherited_entries += len(hits)

        for person_id in hits:
            self.cumulative.inherited_per_person[person_id] = (
                self.cumulative.inherited_per_person.get(person_id, 0) + 1
            )
            self.detail_inherited_per_person[person_id] = (
                self.detail_inherited_per_person.get(person_id, 0) + 1
            )
            if bucket_state is not None:
                bucket_state.inherited_per_person[person_id] = (
                    bucket_state.inherited_per_person.get(person_id, 0) + 1
                )
                bucket_state.detail_inherited_per_person[person_id] = (
                    bucket_state.detail_inherited_per_person.get(person_id, 0) + 1
                )
        self._update_position(chrom, pos)
        self._maybe_rotate_segment()

    def write_mendelian_bad(
        self,
        chrom: str,
        pos: str,
        ref: str,
        alt: str,
        hits: dict[str, HitRecord],
        *,
        bucket: str = "",
    ) -> None:
        if not hits:
            return
        if not self._writers_ready:
            self.set_chrom_mode(CHROM_MODE_AUTOSOMAL)
        if bucket not in self._mendelian_bad:
            raise KeyError(f"Unknown output bucket: {bucket!r}")
        self._mendelian_bad[bucket].append(
            chrom, pos, ref, alt, serialize_payload(hits, short_format=self.short_format)
        )
        self.mendelian_bad_segment_lines += 1
        self.cumulative.mendelian_bad_variants += 1
        self.cumulative.mendelian_bad_entries += len(hits)

        bucket_state = self.bucket_stats.get(bucket)
        if bucket_state is not None:
            bucket_state.mendelian_bad_variants += 1
            bucket_state.mendelian_bad_entries += len(hits)
            bucket_state.detail_mendelian_bad_variants += 1
            bucket_state.detail_mendelian_bad_entries += len(hits)

        for record in hits.values():
            gt_key = _gt_key(record)
            self.cumulative.mendelian_bad_per_gt[gt_key] = (
                self.cumulative.mendelian_bad_per_gt.get(gt_key, 0) + 1
            )
            self.detail_mendelian_bad_per_gt[gt_key] = (
                self.detail_mendelian_bad_per_gt.get(gt_key, 0) + 1
            )
            if bucket_state is not None:
                bucket_state.mendelian_bad_per_gt[gt_key] = (
                    bucket_state.mendelian_bad_per_gt.get(gt_key, 0) + 1
                )
                bucket_state.detail_mendelian_bad_per_gt[gt_key] = (
                    bucket_state.detail_mendelian_bad_per_gt.get(gt_key, 0) + 1
                )
        self._update_position(chrom, pos)
        self._maybe_rotate_segment()

    def write_denovo(
        self,
        chrom: str,
        pos: str,
        ref: str,
        alt: str,
        variant_key: str,
        hits: dict[str, HitRecord],
        *,
        bucket: str = "",
    ) -> None:
        if not hits:
            return
        if not self._writers_ready:
            self.set_chrom_mode(CHROM_MODE_AUTOSOMAL)
        if bucket not in self._denovo:
            raise KeyError(f"Unknown output bucket: {bucket!r}")
        self._denovo[bucket].append(
            chrom, pos, ref, alt, serialize_payload(hits, short_format=self.short_format)
        )
        self.denovo_segment_lines += 1
        self.cumulative.denovo_variants += 1
        self.cumulative.denovo_entries += len(hits)
        self._denovo_per_variant[bucket].append(variant_key, len(hits))

        bucket_state = self.bucket_stats.get(bucket)
        if bucket_state is not None:
            bucket_state.denovo_variants += 1
            bucket_state.denovo_entries += len(hits)
            bucket_state.detail_denovo_variants += 1
            bucket_state.detail_denovo_entries += len(hits)

        for person_id in hits:
            self.cumulative.denovo_per_person[person_id] = (
                self.cumulative.denovo_per_person.get(person_id, 0) + 1
            )
            self.detail_denovo_per_person[person_id] = (
                self.detail_denovo_per_person.get(person_id, 0) + 1
            )
            if bucket_state is not None:
                bucket_state.denovo_per_person[person_id] = (
                    bucket_state.denovo_per_person.get(person_id, 0) + 1
                )
                bucket_state.detail_denovo_per_person[person_id] = (
                    bucket_state.detail_denovo_per_person.get(person_id, 0) + 1
                )
        self._update_position(chrom, pos)
        self._maybe_rotate_segment()

    def _update_position(self, chrom: str, pos: str) -> None:
        self.last_chrom = chrom
        self.last_pos = max(self.last_pos, int(pos))

    def _maybe_rotate_segment(self) -> None:
        if self.shard_mode or self.segment_size <= 0:
            return
        if (
            self.inherited_segment_lines >= self.segment_size
            or self.mendelian_bad_segment_lines >= self.segment_size
            or self.denovo_segment_lines >= self.segment_size
        ):
            self._finish_segment()

    def _checkpoint(self, *, completed: bool) -> Checkpoint:
        last_pos = self.last_pos
        if self.shard_mode and self.shard_end is not None:
            last_pos = max(last_pos, self.shard_end)
        return Checkpoint(
            chrom=self.last_chrom,
            last_pos=last_pos,
            segment_index=self.segment_index,
            cumulative=self.cumulative,
            completed=completed,
            shard_start=self.shard_start if self.shard_mode else None,
            shard_end=self.shard_end if self.shard_mode else None,
        )

    def _finish_segment(self, *, open_next: bool = True) -> None:
        self._close_writers()
        self._append_detail_delta()
        self._write_cumulative_stats()
        save_checkpoint(
            self.output_dir,
            self._checkpoint(completed=False),
            include_details=False,
            file_suffix=self.file_suffix,
        )
        self.segment_index += 1
        if open_next:
            self._open_segment_writers()

    def close(self) -> None:
        if self._writers_ready:
            self._close_writers()

    def finalize(self, *, completed: bool = True) -> None:
        if not self._writers_ready:
            self.set_chrom_mode(self.chrom_mode)
        self._append_detail_delta()
        self._write_cumulative_stats()
        if self.shard_mode or self.segment_size > 0 or completed:
            save_checkpoint(
                self.output_dir,
                self._checkpoint(completed=completed),
                include_details=False,
                file_suffix=self.file_suffix,
            )
        self._merge_inherited_per_variant_files()
        self._merge_denovo_per_variant_files()
        self.save_summary_files()

    def _write_cumulative_stats(self) -> None:
        payload: dict[str, Any] = {
            **self.cumulative.to_dict(include_details=False),
            "last_chrom": self.last_chrom,
            "last_pos": self.last_pos,
            "segment_index": self.segment_index,
            "chrom_mode": self.chrom_mode,
        }
        if self.uses_named_buckets:
            payload["buckets"] = {
                name: {
                    "inherited_entries": stats.inherited_entries,
                    "inherited_variants": stats.inherited_variants,
                    "mendelian_bad_entries": stats.mendelian_bad_entries,
                    "mendelian_bad_variants": stats.mendelian_bad_variants,
                    "denovo_entries": stats.denovo_entries,
                    "denovo_variants": stats.denovo_variants,
                }
                for name, stats in self.bucket_stats.items()
            }
        write_json_atomic(self._output_path("stats_cumulative", ".json"), payload)

    def _append_detail_delta(self) -> None:
        if (
            not self.detail_inherited_per_person
            and not self.detail_denovo_per_person
            and not self.detail_mendelian_bad_per_gt
            and not any(
                stats.detail_inherited_per_person
                or stats.detail_denovo_per_person
                or stats.detail_mendelian_bad_per_gt
                or stats.detail_inherited_variants
                or stats.detail_mendelian_bad_variants
                or stats.detail_denovo_variants
                for stats in self.bucket_stats.values()
            )
        ):
            return
        record: dict[str, Any] = {
            "segment_index": self.segment_index,
            "inherited_per_person": self.detail_inherited_per_person,
            "denovo_per_person": self.detail_denovo_per_person,
            "mendelian_bad_per_gt": self.detail_mendelian_bad_per_gt,
        }
        if self.uses_named_buckets:
            record["buckets"] = {
                name: {
                    "inherited_entries": stats.detail_inherited_entries,
                    "inherited_variants": stats.detail_inherited_variants,
                    "mendelian_bad_entries": stats.detail_mendelian_bad_entries,
                    "mendelian_bad_variants": stats.detail_mendelian_bad_variants,
                    "denovo_entries": stats.detail_denovo_entries,
                    "denovo_variants": stats.detail_denovo_variants,
                    "inherited_per_person": stats.detail_inherited_per_person,
                    "denovo_per_person": stats.detail_denovo_per_person,
                    "mendelian_bad_per_gt": stats.detail_mendelian_bad_per_gt,
                }
                for name, stats in self.bucket_stats.items()
                if (
                    stats.detail_inherited_per_person
                    or stats.detail_denovo_per_person
                    or stats.detail_mendelian_bad_per_gt
                    or stats.detail_inherited_variants
                    or stats.detail_mendelian_bad_variants
                    or stats.detail_denovo_variants
                )
            }
        with self._detail_deltas_path().open("a", encoding="utf-8") as handle:
            json.dump(record, handle, separators=(",", ":"), sort_keys=True)
            handle.write("\n")
        self.detail_inherited_per_person.clear()
        self.detail_denovo_per_person.clear()
        self.detail_mendelian_bad_per_gt.clear()
        for stats in self.bucket_stats.values():
            stats.detail_inherited_entries = 0
            stats.detail_inherited_variants = 0
            stats.detail_mendelian_bad_entries = 0
            stats.detail_mendelian_bad_variants = 0
            stats.detail_denovo_entries = 0
            stats.detail_denovo_variants = 0
            stats.detail_inherited_per_person.clear()
            stats.detail_denovo_per_person.clear()
            stats.detail_mendelian_bad_per_gt.clear()

    def _write_detail_baseline(self, segment_index: int) -> None:
        """Convert a legacy self-contained checkpoint to the delta format."""
        record: dict[str, Any] = {
            "segment_index": segment_index,
            "inherited_per_person": self.cumulative.inherited_per_person,
            "denovo_per_person": self.cumulative.denovo_per_person,
            "mendelian_bad_per_gt": self.cumulative.mendelian_bad_per_gt,
        }
        with self._detail_deltas_path().open("w", encoding="utf-8") as handle:
            json.dump(record, handle, separators=(",", ":"), sort_keys=True)
            handle.write("\n")

    def _load_detail_deltas(self, through_segment: int) -> None:
        path = self._detail_deltas_path()
        if not path.is_file():
            raise FileNotFoundError(
                f"Checkpoint references missing cumulative detail file: {path}"
            )

        self.cumulative.inherited_per_person.clear()
        self.cumulative.denovo_per_person.clear()
        self.cumulative.mendelian_bad_per_gt.clear()
        for stats in self.bucket_stats.values():
            stats.inherited_entries = 0
            stats.inherited_variants = 0
            stats.mendelian_bad_entries = 0
            stats.mendelian_bad_variants = 0
            stats.denovo_entries = 0
            stats.denovo_variants = 0
            stats.inherited_per_person.clear()
            stats.denovo_per_person.clear()
            stats.mendelian_bad_per_gt.clear()

        pending_segment = -1
        pending_record: dict[str, Any] | None = None
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                record = json.loads(line)
                segment_index = int(record["segment_index"])
                if segment_index > through_segment:
                    continue
                if pending_record is not None and segment_index != pending_segment:
                    self._apply_detail_record(pending_record)
                pending_segment = segment_index
                pending_record = record
        if pending_record is not None:
            self._apply_detail_record(pending_record)

    def _apply_detail_record(self, record: dict[str, Any]) -> None:
        for person_id, count in dict(record.get("inherited_per_person", {})).items():
            self.cumulative.inherited_per_person[str(person_id)] = (
                self.cumulative.inherited_per_person.get(str(person_id), 0)
                + int(count)
            )
        for person_id, count in dict(record.get("denovo_per_person", {})).items():
            self.cumulative.denovo_per_person[str(person_id)] = (
                self.cumulative.denovo_per_person.get(str(person_id), 0) + int(count)
            )
        for gt_key, count in dict(record.get("mendelian_bad_per_gt", {})).items():
            self.cumulative.mendelian_bad_per_gt[str(gt_key)] = (
                self.cumulative.mendelian_bad_per_gt.get(str(gt_key), 0)
                + int(count)
            )
        for bucket, payload in dict(record.get("buckets", {})).items():
            stats = self.bucket_stats.setdefault(str(bucket), BucketStats())
            stats.inherited_entries += int(payload.get("inherited_entries", 0))
            stats.inherited_variants += int(payload.get("inherited_variants", 0))
            stats.mendelian_bad_entries += int(payload.get("mendelian_bad_entries", 0))
            stats.mendelian_bad_variants += int(
                payload.get("mendelian_bad_variants", 0)
            )
            stats.denovo_entries += int(payload.get("denovo_entries", 0))
            stats.denovo_variants += int(payload.get("denovo_variants", 0))
            for person_id, count in dict(
                payload.get("inherited_per_person", {})
            ).items():
                stats.inherited_per_person[str(person_id)] = (
                    stats.inherited_per_person.get(str(person_id), 0) + int(count)
                )
            for person_id, count in dict(payload.get("denovo_per_person", {})).items():
                stats.denovo_per_person[str(person_id)] = (
                    stats.denovo_per_person.get(str(person_id), 0) + int(count)
                )
            for gt_key, count in dict(payload.get("mendelian_bad_per_gt", {})).items():
                stats.mendelian_bad_per_gt[str(gt_key)] = (
                    stats.mendelian_bad_per_gt.get(str(gt_key), 0) + int(count)
                )

    def _merge_inherited_per_variant_files(self) -> None:
        self._merge_per_variant_files("inherited")

    def _merge_denovo_per_variant_files(self) -> None:
        self._merge_per_variant_files("denovo")

    def _merge_per_variant_files(self, kind: str) -> None:
        for bucket in self._bucket_names():
            stem = self._kind_stem(f"{kind}_per_variant", bucket)
            output_path = self._output_path(stem, ".json")
            if self.shard_mode:
                pattern = f"{stem}_*_*.json"
            else:
                pattern = f"{stem}_seg*.json"
            if not self.shard_mode and self.segment_size <= 0:
                part_path = self._per_variant_path(kind, bucket, self.segment_index)
                if part_path.is_file():
                    part_path.replace(output_path)
                else:
                    write_json_atomic(output_path, {})
                continue

            merged = JsonObjectStreamWriter(output_path.with_suffix(".json.tmp"))
            try:
                for path in sorted(self.output_dir.glob(pattern)):
                    if path == output_path:
                        continue
                    with path.open(encoding="utf-8") as handle:
                        for line in handle:
                            entry = line.strip()
                            if not entry or entry in ("{", "}", "{}"):
                                continue
                            merged.append_encoded(
                                entry.removeprefix(",").removesuffix(",")
                            )
            finally:
                merged.close()
            merged.path.replace(output_path)

    def save_summary_files(self) -> None:
        if self.uses_named_buckets:
            for bucket, stats in self.bucket_stats.items():
                write_json_atomic(
                    self._output_path(f"inherited_per_person_{bucket}", ".json"),
                    stats.inherited_per_person,
                )
                write_json_atomic(
                    self._output_path(f"denovo_per_person_{bucket}", ".json"),
                    stats.denovo_per_person,
                )
                write_json_atomic(
                    self._output_path(f"mendelian_bad_per_gt_{bucket}", ".json"),
                    stats.mendelian_bad_per_gt,
                )
                write_json_atomic(
                    self._output_path(f"stats_{bucket}", ".json"),
                    {
                        "inherited_entries": stats.inherited_entries,
                        "inherited_variants": stats.inherited_variants,
                        "mendelian_bad_entries": stats.mendelian_bad_entries,
                        "mendelian_bad_variants": stats.mendelian_bad_variants,
                        "denovo_entries": stats.denovo_entries,
                        "denovo_variants": stats.denovo_variants,
                    },
                )
        else:
            write_json_atomic(
                self._output_path("inherited_per_person", ".json"),
                self.cumulative.inherited_per_person,
            )
            write_json_atomic(
                self._output_path("denovo_per_person", ".json"),
                self.cumulative.denovo_per_person,
            )
            write_json_atomic(
                self._output_path("mendelian_bad_per_gt", ".json"),
                self.cumulative.mendelian_bad_per_gt,
            )

        write_json_atomic(
            self._output_path("stats", ".json"),
            {
                "variants_seen": self.cumulative.variants_seen,
                "alleles_tested": self.cumulative.alleles_tested,
                "inherited_entries": self.cumulative.inherited_entries,
                "inherited_variants": self.cumulative.inherited_variants,
                "mendelian_bad_entries": self.cumulative.mendelian_bad_entries,
                "mendelian_bad_variants": self.cumulative.mendelian_bad_variants,
                "denovo_entries": self.cumulative.denovo_entries,
                "denovo_variants": self.cumulative.denovo_variants,
                "chrom_mode": self.chrom_mode,
            },
        )

    @property
    def inherited_entries(self) -> int:
        return self.cumulative.inherited_entries

    @property
    def mendelian_bad_entries(self) -> int:
        return self.cumulative.mendelian_bad_entries

    @property
    def denovo_entries(self) -> int:
        return self.cumulative.denovo_entries

    @property
    def inherited_variants(self) -> int:
        return self.cumulative.inherited_variants

    @property
    def mendelian_bad_variants(self) -> int:
        return self.cumulative.mendelian_bad_variants

    @property
    def denovo_variants(self) -> int:
        return self.cumulative.denovo_variants


def read_result_tsv(
    path: Path,
    *,
    short_format: bool = True,
) -> list[tuple[str, str, str, str, object]]:
    """Parse a result TSV into (chrom, pos, ref, alt, payload) records."""
    records = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt, payload = line.rstrip("\n").split("\t", 5)
            if short_format:
                records.append((chrom, pos, ref, alt, parse_patient_ids(payload)))
            else:
                records.append((chrom, pos, ref, alt, parse_trio_calls(payload)))
    return records


def glob_result_tsvs(output_dir: Path, prefix: str) -> list[Path]:
    segmented = list(output_dir.glob(f"{prefix}_*.tsv"))
    if segmented:
        return sorted(segmented, key=_result_tsv_sort_key)
    single = output_dir / f"{prefix}.tsv"
    return [single] if single.is_file() else []


def _result_tsv_sort_key(path: Path) -> tuple[int, int, int, str]:
    name = path.stem
    parts = name.split("_")
    if len(parts) >= 2 and parts[-1].isdigit() and parts[-2].isdigit():
        return (1, int(parts[-2]), int(parts[-1]), name)
    if parts[-1].isdigit():
        return (0, int(parts[-1]), 0, name)
    return (2, 0, 0, name)
