from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from inherited import __version__
from inherited.af import is_rare, load_af_json
from inherited.checkpoint import load_checkpoint, suffixed_filename, validate_file_suffix
from inherited.constants import (
    DEFAULT_AF_THRESHOLD,
    DEFAULT_BLOCK_SIZE,
    DEFAULT_MEMORY_BLOCK,
    DEFAULT_SEGMENT_SIZE,
)
from inherited.classify import (
    classify_father_son,
    classify_mother_son,
    classify_trio,
)
from inherited.debug import log_memory_if_due
from inherited.families import (
    build_sexed_trio_indices,
    load_family_relations,
    person_labels_for_header,
)
from inherited.genotype import (
    DEFAULT_QUALITY,
    QualityFilters,
    get_good_site,
    is_hom_ref,
    sample_gt_has_alt,
)
from inherited.output import HitRecord, ResultWriter
from inherited.repeats import RepeatIntervalFilter
from inherited.shards import VcfShard
from inherited.xchrom import (
    CHROM_MODE_AUTOSOMAL,
    X_BUCKET_FEMALES,
    X_BUCKET_MALE_NONPAR,
    chrom_mode_for,
    male_x_bucket,
    x_region,
)


def get_nfields(line: str, n: int) -> list[str]:
    end = -1
    for _ in range(n):
        end = line.find("\t", end + 1)
        if end < 0:
            return line.rstrip().split("\t")[:n]
    return line[:end].split("\t")


def get_position(line: str) -> int:
    """Parse the POS column without splitting the full VCF sample row."""
    first_tab = line.find("\t")
    second_tab = line.find("\t", first_tab + 1)
    if first_tab < 0 or second_tab < 0:
        raise ValueError("Malformed VCF record: expected at least three columns")
    return int(line[first_tab + 1 : second_tab])


@dataclass
class AnalysisStats:
    variants_seen: int = 0
    alleles_tested: int = 0
    inherited_entries: int = 0
    mendelian_bad_entries: int = 0
    denovo_entries: int = 0
    inherited_variants: int = 0
    mendelian_bad_variants: int = 0
    denovo_variants: int = 0


def analyze_vcf(
    vcf_path: Path | None,
    af_json_path: Path,
    family_file: Path,
    output_dir: Path,
    *,
    vcf_shards: list[VcfShard] | None = None,
    multiallelic: bool = True,
    af_threshold: float = DEFAULT_AF_THRESHOLD,
    debug: bool = False,
    memory_block: int = DEFAULT_MEMORY_BLOCK,
    block_size: int = DEFAULT_BLOCK_SIZE,
    segment_size: int = DEFAULT_SEGMENT_SIZE,
    short_format: bool = True,
    resume: bool = False,
    repeats_path: Path | None = None,
    family_column_map: dict[str, str] | None = None,
    qc: QualityFilters = DEFAULT_QUALITY,
    file_suffix: str | None = None,
) -> AnalysisStats:
    """Scan a VCF, classify trios, and stream results to segmented TSV files."""
    if (vcf_path is None) == (vcf_shards is None):
        raise ValueError("Specify exactly one of vcf_path or vcf_shards")
    shard_mode = vcf_shards is not None
    if file_suffix is not None:
        file_suffix = validate_file_suffix(file_suffix)
        if shard_mode:
            raise ValueError("--file-suffix cannot be used with --vcf-dir")
        if resume:
            raise ValueError("--file-suffix cannot be used with --resume")
        segment_size = 0
    if resume and segment_size <= 0 and not shard_mode:
        raise ValueError("--resume requires --segment-size > 0")

    if vcf_shards is not None:
        shards = list(vcf_shards)
    else:
        assert vcf_path is not None
        shards = [VcfShard(path=vcf_path, chrom="", start=0, end=0)]

    af_table = load_af_json(af_json_path)
    relations = load_family_relations(family_file, column_map=family_column_map)

    checkpoint = load_checkpoint(output_dir, file_suffix) if resume else None
    if resume and checkpoint is None:
        raise FileNotFoundError(f"No checkpoint found in {output_dir}")
    if resume and checkpoint.completed:
        raise ValueError(f"Checkpoint in {output_dir} is already marked completed")

    if checkpoint is not None:
        writer = ResultWriter.from_checkpoint(
            output_dir,
            checkpoint,
            block_size=block_size,
            segment_size=segment_size,
            short_format=short_format,
            chrom_mode=(
                chrom_mode_for(checkpoint.chrom)
                if checkpoint.chrom
                else CHROM_MODE_AUTOSOMAL
            ),
            shard_mode=shard_mode,
            file_suffix=file_suffix,
        )
        resume_last_pos = checkpoint.last_pos
    else:
        writer = ResultWriter(
            output_dir,
            block_size=block_size,
            segment_size=segment_size,
            short_format=short_format,
            shard_mode=shard_mode,
            file_suffix=file_suffix,
        )
        resume_last_pos = -1

    repeat_filter: RepeatIntervalFilter | None = None
    if repeats_path is not None:
        repeat_filter = RepeatIntervalFilter(repeats_path)
        if resume_last_pos >= 0:
            repeat_filter.advance_past(resume_last_pos)

    try:
        female_trios: list[tuple[int, int, int]] = []
        male_trios: list[tuple[int, int, int]] = []
        all_trios: list[tuple[int, int, int]] = []
        sample_header: list[str] = []
        person_header: list[str] = []
        mode_set = False

        for shard in shards:
            if shard_mode and shard.end <= resume_last_pos:
                continue
            if shard_mode:
                writer.begin_shard(shard.start, shard.end)

            opener = gzip.open if str(shard.path).endswith(".gz") else open
            with opener(shard.path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#CHROM"):
                        header = line.strip().split("\t")[9:]
                        if sample_header:
                            if header != sample_header:
                                raise ValueError(
                                    f"Sample header mismatch in {shard.path}"
                                )
                            continue
                        sample_header = header
                        person_header = person_labels_for_header(
                            sample_header, relations.sample_to_person
                        )
                        female_trios, male_trios = build_sexed_trio_indices(
                            sample_header, relations
                        )
                        all_trios = female_trios + male_trios
                        continue

                    pos = get_position(line)
                    if pos <= resume_last_pos:
                        continue
                    if repeat_filter is not None and repeat_filter.in_repeat(pos):
                        continue

                    if not mode_set:
                        chrom = get_nfields(line, 1)[0]
                        writer.set_chrom_mode(chrom_mode_for(chrom))
                        mode_set = True

                    if multiallelic:
                        _process_multiallelic_line(
                            line,
                            af_table,
                            af_threshold,
                            female_trios,
                            male_trios,
                            all_trios,
                            person_header,
                            writer,
                            qc=qc,
                        )
                    else:
                        _process_biallelic_line(
                            line,
                            af_table,
                            af_threshold,
                            female_trios,
                            male_trios,
                            all_trios,
                            person_header,
                            writer,
                            qc=qc,
                        )

                    log_memory_if_due(
                        writer.cumulative.variants_seen,
                        debug=debug,
                        memory_block=memory_block,
                    )
    finally:
        if repeat_filter is not None:
            repeat_filter.close()
        writer.close()

    writer.finalize(completed=True)
    return _stats_from_writer(writer)


def _stats_from_writer(writer: ResultWriter) -> AnalysisStats:
    return AnalysisStats(
        variants_seen=writer.cumulative.variants_seen,
        alleles_tested=writer.cumulative.alleles_tested,
        inherited_entries=writer.inherited_entries,
        mendelian_bad_entries=writer.mendelian_bad_entries,
        denovo_entries=writer.denovo_entries,
        inherited_variants=writer.inherited_variants,
        mendelian_bad_variants=writer.mendelian_bad_variants,
        denovo_variants=writer.denovo_variants,
    )


def _process_multiallelic_line(
    line: str,
    af_table: dict[str, float],
    af_threshold: float,
    female_trios: list[tuple[int, int, int]],
    male_trios: list[tuple[int, int, int]],
    all_trios: list[tuple[int, int, int]],
    person_header: list[str],
    writer: ResultWriter,
    *,
    qc: QualityFilters,
) -> None:
    chrom, pos, keys, ref, alts = get_nfields(line, 5)
    if len(ref) > 1:
        return

    skeys = keys.split(";")
    salts = alts.split(",")
    writer.cumulative.variants_seen += 1

    alleles_to_process: list[tuple[int, str, str]] = []
    for alt_index, key in enumerate(skeys, start=1):
        if alt_index > len(salts):
            break
        alt = salts[alt_index - 1]
        if len(alt) > 1:
            continue
        if not is_rare(af_table, key, af_threshold):
            continue
        alleles_to_process.append((alt_index, key, alt))

    if not alleles_to_process:
        return

    sample_fields = line.rstrip().split("\t")[9:]
    for alt_index, key, alt in alleles_to_process:
        writer.cumulative.alleles_tested += 1
        _process_allele(
            chrom,
            pos,
            ref,
            alt,
            key,
            alt_index,
            sample_fields,
            person_header,
            female_trios,
            male_trios,
            all_trios,
            writer,
            clean_ad=True,
            qc=qc,
        )


def _process_biallelic_line(
    line: str,
    af_table: dict[str, float],
    af_threshold: float,
    female_trios: list[tuple[int, int, int]],
    male_trios: list[tuple[int, int, int]],
    all_trios: list[tuple[int, int, int]],
    person_header: list[str],
    writer: ResultWriter,
    *,
    qc: QualityFilters,
) -> None:
    chrom, pos, key, ref, alt = get_nfields(line, 5)
    if len(ref) > 1 and len(alt) > 1:
        return
    if not is_rare(af_table, key, af_threshold):
        return

    writer.cumulative.variants_seen += 1
    writer.cumulative.alleles_tested += 1
    sample_fields = line.rstrip().split("\t")[9:]
    _process_allele(
        chrom,
        pos,
        ref,
        alt,
        key,
        1,
        sample_fields,
        person_header,
        female_trios,
        male_trios,
        all_trios,
        writer,
        qc=qc,
    )


def _process_allele(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    variant_key: str,
    alt_index: int,
    sample_fields: list[str],
    person_header: list[str],
    female_trios: list[tuple[int, int, int]],
    male_trios: list[tuple[int, int, int]],
    all_trios: list[tuple[int, int, int]],
    writer: ResultWriter,
    *,
    clean_ad: bool = False,
    qc: QualityFilters,
) -> None:
    if writer.x_mode:
        _process_x_allele(
            chrom,
            pos,
            ref,
            alt,
            variant_key,
            alt_index,
            sample_fields,
            person_header,
            female_trios,
            male_trios,
            writer,
            clean_ad=clean_ad,
            qc=qc,
        )
        return

    if writer.y_mode:
        _process_y_allele(
            chrom,
            pos,
            ref,
            alt,
            variant_key,
            alt_index,
            sample_fields,
            person_header,
            male_trios,
            writer,
            clean_ad=clean_ad,
            qc=qc,
        )
        return

    _process_trios_for_allele(
        chrom,
        pos,
        ref,
        alt,
        variant_key,
        alt_index,
        sample_fields,
        person_header,
        all_trios,
        writer,
        clean_ad=clean_ad,
        bucket="",
        qc=qc,
    )


def _process_x_allele(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    variant_key: str,
    alt_index: int,
    sample_fields: list[str],
    person_header: list[str],
    female_trios: list[tuple[int, int, int]],
    male_trios: list[tuple[int, int, int]],
    writer: ResultWriter,
    *,
    clean_ad: bool = False,
    qc: QualityFilters,
) -> None:
    _process_trios_for_allele(
        chrom,
        pos,
        ref,
        alt,
        variant_key,
        alt_index,
        sample_fields,
        person_header,
        female_trios,
        writer,
        clean_ad=clean_ad,
        bucket=X_BUCKET_FEMALES,
        qc=qc,
    )

    pos_int = int(pos)
    region = x_region(pos_int)
    if region == "nonPar":
        _process_male_nonpar_pairs(
            chrom,
            pos,
            ref,
            alt,
            variant_key,
            alt_index,
            sample_fields,
            person_header,
            male_trios,
            writer,
            clean_ad=clean_ad,
            qc=qc,
        )
        return

    _process_trios_for_allele(
        chrom,
        pos,
        ref,
        alt,
        variant_key,
        alt_index,
        sample_fields,
        person_header,
        male_trios,
        writer,
        clean_ad=clean_ad,
        bucket=male_x_bucket(pos_int),
        qc=qc,
    )


def _process_trios_for_allele(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    variant_key: str,
    alt_index: int,
    sample_fields: list[str],
    person_header: list[str],
    trios_ind: list[tuple[int, int, int]],
    writer: ResultWriter,
    *,
    clean_ad: bool = False,
    bucket: str = "",
    qc: QualityFilters,
) -> None:
    parents_cache: dict[int, list[object]] = {}
    inherited_hits: dict[str, HitRecord] = {}
    bad_hits: dict[str, HitRecord] = {}
    denovo_hits: dict[str, HitRecord] = {}

    for child_idx, mother_idx, father_idx in trios_ind:
        child_sample = sample_fields[child_idx]
        if not sample_gt_has_alt(child_sample, alt_index):
            continue
        ac, child_gt, child_gq = get_good_site(
            child_sample,
            alt_index,
            clean_ad=clean_ad,
            skip_qc_if_no_alt=True,
            qc=qc,
        )
        if ac <= 0:
            continue

        if mother_idx in parents_cache:
            mac, mother_gt, mother_gq = parents_cache[mother_idx]
        else:
            mac, mother_gt, mother_gq = get_good_site(
                sample_fields[mother_idx], alt_index, clean_ad=clean_ad, qc=qc
            )
            parents_cache[mother_idx] = [mac, mother_gt, mother_gq]

        if father_idx in parents_cache:
            fac, father_gt, father_gq = parents_cache[father_idx]
        else:
            fac, father_gt, father_gq = get_good_site(
                sample_fields[father_idx], alt_index, clean_ad=clean_ad, qc=qc
            )
            parents_cache[father_idx] = [fac, father_gt, father_gq]

        if mac < 0 or fac < 0:
            continue

        call_class = classify_trio(
            ac, mac, fac, mother_gt, father_gt, child_gt, alt_index
        )
        if call_class is None:
            continue

        pid = person_header[child_idx]
        record: HitRecord = (mother_gt, father_gt, child_gt, child_gq)

        if call_class == "inherited":
            inherited_hits[pid] = record
        elif call_class == "mendelian_bad":
            bad_hits[pid] = record
        else:
            denovo_hits[pid] = record

    if inherited_hits:
        writer.write_inherited(
            chrom, pos, ref, alt, variant_key, inherited_hits, bucket=bucket
        )
    if bad_hits:
        writer.write_mendelian_bad(
            chrom, pos, ref, alt, bad_hits, bucket=bucket
        )
    if denovo_hits:
        writer.write_denovo(
            chrom, pos, ref, alt, variant_key, denovo_hits, bucket=bucket
        )


def _process_male_nonpar_pairs(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    variant_key: str,
    alt_index: int,
    sample_fields: list[str],
    person_header: list[str],
    male_trios: list[tuple[int, int, int]],
    writer: ResultWriter,
    *,
    clean_ad: bool = False,
    qc: QualityFilters,
) -> None:
    """Classify male nonPAR chrX sites using mother-son pairs only."""
    mothers_cache: dict[int, list[object]] = {}
    inherited_hits: dict[str, HitRecord] = {}
    denovo_hits: dict[str, HitRecord] = {}
    bucket = male_x_bucket(int(pos))

    for child_idx, mother_idx, _father_idx in male_trios:
        child_sample = sample_fields[child_idx]
        if not sample_gt_has_alt(child_sample, alt_index):
            continue
        ac, child_gt, child_gq = get_good_site(
            child_sample,
            alt_index,
            clean_ad=clean_ad,
            haploid=True,
            skip_qc_if_no_alt=True,
            qc=qc,
        )
        if ac <= 0:
            continue

        if mother_idx in mothers_cache:
            mac, mother_gt, mother_gq = mothers_cache[mother_idx]
        else:
            mac, mother_gt, mother_gq = get_good_site(
                sample_fields[mother_idx], alt_index, clean_ad=clean_ad, qc=qc
            )
            mothers_cache[mother_idx] = [mac, mother_gt, mother_gq]

        call_class = classify_mother_son(ac, mac)
        if call_class is None:
            continue
        if call_class == "denovo" and not is_hom_ref(mother_gt):
            continue

        pid = person_header[child_idx]
        record: HitRecord = (mother_gt, child_gt, child_gq)
        if call_class == "inherited":
            inherited_hits[pid] = record
        else:
            denovo_hits[pid] = record

    if inherited_hits:
        writer.write_inherited(
            chrom, pos, ref, alt, variant_key, inherited_hits, bucket=bucket
        )
    if denovo_hits:
        writer.write_denovo(
            chrom, pos, ref, alt, variant_key, denovo_hits, bucket=bucket
        )


def _process_y_allele(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    variant_key: str,
    alt_index: int,
    sample_fields: list[str],
    person_header: list[str],
    male_trios: list[tuple[int, int, int]],
    writer: ResultWriter,
    *,
    clean_ad: bool = False,
    qc: QualityFilters,
) -> None:
    """Classify chrY sites using father-son pairs with haploid QC."""
    fathers_cache: dict[int, list[object]] = {}
    inherited_hits: dict[str, HitRecord] = {}
    denovo_hits: dict[str, HitRecord] = {}
    bucket = X_BUCKET_MALE_NONPAR

    for child_idx, _mother_idx, father_idx in male_trios:
        child_sample = sample_fields[child_idx]
        if not sample_gt_has_alt(child_sample, alt_index):
            continue
        ac, child_gt, child_gq = get_good_site(
            child_sample,
            alt_index,
            clean_ad=clean_ad,
            haploid=True,
            skip_qc_if_no_alt=True,
            qc=qc,
        )
        if ac <= 0:
            continue

        if father_idx in fathers_cache:
            fac, father_gt, father_gq = fathers_cache[father_idx]
        else:
            fac, father_gt, father_gq = get_good_site(
                sample_fields[father_idx],
                alt_index,
                clean_ad=clean_ad,
                haploid=True,
                qc=qc,
            )
            fathers_cache[father_idx] = [fac, father_gt, father_gq]

        call_class = classify_father_son(ac, fac)
        if call_class is None:
            continue
        if call_class == "denovo" and not is_hom_ref(father_gt):
            continue

        pid = person_header[child_idx]
        record: HitRecord = (father_gt, child_gt, child_gq)
        if call_class == "inherited":
            inherited_hits[pid] = record
        else:
            denovo_hits[pid] = record

    if inherited_hits:
        writer.write_inherited(
            chrom, pos, ref, alt, variant_key, inherited_hits, bucket=bucket
        )
    if denovo_hits:
        writer.write_denovo(
            chrom, pos, ref, alt, variant_key, denovo_hits, bucket=bucket
        )


def save_run_params(
    output_dir: Path,
    *,
    vcf_path: Path | None = None,
    vcf_dir: Path | None = None,
    vcf_pattern: str | None = None,
    chrom: str | None = None,
    vcf_files: list[Path] | None = None,
    af_json_path: Path,
    family_file: Path,
    multiallelic: bool,
    af_threshold: float,
    debug: bool = False,
    memory_block: int = DEFAULT_MEMORY_BLOCK,
    block_size: int = DEFAULT_BLOCK_SIZE,
    segment_size: int = DEFAULT_SEGMENT_SIZE,
    short_format: bool = True,
    resume: bool = False,
    repeats_path: Path | None = None,
    family_map_path: Path | None = None,
    qc: QualityFilters = DEFAULT_QUALITY,
    file_suffix: str | None = None,
) -> Path:
    """Write the parameters for this run into the chromosome output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    params_path = output_dir / suffixed_filename("params", ".json", file_suffix)

    payload: dict[str, Any] = {
        "package_version": __version__,
        "run_at": datetime.now(timezone.utc).isoformat(),
        "af_json": str(af_json_path.resolve()),
        "family_file": str(family_file.resolve()),
        "output_dir": str(output_dir.resolve()),
        "multiallelic": multiallelic,
        "af_threshold": af_threshold,
        "debug": debug,
        "memory_block": memory_block,
        "block_size": block_size,
        "segment_size": segment_size,
        "short_format": short_format,
        "resume": resume,
        "remove_repeats": str(repeats_path.resolve()) if repeats_path is not None else None,
        "family_map": str(family_map_path.resolve()) if family_map_path is not None else None,
        "quality_filters": qc.as_params(),
        "file_suffix": file_suffix,
    }
    if vcf_path is not None:
        payload["vcf"] = str(vcf_path.resolve())
    if vcf_dir is not None:
        payload["vcf_dir"] = str(vcf_dir.resolve())
        payload["vcf_pattern"] = vcf_pattern
        payload["chr"] = chrom
        payload["vcf_files"] = [
            str(path.resolve()) for path in (vcf_files or [])
        ]
    with params_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    return params_path
