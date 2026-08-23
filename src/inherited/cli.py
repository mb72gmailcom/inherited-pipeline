from __future__ import annotations

import argparse
import sys
from pathlib import Path

from inherited.analyze import analyze_vcf, save_run_params
from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_AB_HOM,
    DEFAULT_AF_THRESHOLD,
    DEFAULT_BLOCK_SIZE,
    DEFAULT_DP,
    DEFAULT_GQ,
    DEFAULT_HAPLO_AB,
    DEFAULT_HAPLO_DP,
    DEFAULT_MEMORY_BLOCK,
    DEFAULT_SEGMENT_SIZE,
)
from inherited.families import load_family_column_map
from inherited.genotype import QualityFilters
from inherited.shards import discover_vcf_shards


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="inherited",
        description="Classify inherited and mendelian-inconsistent rare variants in family trios",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser(
        "analyze",
        help="Analyze a VCF using gnomAD AF JSON and family relations",
    )
    vcf_input = analyze.add_mutually_exclusive_group(required=True)
    vcf_input.add_argument(
        "--vcf",
        type=Path,
        help="Input VCF (.vcf or .vcf.gz)",
    )
    vcf_input.add_argument(
        "--vcf-dir",
        type=Path,
        help="Directory of split pVCF shards named {pattern}.{chr}_{start}_{end}.vcf.gz",
    )
    analyze.add_argument(
        "--vcf-pattern",
        help=(
            "Filename prefix for --vcf-dir shards named "
            "{pattern}.{chr}_{start}_{end}.vcf.gz. The contig is read from "
            "the filenames; the directory should contain one chromosome."
        ),
    )
    analyze.add_argument(
        "--af-json",
        required=True,
        type=Path,
        help="JSON file with precomputed gnomAD allele frequencies",
    )
    analyze.add_argument(
        "--family-file",
        required=True,
        type=Path,
        help="Tab-separated family relations file",
    )
    analyze.add_argument(
        "--family-map",
        type=Path,
        default=None,
        metavar="FILE",
        help=(
            "JSON object mapping internal family columns "
            "(spid, sfid, father, mother, sex) to headers in --family-file. "
            "Unmapped columns still use built-in aliases "
            "(spid/ind_id, sfid/family_id, father/father_id, mother/mother_id, sex)."
        ),
    )
    analyze.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for segmented inherited/mendelian_bad/denovo TSV output",
    )
    analyze.add_argument(
        "--multiallelic",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Process multiallelic sites per alternate allele (default: True)",
    )
    analyze.add_argument(
        "--af-threshold",
        type=float,
        default=DEFAULT_AF_THRESHOLD,
        help=f"Maximum gnomAD AF to include (default: {DEFAULT_AF_THRESHOLD})",
    )
    analyze.add_argument(
        "--gq-threshold",
        type=int,
        default=DEFAULT_GQ,
        help=f"Minimum genotype quality (default: {DEFAULT_GQ})",
    )
    analyze.add_argument(
        "--dp-threshold",
        type=int,
        default=DEFAULT_DP,
        help=f"Minimum diploid depth (default: {DEFAULT_DP})",
    )
    analyze.add_argument(
        "--dp-haploid-threshold",
        type=int,
        default=DEFAULT_HAPLO_DP,
        help=f"Minimum haploid depth for male nonPAR chrX/chrY (default: {DEFAULT_HAPLO_DP})",
    )
    analyze.add_argument(
        "--ab-threshold",
        type=float,
        default=DEFAULT_AB,
        help=(
            "Diploid heterozygous allele-balance half-band; het AB must fall in "
            f"[value, 1-value] (default: {DEFAULT_AB})"
        ),
    )
    analyze.add_argument(
        "--ab-hom-threshold",
        type=float,
        default=DEFAULT_AB_HOM,
        help=f"Minimum diploid homozygous-alt allele balance (default: {DEFAULT_AB_HOM})",
    )
    analyze.add_argument(
        "--ab-haploid-threshold",
        type=float,
        default=DEFAULT_HAPLO_AB,
        help=f"Minimum haploid allele balance when the alt is present (default: {DEFAULT_HAPLO_AB})",
    )
    analyze.add_argument(
        "--debug",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Print memory usage periodically during analysis (default: False)",
    )
    analyze.add_argument(
        "--memory-block",
        type=int,
        default=DEFAULT_MEMORY_BLOCK,
        help=f"Variant interval for debug memory logging (default: {DEFAULT_MEMORY_BLOCK})",
    )
    analyze.add_argument(
        "--block-size",
        type=int,
        default=DEFAULT_BLOCK_SIZE,
        help=f"Lines per block when streaming TSV output (default: {DEFAULT_BLOCK_SIZE})",
    )
    analyze.add_argument(
        "--segment-size",
        type=int,
        default=DEFAULT_SEGMENT_SIZE,
        help=(
            f"Max result lines per output segment when using --vcf; 0 disables "
            f"segmentation. Ignored with --vcf-dir, where each input shard "
            f"writes files labeled {{start}}_{{end}} "
            f"(default: {DEFAULT_SEGMENT_SIZE})"
        ),
    )
    analyze.add_argument(
        "--short-format",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write only patient IDs in the last TSV column (default: True)",
    )
    analyze.add_argument(
        "--resume",
        action="store_true",
        help="Resume from checkpoint.json in the output directory",
    )
    analyze.add_argument(
        "--remove-repeats",
        type=Path,
        default=None,
        metavar="FILE",
        help="Skip variants inside 0-based BED intervals [start, end) from this file",
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    if args.command == "analyze":
        vcf_shards = None
        if args.vcf_dir is not None:
            if not args.vcf_pattern:
                print("error: --vcf-dir requires --vcf-pattern", file=sys.stderr)
                raise SystemExit(1)
            if not args.vcf_dir.is_dir():
                print(f"error: VCF directory not found: {args.vcf_dir}", file=sys.stderr)
                raise SystemExit(1)
        else:
            if args.vcf_pattern:
                print("error: --vcf-pattern requires --vcf-dir", file=sys.stderr)
                raise SystemExit(1)
            if not args.vcf.is_file():
                print(f"error: VCF not found: {args.vcf}", file=sys.stderr)
                raise SystemExit(1)
        if not args.af_json.is_file():
            print(f"error: AF JSON not found: {args.af_json}", file=sys.stderr)
            raise SystemExit(1)
        if not args.family_file.is_file():
            print(f"error: family file not found: {args.family_file}", file=sys.stderr)
            raise SystemExit(1)
        if args.remove_repeats is not None and not args.remove_repeats.is_file():
            print(f"error: repeat intervals file not found: {args.remove_repeats}", file=sys.stderr)
            raise SystemExit(1)
        if args.family_map is not None and not args.family_map.is_file():
            print(f"error: family-map file not found: {args.family_map}", file=sys.stderr)
            raise SystemExit(1)

        try:
            if args.vcf_dir is not None:
                vcf_shards = discover_vcf_shards(args.vcf_dir, args.vcf_pattern)
            family_column_map = (
                load_family_column_map(args.family_map)
                if args.family_map is not None
                else None
            )
            qc = QualityFilters(
                gq=args.gq_threshold,
                dp=args.dp_threshold,
                ab=args.ab_threshold,
                ab_hom=args.ab_hom_threshold,
                haplo_dp=args.dp_haploid_threshold,
                haplo_ab=args.ab_haploid_threshold,
            )
            stats = analyze_vcf(
                vcf_path=args.vcf,
                af_json_path=args.af_json,
                family_file=args.family_file,
                output_dir=args.output_dir,
                vcf_shards=vcf_shards,
                multiallelic=args.multiallelic,
                af_threshold=args.af_threshold,
                debug=args.debug,
                memory_block=args.memory_block,
                block_size=args.block_size,
                segment_size=args.segment_size,
                short_format=args.short_format,
                resume=args.resume,
                repeats_path=args.remove_repeats,
                family_column_map=family_column_map,
                qc=qc,
            )
        except (FileNotFoundError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            raise SystemExit(1) from exc

        params_path = save_run_params(
            args.output_dir,
            vcf_path=args.vcf,
            vcf_dir=args.vcf_dir,
            vcf_pattern=args.vcf_pattern,
            chrom=vcf_shards[0].chrom if vcf_shards else None,
            vcf_files=[shard.path for shard in vcf_shards] if vcf_shards else None,
            af_json_path=args.af_json,
            family_file=args.family_file,
            multiallelic=args.multiallelic,
            af_threshold=args.af_threshold,
            debug=args.debug,
            memory_block=args.memory_block,
            block_size=args.block_size,
            segment_size=args.segment_size,
            short_format=args.short_format,
            resume=args.resume,
            repeats_path=args.remove_repeats,
            family_map_path=args.family_map,
            qc=qc,
        )
        if vcf_shards is not None:
            output_label = "shard-labeled TSV files"
        elif args.segment_size > 0:
            output_label = "segmented TSV files"
        else:
            output_label = "inherited.tsv / mendelian_bad.tsv / denovo.tsv"
        print(
            f"Wrote {stats.inherited_entries} inherited entries "
            f"({stats.inherited_variants} variants), "
            f"{stats.mendelian_bad_entries} mendelian_bad entries "
            f"({stats.mendelian_bad_variants} variants), and "
            f"{stats.denovo_entries} denovo entries "
            f"({stats.denovo_variants} variants) as {output_label} "
            f"to {args.output_dir}"
        )
        print(f"Wrote parameters to {params_path}")


if __name__ == "__main__":
    main()
