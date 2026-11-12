from __future__ import annotations

import argparse
import sys
from pathlib import Path

from inherited.alleles import validate_allele_mode
from inherited.analyze import analyze_vcf, save_run_params
from inherited.checkpoint import validate_file_suffix
from inherited.filter import PREFIXES, filter_results
from inherited.params import resolved, write_params
from inherited.subtract import subtract_variants
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
from inherited.vcf_stats import (
    aggregate_chromosome_stats,
    collect_vcf_file,
    collect_vcf_stats,
    sum_vcf_stats,
    write_chromosome_stats,
    write_sum_stats,
)


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
        "--trim",
        action="store_true",
        help=(
            "Trim a shared prefix and suffix, then keep a one-base allele. "
            "A padded multi-nucleotide variant such as AT/AG becomes T/G. "
            "Indels and real multi-nucleotide variants are dropped."
        ),
    )
    analyze.add_argument(
        "--keep-indels",
        action="store_true",
        help=(
            "Trim every allele, left-align indels, and keep indels and real "
            "multi-nucleotide variants as well as one-base alleles. Requires --fasta."
        ),
    )
    analyze.add_argument(
        "--fasta",
        type=Path,
        default=None,
        metavar="FILE",
        help=(
            "Uncompressed reference FASTA with a .fai index. "
            "Required with --keep-indels and rejected without it."
        ),
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
        "--ab-hom00-threshold",
        type=float,
        default=None,
        help=(
            "Optional minimum diploid homozygous-ref allele balance. When set, "
            "genotype 0/0 must satisfy ab00 >= value, where ab00 is "
            "AD[REF] / sum(AD) (inclusive). Typical value: 0.9, the same idea "
            "as --ab-hom-threshold. Omitted (default): no 0/0 AB check. "
            "Does not apply to haploid ref or to other ac==0 genotypes such as 0/2"
        ),
    )
    analyze.add_argument(
        "--ab-noncarrier-max",
        type=float,
        default=None,
        help=(
            "Optional maximum share of reads for the queried alt in genotypes that "
            "carry a different alt but not the queried one (e.g. 0/1 or 1/3 when "
            "analysing alt 2; diploid and haploid). When set, such genotypes must "
            "satisfy AD[alt] / sum(AD) <= value (inclusive). Reference-only calls "
            "(0/0, haploid 0) are not affected; see --ab-hom00-threshold. "
            "Omitted (default): no AB check on these genotypes"
        ),
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
    analyze.add_argument(
        "--file-suffix",
        default=None,
        metavar="SUFFIX",
        help=(
            "Label every output file as {stem}_{suffix}.ext. Requires --vcf and "
            "disables output segmentation. Cannot be used with --vcf-dir or --resume."
        ),
    )

    filter_cmd = subparsers.add_parser(
        "filter",
        help="Cap patients per variant and variants per patient in result TSVs",
    )
    filter_cmd.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Parent of chrN, chrX, and chrY result directories",
    )
    filter_cmd.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help=(
            "Directory for histograms and, when a cap is set, filtered VCFs. "
            "Filtered TSVs are written to this path with -sites appended"
        ),
    )
    filter_cmd.add_argument(
        "--prefix",
        required=True,
        choices=PREFIXES,
        help="Result files to read: inherited, denovo, or mendelian_bad",
    )
    filter_cmd.add_argument(
        "--variant-cap",
        type=_nonnegative_int,
        default=None,
        help="Keep a person when their variant count is <= this value",
    )
    filter_cmd.add_argument(
        "--patient-cap",
        type=_nonnegative_int,
        default=None,
        help="Drop a variant when its patient count is > this value",
    )

    subtract_cmd = subparsers.add_parser(
        "subtract",
        help="Keep variants present in one VCF tree and absent from the paired file",
    )
    subtract_cmd.add_argument(
        "--input-dir1",
        required=True,
        type=Path,
        help="Directory of chromosome VCFs to keep",
    )
    subtract_cmd.add_argument(
        "--input-dir2",
        required=True,
        type=Path,
        help="Directory of paired VCFs whose variants are removed",
    )
    subtract_cmd.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for variants_${start}_${end}.vcf files",
    )
    subtract_cmd.add_argument(
        "--no-sites",
        dest="sites",
        action="store_false",
        help=(
            "Subtract VCFs only: do not read {input_dir1}-sites and do not write "
            "{output_dir}-sites TSVs or patients.txt"
        ),
    )

    stats_cmd = subparsers.add_parser(
        "vcf-stats",
        help="Count SNV, MNV, insertion, and deletion alleles in one chromosome of pVCFs",
    )
    stats_input = stats_cmd.add_mutually_exclusive_group(required=True)
    stats_input.add_argument(
        "--input-dir",
        type=Path,
        help="VCF directory, or the parent of chrN stats directories when aggregating.",
    )
    stats_input.add_argument(
        "--vcf",
        type=Path,
        help="One VCF: {prefix}.chrN.vcf.gz or {prefix}.chrN_{start}_{end}.vcf.gz.",
    )
    stats_cmd.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Parent directory for chrN stats. With --sum, the directory that receives the totals.",
    )
    stats_cmd.add_argument(
        "--chr",
        default=None,
        help="Chromosome (1, X, Y, or chr1). With --prefix, scan its VCFs. "
        "Without --prefix, add that chromosome's shard JSON files.",
    )
    stats_cmd.add_argument(
        "--prefix",
        default=None,
        help="Filename prefix before .chrN. With --input-dir and --chr, scan VCFs.",
    )
    stats_cmd.add_argument(
        "--sum",
        action="store_true",
        help="Add unsuffixed chrN/summary.json files under --input-dir. Does not read VCFs.",
    )

    return parser


def _nonnegative_int(text: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer >= 0") from exc
    if value < 0:
        raise argparse.ArgumentTypeError("must be an integer >= 0")
    return value


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    if args.command == "filter":
        _run_filter(args)
        return

    if args.command == "vcf-stats":
        _run_vcf_stats(args)
        return

    if args.command == "subtract":
        _run_subtract(args)
        return

    if args.command == "analyze":
        vcf_shards = None
        file_suffix = args.file_suffix
        if file_suffix is not None:
            if args.vcf_dir is not None:
                print("error: --file-suffix requires --vcf", file=sys.stderr)
                raise SystemExit(1)
            if args.resume:
                print("error: --file-suffix cannot be used with --resume", file=sys.stderr)
                raise SystemExit(1)
            try:
                file_suffix = validate_file_suffix(file_suffix)
            except ValueError as exc:
                print(f"error: {exc}", file=sys.stderr)
                raise SystemExit(1) from exc
        segment_size = 0 if file_suffix is not None else args.segment_size
        try:
            validate_allele_mode(
                trim=args.trim,
                keep_indels=args.keep_indels,
                fasta=args.fasta,
                multiallelic=args.multiallelic,
            )
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            raise SystemExit(1) from exc
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
                ab_hom00=args.ab_hom00_threshold,
                ab_noncarrier_max=args.ab_noncarrier_max,
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
                segment_size=segment_size,
                short_format=args.short_format,
                resume=args.resume,
                repeats_path=args.remove_repeats,
                family_column_map=family_column_map,
                qc=qc,
                file_suffix=file_suffix,
                trim=args.trim,
                keep_indels=args.keep_indels,
                fasta_path=args.fasta,
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
            segment_size=segment_size,
            short_format=args.short_format,
            resume=args.resume,
            repeats_path=args.remove_repeats,
            family_map_path=args.family_map,
            qc=qc,
            file_suffix=file_suffix,
            trim=args.trim,
            keep_indels=args.keep_indels,
            fasta_path=args.fasta,
        )
        if vcf_shards is not None:
            output_label = "shard-labeled TSV files"
        elif file_suffix is not None:
            output_label = (
                f"inherited_{file_suffix}.tsv / mendelian_bad_{file_suffix}.tsv / "
                f"denovo_{file_suffix}.tsv"
            )
        elif segment_size > 0:
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


def _run_filter(args: argparse.Namespace) -> None:
    if not args.input_dir.is_dir():
        print(f"error: input directory not found: {args.input_dir}", file=sys.stderr)
        raise SystemExit(1)
    try:
        sites_dir = filter_results(
            args.input_dir,
            args.output_dir,
            args.prefix,
            variant_cap=args.variant_cap,
            patient_cap=args.patient_cap,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    if sites_dir is None:
        print(f"Wrote histograms to {args.output_dir}")
    else:
        print(
            f"Wrote histograms and filtered VCFs to {args.output_dir} "
            f"and TSVs to {sites_dir}"
        )
    print(f"Wrote parameters to {args.output_dir / 'params.json'}")


def _run_vcf_stats(args: argparse.Namespace) -> None:
    try:
        if args.vcf is not None:
            if args.sum or args.chr or args.prefix:
                raise ValueError("--vcf does not take --sum, --chr, or --prefix")
            mode = "vcf"
            payload, span = collect_vcf_file(args.vcf)
            destination = write_chromosome_stats(args.output_dir, payload, span)
            params_name = "params.json" if span is None else f"params_{span[0]}_{span[1]}.json"
        elif args.sum:
            if args.input_dir is None:
                raise ValueError("--sum requires --input-dir")
            if args.chr or args.prefix:
                raise ValueError("--sum does not take --chr or --prefix")
            mode = "sum"
            payload = sum_vcf_stats(args.input_dir)
            destination = write_sum_stats(args.output_dir, payload)
            params_name = "params.json"
        elif args.prefix:
            if args.input_dir is None or not args.chr:
                raise ValueError("--prefix requires --input-dir and --chr")
            mode = "prefix"
            payload = collect_vcf_stats(args.input_dir, args.chr, args.prefix)
            destination = write_chromosome_stats(args.output_dir, payload)
            params_name = "params.json"
        else:
            if args.input_dir is None or not args.chr:
                raise ValueError("--chr and --input-dir are required to aggregate shard JSON files")
            mode = "aggregate"
            payload = aggregate_chromosome_stats(args.input_dir, args.chr)
            destination = write_chromosome_stats(args.output_dir, payload)
            params_name = "params.json"
        params_path = write_params(
            destination / params_name,
            {
                "command": "vcf-stats",
                "mode": mode,
                "input_dir": resolved(args.input_dir),
                "vcf": resolved(args.vcf),
                "chr": args.chr,
                "prefix": args.prefix,
                "sum": args.sum,
                "output_dir": str(args.output_dir.resolve()),
            },
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(f"Wrote {destination}")
    print(f"Wrote parameters to {params_path}")


def _run_subtract(args: argparse.Namespace) -> None:
    try:
        subtract_variants(
            args.input_dir1, args.input_dir2, args.output_dir, sites=args.sites
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    if args.sites:
        sites_dir = args.output_dir.with_name(args.output_dir.name + "-sites")
        print(f"Wrote subtracted VCFs to {args.output_dir} and TSVs to {sites_dir}")
    else:
        print(f"Wrote subtracted VCFs to {args.output_dir}")
    print(f"Wrote parameters to {args.output_dir / 'params.json'}")


if __name__ == "__main__":
    main()
