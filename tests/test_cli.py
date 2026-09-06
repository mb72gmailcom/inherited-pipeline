import json
from pathlib import Path

import pytest

from inherited.cli import build_parser, main
from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_AB_HOM,
    DEFAULT_DP,
    DEFAULT_GQ,
    DEFAULT_HAPLO_AB,
    DEFAULT_HAPLO_DP,
)

FIXTURES = Path(__file__).parent / "fixtures"

_ANALYZE_MIN = [
    "analyze",
    "--vcf",
    "a.vcf",
    "--af-json",
    "af.json",
    "--family-file",
    "fam.tsv",
    "-o",
    "out",
]


def test_analyze_parser_qc_defaults():
    args = build_parser().parse_args(_ANALYZE_MIN)
    assert args.gq_threshold == DEFAULT_GQ
    assert args.dp_threshold == DEFAULT_DP
    assert args.dp_haploid_threshold == DEFAULT_HAPLO_DP
    assert args.ab_threshold == DEFAULT_AB
    assert args.ab_hom_threshold == DEFAULT_AB_HOM
    assert args.ab_hom00_threshold is None
    assert args.ab_haploid_threshold == DEFAULT_HAPLO_AB
    assert args.vcf_dir is None
    assert args.vcf_pattern is None
    assert args.family_map is None
    assert args.file_suffix is None


def test_analyze_parser_qc_overrides():
    args = build_parser().parse_args(
        [
            *_ANALYZE_MIN,
            "--gq-threshold",
            "30",
            "--dp-threshold",
            "15",
            "--dp-haploid-threshold",
            "8",
            "--ab-threshold",
            "0.3",
            "--ab-hom-threshold",
            "0.95",
            "--ab-hom00-threshold",
            "0.9",
            "--ab-haploid-threshold",
            "0.8",
        ]
    )
    assert args.gq_threshold == 30
    assert args.dp_threshold == 15
    assert args.dp_haploid_threshold == 8
    assert args.ab_threshold == 0.3
    assert args.ab_hom_threshold == 0.95
    assert args.ab_hom00_threshold == 0.9
    assert args.ab_haploid_threshold == 0.8


def test_analyze_parser_vcf_dir_mode():
    args = build_parser().parse_args(
        [
            "analyze",
            "--vcf-dir",
            "shards",
            "--vcf-pattern",
            "SPARK.WGS.2026_08.gatk",
            "--af-json",
            "af.json",
            "--family-file",
            "fam.tsv",
            "-o",
            "out",
        ]
    )
    assert args.vcf is None
    assert args.vcf_dir.as_posix() == "shards"
    assert args.vcf_pattern == "SPARK.WGS.2026_08.gatk"


def test_analyze_parser_family_map():
    args = build_parser().parse_args([*_ANALYZE_MIN, "--family-map", "map.json"])
    assert args.family_map.as_posix() == "map.json"


def test_analyze_parser_rejects_vcf_and_vcf_dir_together():
    with pytest.raises(SystemExit):
        build_parser().parse_args(
            [
                *_ANALYZE_MIN,
                "--vcf-dir",
                "shards",
            ]
        )


def test_analyze_parser_file_suffix():
    args = build_parser().parse_args([*_ANALYZE_MIN, "--file-suffix", "1_2500"])
    assert args.file_suffix == "1_2500"


def test_main_rejects_file_suffix_with_vcf_dir(capsys):
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "analyze",
                "--vcf-dir",
                "shards",
                "--vcf-pattern",
                "callset",
                "--af-json",
                "af.json",
                "--family-file",
                "fam.tsv",
                "-o",
                "out",
                "--file-suffix",
                "1_2500",
            ]
        )
    assert exc.value.code == 1
    assert "--file-suffix requires --vcf" in capsys.readouterr().err


def test_main_rejects_file_suffix_with_resume(capsys):
    with pytest.raises(SystemExit) as exc:
        main([*_ANALYZE_MIN, "--file-suffix", "1_2500", "--resume"])
    assert exc.value.code == 1
    assert "--file-suffix cannot be used with --resume" in capsys.readouterr().err


def test_main_rejects_file_suffix_with_path_separator(capsys):
    with pytest.raises(SystemExit) as exc:
        main([*_ANALYZE_MIN, "--file-suffix", "a/b"])
    assert exc.value.code == 1
    assert "path separator" in capsys.readouterr().err


def test_main_file_suffix_writes_suffixed_params(tmp_path):
    out = tmp_path / "out"
    main(
        [
            "analyze",
            "--vcf",
            str(FIXTURES / "tiny.vcf"),
            "--af-json",
            str(FIXTURES / "tiny_af.json"),
            "--family-file",
            str(FIXTURES / "families.tsv"),
            "-o",
            str(out),
            "--file-suffix",
            "1_2500",
        ]
    )
    params = json.loads((out / "params_1_2500.json").read_text(encoding="utf-8"))
    assert params["file_suffix"] == "1_2500"
    assert params["segment_size"] == 0
    assert not (out / "params.json").exists()
    assert (out / "inherited_1_2500.tsv").is_file()
