import gzip
import json
from pathlib import Path

import pytest

from inherited.cli import main
from inherited.vcf_stats import collect_vcf_stats, sum_vcf_stats, trim_alleles

_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
)
_SAMPLE = "GT\t0/1"


def _write_vcf(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_HEADER + body, encoding="utf-8")


def _site(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t{_SAMPLE}\n"


def test_a_breakend_is_counted_as_symbolic(tmp_path: Path) -> None:
    _write_vcf(tmp_path / "cohort.chr1.vcf", _site("chr1", 10, "A", "A[chr2:100["))

    stats = collect_vcf_stats(tmp_path, "1", "cohort")

    assert stats["symbolic_count"] == 1
    assert stats["bi_ins_count"] == 0
    assert stats["bi_line_count"] == 1


def test_trim_turns_right_padded_snv_into_one_base() -> None:
    assert trim_alleles(20, "AA", "CA") == (20, "A", "C")


def test_whole_chromosome_file_counts_raw_and_trimmed_alleles(tmp_path: Path) -> None:
    body = "".join(
        [
            _site("chr1", 10, "A", "C"),
            _site("chr1", 20, "AA", "CA"),
            _site("chr1", 30, "AT", "A"),
            _site("chr1", 40, "A", "AT"),
            _site("chr1", 50, "AT", "G,AG,TAC"),
            _site("chr1", 60, "A", "G,<NON_REF>"),
            _site("chr1", 70, "AA", "AA"),
            _site("chr2", 80, "A", "G"),
        ]
    )
    _write_vcf(tmp_path / "cohort.chr1.vcf", body)

    stats = collect_vcf_stats(tmp_path, "1", "cohort")

    assert stats["chrom"] == "chr1"
    assert stats["files"] == ["cohort.chr1.vcf"]
    assert stats["bi_line_count"] == 5
    assert stats["multi_line_count"] == 2
    assert stats["bi_snv_count_raw"] == 1
    assert stats["bi_mnv_count_raw"] == 1
    assert stats["bi_snv_count"] == 2
    assert stats["bi_mnv_count"] == 0
    assert stats["bi_del_count"] == 1
    assert stats["bi_ins_count"] == 1
    assert stats["multi_snv_count_raw"] == 1
    assert stats["multi_mnv_count_raw"] == 1
    assert stats["multi_snv_count"] == 2
    assert stats["multi_mnv_count"] == 0
    assert stats["multi_del_count"] == 1
    assert stats["multi_ins_count"] == 1
    assert stats["symbolic_count"] == 1
    assert stats["ref_equals_alt"] == 1
    assert stats["skipped_other_chrom"] == 1
    assert stats["d_indel"] == {"2": 3, "3": 1}
    assert stats["d_dl"] == {"-1": 2, "1": 2}
    assert stats["bi_mnv_count_raw"] - stats["bi_mnv_count"] == (
        stats["bi_snv_count"] - stats["bi_snv_count_raw"]
    )


def _write_gzip_vcf(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(_HEADER + body)


def test_position_shards_are_read_in_start_order(tmp_path: Path) -> None:
    chrom_dir = tmp_path / "chr1"
    _write_gzip_vcf(chrom_dir / "cohort.chr1_20_30.vcf.gz", _site("1", 25, "A", "G"))
    _write_gzip_vcf(chrom_dir / "cohort.chr1_1_10.vcf.gz", _site("chr1", 5, "A", "T"))

    stats = collect_vcf_stats(tmp_path, "chr1", "cohort")

    assert stats["files"] == ["cohort.chr1_1_10.vcf.gz", "cohort.chr1_20_30.vcf.gz"]
    assert stats["bi_snv_count"] == 2


def test_whole_file_and_shards_together_are_an_error(tmp_path: Path) -> None:
    _write_vcf(tmp_path / "cohort.chr1.vcf.gz", _site("chr1", 1, "A", "C"))
    _write_vcf(tmp_path / "chr1" / "cohort.chr1_1_10.vcf", _site("chr1", 2, "A", "G"))
    # whole file is .vcf.gz but was written as plain text. Discovery only checks the name.
    with pytest.raises(ValueError, match="Both a whole-chromosome VCF and position shards"):
        collect_vcf_stats(tmp_path, "1", "cohort")


def test_sum_adds_chromosome_json_files(tmp_path: Path) -> None:
    _write_vcf(tmp_path / "vcfs" / "cohort.chr1.vcf", _site("chr1", 10, "A", "C"))
    _write_vcf(tmp_path / "vcfs" / "cohort.chr2.vcf", _site("chr2", 10, "AA", "CA"))
    first = collect_vcf_stats(tmp_path / "vcfs", "1", "cohort")
    second = collect_vcf_stats(tmp_path / "vcfs", "2", "cohort")
    stats_dir = tmp_path / "stats"
    stats_dir.mkdir()
    (stats_dir / "chr1.json").write_text(json.dumps(first), encoding="utf-8")
    (stats_dir / "chr2.json").write_text(json.dumps(second), encoding="utf-8")
    (stats_dir / "notes.json").write_text("{}", encoding="utf-8")

    total = sum_vcf_stats(stats_dir)

    assert total["sum"] is True
    assert total["chroms"] == ["chr1", "chr2"]
    assert total["bi_line_count"] == 2
    assert total["bi_snv_count"] == 2
    assert total["bi_snv_count_raw"] == 1
    assert total["bi_mnv_count_raw"] == 1


def test_cli_writes_a_chromosome_json(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write_vcf(tmp_path / "cohort.chrX.vcf", _site("X", 10, "A", "G"))
    output = tmp_path / "out" / "chrX.json"
    main(
        [
            "vcf-stats",
            "--input-dir",
            str(tmp_path),
            "--chr",
            "X",
            "--prefix",
            "cohort",
            "--output",
            str(output),
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["bi_snv_count"] == 1
    assert str(output) in capsys.readouterr().out
