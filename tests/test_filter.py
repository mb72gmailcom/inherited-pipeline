from pathlib import Path

import pytest

from inherited.cli import build_parser, main
from inherited.filter import VCF_HEADER, TSV_HEADER, filter_results

SHORT_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(SHORT_HEADER + body, encoding="utf-8")


def _cohort(root: Path) -> None:
    _write(
        root / "chr21" / "inherited_00000.tsv",
        "chr21\t10\t.\tA\tG\tp1;p2;p3\n"
        "chr21\t20\t.\tC\tT\tp1;p4\n"
        "chr21\t30\t.\tG\tA\tp4\n",
    )
    _write(
        root / "chr22" / "inherited_00000.tsv",
        "chr22\t40\t.\tT\tC\tp1\n"
        "chr22\t50\t.\tA\tC\tp2;p4\n",
    )
    _write(
        root / "chrX" / "inherited_males_nonPar_00000.tsv",
        "chrX\t70\t.\tA\tG\tp3\n",
    )
    _write(
        root / "chr22" / "denovo_00000.tsv",
        "chr22\t60\t.\tA\tT\tp1\n",
    )
    _write(
        root / "notes" / "inherited_00000.tsv",
        "chr21\t80\t.\tA\tT\tp9;p9;p9\n",
    )


def test_filter_parser_defaults():
    args = build_parser().parse_args(
        ["filter", "--input-dir", "in", "--output-dir", "out", "--prefix", "denovo"]
    )
    assert args.prefix == "denovo"
    assert args.variant_cap is None
    assert args.patient_cap is None


def test_filter_parser_rejects_negative_cap():
    with pytest.raises(SystemExit):
        build_parser().parse_args(
            [
                "filter",
                "--input-dir",
                "in",
                "--output-dir",
                "out",
                "--prefix",
                "inherited",
                "--variant-cap",
                "-1",
                "--patient-cap",
                "1",
            ]
        )


def test_histograms_without_caps(tmp_path: Path):
    source = tmp_path / "results"
    output = tmp_path / "kept"
    _cohort(source)

    sites = filter_results(source, output, "inherited")

    assert sites is None
    assert not (tmp_path / "kept-sites").exists()
    assert not (output / "chr21").exists()
    assert (output / "inherited_variants_per_person_hist.tsv").read_text(
        encoding="utf-8"
    ) == ("count\tn\n2\t2\n3\t2\n")
    assert (output / "inherited_patients_per_variant_hist.tsv").read_text(
        encoding="utf-8"
    ) == ("count\tn\n1\t3\n2\t2\n3\t1\n")


def test_both_caps_write_vcf_and_sites(tmp_path: Path):
    source = tmp_path / "results"
    output = tmp_path / "kept"
    _cohort(source)

    sites = filter_results(
        source,
        output,
        "inherited",
        variant_cap=2,
        patient_cap=2,
    )

    assert sites == tmp_path / "kept-sites"
    assert (output / "chr21" / "inherited_00000.vcf").read_text(encoding="utf-8") == VCF_HEADER
    assert (sites / "chr21" / "inherited_00000.tsv").read_text(encoding="utf-8") == TSV_HEADER
    assert (output / "chr22" / "inherited_00000.vcf").read_text(encoding="utf-8") == (
        VCF_HEADER + "chr22\t50\t.\tA\tC\n"
    )
    assert (sites / "chr22" / "inherited_00000.tsv").read_text(encoding="utf-8") == (
        TSV_HEADER + "chr22\t50\t.\tA\tC\tp2\n"
    )
    assert (output / "chrX" / "inherited_males_nonPar_00000.vcf").read_text(
        encoding="utf-8"
    ) == (VCF_HEADER + "chrX\t70\t.\tA\tG\n")
    assert (sites / "chrX" / "inherited_males_nonPar_00000.tsv").read_text(
        encoding="utf-8"
    ) == (TSV_HEADER + "chrX\t70\t.\tA\tG\tp3\n")
    assert not (output / "chr22" / "denovo_00000.vcf").exists()
    assert not (sites / "notes").exists()


def test_one_cap_is_an_error(tmp_path: Path, capsys):
    source = tmp_path / "results"
    source.mkdir()
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "filter",
                "--input-dir",
                str(source),
                "--output-dir",
                str(tmp_path / "kept"),
                "--prefix",
                "inherited",
                "--patient-cap",
                "2",
            ]
        )
    assert exc.value.code == 1
    assert "must be given together" in capsys.readouterr().err
    assert not (tmp_path / "kept").exists()


def test_full_format_header_requires_short_format(tmp_path: Path, capsys):
    source = tmp_path / "results" / "chr21"
    source.mkdir(parents=True)
    (source / "inherited_00000.tsv").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tTRIO_CALLS\n"
        "chr21\t10\t.\tA\tG\tp1=0/1|0/0|0/1|30\n",
        encoding="utf-8",
    )
    output = tmp_path / "kept"
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "filter",
                "--input-dir",
                str(tmp_path / "results"),
                "--output-dir",
                str(output),
                "--prefix",
                "inherited",
            ]
        )
    assert exc.value.code == 1
    err = capsys.readouterr().err
    assert err.startswith("error: ")
    assert "inherited_00000.tsv requires input in short format" in err
    assert not output.exists()


def test_full_format_payload_requires_short_format(tmp_path: Path):
    source = tmp_path / "results"
    _write(
        source / "chr21" / "mendelian_bad.tsv",
        "chr21\t10\t.\tA\tG\tp1=0/1|0/0|0/1|30\n",
    )
    with pytest.raises(ValueError, match="requires input in short format"):
        filter_results(source, tmp_path / "kept", "mendelian_bad", variant_cap=5, patient_cap=5)


def test_keeps_patient_order(tmp_path: Path):
    source = tmp_path / "results"
    _write(
        source / "chr2" / "denovo.tsv",
        "chr2\t10\t.\tA\tG\tp3;p1;p2\n",
    )
    sites = filter_results(
        source,
        tmp_path / "kept",
        "denovo",
        variant_cap=1,
        patient_cap=3,
    )
    text = (sites / "chr2" / "denovo.tsv").read_text(encoding="utf-8")
    assert text.endswith("chr2\t10\t.\tA\tG\tp3;p1;p2\n")
