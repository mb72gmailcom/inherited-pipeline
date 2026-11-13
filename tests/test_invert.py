import json
from pathlib import Path

import pytest

from inherited.cli import main
from inherited.invert import invert_sites

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER + body, encoding="utf-8")


def _read(path: Path) -> list[list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "#PATIENT\tN_VARIANTS\tVARIANTS"
    return [line.split("\t") for line in lines[1:]]


def test_invert_both_namings_sorted_and_merged(tmp_path: Path):
    src = tmp_path / "sites"
    # start_end naming; second shard listed first by name order, rows out of order
    _write(
        src / "chr21" / "inherited_12500001_15000000.tsv",
        "chr21\t12500010\t.\tA\tC\tp2\n",
    )
    _write(
        src / "chr21" / "inherited_10000001_12500000.tsv",
        "chr21\t10000030\t.\tG\tGA\tp1;p3\n"
        "chr21\t10000010\t.\tA\tG\tp1;p2\n",
    )
    # old segment naming
    _write(src / "chr2" / "inherited_00000.tsv", "chr2\t500\t.\tT\tA\tp3\n")
    _write(src / "chr2" / "inherited_00001.tsv", "chr2\t900\t.\tC\tT\tp1\n")
    # other prefixes and folders are ignored
    _write(src / "chr2" / "denovo_00000.tsv", "chr2\t1\t.\tA\tT\tp9\n")
    _write(src / "notes" / "inherited_00000.tsv", "chr2\t1\t.\tA\tT\tp9\n")

    out = tmp_path / "out"
    totals = invert_sites(src, out)
    assert totals == {"variants": 5, "pairs": 7, "patients": 3}

    assert _read(out / "chr21" / "patients.tsv") == [
        ["p1", "2", "chr21_10000010_A_G;chr21_10000030_G_GA"],
        ["p2", "2", "chr21_10000010_A_G;chr21_12500010_A_C"],
        ["p3", "1", "chr21_10000030_G_GA"],
    ]
    assert _read(out / "chr2" / "patients.tsv") == [
        ["p1", "1", "chr2_900_C_T"],
        ["p3", "1", "chr2_500_T_A"],
    ]
    assert not (out / "patients.tsv").exists()
    assert (out / "counts.tsv").read_text().splitlines() == [
        "#PATIENT\tTOTAL\tchr2\tchr21",
        "p1\t3\t1\t2",
        "p2\t2\t0\t2",
        "p3\t2\t1\t1",
    ]
    params = json.loads((out / "params.json").read_text())
    assert params["command"] == "invert"
    assert params["chromosomes"] == ["chr2", "chr21"]
    assert params["per_chromosome"]["chr21"]["files"] == 2
    assert params["totals"] == totals
    assert "package_version" in params


def test_invert_duplicate_rows_written_once(tmp_path: Path):
    src = tmp_path / "sites"
    _write(src / "chr1" / "inherited_00000.tsv", "chr1\t10\t.\tA\tG\tp1;p1;p2\n")
    _write(src / "chr1" / "inherited_00001.tsv", "chr1\t10\t.\tA\tG\tp1;p3\n")
    out = tmp_path / "out"
    invert_sites(src, out)
    assert _read(out / "chr1" / "patients.tsv") == [
        ["p1", "1", "chr1_10_A_G"],
        ["p2", "1", "chr1_10_A_G"],
        ["p3", "1", "chr1_10_A_G"],
    ]
    params = json.loads((out / "params.json").read_text())
    assert params["per_chromosome"]["chr1"]["duplicate_rows"] == 1


def test_invert_chrom_selection(tmp_path: Path):
    src = tmp_path / "sites"
    _write(src / "chr1" / "inherited_00000.tsv", "chr1\t10\t.\tA\tG\tp1\n")
    _write(src / "chr22" / "inherited_00000.tsv", "chr22\t10\t.\tA\tG\tp2\n")
    out = tmp_path / "out"
    invert_sites(src, out, chroms=["22"])
    assert not (out / "chr1").exists()
    assert _read(out / "chr22" / "patients.tsv") == [["p2", "1", "chr22_10_A_G"]]
    with pytest.raises(ValueError, match="chr5"):
        invert_sites(src, out, chroms=["chr5"])


@pytest.mark.parametrize(
    "body, message",
    [
        ("chr2\t10\t.\tA\tG\tp1\n", "expected chr1"),
        ("chr1\t10\t.\tA\tG\n", "6 tab-separated"),
        ("chr1\tx\t.\tA\tG\tp1\n", "non-integer POS"),
        ("chr1\t10\t.\tA\tG\tp1=0/1\n", "short format"),
    ],
)
def test_invert_rejects_bad_rows(tmp_path: Path, body: str, message: str):
    src = tmp_path / "sites"
    _write(src / "chr1" / "inherited_00000.tsv", body)
    with pytest.raises(ValueError, match=message):
        invert_sites(src, tmp_path / "out")


def test_invert_missing_files(tmp_path: Path):
    (tmp_path / "sites" / "chr1").mkdir(parents=True)
    with pytest.raises(ValueError, match="no inherited_"):
        invert_sites(tmp_path / "sites", tmp_path / "out")


def test_invert_cli(tmp_path: Path, capsys):
    src = tmp_path / "sites"
    _write(src / "chr1" / "denovo_00000.tsv", "chr1\t10\t.\tA\tG\tp1\n")
    out = tmp_path / "out"
    main(["invert", "--input-dir", str(src), "--output-dir", str(out), "--prefix", "denovo"])
    assert "1 patients, 1 variants, 1 patient-variant pairs" in capsys.readouterr().out
    assert _read(out / "chr1" / "patients.tsv") == [["p1", "1", "chr1_10_A_G"]]

    with pytest.raises(SystemExit):
        main(["invert", "--input-dir", str(tmp_path / "nope"), "--output-dir", str(out)])


def test_invert_refuses_foreign_output_tree(tmp_path: Path):
    src = tmp_path / "sites"
    _write(src / "chr1" / "inherited_00000.tsv", "chr1\t10\t.\tA\tG\tp1\n")
    vcf_tree = tmp_path / "vcfs"
    vcf = vcf_tree / "chr1" / "variants_1_1000.vcf"
    vcf.parent.mkdir(parents=True)
    vcf.write_text("chr1\t10\t.\tA\tG\n")
    with pytest.raises(ValueError, match="not an invert output"):
        invert_sites(src, vcf_tree)
    assert vcf.read_text() == "chr1\t10\t.\tA\tG\n"
    assert not (vcf_tree / "chr1" / "patients.tsv").exists()
    with pytest.raises(ValueError, match="not an invert output"):
        invert_sites(src, src)

    out = tmp_path / "out"
    invert_sites(src, out)
    invert_sites(src, out)  # re-running into its own output is fine
