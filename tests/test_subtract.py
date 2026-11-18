import json
from pathlib import Path

import pytest

from inherited.cli import main
from inherited.subtract import subtract_variants


SITES_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")


def _write_sites(path: Path, body: str) -> None:
    _write(path, SITES_HEADER + body)


def test_subtract_paired_file_only(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(
        dir1 / "chr21" / "variants_10000001_12500000.vcf",
        "chr21\t10000010\t.\tA\tG\n"
        "chr21\t10000020\trs1\tC\tT\n"
        "chr21\t10000030\t.\tG\tA\n",
    )
    _write(
        dir2 / "chr21" / "variants_10000001_12500000.vcf",
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "chr21\t10000020\tother\tC\tT\n",
    )
    _write(
        dir2 / "chr21" / "variants_12500001_15000000.vcf",
        "chr21\t10000030\t.\tG\tA\n",
    )
    _write(
        dir1 / "chr21" / "variants_12500001_15000000.vcf",
        "chr21\t12500010\t.\tA\tC\n",
    )
    _write(dir1 / "chr22" / "variants_1_1000.vcf", "chr22\t10\t.\tT\tA\n")
    _write(dir2 / "chrX" / "variants_1_1000.vcf", "chrX\t10\t.\tA\tG\n")
    _write(dir1 / "chr21" / "variants.vcf", "chr21\t1\t.\tA\tT\n")
    _write(dir1 / "notes" / "variants_1_1000.vcf", "chr21\t1\t.\tA\tT\n")
    _write_sites(
        tmp_path / "a-sites" / "chr21" / "inherited_10000001_12500000.tsv",
        "chr21\t10000010\t.\tA\tG\tp1;p2\n"
        "chr21\t10000020\trs1\tC\tT\tp3\n"
        "chr21\t10000030\t.\tG\tA\tp2;p4\n",
    )
    _write_sites(
        tmp_path / "a-sites" / "chr21" / "inherited_12500001_15000000.tsv",
        "chr21\t12500010\t.\tA\tC\tp5\n",
    )
    _write_sites(
        tmp_path / "a-sites" / "chr22" / "denovo_1_1000.tsv",
        "chr22\t10\t.\tT\tA\tp1\n",
    )
    _write_sites(
        tmp_path / "b-sites" / "chr21" / "inherited_10000001_12500000.tsv",
        "chr21\t10000010\t.\tA\tG\tp9\n",
    )

    output = tmp_path / "out"
    subtract_variants(dir1, dir2, output)

    assert (output / "chr21" / "variants_10000001_12500000.vcf").read_text(
        encoding="utf-8"
    ) == ("chr21\t10000010\t.\tA\tG\n" "chr21\t10000030\t.\tG\tA\n")
    assert (output / "chr21" / "variants_12500001_15000000.vcf").read_text(
        encoding="utf-8"
    ) == "chr21\t12500010\t.\tA\tC\n"
    assert (output / "chr22" / "variants_1_1000.vcf").read_text(encoding="utf-8") == (
        "chr22\t10\t.\tT\tA\n"
    )
    assert not (output / "chrX").exists()
    assert not (output / "chr21" / "variants.vcf").exists()
    assert not (output / "notes").exists()
    sites = tmp_path / "out-sites"
    assert (sites / "chr21" / "inherited_10000001_12500000.tsv").read_text(
        encoding="utf-8"
    ) == (
        SITES_HEADER
        + "chr21\t10000010\t.\tA\tG\tp1;p2\n"
        + "chr21\t10000030\t.\tG\tA\tp2;p4\n"
    )
    assert (sites / "chr22" / "denovo_1_1000.tsv").read_text(encoding="utf-8") == (
        SITES_HEADER + "chr22\t10\t.\tT\tA\tp1\n"
    )
    assert (sites / "patients.txt").read_text(encoding="utf-8") == "p1\np2\np4\np5\n"


def test_subtract_writes_empty_file(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    _write(dir2 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    _write_sites(
        tmp_path / "a-sites" / "chr21" / "inherited_1_1000.tsv",
        "chr21\t10\t.\tA\tG\tp1\n",
    )

    output = tmp_path / "out"
    subtract_variants(dir1, dir2, output)

    assert (output / "chr21" / "variants_1_1000.vcf").read_text(encoding="utf-8") == ""
    assert (tmp_path / "out-sites" / "chr21" / "inherited_1_1000.tsv").read_text(
        encoding="utf-8"
    ) == SITES_HEADER
    assert (tmp_path / "out-sites" / "patients.txt").read_text(encoding="utf-8") == ""


def test_subtract_keeps_different_allele(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    _write(dir2 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tT\n")
    _write_sites(
        tmp_path / "a-sites" / "chr21" / "mendelian_bad_1_1000.tsv",
        "chr21\t10\t.\tA\tG\tp7\n",
    )

    output = tmp_path / "out"
    subtract_variants(dir1, dir2, output)

    assert (output / "chr21" / "variants_1_1000.vcf").read_text(encoding="utf-8") == (
        "chr21\t10\t.\tA\tG\n"
    )
    assert (
        tmp_path / "out-sites" / "chr21" / "mendelian_bad_1_1000.tsv"
    ).read_text(encoding="utf-8") == (SITES_HEADER + "chr21\t10\t.\tA\tG\tp7\n")


def test_subtract_requires_sites_row_for_kept_variant(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    _write_sites(tmp_path / "a-sites" / "chr21" / "inherited_1_1000.tsv", "")
    dir2.mkdir()

    with pytest.raises(ValueError, match="missing from"):
        subtract_variants(dir1, dir2, tmp_path / "out")


def test_subtract_missing_input_dir(tmp_path: Path, capsys):
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "subtract",
                "--input-dir1",
                str(tmp_path / "missing"),
                "--input-dir2",
                str(tmp_path),
                "--output-dir",
                str(tmp_path / "out"),
            ]
        )
    assert exc.value.code == 1
    assert "input directory not found" in capsys.readouterr().err


def test_subtract_no_sites_writes_vcfs_only(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(
        dir1 / "chr21" / "variants_1_1000.vcf",
        "chr21\t10\t.\tA\tG\nchr21\t20\t.\tC\tT\n",
    )
    _write(dir2 / "chr21" / "variants_1_1000.vcf", "chr21\t20\t.\tC\tT\n")

    output = tmp_path / "out"
    main(
        [
            "subtract",
            "--input-dir1",
            str(dir1),
            "--input-dir2",
            str(dir2),
            "--output-dir",
            str(output),
            "--no-sites",
        ]
    )

    assert (output / "chr21" / "variants_1_1000.vcf").read_text(encoding="utf-8") == (
        "chr21\t10\t.\tA\tG\n"
    )
    assert not (tmp_path / "out-sites").exists()
    params = json.loads((output / "params.json").read_text(encoding="utf-8"))
    assert params.pop("package_version")
    assert params.pop("run_at")
    assert params == {
        "command": "subtract",
        "concatenate": False,
        "files": {"chr21": {"files1": 1, "files2": 1}},
        "input_dir1": str(dir1.resolve()),
        "input_dir2": str(dir2.resolve()),
        "matching": "per file: input_dir1 file against the input_dir2 file of the same name",
        "output_dir": str(output.resolve()),
        "sites": False,
    }


def test_subtract_writes_params_with_sites(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    _write_sites(
        tmp_path / "a-sites" / "chr21" / "inherited_1_1000.tsv",
        "chr21\t10\t.\tA\tG\tp1\n",
    )
    dir2.mkdir()

    output = tmp_path / "out"
    subtract_variants(dir1, dir2, output)

    params = json.loads((output / "params.json").read_text(encoding="utf-8"))
    assert params["sites"] is True


def test_subtract_fails_when_no_vcf_files(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "notes.txt", "x\n")
    dir2.mkdir()

    output = tmp_path / "out"
    with pytest.raises(ValueError, match="no variants_<start>_<end>.vcf files"):
        subtract_variants(dir1, dir2, output, sites=False)
    assert not output.exists()


def test_segment_files_need_concatenate(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_00000.vcf", "chr21\t10\t.\tA\tG\n")
    _write(dir2 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")

    output = tmp_path / "out"
    with pytest.raises(ValueError, match="use --concatenate"):
        subtract_variants(dir1, dir2, output, sites=False)
    assert not output.exists()

    # also when only the second tree is segment-numbered
    with pytest.raises(ValueError, match="use --concatenate"):
        subtract_variants(dir2, dir1, output, sites=False)


def _segment_trees(tmp_path: Path) -> tuple[Path, Path]:
    """Two trees of chr21 split by variant count at different places."""
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_00000.vcf", "chr21\t10\t.\tA\tG\nchr21\t20\t.\tC\tT\n")
    _write(dir1 / "chr21" / "variants_00001.vcf", "chr21\t30\t.\tG\tA\nchr21\t40\t.\tT\tC\n")
    _write(dir2 / "chr21" / "variants_00000.vcf", "chr21\t20\t.\tC\tT\n")
    _write(dir2 / "chr21" / "variants_00001.vcf", "chr21\t30\t.\tG\tA\nchr21\t50\t.\tA\tT\n")
    _write_sites(
        tmp_path / "a-sites" / "chr21" / "inherited_00000.tsv",
        "chr21\t10\t.\tA\tG\tp1\nchr21\t20\t.\tC\tT\tp2\n",
    )
    _write_sites(
        tmp_path / "a-sites" / "chr21" / "inherited_00001.tsv",
        "chr21\t30\t.\tG\tA\tp3\nchr21\t40\t.\tT\tC\tp4;p1\n",
    )
    return dir1, dir2


def test_concatenate_matches_across_segment_boundaries(tmp_path: Path):
    dir1, dir2 = _segment_trees(tmp_path)
    output = tmp_path / "out"

    groups = subtract_variants(dir1, dir2, output, concatenate=True)

    assert (output / "chr21" / "variants_00000.vcf").read_text(encoding="utf-8") == (
        "chr21\t10\t.\tA\tG\n"
    )
    assert (output / "chr21" / "variants_00001.vcf").read_text(encoding="utf-8") == (
        "chr21\t40\t.\tT\tC\n"
    )
    sites = tmp_path / "out-sites"
    assert (sites / "chr21" / "inherited_00000.tsv").read_text(encoding="utf-8") == (
        SITES_HEADER + "chr21\t10\t.\tA\tG\tp1\n"
    )
    assert (sites / "chr21" / "inherited_00001.tsv").read_text(encoding="utf-8") == (
        SITES_HEADER + "chr21\t40\t.\tT\tC\tp4;p1\n"
    )
    assert (sites / "patients.txt").read_text(encoding="utf-8") == "p1\np4\n"
    assert groups == {"chr21": {"files1": 2, "files2": 2}}
    params = json.loads((output / "params.json").read_text(encoding="utf-8"))
    assert params["concatenate"] is True
    assert params["files"] == {"chr21": {"files1": 2, "files2": 2}}
    assert params["matching"].startswith("per chromosome")


def test_concatenate_variant_in_a_different_segment(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_00000.vcf", "chr21\t10\t.\tA\tG\n")
    _write(dir1 / "chr21" / "variants_00001.vcf", "chr21\t90\t.\tA\tG\n")
    # tree 2 has chr21:10 in its second segment and chr21:90 in its first
    _write(dir2 / "chr21" / "variants_00000.vcf", "chr21\t90\t.\tA\tG\n")
    _write(dir2 / "chr21" / "variants_00001.vcf", "chr21\t10\t.\tA\tG\n")
    output = tmp_path / "out"

    subtract_variants(dir1, dir2, output, sites=False, concatenate=True)

    assert (output / "chr21" / "variants_00000.vcf").read_text(encoding="utf-8") == ""
    assert (output / "chr21" / "variants_00001.vcf").read_text(encoding="utf-8") == ""


def test_concatenate_old_tree_minus_position_tree(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_00000.vcf", "chr21\t10\t.\tA\tG\nchr21\t5000\t.\tC\tT\n")
    _write(dir2 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    _write(dir2 / "chr21" / "variants_1001_2000.vcf", "chr21\t1500\t.\tA\tG\n")
    output = tmp_path / "out"

    subtract_variants(dir1, dir2, output, sites=False, concatenate=True)

    assert sorted(p.name for p in (output / "chr21").iterdir()) == ["variants_00000.vcf"]
    assert (output / "chr21" / "variants_00000.vcf").read_text(encoding="utf-8") == (
        "chr21\t5000\t.\tC\tT\n"
    )


def test_concatenate_keeps_chrx_groups_apart(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chrX" / "variants_females_00000.vcf", "chrX\t10\t.\tA\tG\n")
    _write(dir1 / "chrX" / "variants_males_nonPar_00000.vcf", "chrX\t10\t.\tA\tG\n")
    _write(dir1 / "chrX" / "variants_males_par1_00000.vcf", "chrX\t20\t.\tA\tG\n")
    _write(dir2 / "chrX" / "variants_females_00003.vcf", "chrX\t10\t.\tA\tG\n")
    output = tmp_path / "out"

    groups = subtract_variants(dir1, dir2, output, sites=False, concatenate=True)

    assert (output / "chrX" / "variants_females_00000.vcf").read_text(encoding="utf-8") == ""
    assert (output / "chrX" / "variants_males_nonPar_00000.vcf").read_text(
        encoding="utf-8"
    ) == "chrX\t10\t.\tA\tG\n"
    assert (output / "chrX" / "variants_males_par1_00000.vcf").read_text(
        encoding="utf-8"
    ) == "chrX\t20\t.\tA\tG\n"
    assert groups == {
        "chrX/females": {"files1": 1, "files2": 1},
        "chrX/males_nonPar": {"files1": 1, "files2": 0},
        "chrX/males_par1": {"files1": 1, "files2": 0},
    }


def test_concatenate_rejects_mixed_naming(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_00000.vcf", "chr21\t10\t.\tA\tG\n")
    _write(dir1 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    dir2.mkdir()

    with pytest.raises(ValueError, match="mixes file naming styles"):
        subtract_variants(dir1, dir2, tmp_path / "out", sites=False, concatenate=True)
    assert not (tmp_path / "out").exists()


def test_by_name_now_includes_chrx_group_files(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chrX" / "variants_females_1_1000.vcf",
           "chrX\t10\t.\tA\tG\nchrX\t20\t.\tC\tT\n")
    _write(dir2 / "chrX" / "variants_females_1_1000.vcf", "chrX\t20\t.\tC\tT\n")
    _write(dir2 / "chrX" / "variants_males_nonPar_1_1000.vcf", "chrX\t10\t.\tA\tG\n")
    output = tmp_path / "out"

    subtract_variants(dir1, dir2, output, sites=False)

    assert (output / "chrX" / "variants_females_1_1000.vcf").read_text(encoding="utf-8") == (
        "chrX\t10\t.\tA\tG\n"
    )


def test_cli_concatenate(tmp_path: Path, capsys):
    dir1, dir2 = _segment_trees(tmp_path)
    output = tmp_path / "out"
    main(["subtract", "--input-dir1", str(dir1), "--input-dir2", str(dir2),
          "--output-dir", str(output), "--concatenate", "--no-sites"])
    assert (output / "chr21" / "variants_00001.vcf").read_text(encoding="utf-8") == (
        "chr21\t40\t.\tT\tC\n"
    )
    assert json.loads((output / "params.json").read_text(encoding="utf-8"))["concatenate"] is True
