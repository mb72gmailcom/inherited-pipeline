from pathlib import Path

import pytest

from inherited.cli import main
from inherited.subtract import subtract_variants


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")


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


def test_subtract_writes_empty_file(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    _write(dir2 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")

    output = tmp_path / "out"
    subtract_variants(dir1, dir2, output)

    assert (output / "chr21" / "variants_1_1000.vcf").read_text(encoding="utf-8") == ""


def test_subtract_keeps_different_allele(tmp_path: Path):
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    _write(dir1 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tG\n")
    _write(dir2 / "chr21" / "variants_1_1000.vcf", "chr21\t10\t.\tA\tT\n")

    output = tmp_path / "out"
    subtract_variants(dir1, dir2, output)

    assert (output / "chr21" / "variants_1_1000.vcf").read_text(encoding="utf-8") == (
        "chr21\t10\t.\tA\tG\n"
    )


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
