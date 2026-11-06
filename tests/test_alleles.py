from pathlib import Path

import pytest

from inherited.alleles import (
    ReferenceFasta,
    allele_key,
    is_symbolic,
    prepare_allele,
    realign_left,
    trim_alleles,
    validate_allele_mode,
)

# 1:G 2:C 3-12:T 13:G 14-23:A
_SEQUENCE = "GC" + "T" * 10 + "G" + "A" * 10


def _base_at(pos: int) -> str | None:
    if pos < 1 or pos > len(_SEQUENCE):
        return None
    return _SEQUENCE[pos - 1]


def test_trim_keeps_a_one_base_change_and_a_real_mnv() -> None:
    assert trim_alleles(14, "AA", "AG") == (15, "A", "G")
    assert trim_alleles(14, "AA", "GG") == (14, "AA", "GG")
    assert trim_alleles(14, "AA", "AA") is None


def test_realign_left_moves_a_homopolymer_indel() -> None:
    assert realign_left(11, "TT", "T", _base_at) == (2, "CT", "C")
    assert realign_left(12, "T", "TT", _base_at) == (2, "C", "CT")


def test_realign_left_leaves_substitutions_in_place() -> None:
    assert realign_left(14, "AA", "AG", _base_at) == (15, "A", "G")
    assert realign_left(14, "AA", "GG", _base_at) == (14, "AA", "GG")
    assert realign_left(13, "GA", "G", _base_at) == (13, "GA", "G")


def test_prepare_allele_keeps_indels_only_when_asked() -> None:
    assert prepare_allele(14, "AA", "AG", keep_indels=False) == (15, "A", "G")
    assert prepare_allele(14, "AA", "GG", keep_indels=False) is None
    assert prepare_allele(13, "GA", "G", keep_indels=False) is None
    assert prepare_allele(14, "AA", "GG", keep_indels=True, base_at=_base_at) == (14, "AA", "GG")
    assert prepare_allele(11, "TT", "T", keep_indels=True, base_at=_base_at) == (2, "CT", "C")
    assert prepare_allele(14, "<DEL>", "A", keep_indels=True, base_at=_base_at) is None
    assert is_symbolic("A[chr2:1[")


def test_allele_key_uses_a_chr_prefix() -> None:
    assert allele_key("22", 15, "A", "G") == "chr22_15_A_G"
    assert allele_key("chr22", 15, "A", "G") == "chr22_15_A_G"
    assert allele_key("MT", 5, "A", "G") == "chrM_5_A_G"


def test_validate_allele_mode_rejects_illegal_combinations() -> None:
    with pytest.raises(ValueError, match="--no-multiallelic"):
        validate_allele_mode(trim=True, keep_indels=False, fasta=None, multiallelic=False)
    with pytest.raises(ValueError, match="requires --fasta"):
        validate_allele_mode(trim=False, keep_indels=True, fasta=None, multiallelic=True)
    with pytest.raises(ValueError, match="requires --keep-indels"):
        validate_allele_mode(
            trim=False, keep_indels=False, fasta=Path("ref.fa"), multiallelic=True
        )


def _write_fasta(path: Path, name: str, sequence: str, width: int = 10) -> None:
    blob = bytearray(f">{name}\n".encode())
    offset = len(blob)
    for start in range(0, len(sequence), width):
        blob += sequence[start : start + width].encode() + b"\n"
    path.write_bytes(blob)
    Path(str(path) + ".fai").write_text(
        f"{name}\t{len(sequence)}\t{offset}\t{width}\t{width + 1}\n",
        encoding="utf-8",
    )


def test_reference_fasta_reads_across_wrapped_lines(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fa"
    _write_fasta(fasta, "chr22", _SEQUENCE, width=10)
    genome = ReferenceFasta(fasta)
    try:
        assert genome.base("22", 1) == "G"
        assert genome.base("22", 11) == "T"
        assert genome.base("chr22", 13) == "G"
        assert realign_left(11, "TT", "T", lambda pos: genome.base("22", pos)) == (2, "CT", "C")
    finally:
        genome.close()


def test_reference_fasta_requires_an_uncompressed_index(tmp_path: Path) -> None:
    gzip_path = tmp_path / "ref.fa.gz"
    gzip_path.write_bytes(b"")
    with pytest.raises(ValueError, match="uncompressed"):
        ReferenceFasta(gzip_path)
    plain = tmp_path / "ref.fa"
    plain.write_text(">chr22\nAC\n", encoding="utf-8")
    with pytest.raises(FileNotFoundError, match="index"):
        ReferenceFasta(plain)
