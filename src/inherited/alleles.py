"""Trim alleles and left-align indels.

Trimming removes a shared prefix and suffix and does not use the reference.
Left-alignment follows the bcftools ``norm`` realignment: a shared suffix may
be removed until an allele is empty, the reference base before the variant is
then prepended, and a shared prefix is removed while both alleles stay longer
than one base. Equal-length alleles never need the reference.
"""

from __future__ import annotations

from pathlib import Path


def validate_allele_mode(
    *,
    trim: bool,
    keep_indels: bool,
    fasta: Path | None,
    multiallelic: bool,
) -> None:
    if (trim or keep_indels) and not multiallelic:
        raise ValueError("--trim and --keep-indels cannot be used with --no-multiallelic")
    if keep_indels and fasta is None:
        raise ValueError("--keep-indels requires --fasta")
    if fasta is not None and not keep_indels:
        raise ValueError("--fasta requires --keep-indels")


def is_symbolic(allele: str) -> bool:
    """True for ``*``, ``<DEL>``, and breakend alleles."""
    return allele == "*" or allele.startswith("<") or "[" in allele or "]" in allele


def trim_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str] | None:
    """Trim a shared prefix and suffix. A prefix trim shifts pos.

    Returns None when the remaining REF and ALT are identical.
    """
    ref_bases = ref.upper()
    alt_bases = alt.upper()
    while len(ref_bases) > 1 and len(alt_bases) > 1 and ref_bases[-1] == alt_bases[-1]:
        ref_bases = ref_bases[:-1]
        alt_bases = alt_bases[:-1]
    while len(ref_bases) > 1 and len(alt_bases) > 1 and ref_bases[0] == alt_bases[0]:
        ref_bases = ref_bases[1:]
        alt_bases = alt_bases[1:]
        pos += 1
    if ref_bases == alt_bases:
        return None
    return pos, ref_bases, alt_bases


def realign_left(
    pos: int,
    ref: str,
    alt: str,
    base_at,
) -> tuple[int, str, str] | None:
    """Return the left-aligned allele, or None when it matches the reference.

    ``base_at(position)`` returns the uppercase reference base at a 1-based
    position, or None when that position cannot be read.
    """
    ref_bases = ref.upper()
    alt_bases = alt.upper()
    steps = pos + len(ref_bases) + len(alt_bases) + 2
    for _ in range(steps):
        if ref_bases == alt_bases:
            return None
        if not ref_bases or not alt_bases or ref_bases[-1] != alt_bases[-1]:
            break
        if min(len(ref_bases), len(alt_bases)) <= 1 and pos == 1:
            break
        ref_bases = ref_bases[:-1]
        alt_bases = alt_bases[:-1]
        if ref_bases and alt_bases:
            continue
        if pos <= 1:
            return None
        pos -= 1
        base = base_at(pos)
        if not base:
            return None
        ref_bases = base + ref_bases
        alt_bases = base + alt_bases
    else:
        raise ValueError(f"Left-alignment did not finish at position {pos}")

    while len(ref_bases) > 1 and len(alt_bases) > 1 and ref_bases[0] == alt_bases[0]:
        ref_bases = ref_bases[1:]
        alt_bases = alt_bases[1:]
        pos += 1
    if not ref_bases or not alt_bases or ref_bases == alt_bases:
        return None
    return pos, ref_bases, alt_bases


def prepare_allele(
    pos: int,
    ref: str,
    alt: str,
    *,
    keep_indels: bool,
    base_at=None,
) -> tuple[int, str, str] | None:
    """Trim, and left-align when ``keep_indels`` is set.

    ``--trim`` keeps a one-base allele. ``--keep-indels`` also keeps an indel
    and a real multi-nucleotide variant.
    """
    if is_symbolic(ref) or is_symbolic(alt):
        return None
    if keep_indels:
        normalized = realign_left(pos, ref, alt, base_at)
    else:
        normalized = trim_alleles(pos, ref, alt)
    if normalized is None:
        return None
    new_pos, new_ref, new_alt = normalized
    if len(new_ref) == 1 and len(new_alt) == 1:
        return normalized
    if keep_indels:
        return normalized
    return None


def allele_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Frequency key ``chr22_12345_A_G``. ``22`` and ``chr22`` both match."""
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.upper() == "MT":
        name = "M"
    return f"chr{name}_{pos}_{ref}_{alt}"


class ReferenceFasta:
    """Random access to an uncompressed FASTA via its ``.fai`` index."""

    def __init__(self, path: Path):
        if not path.is_file():
            raise FileNotFoundError(f"FASTA not found: {path}")
        if path.name.endswith(".gz"):
            raise ValueError(
                "--fasta must be an uncompressed FASTA indexed with samtools faidx"
            )
        index_path = Path(str(path) + ".fai")
        if not index_path.is_file():
            raise FileNotFoundError(f"FASTA index not found: {index_path}")
        self.path = path
        self._contigs: dict[str, tuple[int, int, int, int]] = {}
        for line in index_path.read_text(encoding="utf-8").splitlines():
            if not line:
                continue
            name, length, offset, bases, width = line.split("\t")
            self._contigs[name] = (int(length), int(offset), int(bases), int(width))
        self._handle = path.open("rb")
        self._window_chrom = ""
        self._window_start = 0
        self._window = ""

    def close(self) -> None:
        self._handle.close()

    def base(self, chrom: str, pos: int) -> str | None:
        """Return the uppercase base at a 1-based position."""
        if pos < 1:
            return None
        name = self._resolve(chrom)
        length = self._contigs[name][0]
        if pos > length:
            return None
        if (
            name == self._window_chrom
            and self._window_start <= pos < self._window_start + len(self._window)
        ):
            return self._window[pos - self._window_start]
        begin = max(1, pos - 1023)
        self._window = self._fetch(name, begin, pos)
        self._window_chrom = name
        self._window_start = begin
        return self._window[-1]

    def _resolve(self, chrom: str) -> str:
        if chrom in self._contigs:
            return chrom
        if chrom.lower().startswith("chr"):
            bare = chrom[3:]
            if bare in self._contigs:
                return bare
        else:
            prefixed = "chr" + chrom
            if prefixed in self._contigs:
                return prefixed
            if chrom.upper() == "M" and "chrMT" in self._contigs:
                return "chrMT"
            if chrom.upper() == "MT" and "chrM" in self._contigs:
                return "chrM"
        raise ValueError(f"Chromosome {chrom} is not in {self.path}")

    def _fetch(self, name: str, start: int, end: int) -> str:
        """Return bases at 1-based positions ``start`` through ``end``."""
        _length, offset, bases_per_line, bytes_per_line = self._contigs[name]
        chars: list[str] = []
        pos0 = start - 1
        last = end
        while pos0 < last:
            line = pos0 // bases_per_line
            within = pos0 % bases_per_line
            take = min(bases_per_line - within, last - pos0)
            self._handle.seek(offset + line * bytes_per_line + within)
            chars.append(self._handle.read(take).decode("ascii"))
            pos0 += take
        return "".join(chars).upper()
