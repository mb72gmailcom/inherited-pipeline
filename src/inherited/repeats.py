from __future__ import annotations

from pathlib import Path


class RepeatIntervalFilter:
    """Stream repeat intervals and test whether a position falls inside one.

    Interval files are 0-based BED-style ``chrom start end`` rows with
    half-open intervals ``[start, end)``. ``pos`` arguments are 1-based VCF
    coordinates.
    """

    def __init__(self, path: Path) -> None:
        self._handle = path.open(encoding="utf-8")
        self._start = 0
        self._end = 0
        self._exhausted = False
        self._load_next()

    def close(self) -> None:
        self._handle.close()

    def _load_next(self) -> None:
        for line in self._handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            self._start = int(parts[1])
            self._end = int(parts[2])
            return
        self._exhausted = True
        self._start = self._end = 2**63

    def advance_past(self, pos: int) -> None:
        """Drop intervals that end at or before 1-based ``pos``."""
        zero_based = pos - 1
        while not self._exhausted and zero_based >= self._end:
            self._load_next()

    def in_repeat(self, pos: int) -> bool:
        """Return whether 1-based VCF ``pos`` falls in a 0-based BED interval."""
        self.advance_past(pos)
        if self._exhausted:
            return False
        zero_based = pos - 1
        return self._start <= zero_based < self._end
