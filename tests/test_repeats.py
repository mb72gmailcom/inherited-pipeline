from pathlib import Path

from inherited.repeats import RepeatIntervalFilter

FIXTURES = Path(__file__).parent / "fixtures"


def test_in_repeat_inside_interval():
    filt = RepeatIntervalFilter(FIXTURES / "tiny_repeats.bed")
    try:
        assert filt.in_repeat(2501) is True
        assert filt.in_repeat(3000) is True
        assert filt.in_repeat(3499) is True
        assert filt.in_repeat(3500) is True
    finally:
        filt.close()


def test_in_repeat_outside_interval():
    filt = RepeatIntervalFilter(FIXTURES / "tiny_repeats.bed")
    try:
        assert filt.in_repeat(2499) is False
        assert filt.in_repeat(2500) is False
        assert filt.in_repeat(3501) is False
        assert filt.in_repeat(4000) is False
    finally:
        filt.close()


def test_advance_past_skips_intervals_before_position(tmp_path):
    path = tmp_path / "repeats.bed"
    path.write_text("22\t0\t1000\n22\t2000\t3000\n22\t5000\t6000\n", encoding="utf-8")
    filt = RepeatIntervalFilter(path)
    try:
        filt.advance_past(4500)
        assert filt.in_repeat(2500) is False
        assert filt.in_repeat(4500) is False
        assert filt.in_repeat(5500) is True
    finally:
        filt.close()


def test_empty_file_never_matches(tmp_path):
    path = tmp_path / "empty.bed"
    path.write_text("", encoding="utf-8")
    filt = RepeatIntervalFilter(path)
    try:
        assert filt.in_repeat(1000) is False
    finally:
        filt.close()
