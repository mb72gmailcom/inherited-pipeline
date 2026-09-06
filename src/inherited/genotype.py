from __future__ import annotations

from dataclasses import asdict, dataclass

from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_AB_HOM,
    DEFAULT_DP,
    DEFAULT_GQ,
    DEFAULT_HAPLO_AB,
    DEFAULT_HAPLO_DP,
)


@dataclass(frozen=True)
class QualityFilters:
    """Per-run genotype QC thresholds. Defaults match ``constants.py``.

    ``ab_hom00`` is optional; ``None`` means diploid ``0/0`` skips AB.
    """

    gq: int = DEFAULT_GQ
    dp: int = DEFAULT_DP
    ab: float = DEFAULT_AB
    ab_hom: float = DEFAULT_AB_HOM
    ab_hom00: float | None = None
    haplo_dp: int = DEFAULT_HAPLO_DP
    haplo_ab: float = DEFAULT_HAPLO_AB

    def as_params(self) -> dict[str, int | float]:
        return asdict(self)


DEFAULT_QUALITY = QualityFilters()


def _parse_int_field(value: str) -> int | None:
    value = value.strip()
    if value in ("", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _parse_ad_field(ad: str, clean_missing_as_zero: bool = False) -> list[int] | None:
    if not ad or ad.strip() == ".":
        return None

    values: list[int] = []
    for part in ad.split(","):
        part = part.strip()
        if part in ("", "."):
            if clean_missing_as_zero:
                values.append(0)
                continue
            return None
        try:
            values.append(int(part))
        except ValueError:
            return None
    return values or None


def _gt_alleles(gt: str) -> list[str]:
    return gt.replace("|", "/").split("/")


def _diploid_alleles(gt: str) -> tuple[str, str] | None:
    """Parse a diploid GT into two allele strings, or None if unusable."""
    # Common unphased/phased single-digit forms: "0/1", "0|1".
    if len(gt) == 3 and gt[1] in "/|":
        left, right = gt[0], gt[2]
        if left == "." or right == ".":
            return None
        return left, right

    alleles = _gt_alleles(gt)
    if len(alleles) != 2 or alleles[0] == "." or alleles[1] == ".":
        return None
    return alleles[0], alleles[1]


def is_hom_ref(gt: str) -> bool:
    """Return True when every allele in ``gt`` is reference (``0``)."""
    if gt in ("0/0", "0|0", "0"):
        return True
    alleles = _gt_alleles(gt)
    return bool(alleles) and all(allele == "0" for allele in alleles)


def alleles_only_ref_and_alt(gt: str, alt_token: str) -> bool:
    """True when GT uses only reference and the queried alt (no other alts)."""
    alleles = _diploid_alleles(gt)
    if alleles is None:
        return False
    return all(allele == "0" or allele == alt_token for allele in alleles)


def is_mendelian_diploid(m_gt: str, f_gt: str, c_gt: str) -> bool:
    """Return True if the child diploid genotype is Mendelian given both parents.

    Checks whether there exist one maternal and one paternal allele that together
    match the child's genotype (order-independent).
    """
    mother = _diploid_alleles(m_gt)
    father = _diploid_alleles(f_gt)
    child = _diploid_alleles(c_gt)
    if mother is None or father is None or child is None:
        return False

    c0, c1 = child
    for maternal in mother:
        for paternal in father:
            if (maternal == c0 and paternal == c1) or (
                maternal == c1 and paternal == c0
            ):
                return True
    return False


def is_good(
    gt: str,
    dp: str,
    ad: str,
    sb: str,
    gq: str,
    alt_index: int = 1,
    *,
    clean_missing_ad_as_zero: bool = False,
    dp_min: int = DEFAULT_DP,
    ab_min: float = DEFAULT_AB,
    ab_hom_min: float = DEFAULT_AB_HOM,
    ab_hom00_min: float | None = None,
    gq_min: int = DEFAULT_GQ,
    alleles: list[str] | None = None,
    ac: int | None = None,
    haploid: bool = False,
) -> bool:
    if "." in gt:
        return False

    dp_value = _parse_int_field(dp)
    if dp_value is None or dp_value < dp_min:
        return False

    gq_value = _parse_int_field(gq)
    if gq_value is None or gq_value < gq_min:
        return False

    ads = _parse_ad_field(ad, clean_missing_as_zero=clean_missing_ad_as_zero)
    if ads is None or len(ads) <= alt_index:
        return False
    total = sum(ads)
    if total == 0:
        return False

    parsed = alleles if alleles is not None else _gt_alleles(gt)
    if ac is None:
        ac = sum(allele.isdigit() and int(allele) == alt_index for allele in parsed)
    if ac <= 0:
        if (
            ab_hom00_min is not None
            and not haploid
            and parsed
            and all(allele == "0" for allele in parsed)
        ):
            return ads[0] / total >= ab_hom00_min
        return True

    ab = ads[alt_index] / total
    if haploid:
        return ab >= ab_min
    if ac == 1:
        return ab_min <= ab <= 1.0 - ab_min
    return ab >= ab_hom_min


def sample_gt_has_alt(sample_field: str, alt_index: int) -> bool:
    """Return True if the GT field carries ``alt_index``. Ignores DP/AD/GQ.

    Used to skip non-carrier children before ``get_good_site``. The common
    homozygous-reference forms are rejected without splitting FORMAT.
    """
    if sample_field.startswith(("0/0:", "0|0:", "0:")):
        return False

    colon = sample_field.find(":")
    gt = sample_field if colon < 0 else sample_field[:colon]
    if not gt or "." in gt or gt in ("0/0", "0|0", "0"):
        return False
    if len(gt) == 3 and gt[1] in "/|":
        if 0 <= alt_index <= 9:
            token = "0123456789"[alt_index]
            return gt[0] == token or gt[2] == token
        return False
    if len(gt) == 1:
        return gt.isdigit() and int(gt) == alt_index
    return any(
        allele.isdigit() and int(allele) == alt_index for allele in _gt_alleles(gt)
    )


def get_good_site(
    sample_field: str,
    alt_index: int = 1,
    *,
    clean_ad: bool = False,
    haploid: bool = False,
    skip_qc_if_no_alt: bool = False,
    qc: QualityFilters = DEFAULT_QUALITY,
) -> tuple[int, str, str]:
    """Return allele count, GT, and GQ for one alternate allele.

    Returns ``-1`` when the site fails quality filters or FORMAT is incomplete.
    Returns ``0`` for a good-quality homozygous-reference call.

    When ``clean_ad`` is True, missing AD values (``.``) are treated as zero depth
    for that allele. Used on the multiallelic path so one missing ALT depth does
    not reject all alleles at the site.

    When ``haploid`` is True, use haploid depth/AB thresholds while still counting
    alternate alleles from the GT field. Diploid sites use a het AB band
    ``[ab_min, 1 - ab_min]`` and a separate homozygous-alt floor ``ab_hom_min``.
    When ``qc.ab_hom00`` is set, diploid ``0/0`` also requires
    ``AD[0] / sum(AD) >= ab_hom00``. Haploid ref and other ``ac==0`` genotypes
    skip that check.

    When ``skip_qc_if_no_alt`` is True and the genotype carries no copies of the
    queried alt, return ``(0, gt, gq)`` without DP/AD/GQ checks. Used for child
    lookups where non-carriers are discarded immediately.
    """
    parts = sample_field.split(":")
    if len(parts) < 6:
        return -1, ".", "0"

    gt, dp, ad, sb, gq = parts[0], parts[1], parts[2], parts[3], parts[4]
    if "." in gt:
        return -1, ".", "0"

    alleles = _gt_alleles(gt)
    ac = sum(allele.isdigit() and int(allele) == alt_index for allele in alleles)
    if skip_qc_if_no_alt and ac == 0:
        return 0, gt, gq

    dp_min = qc.haplo_dp if haploid else qc.dp
    ab_min = qc.haplo_ab if haploid else qc.ab
    if is_good(
        gt,
        dp,
        ad,
        sb,
        gq,
        alt_index,
        clean_missing_ad_as_zero=clean_ad,
        dp_min=dp_min,
        ab_min=ab_min,
        ab_hom_min=qc.ab_hom,
        ab_hom00_min=qc.ab_hom00,
        gq_min=qc.gq,
        alleles=alleles,
        ac=ac,
        haploid=haploid,
    ):
        return ac, gt, gq
    return -1, ".", "0"
