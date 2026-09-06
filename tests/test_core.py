from pathlib import Path

import pytest

from inherited.af import is_rare, load_af_json
from inherited.analyze import get_position
from inherited.classify import (
    classify_father_son,
    classify_mother_son,
    classify_trio,
)
from inherited.families import (
    has_parent_id,
    load_family_column_map,
    load_family_relations,
    normalize_sex,
    person_labels_for_header,
    resolve_family_columns,
)
from inherited.genotype import (
    QualityFilters,
    get_good_site,
    is_good,
    is_hom_ref,
    is_mendelian_diploid,
    sample_gt_has_alt,
)
from inherited.xchrom import chrom_mode_for, is_x_chrom, is_y_chrom, male_x_bucket, x_region


FIXTURES = Path(__file__).parent / "fixtures"


def test_get_position_parses_without_sample_fields():
    assert get_position("chr22\t12345\tid\tA\tG\t.\t.\t.\tGT\t0/1\n") == 12345


def test_get_position_rejects_malformed_record():
    with pytest.raises(ValueError, match="Malformed VCF"):
        get_position("chr22\t12345")


def test_load_af_json_scalar_and_object(tmp_path):
    table = load_af_json(FIXTURES / "tiny_af.json")
    assert table["var_rare"] == 0.001
    assert is_rare(table, "var_common") is False
    assert is_rare(table, "var_rare") is True


def test_load_af_json_object_value(tmp_path):
    path = tmp_path / "af.json"
    path.write_text('{"k1": {"AF": 0.2, "AF_EUR": 0.01}, "k2": {"AF": 0.005}}')
    table = load_af_json(path)
    assert table["k1"] == 0.01
    assert table["k2"] == 0.005


def test_load_family_relations():
    rel = load_family_relations(FIXTURES / "families.tsv")
    assert rel.trio_cl["child1"] == ("ma1", "fa1")
    assert rel.trios_ids == [["child1", "fa1", "ma1"]]
    assert rel.family_size["fam1"] == 3
    assert "child1" in rel.female_children
    assert not rel.male_children
    assert rel.sample_to_person == {}


def test_load_family_relations_accepts_spark_column_aliases(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text(
        "ind_id\tfamily_id\tfather_id\tmother_id\tsex\n"
        "child1\tfam1\tfa1\tma1\tFemale\n"
        "fa1\tfam1\t0\t0\tMale\n"
        "ma1\tfam1\t0\t0\tFemale\n"
    )
    rel = load_family_relations(path)
    assert rel.trio_cl["child1"] == ("ma1", "fa1")
    assert rel.trios_ids == [["child1", "fa1", "ma1"]]
    assert rel.family_size["fam1"] == 3


def test_resolve_family_columns_rejects_ambiguous_aliases():
    with pytest.raises(ValueError, match="ambiguous"):
        resolve_family_columns(["spid", "ind_id", "sfid", "father", "mother", "sex"])


def test_family_map_selects_among_ambiguous_columns(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text(
        "spid\tind_id\tsfid\tfather\tmother\tsex\n"
        "wrong\tchild1\tfam1\tfa1\tma1\tFemale\n"
        "x\tfa1\tfam1\t0\t0\tMale\n"
        "y\tma1\tfam1\t0\t0\tFemale\n"
    )
    mapping = tmp_path / "map.json"
    mapping.write_text('{"spid": "ind_id"}')
    rel = load_family_relations(path, column_map=load_family_column_map(mapping))
    assert "child1" in rel.trio_cl
    assert "wrong" not in rel.trio_cl


def test_family_map_supports_custom_headers(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text(
        "IID\tFID\tPAT\tMAT\tsex\n"
        "child1\tfam1\tfa1\tma1\tFemale\n"
        "fa1\tfam1\t0\t0\tMale\n"
        "ma1\tfam1\t0\t0\tFemale\n"
    )
    mapping = tmp_path / "map.json"
    mapping.write_text(
        '{"spid": "IID", "sfid": "FID", "father": "PAT", "mother": "MAT"}'
    )
    rel = load_family_relations(path, column_map=load_family_column_map(mapping))
    assert rel.trio_cl["child1"] == ("ma1", "fa1")


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("fa1", True),
        ("0", False),
        ("", False),
        ("   ", False),
        ("False", False),
        ("false", False),
        ("FALSE", False),
    ],
)
def test_has_parent_id(value, expected):
    assert has_parent_id(value) is expected


def test_load_family_relations_maps_ind_id_to_sample_id(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text(
        "ind_id\tsample_id\tfamily_id\tfather_id\tmother_id\tsex\n"
        "SP_child\tSDSM-child\tfam1\tSP_fa\tSP_ma\tFemale\n"
        "SP_fa\tSDSM-fa\tfam1\t\t\tMale\n"
        "SP_ma\tSDSM-ma\tfam1\t\t\tFemale\n"
    )
    rel = load_family_relations(path)
    assert rel.trio_cl["SDSM-child"] == ("SDSM-ma", "SDSM-fa")
    assert rel.trios_ids == [["SDSM-child", "SDSM-fa", "SDSM-ma"]]
    assert "SDSM-child" in rel.female_children
    assert "SP_child" not in rel.trio_cl
    assert rel.trio["SP_child"] == ("SP_ma", "SP_fa")
    assert rel.family_size["fam1"] == 3
    assert rel.sample_to_person == {
        "SDSM-child": "SP_child",
        "SDSM-fa": "SP_fa",
        "SDSM-ma": "SP_ma",
    }


def test_load_family_relations_excludes_trio_without_sample_id(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text(
        "ind_id\tsample_id\tfamily_id\tfather_id\tmother_id\tsex\n"
        "SP_child\tSDSM-child\tfam1\tSP_fa\tSP_ma\tFemale\n"
        "SP_fa\t\tfam1\t\t\tMale\n"
        "SP_ma\tSDSM-ma\tfam1\t\t\tFemale\n"
        "SP_ok\tSDSM-ok\tfam2\tSP_fa2\tSP_ma2\tMale\n"
        "SP_fa2\tSDSM-fa2\tfam2\t\t\tMale\n"
        "SP_ma2\tSDSM-ma2\tfam2\t\t\tFemale\n"
    )
    rel = load_family_relations(path)
    assert set(rel.trio_cl) == {"SDSM-ok"}
    assert rel.trio_cl["SDSM-ok"] == ("SDSM-ma2", "SDSM-fa2")
    assert "SP_child" not in rel.trio_cl
    assert "SDSM-child" not in rel.trio_cl
    assert "SDSM-ok" in rel.male_children


def test_load_family_relations_rejects_conflicting_sample_ids(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text(
        "ind_id\tsample_id\tfamily_id\tfather_id\tmother_id\tsex\n"
        "SP_child\tSDSM-a\tfam1\tSP_fa\tSP_ma\tFemale\n"
        "SP_child\tSDSM-b\tfam1\tSP_fa\tSP_ma\tFemale\n"
        "SP_fa\tSDSM-fa\tfam1\t\t\tMale\n"
        "SP_ma\tSDSM-ma\tfam1\t\t\tFemale\n"
    )
    with pytest.raises(ValueError, match="Conflicting sample_id"):
        load_family_relations(path)


def test_load_family_relations_rejects_conflicting_person_ids(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text(
        "ind_id\tsample_id\tfamily_id\tfather_id\tmother_id\tsex\n"
        "SP_a\tSDSM-same\tfam1\tSP_fa\tSP_ma\tFemale\n"
        "SP_b\tSDSM-same\tfam1\tSP_fa\tSP_ma\tMale\n"
        "SP_fa\tSDSM-fa\tfam1\t\t\tMale\n"
        "SP_ma\tSDSM-ma\tfam1\t\t\tFemale\n"
    )
    with pytest.raises(ValueError, match="Conflicting person_id"):
        load_family_relations(path)


def test_person_labels_for_header_identity_when_map_empty():
    header = ["child1", "fa1", "ma1"]
    assert person_labels_for_header(header, {}) is header


def test_person_labels_for_header_maps_sample_ids():
    header = ["SDSM-child", "SDSM-fa", "SDSM-ma", "other"]
    labels = person_labels_for_header(
        header,
        {"SDSM-child": "SP_child", "SDSM-fa": "SP_fa", "SDSM-ma": "SP_ma"},
    )
    assert labels == ["SP_child", "SP_fa", "SP_ma", "other"]


def test_load_family_relations_skips_blank_or_false_parents(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text(
        "ind_id\tfamily_id\tfather_id\tmother_id\tsex\n"
        "child_ok\tfam1\tfa1\tma1\tFemale\n"
        "child_blank\tfam1\t\t\tMale\n"
        "child_false\tfam1\tFalse\tFalse\tFemale\n"
        "child_one\tfam1\t\tma1\tMale\n"
        "fa1\tfam1\t0\t0\tMale\n"
        "ma1\tfam1\t0\t0\tFemale\n"
    )
    rel = load_family_relations(path)
    assert set(rel.trio_cl) == {"child_ok"}
    assert "child_blank" not in rel.trio
    assert "child_false" not in rel.trio
    assert rel.trio["child_one"] == ("ma1", "")


def test_load_family_relations_splits_sex_and_skips_unknown():
    rel = load_family_relations(FIXTURES / "families_x.tsv")
    assert rel.female_children == {"girl1"}
    assert rel.male_children == {"boy1"}
    assert "unknown1" not in rel.trio_cl


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("1", "male"),
        ("Male", "male"),
        ("male", "male"),
        ("2", "female"),
        ("Female", "female"),
        ("female", "female"),
        (".", None),
        ("", None),
        ("other", None),
    ],
)
def test_normalize_sex(value, expected):
    assert normalize_sex(value) is expected


def test_is_good_rejects_low_dp():
    assert not is_good("0/1", "5", "2,3", "0,0,0,0", "30", 1)


def test_is_good_rejects_missing_ad():
    assert not is_good("0/1", "30", ".", "0,0,0,0", "30", 1)
    assert not is_good("0/1", "30", "15,.", "0,0,0,0", "30", 1)


def test_is_good_allows_homozygous_reference():
    assert is_good("0/0", "30", "30,0", "0,0,0,0", "30", 1)
    assert is_good("0/0", "30", "80,5,15", "0,0,0,0", "30", 1)


def test_is_good_diploid_hom00_ab_floor():
    assert is_good("0/0", "30", "27,3", "0,0,0,0", "30", 1, ab_hom00_min=0.9)
    assert not is_good("0/0", "30", "24,6", "0,0,0,0", "30", 1, ab_hom00_min=0.9)
    assert not is_good("0/0", "30", "80,5,15", "0,0,0,0", "30", 1, ab_hom00_min=0.9)
    assert is_good("0/2", "30", "80,5,15", "0,0,0,0", "30", 1, ab_hom00_min=0.9)
    assert is_good("0", "10", "8,2", "0,0,0,0", "30", 1, haploid=True, ab_hom00_min=0.9)


def test_is_good_diploid_het_ab_band():
    assert is_good("0/1", "30", "15,15", "0,0,0,0", "30", 1)
    assert not is_good("0/1", "30", "18,2", "0,0,0,0", "30", 1)
    assert not is_good("0/1", "30", "2,18", "0,0,0,0", "30", 1)


def test_is_good_diploid_hom_alt_ab_floor():
    assert is_good("1/1", "30", "3,27", "0,0,0,0", "30", 1)
    assert not is_good("1/1", "30", "6,24", "0,0,0,0", "30", 1)


def test_is_good_multiallelic_cleans_missing_ad_as_zero():
    assert not is_good("0/1", "30", "10,4,.,6", "0,0,0,0", "30", 2, clean_missing_ad_as_zero=True)
    assert not is_good("0/1", "30", "10,4,.,6", "0,0,0,0", "30", 1, clean_missing_ad_as_zero=True)
    assert is_good("0/3", "30", "10,4,.,6", "0,0,0,0", "30", 3, clean_missing_ad_as_zero=True)


def test_is_good_multiallelic_strict_without_clean():
    assert not is_good("0/3", "30", "10,4,.,6", "0,0,0,0", "30", 3)


def test_get_good_site_handles_missing_ad():
    sample = "0/1:30:.:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == -1


def test_get_good_site_homozygous_reference_returns_zero():
    sample = "0/0:30:30,0:0,0,0,0:30:0,30,30:."
    assert get_good_site(sample, 1) == (0, "0/0", "30")
    leaky = "0/0:30:24,6:0,0,0,0:30:0,30,30:."
    assert get_good_site(leaky, 1) == (0, "0/0", "30")
    assert get_good_site(leaky, 1, qc=QualityFilters(ab_hom00=0.9)) == (-1, ".", "0")
    assert get_good_site(sample, 1, qc=QualityFilters(ab_hom00=0.9)) == (0, "0/0", "30")


def test_get_good_site_skip_qc_if_no_alt_bypasses_depth():
    # Low DP would fail normal QC, but non-carrier children skip QC entirely.
    sample = "0/0:1:1,0:0,0,0,0:1:0,1,1:."
    assert get_good_site(sample, 1) == (-1, ".", "0")
    assert get_good_site(sample, 1, skip_qc_if_no_alt=True) == (0, "0/0", "1")


def test_get_good_site_skip_qc_if_no_alt_still_filters_carriers():
    sample = "0/1:5:2,3:0,0,0,0:30:0,30,30:."
    assert get_good_site(sample, 1, skip_qc_if_no_alt=True) == (-1, ".", "0")


def test_get_good_site_haploid_thresholds():
    low_ab = "1/1:10:5,5:0,0,0,0:30:0,30,30:."
    assert get_good_site(low_ab, 1, haploid=True) == (-1, ".", "0")
    high_ab = "1/1:10:1,9:0,0,0,0:30:0,30,30:."
    assert get_good_site(high_ab, 1, haploid=True) == (2, "1/1", "30")


def test_get_good_site_respects_custom_qc():
    sample = "0/1:30:15,15:0,0,0,0:25:0,30,30:."
    assert get_good_site(sample, 1) == (1, "0/1", "25")
    assert get_good_site(sample, 1, qc=QualityFilters(gq=30)) == (-1, ".", "0")
    unbalanced = "0/1:30:18,2:0,0,0,0:30:0,30,30:."
    assert get_good_site(unbalanced, 1) == (-1, ".", "0")
    assert get_good_site(unbalanced, 1, qc=QualityFilters(ab=0.1)) == (1, "0/1", "30")


def test_get_good_site_multiallelic_clean_ad():
    sample = "0/3:30:10,4,.,6:0,0,0,0:30:0,30,30:."
    assert get_good_site(sample, 2, clean_ad=True) == (-1, ".", "0")
    assert get_good_site(sample, 3, clean_ad=True) == (1, "0/3", "30")
    assert get_good_site(sample, 3, clean_ad=False) == (-1, ".", "0")


def test_sample_gt_has_alt_skips_hom_ref_without_format_fields():
    assert sample_gt_has_alt("0/0:1:1,0:0,0,0,0:1:0,1,1:.", 1) is False
    assert sample_gt_has_alt("0|0:30:30,0:0,0,0,0:30:0,30,30:.", 1) is False
    assert sample_gt_has_alt("0:5:5,0:0,0,0,0:20:0,20:.", 1) is False
    assert sample_gt_has_alt("0/0", 1) is False
    assert sample_gt_has_alt("./.:30:15,15:0,0,0,0:30:0,30,30:.", 1) is False


def test_sample_gt_has_alt_detects_carriers():
    assert sample_gt_has_alt("0/1:5:2,3:0,0,0,0:30:0,30,30:.", 1) is True
    assert sample_gt_has_alt("1|0:30:15,15:0,0,0,0:30:0,30,30:.", 1) is True
    assert sample_gt_has_alt("1/1:10:1,9:0,0,0,0:30:0,30,30:.", 1) is True
    assert sample_gt_has_alt("0/2:30:10,0,20:0,0,0,0:30:0,30,30:.", 1) is False
    assert sample_gt_has_alt("0/2:30:10,0,20:0,0,0,0:30:0,30,30:.", 2) is True
    assert sample_gt_has_alt("1/2:30:0,15,15:0,0,0,0:30:0,30,30:.", 1) is True
    assert sample_gt_has_alt("0/10:30:10,0,0,0,0,0,0,0,0,0,20:0,0,0,0:30:.", 1) is False
    assert sample_gt_has_alt("0/10:30:10,0,0,0,0,0,0,0,0,0,20:0,0,0,0:30:.", 10) is True


def test_get_good_site_counts_alt():
    sample = "0/1:30:15,15:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == 1
    assert gt == "0/1"
    assert gq == "30"


def test_get_good_site_returns_negative_one_when_not_good():
    sample = "0/1:5:2,3:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == -1
    assert gt == "."
    assert gq == "0"


@pytest.mark.parametrize(
    ("ac", "mac", "fac", "m_gt", "f_gt", "c_gt", "alt_index", "expected"),
    [
        (1, 2, 2, "1/1", "1/1", "0/1", 1, "mendelian_bad"),
        (2, 2, 2, "1/1", "1/1", "1/1", 1, "inherited"),
        (2, 1, 0, "0/1", "0/0", "1/1", 1, "mendelian_bad"),
        (1, 1, 0, "0/1", "0/0", "0/1", 1, "inherited"),
        (2, 0, 1, "0/0", "0/1", "1/1", 1, "mendelian_bad"),
        (1, 0, 1, "0/0", "0/1", "0/1", 1, "inherited"),
        (1, 0, 0, "0/0", "0/0", "0/1", 1, "denovo"),
        (2, 0, 0, "0/0", "0/0", "1/1", 1, "denovo"),
        (1, 0, 0, "0/0", "0/0", "1/2", 2, "mendelian_bad"),
        (1, 0, 0, "0/0", "0/0", "1/2", 1, "mendelian_bad"),
        (1, -1, 0, "0/1", "0/0", "0/1", 1, None),
        (1, 0, -1, "0/0", "0/1", "0/1", 1, None),
        (0, 0, 0, "0/0", "0/0", "0/0", 1, None),
        # Multiallelic: allele-2 counts look inherited, but child has unexplained allele 1.
        (1, 1, 1, "0/2", "0/2", "1/2", 2, "mendelian_bad"),
        (1, 1, 1, "0/2", "0/2", "0/2", 2, "inherited"),
        (1, 0, 0, "0/2", "0/2", "0/3", 3, None),
        (1, 1, 1, "0/2", "0/1", "1/2", 2, "inherited"),
    ],
)
def test_classify_trio(ac, mac, fac, m_gt, f_gt, c_gt, alt_index, expected):
    assert classify_trio(ac, mac, fac, m_gt, f_gt, c_gt, alt_index) == expected


@pytest.mark.parametrize(
    ("ac", "mac", "expected"),
    [
        (1, 1, "inherited"),
        (1, 0, "denovo"),
        (1, -1, None),
        (0, 0, None),
    ],
)
def test_classify_mother_son(ac, mac, expected):
    assert classify_mother_son(ac, mac) == expected


@pytest.mark.parametrize(
    ("ac", "fac", "expected"),
    [
        (1, 1, "inherited"),
        (1, 0, "denovo"),
        (1, -1, None),
        (0, 0, None),
    ],
)
def test_classify_father_son(ac, fac, expected):
    assert classify_father_son(ac, fac) == expected


@pytest.mark.parametrize(
    ("gt", "expected"),
    [
        ("0/0", True),
        ("0|0", True),
        ("0", True),
        ("0/1", False),
        ("0/2", False),
        ("1/1", False),
        (".", False),
        ("0/.", False),
    ],
)
def test_is_hom_ref(gt, expected):
    assert is_hom_ref(gt) is expected


@pytest.mark.parametrize(
    ("m_gt", "f_gt", "c_gt", "expected"),
    [
        ("0/0", "0/1", "0/1", True),
        ("0/1", "0/0", "0/1", True),
        ("0/1", "0/1", "0/1", True),
        ("0/1", "0/1", "1/1", True),
        ("1/1", "1/1", "1/1", True),
        ("1/1", "1/1", "0/1", False),
        ("0/0", "0/0", "0/1", False),
        ("0/2", "0/2", "0/2", True),
        ("0/2", "0/2", "1/2", False),
        ("0/2", "0/1", "1/2", True),
        ("0/1", ".", "0/1", False),
        ("0", "0/1", "0/1", False),
    ],
)
def test_is_mendelian_diploid(m_gt, f_gt, c_gt, expected):
    assert is_mendelian_diploid(m_gt, f_gt, c_gt) is expected


def test_x_region_and_buckets():
    assert is_x_chrom("X")
    assert is_x_chrom("chrX")
    assert not is_x_chrom("22")
    assert is_y_chrom("Y")
    assert is_y_chrom("chrY")
    assert chrom_mode_for("chrX") == "x"
    assert chrom_mode_for("chrY") == "y"
    assert chrom_mode_for("22") == "autosomal"
    assert x_region(10001) == "par1"
    assert x_region(2781479) == "par1"
    assert x_region(2781480) == "nonPar"
    assert x_region(155701383) == "par2"
    assert male_x_bucket(15000) == "males_par1"
    assert male_x_bucket(5_000_000) == "males_nonPar"
    assert male_x_bucket(155800000) == "males_par2"
