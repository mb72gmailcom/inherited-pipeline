import json
from pathlib import Path

import pytest

from inherited.analyze import analyze_vcf
from inherited.checkpoint import Checkpoint, CumulativeStats, load_checkpoint, save_checkpoint
from inherited.output import glob_result_tsvs, read_result_tsv, serialize_payload
from inherited.shards import VcfShard

FIXTURES = Path(__file__).parent / "fixtures"
VCF_HEADER = "".join(
    line
    for line in (FIXTURES / "tiny.vcf").read_text(encoding="utf-8").splitlines(True)
    if line.startswith("#")
)


def _write_shard(path: Path, record_lines: list[str]) -> None:
    path.write_text(VCF_HEADER + "".join(record_lines), encoding="utf-8")


def _tiny_record_lines() -> list[str]:
    return [
        line
        for line in (FIXTURES / "tiny.vcf").read_text(encoding="utf-8").splitlines(True)
        if not line.startswith("#")
    ]


def test_analyze_vcf_writes_short_format_single_file(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        multiallelic=True,
        block_size=1,
        segment_size=0,
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=True)
    bad_records = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv", short_format=True)
    denovo_records = read_result_tsv(tmp_path / "out" / "denovo.tsv", short_format=True)

    assert len(inherited_records) == 1
    chrom, pos, ref, alt, patient_ids = inherited_records[0]
    assert chrom == "22" and pos == "3000" and patient_ids == ["child1"]
    assert len(bad_records) == 1
    assert len(denovo_records) == 1
    assert denovo_records[0][1] == "1000"
    assert denovo_records[0][4] == ["child1"]
    assert stats.inherited_variants == 1
    assert stats.denovo_variants == 1
    assert json.loads((tmp_path / "out" / "inherited_per_variant.json").read_text())["var_inh"] == 1
    assert json.loads((tmp_path / "out" / "denovo_per_variant.json").read_text())["var_rare"] == 1
    assert json.loads((tmp_path / "out" / "denovo_per_person.json").read_text())["child1"] == 1


def test_analyze_vcf_file_suffix_labels_outputs_and_skips_segmentation(tmp_path):
    out = tmp_path / "out"
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        multiallelic=True,
        block_size=1,
        file_suffix="1_2500",
    )

    suffix = "1_2500"
    inherited_records = read_result_tsv(out / f"inherited_{suffix}.tsv", short_format=True)
    denovo_records = read_result_tsv(out / f"denovo_{suffix}.tsv", short_format=True)
    assert inherited_records[0][1] == "3000"
    assert denovo_records[0][1] == "1000"
    assert stats.inherited_variants == 1
    assert stats.denovo_variants == 1
    assert not (out / "inherited.tsv").exists()
    assert not (out / "inherited_00000.tsv").exists()
    assert not (out / "stats.json").exists()
    assert json.loads((out / f"inherited_per_variant_{suffix}.json").read_text())["var_inh"] == 1
    assert json.loads((out / f"stats_{suffix}.json").read_text())["inherited_variants"] == 1


def test_analyze_vcf_file_suffixes_do_not_overwrite_in_shared_dir(tmp_path):
    out = tmp_path / "out"
    analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        segment_size=0,
        file_suffix="1_2500",
    )
    analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        segment_size=0,
        file_suffix="2501_5000",
    )

    assert (out / "inherited_1_2500.tsv").is_file()
    assert (out / "inherited_2501_5000.tsv").is_file()
    assert (out / "stats_1_2500.json").is_file()
    assert (out / "stats_2501_5000.json").is_file()
    assert not (out / "stats.json").exists()


def test_analyze_vcf_rejects_file_suffix_with_shards(tmp_path):
    shard = VcfShard(path=FIXTURES / "tiny.vcf", chrom="22", start=1, end=5000)
    with pytest.raises(ValueError, match="--file-suffix cannot be used with --vcf-dir"):
        analyze_vcf(
            vcf_path=None,
            af_json_path=FIXTURES / "tiny_af.json",
            family_file=FIXTURES / "families.tsv",
            output_dir=tmp_path / "out",
            vcf_shards=[shard],
            file_suffix="1_5000",
        )


def test_analyze_vcf_labels_hits_with_person_id(tmp_path):
    family = tmp_path / "families.tsv"
    family.write_text(
        "ind_id\tsample_id\tfamily_id\tfather_id\tmother_id\tsex\n"
        "SP_child\tchild1\tfam1\tSP_fa\tSP_ma\tFemale\n"
        "SP_fa\tfa1\tfam1\t0\t0\tMale\n"
        "SP_ma\tma1\tfam1\t0\t0\tFemale\n"
    )
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=family,
        output_dir=tmp_path / "out",
        multiallelic=True,
        block_size=1,
        segment_size=0,
    )
    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=True)
    denovo_records = read_result_tsv(tmp_path / "out" / "denovo.tsv", short_format=True)
    assert inherited_records[0][4] == ["SP_child"]
    assert denovo_records[0][4] == ["SP_child"]
    assert json.loads((tmp_path / "out" / "denovo_per_person.json").read_text())["SP_child"] == 1
    assert "child1" not in json.loads((tmp_path / "out" / "denovo_per_person.json").read_text())
    assert stats.inherited_variants == 1
    assert stats.denovo_variants == 1


def test_analyze_vcf_writes_segmented_output(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        segment_size=1,
        block_size=1,
    )

    inherited_files = glob_result_tsvs(tmp_path / "out", "inherited")
    bad_files = glob_result_tsvs(tmp_path / "out", "mendelian_bad")
    assert len(inherited_files) >= 1
    assert len(bad_files) >= 1
    assert inherited_files[0].name == "inherited_00000.tsv"
    assert (tmp_path / "out" / "checkpoint.json").is_file()
    assert (tmp_path / "out" / "stats_cumulative.json").is_file()
    assert stats.inherited_variants == 1
    checkpoint_data = json.loads((tmp_path / "out" / "checkpoint.json").read_text())
    assert checkpoint_data["details_external"] is True
    assert "inherited_per_person" not in checkpoint_data["cumulative"]
    cumulative_data = json.loads(
        (tmp_path / "out" / "stats_cumulative.json").read_text()
    )
    assert "inherited_per_person" not in cumulative_data


def test_analyze_vcf_skips_variants_in_repeat_intervals(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        segment_size=0,
        repeats_path=FIXTURES / "tiny_repeats.bed",
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=True)
    bad_records = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv", short_format=True)
    denovo_records = read_result_tsv(tmp_path / "out" / "denovo.tsv", short_format=True)

    assert len(inherited_records) == 0
    assert len(bad_records) == 1
    assert bad_records[0][1] == "4000"
    assert len(denovo_records) == 1
    assert denovo_records[0][1] == "1000"
    assert stats.inherited_variants == 0
    assert stats.mendelian_bad_variants == 1
    assert stats.denovo_variants == 1


def test_resume_continues_from_checkpoint(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    save_checkpoint(
        out,
        Checkpoint(
            chrom="22",
            last_pos=3000,
            segment_index=0,
            cumulative=CumulativeStats(
                variants_seen=3,
                alleles_tested=1,
                inherited_entries=1,
                inherited_variants=1,
                inherited_per_person={"child1": 1},
            ),
            completed=False,
        ),
    )

    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        segment_size=1000,
        block_size=1,
        resume=True,
    )

    assert stats.inherited_variants == 1
    assert stats.mendelian_bad_variants == 1
    assert stats.denovo_variants == 0
    assert stats.variants_seen == 4
    assert len(read_result_tsv(glob_result_tsvs(out, "mendelian_bad")[0])) == 1


def test_resume_rejects_completed_checkpoint(tmp_path):
    out = tmp_path / "out"
    analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        segment_size=1000,
    )
    assert load_checkpoint(out).completed is True

    with pytest.raises(ValueError, match="completed"):
        analyze_vcf(
            vcf_path=FIXTURES / "tiny.vcf",
            af_json_path=FIXTURES / "tiny_af.json",
            family_file=FIXTURES / "families.tsv",
            output_dir=out,
            segment_size=1000,
            resume=True,
        )


def test_serialize_payload_short_and_full():
    hits = {"child1": ("0/0", "0/1", "0/1", "30")}
    assert serialize_payload(hits, short_format=True) == "child1"
    assert serialize_payload(hits, short_format=False) == "child1=0/0|0/1|0/1|30"
    pair_hits = {"boy1": ("0/1", "1/1", "30")}
    assert serialize_payload(pair_hits, short_format=False) == "boy1=0/1|1/1|30"


def test_analyze_chrx_splits_sex_and_region_outputs(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny_x.vcf",
        af_json_path=FIXTURES / "tiny_x_af.json",
        family_file=FIXTURES / "families_x.tsv",
        output_dir=tmp_path / "out",
        segment_size=0,
        block_size=1,
        short_format=False,
    )

    out = tmp_path / "out"
    females = read_result_tsv(out / "inherited_females.tsv", short_format=False)
    male_par1 = read_result_tsv(out / "inherited_males_par1.tsv", short_format=False)
    male_nonpar = read_result_tsv(out / "inherited_males_nonPar.tsv", short_format=False)
    male_nonpar_denovo = read_result_tsv(
        out / "denovo_males_nonPar.tsv", short_format=False
    )
    male_par2 = read_result_tsv(out / "inherited_males_par2.tsv", short_format=False)

    assert {rec[1] for rec in females} == {"15000", "5000000", "155800000"}
    assert all("girl1" in rec[4] for rec in females)
    assert all("unknown1" not in rec[4] for rec in females)

    assert len(male_par1) == 1 and male_par1[0][1] == "15000"
    assert male_par1[0][4] == {"boy1": ("0/0", "0/1", "0/1", "30")}

    assert len(male_nonpar) == 1 and male_nonpar[0][1] == "5000000"
    assert male_nonpar[0][4] == {"boy1": ("0/1", "1/1", "30")}

    assert len(male_nonpar_denovo) == 1 and male_nonpar_denovo[0][1] == "5000100"
    assert male_nonpar_denovo[0][4] == {"boy1": ("0/0", "1/1", "30")}

    assert len(male_par2) == 1 and male_par2[0][1] == "155800000"
    assert male_par2[0][4] == {"boy1": ("0/1", "0/0", "0/1", "30")}
    assert stats.inherited_variants >= 1
    assert stats.denovo_variants >= 1
    assert (out / "inherited_per_person_females.json").is_file()
    assert (out / "denovo_per_person_males_nonPar.json").is_file()
    assert (out / "denovo_per_variant_males_nonPar.json").is_file()
    assert (out / "inherited_per_person_males_nonPar.json").is_file()
    assert (out / "stats_males_par1.json").is_file()


def test_analyze_chry_uses_only_male_nonpar_bucket(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny_y.vcf",
        af_json_path=FIXTURES / "tiny_y_af.json",
        family_file=FIXTURES / "families_x.tsv",
        output_dir=tmp_path / "out",
        segment_size=0,
        block_size=1,
        short_format=False,
    )

    out = tmp_path / "out"
    inherited = read_result_tsv(out / "inherited_males_nonPar.tsv", short_format=False)
    denovo = read_result_tsv(out / "denovo_males_nonPar.tsv", short_format=False)

    assert len(inherited) == 1 and inherited[0][1] == "100000"
    assert inherited[0][4] == {"boy1": ("1/1", "1/1", "30")}
    assert len(denovo) == 1 and denovo[0][1] == "200000"
    assert denovo[0][4] == {"boy1": ("0/0", "1/1", "30")}

    assert not (out / "inherited_females.tsv").exists()
    assert not (out / "inherited_males_par1.tsv").exists()
    assert not (out / "inherited_males_par2.tsv").exists()
    assert (out / "inherited_per_person_males_nonPar.json").is_file()
    assert (out / "denovo_per_person_males_nonPar.json").is_file()
    assert (out / "denovo_per_variant_males_nonPar.json").is_file()
    assert (out / "stats_males_nonPar.json").is_file()
    assert stats.inherited_variants == 1
    assert stats.denovo_variants == 1
    assert stats.mendelian_bad_variants == 0


def test_multiallelic_denovo_requires_parental_hom_ref(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny_multi.vcf",
        af_json_path=FIXTURES / "tiny_multi_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        multiallelic=True,
        segment_size=0,
        block_size=1,
        short_format=False,
    )

    denovo = read_result_tsv(tmp_path / "out" / "denovo.tsv", short_format=False)
    bad = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv", short_format=False)
    inherited = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=False)

    # Allele 3 with parents 0/2 is skipped; allele 2 with parents 0/0 is kept.
    assert len(denovo) == 1
    assert denovo[0][1] == "6000"
    assert denovo[0][3] == "T"
    assert denovo[0][4] == {"child1": ("0/0", "0/0", "0/2", "30")}

    bad_by_pos = {(rec[1], rec[3]): rec[4] for rec in bad}
    # Parents 0/2 x 0/2 with child 1/2 is Mendelian-inconsistent for allele 2.
    assert bad_by_pos[("5500", "G")] == {"child1": ("0/2", "0/2", "1/2", "30")}
    # Parents 0/0 x 0/0 with child 1/2 is not a clean denovo of either alt.
    assert bad_by_pos[("6500", "C")] == {"child1": ("0/0", "0/0", "1/2", "30")}
    assert bad_by_pos[("6500", "G")] == {"child1": ("0/0", "0/0", "1/2", "30")}
    assert len(bad) == 3

    assert len(inherited) == 0
    assert stats.denovo_variants == 1
    assert stats.mendelian_bad_variants == 3
    assert stats.inherited_variants == 0


def _two_tiny_shards(directory: Path) -> list[VcfShard]:
    directory.mkdir(parents=True, exist_ok=True)
    records = _tiny_record_lines()
    first = directory / "callset.chr22_1_2500.vcf"
    second = directory / "callset.chr22_2501_5000.vcf"
    _write_shard(first, records[:2])
    _write_shard(second, records[2:])
    return [
        VcfShard(path=first, chrom="chr22", start=1, end=2500),
        VcfShard(path=second, chrom="chr22", start=2501, end=5000),
    ]


def test_analyze_vcf_shards_label_output_by_coordinates(tmp_path):
    shards = _two_tiny_shards(tmp_path / "vcf")
    out = tmp_path / "out"
    stats = analyze_vcf(
        vcf_path=None,
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        vcf_shards=shards,
        segment_size=1,
        block_size=1,
    )

    inherited_files = glob_result_tsvs(out, "inherited")
    denovo_files = glob_result_tsvs(out, "denovo")
    bad_files = glob_result_tsvs(out, "mendelian_bad")
    assert [path.name for path in inherited_files] == [
        "inherited_1_2500.tsv",
        "inherited_2501_5000.tsv",
    ]
    assert [path.name for path in denovo_files] == [
        "denovo_1_2500.tsv",
        "denovo_2501_5000.tsv",
    ]
    assert [path.name for path in bad_files] == [
        "mendelian_bad_1_2500.tsv",
        "mendelian_bad_2501_5000.tsv",
    ]
    assert not (out / "inherited_00000.tsv").exists()

    assert len(read_result_tsv(out / "denovo_1_2500.tsv")) == 1
    assert read_result_tsv(out / "denovo_1_2500.tsv")[0][1] == "1000"
    assert len(read_result_tsv(out / "inherited_2501_5000.tsv")) == 1
    assert read_result_tsv(out / "inherited_2501_5000.tsv")[0][1] == "3000"
    assert len(read_result_tsv(out / "mendelian_bad_2501_5000.tsv")) == 1
    assert json.loads((out / "inherited_per_variant.json").read_text())["var_inh"] == 1
    assert json.loads((out / "denovo_per_variant.json").read_text())["var_rare"] == 1
    assert stats.inherited_variants == 1
    assert stats.denovo_variants == 1
    assert stats.mendelian_bad_variants == 1
    checkpoint = json.loads((out / "checkpoint.json").read_text())
    assert checkpoint["completed"] is True
    assert checkpoint["shard_end"] == 5000
    assert checkpoint["last_pos"] >= 5000


def test_resume_skips_completed_vcf_shards(tmp_path):
    shards = _two_tiny_shards(tmp_path / "vcf")
    out = tmp_path / "out"
    out.mkdir()
    save_checkpoint(
        out,
        Checkpoint(
            chrom="22",
            last_pos=2500,
            segment_index=0,
            cumulative=CumulativeStats(
                variants_seen=2,
                alleles_tested=1,
                denovo_entries=1,
                denovo_variants=1,
                denovo_per_person={"child1": 1},
            ),
            completed=False,
            shard_start=1,
            shard_end=2500,
        ),
    )

    stats = analyze_vcf(
        vcf_path=None,
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=out,
        vcf_shards=shards,
        block_size=1,
        resume=True,
    )

    assert stats.variants_seen == 4
    assert stats.denovo_variants == 1
    assert stats.inherited_variants == 1
    assert stats.mendelian_bad_variants == 1
    assert (out / "inherited_2501_5000.tsv").is_file()
    assert not (out / "inherited_1_2500.tsv").exists()
    assert len(read_result_tsv(out / "inherited_2501_5000.tsv")) == 1
