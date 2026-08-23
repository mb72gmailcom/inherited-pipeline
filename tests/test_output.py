import json

from inherited.checkpoint import load_checkpoint
from inherited.output import BlockWriter, ResultWriter, serialize_patient_ids


def test_block_writer_flushes_in_blocks(tmp_path):
    path = tmp_path / "out.tsv"
    writer = BlockWriter(path, block_size=2)
    hits = {"p1": ("0/1", "0/0", "0/1", "30")}

    writer.append("22", "100", "A", "G", serialize_patient_ids(hits))
    writer.append("22", "200", "C", "T", serialize_patient_ids(hits))
    assert writer.lines_written == 2
    assert writer.block == []

    writer.append("22", "300", "G", "A", serialize_patient_ids(hits))
    assert writer.lines_written == 2
    assert len(writer.block) == 1

    writer.close()
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0].startswith("#CHROM")
    assert len(lines) == 4


def test_result_writer_streams_and_merges_per_variant_counts(tmp_path):
    writer = ResultWriter(tmp_path, block_size=1, segment_size=1)
    writer.write_inherited(
        "22", "100", "A", "G", "variant_a", {"p1": ("0/0", "0/1", "0/1", "30")}
    )
    writer.write_inherited(
        "22",
        "200",
        "C",
        "T",
        "variant_b",
        {
            "p1": ("0/0", "0/1", "0/1", "30"),
            "p2": ("0/1", "0/0", "0/1", "30"),
        },
    )
    writer.write_denovo(
        "22", "300", "A", "T", "variant_c", {"p1": ("0/0", "0/0", "0/1", "30")}
    )
    writer.close()
    writer.finalize()

    assert json.loads((tmp_path / "inherited_per_variant.json").read_text()) == {
        "variant_a": 1,
        "variant_b": 2,
    }
    assert json.loads((tmp_path / "denovo_per_variant.json").read_text()) == {
        "variant_c": 1
    }
    assert json.loads((tmp_path / "denovo_per_person.json").read_text()) == {"p1": 1}


def test_result_writer_streams_unsegmented_per_variant_counts(tmp_path):
    writer = ResultWriter(tmp_path, block_size=1, segment_size=0)
    writer.write_inherited(
        "22", "100", "A", "G", "variant_a", {"p1": ("0/0", "0/1", "0/1", "30")}
    )
    writer.close()
    writer.finalize()

    assert json.loads((tmp_path / "inherited_per_variant.json").read_text()) == {
        "variant_a": 1
    }


def test_result_writer_restores_details_from_compact_checkpoint(tmp_path):
    writer = ResultWriter(tmp_path, block_size=1, segment_size=1)
    writer.write_inherited(
        "22", "100", "A", "G", "variant_a", {"p1": ("0/0", "0/1", "0/1", "30")}
    )
    writer.close()

    checkpoint = load_checkpoint(tmp_path)
    assert checkpoint is not None
    assert checkpoint.details_external is True

    resumed = ResultWriter.from_checkpoint(
        tmp_path, checkpoint, block_size=1, segment_size=1
    )
    assert resumed.cumulative.inherited_per_person == {"p1": 1}
    resumed.close()


def test_result_writer_labels_shard_files_by_coordinates(tmp_path):
    writer = ResultWriter(tmp_path, block_size=1, segment_size=1, shard_mode=True)
    writer.begin_shard(1, 2500)
    writer.write_denovo(
        "22", "1000", "A", "G", "variant_a", {"p1": ("0/0", "0/0", "0/1", "30")}
    )
    writer.begin_shard(2501, 5000)
    writer.write_inherited(
        "22", "3000", "A", "G", "variant_b", {"p1": ("0/0", "0/1", "0/1", "30")}
    )
    writer.close()
    writer.finalize()

    assert (tmp_path / "denovo_1_2500.tsv").is_file()
    assert (tmp_path / "inherited_2501_5000.tsv").is_file()
    assert not (tmp_path / "inherited_00000.tsv").exists()
    assert json.loads((tmp_path / "inherited_per_variant.json").read_text()) == {
        "variant_b": 1
    }
    assert json.loads((tmp_path / "denovo_per_variant.json").read_text()) == {
        "variant_a": 1
    }
    checkpoint = load_checkpoint(tmp_path)
    assert checkpoint is not None
    assert checkpoint.shard_end == 5000


def test_result_writer_suffixes_all_output_files(tmp_path):
    writer = ResultWriter(tmp_path, block_size=1, segment_size=1, file_suffix="1_2500")
    writer.write_inherited(
        "22", "100", "A", "G", "variant_a", {"p1": ("0/0", "0/1", "0/1", "30")}
    )
    writer.write_denovo(
        "22", "300", "A", "T", "variant_c", {"p1": ("0/0", "0/0", "0/1", "30")}
    )
    writer.close()
    writer.finalize()

    suffix = "1_2500"
    assert (tmp_path / f"inherited_{suffix}.tsv").is_file()
    assert (tmp_path / f"denovo_{suffix}.tsv").is_file()
    assert (tmp_path / f"mendelian_bad_{suffix}.tsv").is_file()
    assert not (tmp_path / "inherited.tsv").exists()
    assert not (tmp_path / "inherited_00000.tsv").exists()
    assert json.loads((tmp_path / f"inherited_per_variant_{suffix}.json").read_text()) == {
        "variant_a": 1
    }
    assert json.loads((tmp_path / f"denovo_per_variant_{suffix}.json").read_text()) == {
        "variant_c": 1
    }
    assert json.loads((tmp_path / f"denovo_per_person_{suffix}.json").read_text()) == {"p1": 1}
    assert json.loads((tmp_path / f"stats_{suffix}.json").read_text())["inherited_variants"] == 1
    assert (tmp_path / f"stats_cumulative_{suffix}.json").is_file()
    assert (tmp_path / f"checkpoint_{suffix}.json").is_file()
    assert (tmp_path / f"cumulative_detail_deltas_{suffix}.jsonl").is_file()
    assert not (tmp_path / "stats.json").exists()
    assert not (tmp_path / "checkpoint.json").exists()
    checkpoint = load_checkpoint(tmp_path, suffix)
    assert checkpoint is not None
    assert checkpoint.completed is True
