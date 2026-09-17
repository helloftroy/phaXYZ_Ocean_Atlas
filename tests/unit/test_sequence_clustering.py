import csv

from phaatlas.pipeline import sequence_clustering as sc


def test_collect_wanted_target_ids_across_multiple_files(tmp_path):
    p1 = tmp_path / "phaC_unique_targets_with_metadata.tsv"
    with open(p1, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["target_id", "genome"])
        writer.writerow(["T1", "G1"])
        writer.writerow(["T2", "G2"])
        writer.writerow(["T1", "G3"])  # same target_id, different genome row -- must dedupe

    p2 = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    with open(p2, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["target_id", "genome"])
        writer.writerow(["T3", "G4"])

    ids = sc.collect_wanted_target_ids([p1, p2])
    assert ids == {"T1", "T2", "T3"}


def test_write_wanted_ids_file_is_sorted_one_per_line(tmp_path):
    out_path = tmp_path / "wanted_ids.txt"
    n = sc.write_wanted_ids_file({"T3", "T1", "T2"}, out_path)
    assert n == 3
    assert out_path.read_text() == "T1\nT2\nT3\n"


def test_load_cluster_assignments_standard_mmseqs_format(tmp_path):
    """mmseqs cluster createtsv format: representative \\t member, one row
    per member -- including the representative's own row (representative
    is a member of its own cluster)."""
    cluster_tsv = tmp_path / "phaC_cluster.tsv"
    cluster_tsv.write_text(
        "REP1\tREP1\n"
        "REP1\tMEM_A\n"
        "REP1\tMEM_B\n"
        "REP2\tREP2\n"  # singleton cluster
    )

    assignments = sc.load_cluster_assignments(cluster_tsv)
    assert assignments == {
        "REP1": "REP1",
        "MEM_A": "REP1",
        "MEM_B": "REP1",
        "REP2": "REP2",
    }
