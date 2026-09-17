import csv

import numpy as np

from phaatlas.pipeline import sequence_embedding as se


def test_read_fasta_parses_multiline_records(tmp_path):
    path = tmp_path / "seqs.faa"
    path.write_text(">T1 some description\nMKV\nLNR\n>T2\nAAAA\n")
    records = se.read_fasta(path)
    assert records == [("T1", "MKVLNR"), ("T2", "AAAA")]


def test_kmer_composition_vector_normalizes_and_skips_unknown_chars():
    index = se.kmer_index(k=1, alphabet="AC")
    # "AACX" -- 4 windows of size 1: A, A, C, X. X is outside the alphabet -> skipped.
    vec = se.kmer_composition_vector("AACX", k=1, index=index)
    assert vec[index["A"]] == 2 / 3
    assert vec[index["C"]] == 1 / 3


def test_kmer_composition_vector_empty_sequence_is_all_zero():
    index = se.kmer_index(k=2, alphabet="AC")
    vec = se.kmer_composition_vector("", k=2, index=index)
    assert vec.sum() == 0


def test_build_feature_matrix_shape_matches_ids(tmp_path):
    path = tmp_path / "seqs.faa"
    path.write_text(">T1\nMKVLNR\n>T2\nAAAAAA\n")
    ids, matrix = se.build_feature_matrix(path, k=2)
    assert ids == ["T1", "T2"]
    assert matrix.shape[0] == 2
    assert matrix.shape[1] == len(se.kmer_index(2))


def test_pca_2d_handles_small_and_degenerate_inputs():
    assert se.pca_2d(np.zeros((0, 5))).shape == (0, 2)
    assert se.pca_2d(np.zeros((1, 5))).shape == (1, 2)
    coords = se.pca_2d(np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]))
    assert coords.shape == (3, 2)


def test_umap_2d_returns_none_when_not_installed_or_too_few_points():
    # umap-learn is an optional extra -- not installed in the base test env,
    # so this exercises the graceful-fallback path either way.
    result = se.umap_2d(np.random.rand(2, 5))
    assert result is None


def _write_depth_metadata(path, rows):
    fieldnames = ["target_id", "genome", "gtdb_genus", "gtdb_phylum", "study_id", "depth_zone"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow({fn: r.get(fn, "") for fn in fieldnames})


def test_collect_target_annotations_uses_first_row_per_target(tmp_path):
    path = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    _write_depth_metadata(path, [
        {"target_id": "T1", "genome": "G1", "gtdb_genus": "Moritella", "depth_zone": ">4000m"},
        {"target_id": "T1", "genome": "G2", "gtdb_genus": "OtherGenus", "depth_zone": "0-50m"},  # ignored, T1 already seen
        {"target_id": "T2", "genome": "G3", "gtdb_genus": "Planktomarina", "depth_zone": "0-50m"},
    ])
    cluster_assignments = {"T1": "REP1", "T2": "REP2"}
    genome_architecture = {"G1": "ABC", "G3": "AB"}

    annotations = se.collect_target_annotations(path, cluster_assignments, genome_architecture)
    assert annotations["T1"]["genus"] == "Moritella"
    assert annotations["T1"]["cluster_id"] == "REP1"
    assert annotations["T1"]["architecture"] == "ABC"
    assert annotations["T2"]["genus"] == "Planktomarina"


def test_write_embedding_tsv(tmp_path):
    ids = ["T1", "T2"]
    pca_coords = np.array([[1.0, 2.0], [3.0, 4.0]])
    annotations = {"T1": {"cluster_id": "REP1", "genus": "Moritella"}}
    out_path = tmp_path / "embedding.tsv"
    n = se.write_embedding_tsv(ids, pca_coords, None, annotations, out_path)
    assert n == 2
    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows[0]["target_id"] == "T1"
    assert rows[0]["pca_x"] == "1.000000"
    assert rows[0]["cluster_id"] == "REP1"
    assert rows[1]["cluster_id"] == ""  # T2 has no annotation entry -- blank, not an error
