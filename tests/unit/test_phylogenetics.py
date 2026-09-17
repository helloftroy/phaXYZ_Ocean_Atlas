import csv
import shutil

import pytest

from phaatlas.pipeline import phylogenetics as phylo

MAFFT_AVAILABLE = shutil.which("mafft") is not None
FASTTREE_AVAILABLE = shutil.which("FastTree") is not None


def test_select_representatives_ranks_by_genome_count_and_appends_extras(tmp_path):
    path = tmp_path / "ecology.tsv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["cluster_id", "n_genomes"])
        writer.writerow(["C_SMALL", "5"])
        writer.writerow(["C_BIG", "100"])
        writer.writerow(["C_MED", "20"])

    result = phylo.select_representatives(path, top_n=2, extra_cluster_ids=["C_SMALL", "C_MED"])
    assert result == ["C_BIG", "C_MED", "C_SMALL"]  # top 2 by count, then C_SMALL appended (C_MED already present)


def test_select_representatives_no_extras(tmp_path):
    path = tmp_path / "ecology.tsv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["cluster_id", "n_genomes"])
        writer.writerow(["A", "1"])
        writer.writerow(["B", "2"])

    result = phylo.select_representatives(path, top_n=1)
    assert result == ["B"]


def test_extract_sequences_subset_filters_and_preserves_order(tmp_path):
    fasta_path = tmp_path / "seqs.faa"
    fasta_path.write_text(">T1\nMKVL\n>T2\nAAAA\n>T3\nGGGG\n")
    out_path = tmp_path / "subset.faa"

    n = phylo.extract_sequences_subset(fasta_path, ["T3", "T1", "T_MISSING"], out_path)
    assert n == 2  # T_MISSING silently skipped, not an error
    text = out_path.read_text()
    assert text.startswith(">T3")  # wanted-list order preserved
    assert ">T1" in text
    assert ">T2" not in text
    assert ">T_MISSING" not in text


def _write_annotation_metadata(path, rows):
    fieldnames = ["target_id", "genome", "gtdb_genus", "gtdb_phylum", "depth_zone", "depth_m",
                  "study_id", "latitude_degN", "longitude_degE"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow({fn: r.get(fn, "") for fn in fieldnames})


def test_annotate_clusters_takes_modal_values_per_cluster(tmp_path):
    path = tmp_path / "metadata_depth.tsv"
    _write_annotation_metadata(path, [
        # cluster REP1: 2 genomes, both Moritella/Pseudomonadota, one deep one unknown depth
        {"target_id": "T1", "genome": "G1", "gtdb_genus": "Moritella", "gtdb_phylum": "Pseudomonadota",
         "depth_zone": ">4000m", "depth_m": "4500", "study_id": "S1", "latitude_degN": "-70.0", "longitude_degE": "-100.0"},
        {"target_id": "T2", "genome": "G2", "gtdb_genus": "Moritella", "gtdb_phylum": "Pseudomonadota",
         "study_id": "S1", "latitude_degN": "-72.0", "longitude_degE": "-105.0"},
        # cluster REP2: not in wanted set -- must be ignored entirely
        {"target_id": "T3", "genome": "G3", "gtdb_genus": "Vibrio", "gtdb_phylum": "Pseudomonadota"},
    ])
    cluster_assignments = {"T1": "REP1", "T2": "REP1", "T3": "REP2"}
    genome_architecture = {"G1": "ABC", "G2": "ABC"}

    result = phylo.annotate_clusters(path, cluster_assignments, genome_architecture, wanted_cluster_ids=["REP1"])

    assert set(result) == {"REP1"}
    a = result["REP1"]
    assert a.n_genomes == 2
    assert a.dominant_genus == "Moritella"
    assert a.dominant_phylum == "Pseudomonadota"
    assert a.median_depth_m == 4500.0  # only one genome had a resolved depth
    assert a.dominant_architecture == "ABC"
    assert a.dominant_ocean_basin == "Southern Ocean"  # both genomes' lat <= -60


def test_annotate_clusters_handles_missing_data_gracefully(tmp_path):
    path = tmp_path / "metadata_depth.tsv"
    _write_annotation_metadata(path, [
        {"target_id": "T1", "genome": "G1"},  # no genus/phylum/depth/basin/study at all
    ])
    cluster_assignments = {"T1": "REP1"}
    result = phylo.annotate_clusters(path, cluster_assignments, {}, wanted_cluster_ids=["REP1"])
    a = result["REP1"]
    assert a.n_genomes == 1
    assert a.dominant_genus == ""
    assert a.median_depth_m is None
    assert a.dominant_architecture == ""


def test_write_cluster_annotations_round_trip(tmp_path):
    ann = {"REP1": phylo.ClusterAnnotation(
        cluster_id="REP1", n_genomes=5, dominant_genus="Moritella", dominant_phylum="Pseudomonadota",
        median_depth_m=4500.0, dominant_depth_zone=">4000m", dominant_ocean_basin="Southern Ocean",
        dominant_architecture="ABC", top_study="S1",
    )}
    out_path = tmp_path / "annotations.tsv"
    phylo.write_cluster_annotations(ann, out_path)
    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows[0]["cluster_id"] == "REP1"
    assert rows[0]["median_depth_m"] == "4500.00"
    assert rows[0]["dominant_genus"] == "Moritella"


def test_compute_tree_layout_matches_hand_computed_positions(tmp_path):
    # ((A:1,B:2)0.9:0.5,C:3);
    # tip order left-to-right: A, B, C -> y = 0, 1, 2
    # internal node x = 0.5, y = mean(0,1) = 0.5
    # A x = 0.5+1=1.5, B x = 0.5+2=2.5, C x = 0+3=3.0
    # root x = 0, y = mean(internal.y=0.5, C.y=2) = 1.25
    newick_path = tmp_path / "test.nwk"
    newick_path.write_text("((A:1,B:2)0.9:0.5,C:3);\n")

    nodes = phylo.compute_tree_layout(newick_path)
    by_name = {n.name: n for n in nodes if n.is_tip}

    assert by_name["A"].x == pytest.approx(1.5)
    assert by_name["A"].y == pytest.approx(0.0)
    assert by_name["B"].x == pytest.approx(2.5)
    assert by_name["B"].y == pytest.approx(1.0)
    assert by_name["C"].x == pytest.approx(3.0)
    assert by_name["C"].y == pytest.approx(2.0)

    internal = [n for n in nodes if not n.is_tip and n.parent_id is not None][0]
    assert internal.x == pytest.approx(0.5)
    assert internal.y == pytest.approx(0.5)
    assert internal.support == pytest.approx(0.9)

    root = [n for n in nodes if n.parent_id is None][0]
    assert root.x == pytest.approx(0.0)
    assert root.y == pytest.approx(1.25)

    # every non-root node's parent_id must point at a real, earlier-or-equal node in the flat list
    ids = {n.id for n in nodes}
    for n in nodes:
        if n.parent_id is not None:
            assert n.parent_id in ids


def test_write_tree_layout_round_trip(tmp_path):
    nodes = [
        phylo.TreeLayoutNode(id=0, name="", x=0.0, y=1.0, parent_id=None, support=None, is_tip=False),
        phylo.TreeLayoutNode(id=1, name="A", x=1.5, y=0.0, parent_id=0, support=0.9, is_tip=True),
    ]
    out_path = tmp_path / "layout.tsv"
    phylo.write_tree_layout(nodes, out_path)
    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows[1]["name"] == "A"
    assert rows[1]["parent_id"] == "0"
    assert rows[0]["parent_id"] == ""


@pytest.mark.skipif(not MAFFT_AVAILABLE, reason="mafft binary not on PATH")
def test_run_mafft_produces_aligned_fasta_same_length_sequences(tmp_path):
    input_fasta = tmp_path / "in.faa"
    input_fasta.write_text(
        ">A\nMKVLNRQAVASLKELQASAAAINSNPFAAAKPAEIQGLARFVQAAKADPAGAF\n"
        ">B\nMKVLNRQAVASLKELQASAAAINSNPFAAAKPAEIQGLARFVQAAKADPAGA\n"
        ">C\nMKVLNRQAVASLKELQASAAAINSNPFAAAKPAEIQGLARFVQAAKADPAGAFA\n"
    )
    output_fasta = tmp_path / "out.faa"
    phylo.run_mafft(input_fasta, output_fasta)
    records = phylo.read_fasta(output_fasta)
    assert len(records) == 3
    lengths = {len(seq) for _, seq in records}
    assert len(lengths) == 1  # aligned -- every sequence padded to the same length with gaps


@pytest.mark.skipif(not (MAFFT_AVAILABLE and FASTTREE_AVAILABLE), reason="mafft/FastTree not both on PATH")
def test_run_fasttree_produces_valid_newick(tmp_path):
    input_fasta = tmp_path / "in.faa"
    input_fasta.write_text(
        ">A\nMKVLNRQAVASLKELQASAAAINSNPFAAAKPAEIQGLARFVQAAKADPAGAF\n"
        ">B\nMKVLNRQAVASLKELQASAAAINSNPFAAAKPAEIQGLARFVQAAKADPAGA\n"
        ">C\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
    )
    aligned = tmp_path / "aligned.faa"
    phylo.run_mafft(input_fasta, aligned)
    tree_path = tmp_path / "tree.nwk"
    phylo.run_fasttree(aligned, tree_path)
    newick = tree_path.read_text()
    assert newick.strip().endswith(";")
    for tip in ["A", "B", "C"]:
        assert tip in newick
