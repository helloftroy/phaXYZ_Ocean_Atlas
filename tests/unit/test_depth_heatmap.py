import csv

from phaatlas.pipeline import depth_heatmap as dh


def _write_depth_metadata(path, rows):
    """rows: list of (target_id, genome, best_query, depth_zone)."""
    fieldnames = ["target_id", "genome", "best_query", "depth_zone"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for target_id, genome, best_query, zone in rows:
            writer.writerow({"target_id": target_id, "genome": genome, "best_query": best_query, "depth_zone": zone})


def test_build_clade_depth_matrix_counts_distinct_genomes_and_targets(tmp_path):
    path = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    _write_depth_metadata(path, [
        ("T1", "G1", "Q_deep", "1000-4000m"),
        ("T2", "G1", "Q_deep", "1000-4000m"),  # same genome, second target -> n_distinct_target_ids=2, n_distinct_genomes=1
        ("T3", "G2", "Q_deep", "1000-4000m"),
        ("T4", "G3", "Q_shallow", "0-50m"),
        ("T5", "G4", "Q_shallow", "0-50m"),
        ("T6", "G5", "Q_shallow", ""),  # missing depth -> unknown bucket, not dropped
    ])

    clades_ranked, genome_counts, target_id_counts = dh.build_clade_depth_matrix(path)

    assert clades_ranked[0] == "Q_shallow"  # 3 total genomes vs Q_deep's 2 -> ranked first
    assert genome_counts["Q_deep"]["1000-4000m"] == 2
    assert target_id_counts["Q_deep"]["1000-4000m"] == 3
    assert genome_counts["Q_shallow"]["0-50m"] == 2
    assert genome_counts["Q_shallow"][dh.UNKNOWN_DEPTH_LABEL] == 1


def test_write_long_format_omits_empty_cells_and_orders_zones(tmp_path):
    path = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    _write_depth_metadata(path, [
        ("T1", "G1", "Q1", ">4000m"),
        ("T2", "G2", "Q1", "0-50m"),
    ])
    clades_ranked, genome_counts, target_id_counts = dh.build_clade_depth_matrix(path)

    out_path = tmp_path / "long.tsv"
    n = dh.write_long_format(clades_ranked, genome_counts, target_id_counts, out_path)
    assert n == 2
    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    zones_in_order = [r["depth_zone"] for r in rows]
    # "0-50m" must appear before ">4000m" despite reverse alphabetical order
    assert zones_in_order.index("0-50m") < zones_in_order.index(">4000m")


def test_write_wide_matrix_includes_all_zone_columns_even_when_zero(tmp_path):
    path = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    _write_depth_metadata(path, [("T1", "G1", "Q1", "0-50m")])
    clades_ranked, genome_counts, _ = dh.build_clade_depth_matrix(path)

    out_path = tmp_path / "wide.tsv"
    dh.write_wide_matrix(clades_ranked, genome_counts, out_path)
    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows[0]["clade"] == "Q1"
    assert rows[0]["0-50m"] == "1"
    assert rows[0][">4000m"] == "0"
    assert rows[0][dh.UNKNOWN_DEPTH_LABEL] == "0"
