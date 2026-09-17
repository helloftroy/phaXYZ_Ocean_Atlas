import csv

import pytest

from phaatlas.pipeline import cluster_ecology as ce


def test_haversine_km_known_distance():
    # 1 degree of latitude along a meridian is ~111.19 km
    d = ce.haversine_km(0.0, 0.0, 1.0, 0.0)
    assert 110.0 < d < 112.5


def test_spherical_dispersion_identical_points_gives_r_equal_1():
    geo = ce.spherical_dispersion([(10.0, 20.0)] * 5)
    assert geo.n == 5
    assert geo.mean_resultant_length == pytest.approx(1.0)
    assert geo.max_pairwise_km == pytest.approx(0.0)


def test_spherical_dispersion_antipodal_points_gives_r_near_0():
    geo = ce.spherical_dispersion([(0.0, 0.0), (0.0, 180.0)])
    assert geo.mean_resultant_length < 0.01


def test_spherical_dispersion_tight_cluster_gives_high_r():
    # small variation around one point -- should be much closer to 1 than spread-out points
    tight = ce.spherical_dispersion([(10.0, 20.0), (10.1, 20.1), (9.9, 19.9), (10.05, 19.95)])
    spread = ce.spherical_dispersion([(10.0, 20.0), (-40.0, 100.0), (60.0, -150.0), (0.0, 0.0)])
    assert tight.mean_resultant_length > 0.999
    assert spread.mean_resultant_length < tight.mean_resultant_length


def test_spherical_dispersion_empty_returns_none_fields():
    geo = ce.spherical_dispersion([])
    assert geo.n == 0
    assert geo.mean_resultant_length is None
    assert geo.max_pairwise_km is None


def _write_depth_metadata(path, rows):
    """rows: list of dicts with target_id, genome, latitude_degN,
    longitude_degE, depth_m, gtdb_genus, gtdb_phylum, study_id."""
    fieldnames = ["target_id", "genome", "latitude_degN", "longitude_degE", "depth_m", "gtdb_genus", "gtdb_phylum", "study_id"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow({fn: r.get(fn, "") for fn in fieldnames})


def test_summarize_cluster_ecology_aggregates_correctly(tmp_path):
    path = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    _write_depth_metadata(path, [
        {"target_id": "T1", "genome": "G1", "latitude_degN": "10.0", "longitude_degE": "20.0",
         "depth_m": "10", "gtdb_genus": "Moritella", "gtdb_phylum": "Pseudomonadota", "study_id": "S1"},
        {"target_id": "T2", "genome": "G2", "latitude_degN": "10.1", "longitude_degE": "20.1",
         "depth_m": "30", "gtdb_genus": "Moritella", "gtdb_phylum": "Pseudomonadota", "study_id": "S1"},
        # same genome as G1 appearing again under a different target_id in the SAME cluster -- must not double count
        {"target_id": "T3", "genome": "G1", "latitude_degN": "10.0", "longitude_degE": "20.0",
         "depth_m": "10", "gtdb_genus": "Moritella", "gtdb_phylum": "Pseudomonadota", "study_id": "S1"},
        # a different cluster entirely
        {"target_id": "T4", "genome": "G3", "latitude_degN": "-40.0", "longitude_degE": "100.0",
         "depth_m": "", "gtdb_genus": "Planktomarina", "gtdb_phylum": "Pseudomonadota", "study_id": "S2"},
    ])

    cluster_assignments = {"T1": "REP1", "T2": "REP1", "T3": "REP1", "T4": "REP2"}
    stats = ce.summarize_cluster_ecology(path, cluster_assignments)
    by_cluster = {s.cluster_id: s for s in stats}

    rep1 = by_cluster["REP1"]
    assert rep1.n_target_ids == 3   # T1, T2, T3
    assert rep1.n_genomes == 2       # G1, G2 (G1 not double counted)
    assert rep1.geo.n == 2
    assert rep1.median_depth_m == 20.0  # median of [10, 30]
    assert rep1.n_genomes_with_depth == 2
    assert rep1.top_genera == [("Moritella", 2)]

    rep2 = by_cluster["REP2"]
    assert rep2.n_target_ids == 1
    assert rep2.n_genomes == 1
    assert rep2.median_depth_m is None  # no depth reported for G3
    assert rep2.n_genomes_with_depth == 0
    assert rep2.top_genera == [("Planktomarina", 1)]

    # sorted most genomes first
    assert stats[0].cluster_id == "REP1"


def test_summarize_cluster_ecology_skips_targets_with_no_cluster_assignment(tmp_path):
    path = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    _write_depth_metadata(path, [
        {"target_id": "T1", "genome": "G1", "latitude_degN": "10.0", "longitude_degE": "20.0"},
    ])
    stats = ce.summarize_cluster_ecology(path, cluster_assignments={})  # T1 not clustered
    assert stats == []


def test_write_cluster_ecology_round_trip(tmp_path):
    path = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    _write_depth_metadata(path, [
        {"target_id": "T1", "genome": "G1", "latitude_degN": "10.0", "longitude_degE": "20.0",
         "depth_m": "500", "gtdb_genus": "Moritella", "gtdb_phylum": "Pseudomonadota", "study_id": "S1"},
    ])
    stats = ce.summarize_cluster_ecology(path, cluster_assignments={"T1": "REP1"})
    out_path = tmp_path / "cluster_ecology.tsv"
    ce.write_cluster_ecology(stats, out_path)

    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows[0]["cluster_id"] == "REP1"
    assert rows[0]["median_depth_m"] == "500.0000"
    assert "Moritella" in rows[0]["top_genera"]
