import csv

import numpy as np
import pytest
from scipy.spatial.distance import pdist, squareform

from phaatlas.pipeline import community_ordination as co


def test_classify_ocean_basin_known_points():
    assert co.classify_ocean_basin(-70.0, 0.0) == "Southern Ocean"
    assert co.classify_ocean_basin(75.0, 0.0) == "Arctic Ocean"
    assert co.classify_ocean_basin(38.0, 15.0) == "Mediterranean/Black Sea"
    assert co.classify_ocean_basin(0.0, -30.0) == "Atlantic Ocean"
    assert co.classify_ocean_basin(0.0, 70.0) == "Indian Ocean"
    assert co.classify_ocean_basin(0.0, 150.0) == "Pacific Ocean"
    assert co.classify_ocean_basin(0.0, -150.0) == "Pacific Ocean"


def _write_site_metadata(path, rows):
    fieldnames = ["target_id", "sample_id", "latitude_degN", "longitude_degE",
                  "depth_zone", "ecosystem_type", "ecosystem_name", "study_id"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow({fn: r.get(fn, "") for fn in fieldnames})


def test_build_site_cluster_matrix_filters_richness_and_missing_latlon(tmp_path):
    path = tmp_path / "metadata_depth.tsv"
    _write_site_metadata(path, [
        # SITE_A: 3 distinct clusters -- meets min_richness=3
        {"target_id": "T1", "sample_id": "SITE_A", "latitude_degN": "10.0", "longitude_degE": "70.0", "study_id": "S1"},
        {"target_id": "T2", "sample_id": "SITE_A", "latitude_degN": "10.0", "longitude_degE": "70.0", "study_id": "S1"},
        {"target_id": "T3", "sample_id": "SITE_A", "latitude_degN": "10.0", "longitude_degE": "70.0", "study_id": "S1"},
        # SITE_A again, second hit of the SAME cluster as T1 -- must sum to count 2, not overwrite
        {"target_id": "T1b", "sample_id": "SITE_A", "latitude_degN": "10.0", "longitude_degE": "70.0", "study_id": "S1"},
        # SITE_B: only 1 distinct cluster -- below min_richness=3, must be dropped
        {"target_id": "T4", "sample_id": "SITE_B", "latitude_degN": "-40.0", "longitude_degE": "100.0", "study_id": "S2"},
        # SITE_C: 3 distinct clusters, but no lat/lon -- must be dropped
        {"target_id": "T1", "sample_id": "SITE_C", "latitude_degN": "", "longitude_degE": "", "study_id": "S3"},
        {"target_id": "T2", "sample_id": "SITE_C", "latitude_degN": "", "longitude_degE": "", "study_id": "S3"},
        {"target_id": "T3", "sample_id": "SITE_C", "latitude_degN": "", "longitude_degE": "", "study_id": "S3"},
    ])
    cluster_assignments = {"T1": "CLUST_1", "T1b": "CLUST_1", "T2": "CLUST_2", "T3": "CLUST_3", "T4": "CLUST_1"}

    result = co.build_site_cluster_matrix(path, cluster_assignments, min_richness=3, max_sites=None)

    assert result.site_ids == ["SITE_A"]  # SITE_B (too sparse) and SITE_C (no lat/lon) both excluded
    assert set(result.cluster_ids) == {"CLUST_1", "CLUST_2", "CLUST_3"}
    row = result.matrix[result.site_ids.index("SITE_A")]
    c1_idx = result.cluster_ids.index("CLUST_1")
    assert row[c1_idx] == 2  # T1 and T1b both mapped to CLUST_1 at the same site
    assert result.site_meta["SITE_A"]["lat"] == 10.0
    assert result.site_meta["SITE_A"]["ocean_basin"] == "Indian Ocean"


def test_build_site_cluster_matrix_caps_at_max_sites_keeping_richest(tmp_path):
    path = tmp_path / "metadata_depth.tsv"
    rows = []
    # SITE_RICH: 5 distinct clusters
    for i in range(5):
        rows.append({"target_id": f"R{i}", "sample_id": "SITE_RICH", "latitude_degN": "1.0", "longitude_degE": "1.0"})
    # SITE_POOR: 3 distinct clusters (still meets min_richness=3, but is the poorer of the two)
    for i in range(3):
        rows.append({"target_id": f"P{i}", "sample_id": "SITE_POOR", "latitude_degN": "2.0", "longitude_degE": "2.0"})
    _write_site_metadata(path, rows)
    cluster_assignments = {f"R{i}": f"CR{i}" for i in range(5)} | {f"P{i}": f"CP{i}" for i in range(3)}

    result = co.build_site_cluster_matrix(path, cluster_assignments, min_richness=3, max_sites=1)
    assert result.site_ids == ["SITE_RICH"]  # only the richer site kept


def test_bray_curtis_distance_matrix_identical_rows_zero_disjoint_rows_one():
    matrix = np.array([
        [1.0, 2.0, 0.0],
        [1.0, 2.0, 0.0],  # identical to row 0
        [0.0, 0.0, 5.0],  # disjoint from row 0
    ])
    d = co.bray_curtis_distance_matrix(matrix)
    assert d[0, 1] == pytest.approx(0.0, abs=1e-9)
    assert d[0, 2] == pytest.approx(1.0, abs=1e-9)


def test_jaccard_distance_matrix_binarizes_presence_absence():
    matrix = np.array([
        [5.0, 0.0, 1.0],
        [1.0, 0.0, 99.0],  # same presence/absence pattern as row 0 despite different counts
        [0.0, 3.0, 0.0],   # completely disjoint presence pattern
    ])
    d = co.jaccard_distance_matrix(matrix)
    assert d[0, 1] == pytest.approx(0.0, abs=1e-9)
    assert d[0, 2] == pytest.approx(1.0, abs=1e-9)


def test_geographic_distance_matrix_matches_haversine_and_is_symmetric():
    site_ids = ["A", "B", "C"]
    site_meta = {
        "A": {"lat": 0.0, "lon": 0.0},
        "B": {"lat": 1.0, "lon": 0.0},
        "C": {"lat": 0.0, "lon": 1.0},
    }
    d = co.geographic_distance_matrix(site_ids, site_meta)
    assert d[0, 0] == 0.0
    assert d[0, 1] == pytest.approx(co.haversine_km(0.0, 0.0, 1.0, 0.0))
    assert d[0, 1] == d[1, 0]


def test_pcoa_recovers_known_euclidean_distances():
    points = np.array([[0.0, 0.0], [3.0, 0.0], [0.0, 4.0], [3.0, 4.0]])
    d = squareform(pdist(points))
    coords, pct_variance = co.pcoa(d, n_components=2)
    recovered_d = squareform(pdist(coords))
    assert np.allclose(d, recovered_d, atol=1e-6)
    assert pct_variance.sum() == pytest.approx(1.0, abs=1e-6)


def test_nmds_returns_correct_shape_and_positive_stress():
    rng = np.random.default_rng(0)
    points = rng.random((15, 3))
    d = squareform(pdist(points))
    coords, stress = co.nmds(d, n_components=2, random_state=0, n_init=1, max_iter=100)
    assert coords.shape == (15, 2)
    assert stress >= 0


def test_mantel_test_detects_real_correlation():
    rng = np.random.default_rng(0)
    n = 30
    base = rng.random((n, n))
    dist_a = (base + base.T) / 2
    np.fill_diagonal(dist_a, 0)
    noise = rng.normal(0, 0.01, (n, n))
    noise = (noise + noise.T) / 2
    dist_b = dist_a + noise
    np.fill_diagonal(dist_b, 0)

    r, p = co.mantel_test(dist_a, dist_b, n_permutations=199, random_state=1)
    assert r > 0.9
    assert p < 0.05


def test_mantel_test_no_signal_for_unrelated_matrices():
    rng = np.random.default_rng(0)
    n = 30
    a = rng.random((n, n)); dist_a = (a + a.T) / 2; np.fill_diagonal(dist_a, 0)
    b = rng.random((n, n)); dist_b = (b + b.T) / 2; np.fill_diagonal(dist_b, 0)

    r, p = co.mantel_test(dist_a, dist_b, n_permutations=199, random_state=1)
    assert abs(r) < 0.35
    assert p > 0.05


def test_mantel_test_p_value_never_exactly_zero():
    rng = np.random.default_rng(0)
    n = 10
    a = rng.random((n, n)); dist_a = (a + a.T) / 2; np.fill_diagonal(dist_a, 0)
    r, p = co.mantel_test(dist_a, dist_a.copy(), n_permutations=50, random_state=1)
    assert r == pytest.approx(1.0)
    assert p > 0  # add-one correction -- never literally zero
