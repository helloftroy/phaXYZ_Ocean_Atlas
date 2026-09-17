"""Site x gene-cluster community ecology: treats each sampling site (an
OMDB SAMPLE -- the same atomic geographic unit pipeline/ncbi_depth.py
enriches with depth) as an ecological "sample" and each protein sequence
cluster (see pipeline/sequence_clustering.py) as a "species", then asks
the same questions a community ecologist would ask of any
species-abundance matrix:

1. Do sites separate in ordination space (NMDS/PCoA over Bray-Curtis or
   Jaccard dissimilarity) by ocean basin, depth, or biome?
2. Is there a measurable isolation-by-distance signal -- do
   compositionally similar sites tend to be geographically close (Mantel
   test, community dissimilarity vs. geographic distance) -- or does
   gene-cluster composition look more like well-mixed, cosmopolitan
   dispersal with no geographic signal at all?

Two honest approximations, both surfaced in the summary this module
returns rather than left implicit:
- classify_ocean_basin is a simple longitude/latitude bounding-box
  heuristic, not real bathymetric/coastline-aware basin boundaries
  (which follow coastlines, not straight lines) -- good enough to
  separate "Pacific-ish" from "Atlantic-ish" for a first exploratory
  pass, not authoritative.
- Sites are capped (richest-first) for ordination/Mantel runtime
  tractability -- confirmed live that NMDS over ~6,300 real sites took
  82s with a single initialization; 2,000 keeps this interactive
  (~15s) without discarding the sparsest, least-informative sites.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS

from phaatlas.pipeline.cluster_ecology import haversine_km


@dataclass
class SiteCommunityMatrix:
    site_ids: list[str]
    cluster_ids: list[str]
    matrix: np.ndarray  # (n_sites, n_clusters) counts -- rows/cols in the order of site_ids/cluster_ids
    site_meta: dict[str, dict]  # site_id -> {lat, lon, depth_zone, ecosystem_type, ecosystem_name, study_id, ocean_basin}


def classify_ocean_basin(lat: float, lon: float) -> str:
    """Deliberately simple lon/lat bounding-box heuristic -- see module
    docstring for why this is an approximation, not authoritative."""
    if lat <= -60:
        return "Southern Ocean"
    if lat >= 66:
        return "Arctic Ocean"
    if 30 <= lat <= 46 and -6 <= lon <= 42:
        return "Mediterranean/Black Sea"
    if -70 <= lon <= 20:
        return "Atlantic Ocean"
    if 20 < lon <= 120:
        return "Indian Ocean"
    return "Pacific Ocean"  # lon > 120 or lon < -70 (wraps the antimeridian)


def build_site_cluster_matrix(
    metadata_depth_path: Path,
    cluster_assignments: dict[str, str],
    min_richness: int = 3,
    max_sites: int | None = 2000,
) -> SiteCommunityMatrix:
    """Site = sample_id. Drops sites with no resolvable lat/lon (nothing
    to test geographically) and sites with fewer than min_richness
    distinct clusters observed (too sparse for a meaningful community
    comparison). If more than max_sites qualify, keeps the max_sites
    richest -- a documented tractability tradeoff, not a silent one (see
    module docstring)."""
    import csv

    site_cluster_counts: dict[str, dict[str, int]] = {}
    site_meta: dict[str, dict] = {}

    with open(metadata_depth_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            site = row.get("sample_id", "")
            lat_raw, lon_raw = row.get("latitude_degN", ""), row.get("longitude_degE", "")
            if not site or not lat_raw or not lon_raw:
                continue
            cluster_id = cluster_assignments.get(row.get("target_id", ""))
            if cluster_id is None:
                continue

            counts = site_cluster_counts.setdefault(site, {})
            counts[cluster_id] = counts.get(cluster_id, 0) + 1

            if site not in site_meta:
                lat, lon = float(lat_raw), float(lon_raw)
                site_meta[site] = {
                    "lat": lat, "lon": lon,
                    "depth_zone": row.get("depth_zone", ""),
                    "ecosystem_type": row.get("ecosystem_type", ""),
                    "ecosystem_name": row.get("ecosystem_name", ""),
                    "study_id": row.get("study_id", ""),
                    "ocean_basin": classify_ocean_basin(lat, lon),
                }

    qualifying = [s for s, counts in site_cluster_counts.items() if len(counts) >= min_richness]
    qualifying.sort(key=lambda s: -len(site_cluster_counts[s]))
    if max_sites is not None and len(qualifying) > max_sites:
        qualifying = qualifying[:max_sites]
    site_ids = sorted(qualifying)  # alphabetical for a stable, reproducible matrix row order

    cluster_id_set: set[str] = set()
    for s in site_ids:
        cluster_id_set.update(site_cluster_counts[s].keys())
    cluster_ids = sorted(cluster_id_set)
    cluster_index = {c: i for i, c in enumerate(cluster_ids)}

    matrix = np.zeros((len(site_ids), len(cluster_ids)), dtype=np.float64)
    for i, s in enumerate(site_ids):
        for c, n in site_cluster_counts[s].items():
            matrix[i, cluster_index[c]] = n

    return SiteCommunityMatrix(
        site_ids=site_ids, cluster_ids=cluster_ids, matrix=matrix,
        site_meta={s: site_meta[s] for s in site_ids},
    )


def bray_curtis_distance_matrix(matrix: np.ndarray) -> np.ndarray:
    return squareform(pdist(matrix, metric="braycurtis"))


def jaccard_distance_matrix(matrix: np.ndarray) -> np.ndarray:
    return squareform(pdist((matrix > 0).astype(np.float64), metric="jaccard"))


def geographic_distance_matrix(site_ids: list[str], site_meta: dict[str, dict]) -> np.ndarray:
    n = len(site_ids)
    coords = [(site_meta[s]["lat"], site_meta[s]["lon"]) for s in site_ids]
    d = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            km = haversine_km(coords[i][0], coords[i][1], coords[j][0], coords[j][1])
            d[i, j] = d[j, i] = km
    return d


def nmds(
    distance_matrix: np.ndarray,
    n_components: int = 2,
    random_state: int = 0,
    n_init: int = 2,
    max_iter: int = 150,
) -> tuple[np.ndarray, float]:
    """Non-metric MDS (rank-order preserving) via scikit-learn's SMACOF
    implementation -- the standard ordination for community-composition
    dissimilarities like Bray-Curtis, which are not a true Euclidean
    metric. Returns (coords, stress) where stress is normalized (Stress-1,
    scikit-learn's own default for non-metric MDS): 0 perfect, 0.025
    excellent, 0.05 good, 0.1 fair, 0.2 poor -- NOT raw stress, whose
    magnitude scales with n_sites and is not comparable across runs."""
    model = MDS(
        n_components=n_components, metric=False, dissimilarity="precomputed",
        random_state=random_state, n_init=n_init, max_iter=max_iter, normalized_stress=True,
    )
    coords = model.fit_transform(distance_matrix)
    return coords, model.stress_


def pcoa(distance_matrix: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Classical (metric) MDS / Principal Coordinates Analysis via
    eigendecomposition of the double-centered Gower matrix -- the
    standard closed-form alternative to NMDS's iterative stress
    minimization. Returns (coords, pct_variance_explained) for the
    requested number of components; negative eigenvalues (common with a
    non-Euclidean dissimilarity like Bray-Curtis) are dropped rather than
    producing meaningless imaginary coordinates."""
    n = distance_matrix.shape[0]
    d2 = distance_matrix ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigenvalues, eigenvectors = np.linalg.eigh(b)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]

    positive = eigenvalues > 1e-10
    total_positive = eigenvalues[positive].sum()

    k = min(n_components, positive.sum())
    coords = eigenvectors[:, :k] * np.sqrt(eigenvalues[:k])
    if k < n_components:
        coords = np.pad(coords, ((0, 0), (0, n_components - k)))
    pct_variance = eigenvalues[:n_components].clip(min=0) / total_positive if total_positive > 0 else np.zeros(n_components)
    return coords, pct_variance


def mantel_test(
    dist_a: np.ndarray,
    dist_b: np.ndarray,
    n_permutations: int = 999,
    random_state: int = 0,
) -> tuple[float, float]:
    """Pearson correlation between the condensed (upper-triangle,
    excluding diagonal) forms of two square distance matrices over the
    SAME entities in the SAME order, with significance from a permutation
    test: dist_b's rows/columns are jointly shuffled (the correct Mantel
    procedure -- shuffling the condensed vector directly instead would
    silently break the underlying distance-matrix structure) and the
    correlation recomputed each time. Returns (r, p_value); p_value uses
    the standard "add-one" correction so it is never reported as exactly
    zero regardless of n_permutations."""
    n = dist_a.shape[0]
    iu = np.triu_indices(n, k=1)
    a_cond = dist_a[iu]
    b_cond = dist_b[iu]
    r_obs = float(np.corrcoef(a_cond, b_cond)[0, 1])

    rng = np.random.default_rng(random_state)
    count_as_extreme = 0
    for _ in range(n_permutations):
        perm = rng.permutation(n)
        b_perm_cond = dist_b[np.ix_(perm, perm)][iu]
        r_perm = np.corrcoef(a_cond, b_perm_cond)[0, 1]
        if abs(r_perm) >= abs(r_obs):
            count_as_extreme += 1
    p_value = (count_as_extreme + 1) / (n_permutations + 1)
    return r_obs, p_value
