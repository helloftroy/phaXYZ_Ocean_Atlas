"""For each protein sequence cluster (see pipeline/sequence_clustering.py
for how clusters are built), summarizes whether its host genomes'
sample locations are geographically tight or scattered, what its
median sample depth is (where resolvable), and which taxa carry it --
directly answering "do individual protein clusters occupy distinct
environmental niches".

Geographic spread uses proper spherical statistics, not a flat lat/lon
standard deviation (which badly distorts near the poles and across the
antimeridian): each sample location becomes a unit vector on the sphere,
and the length of their mean vector -- the "mean resultant length", R,
a standard circular/directional-statistics concentration measure -- is 1
when every point coincides and approaches 0 as points spread out over
the whole globe. Reported alongside the more intuitive max pairwise
great-circle distance (km) between locations in the cluster.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

EARTH_RADIUS_KM = 6371.0088


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(a))


def _to_unit_vector(lat_deg: float, lon_deg: float) -> tuple[float, float, float]:
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    return (math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat))


def _subsample_deterministic(points: list, cap: int) -> list:
    """Evenly-strided, not random -- keeps max-pairwise-distance estimates
    reproducible across runs for large clusters, at the cost of being an
    estimate rather than the exact max once a cluster exceeds cap points."""
    if len(points) <= cap:
        return points
    step = len(points) / cap
    return [points[int(i * step)] for i in range(cap)]


@dataclass
class GeoDispersion:
    n: int
    mean_resultant_length: float | None  # 1.0 = all identical location, ->0 = spread over whole globe
    centroid_lat: float | None
    centroid_lon: float | None
    max_pairwise_km: float | None  # exact if n <= max_pairwise_sample_cap, else an estimate from a deterministic subsample


def spherical_dispersion(points: list[tuple[float, float]], max_pairwise_sample_cap: int = 200) -> GeoDispersion:
    if not points:
        return GeoDispersion(n=0, mean_resultant_length=None, centroid_lat=None, centroid_lon=None, max_pairwise_km=None)

    vectors = [_to_unit_vector(lat, lon) for lat, lon in points]
    n = len(vectors)
    mx = sum(v[0] for v in vectors) / n
    my = sum(v[1] for v in vectors) / n
    mz = sum(v[2] for v in vectors) / n
    r = math.sqrt(mx * mx + my * my + mz * mz)
    centroid_lat = math.degrees(math.atan2(mz, math.sqrt(mx * mx + my * my)))
    centroid_lon = math.degrees(math.atan2(my, mx))

    sample = _subsample_deterministic(points, max_pairwise_sample_cap)
    max_km = 0.0
    for i in range(len(sample)):
        for j in range(i + 1, len(sample)):
            d = haversine_km(sample[i][0], sample[i][1], sample[j][0], sample[j][1])
            if d > max_km:
                max_km = d

    return GeoDispersion(n=n, mean_resultant_length=r, centroid_lat=centroid_lat, centroid_lon=centroid_lon, max_pairwise_km=max_km)


@dataclass
class ClusterEcologyStats:
    cluster_id: str
    n_target_ids: int
    n_genomes: int
    geo: GeoDispersion
    median_depth_m: float | None
    n_genomes_with_depth: int
    top_genera: list[tuple[str, int]]
    top_phyla: list[tuple[str, int]]
    top_studies: list[tuple[str, int]]


@dataclass
class _ClusterAccumulator:
    target_ids: set[str] = field(default_factory=set)
    genomes: dict[str, dict] = field(default_factory=dict)  # genome -> {lat, lon, depth_m, genus, phylum, study}


def summarize_cluster_ecology(
    metadata_depth_path: Path,
    cluster_assignments: dict[str, str],
    top_n: int = 5,
) -> list[ClusterEcologyStats]:
    """metadata_depth_path: a <family>_unique_targets_with_metadata_depth.tsv
    (depth columns optional -- median_depth_m/n_genomes_with_depth are
    just 0/None if absent). cluster_assignments: {target_id: cluster_id},
    e.g. from sequence_clustering.load_cluster_assignments. A target_id
    with no entry in cluster_assignments is skipped (its cluster wasn't
    built -- see that module for why this can legitimately happen for a
    sequence mmseqs' own clustering treated differently)."""
    accs: dict[str, _ClusterAccumulator] = {}

    with open(metadata_depth_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            target_id = row.get("target_id", "")
            cluster_id = cluster_assignments.get(target_id)
            if cluster_id is None:
                continue
            genome = row.get("genome", "")
            if not genome:
                continue

            acc = accs.setdefault(cluster_id, _ClusterAccumulator())
            acc.target_ids.add(target_id)
            if genome not in acc.genomes:
                lat = _to_float(row.get("latitude_degN"))
                lon = _to_float(row.get("longitude_degE"))
                depth_m = _to_float(row.get("depth_m"))
                acc.genomes[genome] = {
                    "lat": lat, "lon": lon, "depth_m": depth_m,
                    "genus": row.get("gtdb_genus", ""), "phylum": row.get("gtdb_phylum", ""),
                    "study": row.get("study_id", ""),
                }

    results: list[ClusterEcologyStats] = []
    for cluster_id, acc in accs.items():
        points = [(g["lat"], g["lon"]) for g in acc.genomes.values() if g["lat"] is not None and g["lon"] is not None]
        depths = [g["depth_m"] for g in acc.genomes.values() if g["depth_m"] is not None]
        genera = Counter(g["genus"] for g in acc.genomes.values() if g["genus"])
        phyla = Counter(g["phylum"] for g in acc.genomes.values() if g["phylum"])
        studies = Counter(g["study"] for g in acc.genomes.values() if g["study"])

        results.append(ClusterEcologyStats(
            cluster_id=cluster_id,
            n_target_ids=len(acc.target_ids),
            n_genomes=len(acc.genomes),
            geo=spherical_dispersion(points),
            median_depth_m=median(depths) if depths else None,
            n_genomes_with_depth=len(depths),
            top_genera=genera.most_common(top_n),
            top_phyla=phyla.most_common(top_n),
            top_studies=studies.most_common(top_n),
        ))

    results.sort(key=lambda s: -s.n_genomes)
    return results


def _to_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def write_cluster_ecology(stats: list[ClusterEcologyStats], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "cluster_id", "n_target_ids", "n_genomes",
            "geo_mean_resultant_length", "geo_centroid_lat", "geo_centroid_lon", "geo_max_pairwise_km", "geo_n",
            "median_depth_m", "n_genomes_with_depth",
            "top_genera", "top_phyla", "top_studies",
        ])
        for s in stats:
            g = s.geo
            writer.writerow([
                s.cluster_id, s.n_target_ids, s.n_genomes,
                _fmt(g.mean_resultant_length), _fmt(g.centroid_lat), _fmt(g.centroid_lon), _fmt(g.max_pairwise_km), g.n,
                _fmt(s.median_depth_m), s.n_genomes_with_depth,
                "; ".join(f"{k} ({v})" for k, v in s.top_genera),
                "; ".join(f"{k} ({v})" for k, v in s.top_phyla),
                "; ".join(f"{k} ({v})" for k, v in s.top_studies),
            ])


def _fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.4f}"
