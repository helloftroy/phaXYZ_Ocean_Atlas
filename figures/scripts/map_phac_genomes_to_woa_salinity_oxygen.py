"""Map phaC genome locations/depths to WOA23 annual salinity and dissolved
oxygen -- the same join already done for temperature
(map_phac_genomes_to_woa_temperature.py), extended to the other two
variables. Self-contained (duplicates that script's matching logic
rather than importing it) rather than risk touching already-working,
separately-maintained code for a one-off extension.

Inputs:
  - PHA_bioprospecting/omdb_search/results/phaC_unique_targets_with_metadata_depth.tsv
  - physiochem_data_sources/WOA23/woa23_decav_s00an01.csv.gz (salinity, PSU)
  - physiochem_data_sources/WOA23/woa23_all_o00an01.csv.gz (dissolved O2, umol/kg --
    this variable uses the "all" multi-decade statistical period, not "decav",
    since less historical oxygen data exists than temperature/salinity; confirmed
    live against the NCEI server, "decav" 404s for oxygen, "all" is the real path)

Output:
  - PHA_bioprospecting/omdb_search/results/phaC_genomes_woa23_salinity_oxygen.tsv

What WOA23 actually gives per grid cell/depth, for anyone reusing this: a
single OBJECTIVELY ANALYZED CLIMATOLOGICAL MEAN (the "an" product) -- not
a raw in-situ measurement, and not a min/max range. A companion STANDARD
DEVIATION product ("sd" in place of "an" in the filename) exists for
every WOA23 variable if a breadth/spread proxy is ever wanted instead of
or alongside the mean -- not fetched here, since the direct ask was for
the mean values themselves.

QC: full load_bad_targets() (BAD_QUERIES-derived target exclusions +
BAD_TARGET_IDS), not just the best_query-only check the original
temperature script used -- matches this project's current standard
convention. This matters here specifically: the phaC reference-query fix
(2026-09-22, PHA_CLEAN_RESULTS.md section 2) more than halved the
QC-passing phaC-positive genome count (68,424 -> 31,464), and only
~34% of even that corrected population has lat/lon/depth populated
at all (confirmed live below) -- both real constraints on what this
join can possibly cover, not a bug in this script.
"""
import csv
import gzip
import math
import sys
from bisect import bisect_left
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "figures" / "scripts"))
import _phac_qc

FA = ROOT / "PHA_bioprospecting" / "omdb_search" / "results"
PHAC = FA / "phaC_unique_targets_with_metadata_depth.tsv"
WOA_DIR = Path("/Users/hellpark/multimodal_seusmbol/physiochem_data_sources/WOA23")
WOA_SALINITY = WOA_DIR / "woa23_decav_s00an01.csv.gz"
WOA_OXYGEN = WOA_DIR / "woa23_all_o00an01.csv.gz"
OUT = FA / "phaC_genomes_woa23_salinity_oxygen.tsv"


def parse_depths_from_header(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#COMMA SEPARATED"):
                marker = "DEPTHS (M):"
                return [float(x) for x in line.split(marker, 1)[1].split(",")]
    raise ValueError(f"Could not find WOA depth header in {path}")


def nearest_center(value):
    return math.floor(value) + 0.5


def norm_lon(lon):
    return ((lon + 180.0) % 360.0) - 180.0


def load_woa(path):
    depths = parse_depths_from_header(path)
    grid = {}
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.reader(row for row in f if not row.startswith("#"))
        for row in reader:
            if len(row) < 3:
                continue
            lat = round(float(row[0]), 3)
            lon = round(float(row[1]), 3)
            vals = []
            for raw in row[2:]:
                raw = raw.strip()
                vals.append(float(raw) if raw else None)
            if len(vals) < len(depths):
                vals.extend([None] * (len(depths) - len(vals)))
            grid[(lat, lon)] = vals[:len(depths)]
    return depths, grid


def horizontal_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(((lon2 - lon1 + 180.0) % 360.0) - 180.0)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def interpolate_at_depth(depths, vals, target_depth):
    valid = [(d, v) for d, v in zip(depths, vals) if v is not None]
    if not valid:
        return None
    vd = [d for d, _ in valid]
    i = bisect_left(vd, target_depth)
    if i < len(valid) and valid[i][0] == target_depth:
        d, v = valid[i]
        return v, "exact_depth"
    if 0 < i < len(valid):
        d0, v0 = valid[i - 1]
        d1, v1 = valid[i]
        frac = (target_depth - d0) / (d1 - d0)
        return v0 + frac * (v1 - v0), "depth_interpolated"
    d, v = valid[0] if i == 0 else valid[-1]
    return v, "nearest_depth_outside_available_range"


def candidate_cells(lat0, lon0, radius_deg):
    lat_start = math.floor(lat0 - radius_deg) + 0.5
    lat_end = math.floor(lat0 + radius_deg) + 0.5
    lat = lat_start
    while lat <= lat_end + 1e-9:
        lon_start = math.floor(lon0 - radius_deg) + 0.5
        lon_end = math.floor(lon0 + radius_deg) + 0.5
        lon = lon_start
        while lon <= lon_end + 1e-9:
            yield round(lat, 3), round(norm_lon(lon), 3)
            lon += 1.0
        lat += 1.0


def match_woa(depths, grid, lat, lon, depth_m):
    lat0 = round(nearest_center(lat), 3)
    lon0 = round(nearest_center(norm_lon(lon)), 3)
    best = None
    for radius in (0, 1, 2, 3):
        for cell in candidate_cells(lat0, lon0, radius):
            vals = grid.get(cell)
            if vals is None:
                continue
            interp = interpolate_at_depth(depths, vals, depth_m)
            if interp is None:
                continue
            value, depth_method = interp
            dist = horizontal_km(lat, lon, cell[0], cell[1])
            record = (dist, radius, cell, value, depth_method)
            if best is None or record[0] < best[0]:
                best = record
        if best is not None:
            break
    if best is None:
        return None
    dist, radius, cell, value, depth_method = best
    status = "matched_exact_1deg_cell" if radius == 0 else f"matched_neighbor_within_{radius}deg"
    return value, dist, depth_method, status


def unique_genome_rows(path, bad_targets):
    wanted = [
        "genome", "sample_id", "study_id", "biosample", "bioproject",
        "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
        "gtdb_family", "gtdb_genus", "gtdb_species",
        "latitude_degN", "longitude_degE", "depth_raw", "depth_m", "depth_zone",
    ]
    seen = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("target_id", "") in bad_targets:
                continue
            genome = row.get("genome", "")
            if not genome or genome in seen:
                continue
            seen[genome] = {k: row.get(k, "") for k in wanted}
    return seen.values()


def main():
    bad_targets = _phac_qc.load_bad_targets()
    print(f"loading WOA23 salinity grid ({WOA_SALINITY.name}) ...")
    s_depths, s_grid = load_woa(WOA_SALINITY)
    print(f"loading WOA23 oxygen grid ({WOA_OXYGEN.name}) ...")
    o_depths, o_grid = load_woa(WOA_OXYGEN)

    rows = list(unique_genome_rows(PHAC, bad_targets))
    n_with_latlondepth = sum(1 for r in rows if r["latitude_degN"] and r["longitude_degE"] and r["depth_m"])
    print(f"{len(rows):,} QC-passing distinct phaC-positive genomes in {PHAC.name}")
    print(f"{n_with_latlondepth:,}/{len(rows):,} ({100*n_with_latlondepth/len(rows):.1f}%) have lat/lon/depth populated "
          f"-- WOA matching can only ever cover this subset")

    fieldnames = [
        "genome", "sample_id", "study_id", "biosample", "bioproject",
        "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
        "gtdb_family", "gtdb_genus", "gtdb_species",
        "latitude_degN", "longitude_degE", "depth_raw", "depth_m", "depth_zone",
        "woa23_salinity_psu", "woa23_salinity_distance_km", "woa23_salinity_depth_method", "woa23_salinity_match_status",
        "woa23_oxygen_umol_kg", "woa23_oxygen_distance_km", "woa23_oxygen_depth_method", "woa23_oxygen_match_status",
    ]
    counts = {"attempted": 0, "both_matched": 0, "missing_input": 0}
    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = {k: row.get(k, "") for k in fieldnames}
            lat_raw, lon_raw, depth_raw = row["latitude_degN"], row["longitude_degE"], row["depth_m"]
            if not (lat_raw and lon_raw and depth_raw):
                out["woa23_salinity_match_status"] = "missing_lat_lon_or_depth"
                out["woa23_oxygen_match_status"] = "missing_lat_lon_or_depth"
                counts["missing_input"] += 1
                writer.writerow(out)
                continue
            counts["attempted"] += 1
            lat, lon, depth = float(lat_raw), float(lon_raw), float(depth_raw)

            s_match = match_woa(s_depths, s_grid, lat, lon, depth)
            if s_match:
                value, dist, method, status = s_match
                out["woa23_salinity_psu"] = f"{value:.4f}"
                out["woa23_salinity_distance_km"] = f"{dist:.4f}"
                out["woa23_salinity_depth_method"] = method
                out["woa23_salinity_match_status"] = status
            else:
                out["woa23_salinity_match_status"] = "no_woa_value_within_3deg"

            o_match = match_woa(o_depths, o_grid, lat, lon, depth)
            if o_match:
                value, dist, method, status = o_match
                out["woa23_oxygen_umol_kg"] = f"{value:.4f}"
                out["woa23_oxygen_distance_km"] = f"{dist:.4f}"
                out["woa23_oxygen_depth_method"] = method
                out["woa23_oxygen_match_status"] = status
            else:
                out["woa23_oxygen_match_status"] = "no_woa_value_within_3deg"

            if s_match and o_match:
                counts["both_matched"] += 1
            writer.writerow(out)

    print(f"wrote {OUT}")
    print(counts)


if __name__ == "__main__":
    main()
