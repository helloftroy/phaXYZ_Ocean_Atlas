"""Map phaC genome locations/depths to WOA23 annual temperature.

Inputs:
  - fair_ocean_agent/phaC_unique_targets_with_metadata_depth.tsv
  - physiochem_data_sources/WOA23/woa23_decav_t00an01.csv.gz

Output:
  - fair_ocean_agent/phaC_genomes_woa23_annual_temperature.tsv

The WOA CSV is the NOAA/NCEI WOA23 1-degree annual temperature
climatological mean product. This script uses nearest horizontal 1-degree
grid-cell center and linear interpolation between WOA standard depths when
both bracket values are present. If the exact nearest horizontal cell lacks
data at the requested depth, it searches nearby cells out to 3 degrees and
records the match status.
"""
import csv
import gzip
import math
from bisect import bisect_left
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent.parent
WORKSPACE = ROOT.parent
PHAC = WORKSPACE / "fair_ocean_agent" / "phaC_unique_targets_with_metadata_depth.tsv"
WOA = WORKSPACE / "physiochem_data_sources" / "WOA23" / "woa23_decav_t00an01.csv.gz"
OUT = WORKSPACE / "fair_ocean_agent" / "phaC_genomes_woa23_annual_temperature.tsv"


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
    lon = ((lon + 180.0) % 360.0) - 180.0
    # Keep +180 mapped to -180 convention, then center with floor()+0.5.
    return lon


def load_woa(path):
    depths = parse_depths_from_header(path)
    grid = {}
    by_latlon_index = {}
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
            by_latlon_index.setdefault(lat, set()).add(lon)
    return depths, grid, by_latlon_index


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
        return v, d, d, v, v, "exact_depth"
    if 0 < i < len(valid):
        d0, v0 = valid[i - 1]
        d1, v1 = valid[i]
        frac = (target_depth - d0) / (d1 - d0)
        return v0 + frac * (v1 - v0), d0, d1, v0, v1, "depth_interpolated"
    # Outside available depth range: nearest valid depth only.
    if i == 0:
        d, v = valid[0]
    else:
        d, v = valid[-1]
    return v, d, d, v, v, "nearest_depth_outside_available_range"


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
            temp, d0, d1, v0, v1, depth_method = interp
            dist = horizontal_km(lat, lon, cell[0], cell[1])
            record = (dist, radius, cell, temp, d0, d1, v0, v1, depth_method)
            if best is None or record[0] < best[0]:
                best = record
        if best is not None:
            break
    if best is None:
        return None
    dist, radius, cell, temp, d0, d1, v0, v1, depth_method = best
    status = "matched_exact_1deg_cell" if radius == 0 else f"matched_neighbor_within_{radius}deg"
    return {
        "woa23_temp_annual_degC": temp,
        "woa23_lat": cell[0],
        "woa23_lon": cell[1],
        "woa23_horizontal_distance_km": dist,
        "woa23_depth_lower_m": d0,
        "woa23_depth_upper_m": d1,
        "woa23_temp_lower_degC": v0,
        "woa23_temp_upper_degC": v1,
        "woa23_depth_method": depth_method,
        "woa23_match_status": status,
    }


def unique_genome_rows(path):
    wanted = [
        "genome", "sample_id", "study_id", "biosample", "bioproject",
        "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
        "gtdb_family", "gtdb_genus", "gtdb_species",
        "latitude_degN", "longitude_degE", "depth_raw", "depth_m", "depth_zone",
    ]
    seen = {}
    target_counts = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genome = row.get("genome", "")
            if not genome:
                continue
            target_counts[genome] = target_counts.get(genome, 0) + 1
            if genome not in seen:
                seen[genome] = {k: row.get(k, "") for k in wanted}
    for genome, row in seen.items():
        row["n_phac_target_rows"] = str(target_counts.get(genome, 0))
        yield row


def main():
    depths, grid, _ = load_woa(WOA)
    rows = list(unique_genome_rows(PHAC))
    fieldnames = [
        "genome", "sample_id", "study_id", "biosample", "bioproject",
        "n_phac_target_rows", "gtdb_domain", "gtdb_phylum", "gtdb_class",
        "gtdb_order", "gtdb_family", "gtdb_genus", "gtdb_species",
        "latitude_degN", "longitude_degE", "depth_raw", "depth_m", "depth_zone",
        "woa23_temp_annual_degC", "woa23_lat", "woa23_lon",
        "woa23_horizontal_distance_km", "woa23_depth_lower_m",
        "woa23_depth_upper_m", "woa23_temp_lower_degC", "woa23_temp_upper_degC",
        "woa23_depth_method", "woa23_match_status",
    ]

    counts = {"rows": 0, "attempted": 0, "matched": 0, "missing_input": 0, "unmatched": 0}
    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            counts["rows"] += 1
            lat_raw, lon_raw, depth_raw = row.get("latitude_degN", ""), row.get("longitude_degE", ""), row.get("depth_m", "")
            out = {k: row.get(k, "") for k in fieldnames}
            if not (lat_raw and lon_raw and depth_raw):
                out["woa23_match_status"] = "missing_lat_lon_or_depth"
                counts["missing_input"] += 1
            else:
                counts["attempted"] += 1
                match = match_woa(depths, grid, float(lat_raw), float(lon_raw), float(depth_raw))
                if match is None:
                    out["woa23_match_status"] = "no_woa_value_within_3deg"
                    counts["unmatched"] += 1
                else:
                    counts["matched"] += 1
                    out.update({
                        k: (f"{v:.4f}" if isinstance(v, float) else v)
                        for k, v in match.items()
                    })
            writer.writerow(out)

    print(f"wrote {OUT}")
    print(counts)


if __name__ == "__main__":
    main()
