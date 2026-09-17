"""Export all OMDB genomes with sample/study/location metadata.

This builds the denominator table for "% of all genomes with phaC".
OMDBv2.0_data.tsv.gz lists every genome and its SAMPLE/STUDY, but not
lat/lon; lat/lon are fetched from OMDB's sample metadata API.

Default output:
  data/all_genomes/omdb_all_genomes_with_locations.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import time
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent.parent


DEFAULT_CATALOG = ROOT / "PHA_bioprospecting" / "databases" / "OMDBv2" / "OMDBv2.0_data.tsv.gz"
DEFAULT_OUT = ROOT / "data" / "all_genomes" / "omdb_all_genomes_with_locations.tsv"
OMDB_API_BASE = "https://omdb.microbiomics.io"
SAMPLE_ROW_COLUMNS = [
    "sample_metadata_found", "environment", "ecosystem", "ecosystem_type", "ecosystem_name",
    "ecosystem_compartment", "sample_source", "latitude_degN", "longitude_degE",
    "biosample", "bioproject", "public_sample_link",
]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def read_catalog(path: Path) -> list[dict[str, str]]:
    rows = []
    with open_text(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append({
                "genome": row.get("GENOME", ""),
                "sample_id": row.get("SAMPLE", ""),
                "study_id": row.get("STUDY", ""),
                "genome_file": row.get("GENOME_FILE", ""),
                "genes_aa_file": row.get("GENES_AA_FILE", ""),
            })
    return rows


def chunks(items: list[str], size: int):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def anchored_regex(ids: list[str]) -> str:
    return "^(?:" + "|".join(re.escape(i) for i in ids) + ")$"


def post_json(url: str, body: dict, retries: int = 5) -> dict:
    payload = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "PHA-Ocean-Atlas/0.1 denominator-export",
        },
        method="POST",
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(min(30, 2 ** attempt))
    raise RuntimeError("unreachable")


def fetch_sample_metadata(
    sample_ids: list[str],
    batch_size: int,
    api_base: str,
    sleep_seconds: float,
) -> dict[str, dict]:
    out = {}
    for batch in chunks(sorted(set(sample_ids)), batch_size):
        body = {
            "repo_id": "ocean",
            "global_search_value": "",
            "search_value_pair": [
                {"key": "sushisample", "value": anchored_regex(batch), "regex": True}
            ],
            "start": 0,
            "length": len(batch),
            "window": 0,
            "fetchall": 1,
            "order": [],
            "order_cols": [],
        }
        data = post_json(f"{api_base}/api/sample-cols", body)
        for row in data.get("data", []):
            sample_id = row.get("sushisample")
            if sample_id:
                out[sample_id] = row
        if sleep_seconds:
            time.sleep(sleep_seconds)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--sample-batch-size", type=int, default=100)
    ap.add_argument("--sleep-seconds", type=float, default=0.3)
    ap.add_argument("--api-base", default=OMDB_API_BASE)
    args = ap.parse_args()

    if not args.catalog.exists():
        raise SystemExit(
            f"OMDB catalog not found: {args.catalog}\n"
            "Download it on the cluster with:\n"
            "  cd PHA_bioprospecting/scripts\n"
            "  ./download_omdb.sh data"
        )

    genomes = read_catalog(args.catalog)
    sample_ids = sorted({r["sample_id"] for r in genomes if r["sample_id"]})
    print(f"catalog genomes: {len(genomes):,}")
    print(f"distinct samples to resolve: {len(sample_ids):,}")

    sample_meta = fetch_sample_metadata(
        sample_ids,
        batch_size=args.sample_batch_size,
        api_base=args.api_base,
        sleep_seconds=args.sleep_seconds,
    )
    print(f"samples resolved by OMDB API: {len(sample_meta):,}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "genome", "sample_id", "study_id", "genome_file", "genes_aa_file",
    ] + SAMPLE_ROW_COLUMNS
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in genomes:
            meta = sample_meta.get(row["sample_id"], {})
            out = dict(row)
            out.update({k: meta.get(k, "") for k in SAMPLE_ROW_COLUMNS})
            out["sample_metadata_found"] = "yes" if meta else "no"
            writer.writerow(out)

    n_with_location = sum(
        1 for r in genomes
        if sample_meta.get(r["sample_id"], {}).get("latitude_degN")
        and sample_meta.get(r["sample_id"], {}).get("longitude_degE")
    )
    print(f"genomes with lat/lon: {n_with_location:,}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
