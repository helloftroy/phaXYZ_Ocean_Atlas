"""Joins OMDB search hits (target_id = NR100 cluster representative) back
to genome/sample/study identity, GTDB taxonomy, and location/ecosystem
metadata -- entirely from data OMDB itself publishes, but not as a single
ready-made table. Two separate sources are needed:

1. OMDBv2.0_AA_G_NR100.cluster.tsv.gz (~4.2GB gzipped, ~249.5M rows -- one
   per NR100 cluster) -- maps each NR100 representative ("CLUSTER" column,
   our own target_id) to every member gene ID ("MEMBERS" column,
   semicolon-separated). Confirmed live via an HTTP range fetch of the
   real file: each row already carries the FULL member list, so a single
   streaming pass (no cross-row accumulation needed) resolves every
   wanted target_id's genome(s). A member gene ID has the form
   <GENOME>-scaffold_<N>_<gene#> (confirmed live by downloading one real
   genome's own .genes.faa.gz and inspecting its headers), so GENOME is
   recovered by string-splitting alone, no extra lookup.

2. omdb.microbiomics.io's own public JSON API (/api/genome-cols,
   /api/sample-cols) -- the same endpoints their own web genome browser
   calls, discovered by reading their bundled JS (genome_collection.js /
   tablerenderer.js), not documented anywhere and not linked from the bulk
   Downloads page. No auth required. Confirmed live: genome-cols returns
   full GTDB taxonomy (domain..species) per GENOME plus its SAMPLE/STUDY;
   sample-cols returns latitude_degN/longitude_degE plus
   ecosystem_type/ecosystem_name/ecosystem_compartment/sample_source.
   OMDB has no single raw "depth in meters" field anywhere (confirmed by
   reading their own column definitions in genome_collection.js) -- these
   categorical ecosystem fields are the closest available placement/
   depth-zone signal, not a substitute for one.

   The API supports batching many IDs into one request via a regex OR
   search_value_pair entry (confirmed live). IMPORTANT, also confirmed
   live: the API's regex match is NOT anchored on its own -- a truncated
   PREFIX of a real ID still matched as a substring. Every ID in a batch
   is therefore re.escape()'d and the whole alternation is wrapped in
   ^(?:...)$ here; skipping this would risk a batch silently pulling in
   wrong rows whenever one ID happens to be a substring of another.

No local database write -- this reads a family's existing
<family>_unique_targets.tsv and writes a sibling
<family>_unique_targets_with_metadata.tsv, one row per (target_id,
genome). A 100%-identity NR100 cluster can span more than one
genome/sample (the same exact protein sequence found in multiple
organisms/samples) -- that breadth is itself a meaningful signal, so it
is kept as separate rows rather than collapsed into one.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path

import httpx
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

OMDB_API_BASE = "https://omdb.microbiomics.io"
OMDB_REPO_ID = "ocean"
DEFAULT_CONTACT = "REPLACE_WITH_CONTACT_EMAIL@example.org"

CLUSTER_SUMMARY_COLUMNS = ["n_members_in_cluster", "n_genomes_in_cluster_total", "n_genomes_in_cluster_shown"]
GENOME_ROW_COLUMNS = [
    "genome_index", "genome", "genome_metadata_found",
    "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order", "gtdb_family", "gtdb_genus", "gtdb_species",
    "gtdb_version", "completeness", "contamination", "sample_id", "study_id",
]
SAMPLE_ROW_COLUMNS = [
    "sample_metadata_found", "environment", "ecosystem", "ecosystem_type", "ecosystem_name",
    "ecosystem_compartment", "sample_source", "latitude_degN", "longitude_degE",
    "biosample", "bioproject", "public_sample_link",
]


def _user_agent() -> str:
    contact = os.environ.get("PHA_REFERENCE_CONTACT_EMAIL", DEFAULT_CONTACT)
    return f"PHA-Ocean-Atlas-pha-reference/0.1 (research pipeline; contact: {contact})"


def _client() -> httpx.Client:
    return httpx.Client(headers={"User-Agent": _user_agent()}, timeout=60.0)


@retry(
    reraise=True,
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((httpx.TransportError, httpx.HTTPStatusError)),
)
def _post(client: httpx.Client, url: str, json_body: dict) -> httpx.Response:
    resp = client.post(url, json=json_body)
    if resp.status_code == 429:
        wait_s = float(resp.headers.get("Retry-After", "5"))
        time.sleep(wait_s)
        resp.raise_for_status()
    resp.raise_for_status()
    return resp


def parse_genome_from_member(member_id: str) -> str:
    """<GENOME>-scaffold_<N>_<gene#> -> <GENOME>, confirmed against a real
    downloaded .genes.faa.gz's own headers."""
    marker = "-scaffold_"
    idx = member_id.find(marker)
    if idx == -1:
        raise ValueError(f"Member id {member_id!r} does not contain the expected {marker!r} marker")
    return member_id[:idx]


@dataclass
class TargetClusterInfo:
    target_id: str
    n_members: int
    n_distinct_genomes: int  # true count, NOT capped
    genomes: list[str]  # capped at max_genomes_per_target


def stream_target_cluster_genomes(
    cluster_tsv_gz_path: Path,
    wanted_target_ids: set[str],
    max_genomes_per_target: int = 5,
) -> dict[str, TargetClusterInfo]:
    """Single streaming pass over OMDBv2.0_AA_G_NR100.cluster.tsv.gz
    (~249.5M rows). Every row already carries its cluster's full member
    list (see module docstring), so no cross-row accumulation is needed --
    only rows whose CLUSTER is in wanted_target_ids are even split, via a
    cheap substring check on the first tab-delimited field before paying
    for a full split, since the overwhelming majority of rows will not
    match for any one family's wanted set.
    """
    result: dict[str, TargetClusterInfo] = {}
    with gzip.open(cluster_tsv_gz_path, "rt", newline="") as f:
        f.readline()  # header: CLUSTER LENGTH SIZE REPRESENTATIVE MEMBERS
        for line in f:
            tab1 = line.find("\t")
            cluster_id = line[:tab1]
            if cluster_id not in wanted_target_ids:
                continue
            fields = line.rstrip("\n").split("\t")
            members = fields[4].split(";")
            distinct_genomes = sorted({parse_genome_from_member(m) for m in members})
            result[cluster_id] = TargetClusterInfo(
                target_id=cluster_id,
                n_members=int(fields[2]),
                n_distinct_genomes=len(distinct_genomes),
                genomes=distinct_genomes[:max_genomes_per_target],
            )
    return result


def _anchor_regex_batch(ids: list[str]) -> str:
    """Confirmed live against the real API: its regex search is NOT
    anchored on its own (a truncated prefix of a real ID still matched),
    so every alternative here is escaped and the whole pattern is wrapped
    in ^(?:...)$ to guarantee exact, unambiguous matches only."""
    escaped = "|".join(re.escape(i) for i in ids)
    return f"^(?:{escaped})$"


def _chunks(items: list[str], size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _query_omdb_api(
    client: httpx.Client,
    endpoint: str,
    key_field: str,
    ids: list[str],
    batch_size: int = 100,
    api_base: str = OMDB_API_BASE,
    sleep_seconds: float = 0.3,
) -> dict[str, dict]:
    """POSTs to <api_base>/api/<endpoint> in batches of batch_size ids,
    matching key_field via an anchored regex OR pattern (see
    _anchor_regex_batch). Returns {id: row_dict} for every id that had a
    match; ids with no row in the response are simply absent from the
    result -- callers report those as unresolved rather than guessing."""
    out: dict[str, dict] = {}
    unique_ids = sorted(set(ids))
    for batch in _chunks(unique_ids, batch_size):
        body = {
            "repo_id": OMDB_REPO_ID,
            "global_search_value": "",
            "search_value_pair": [{"key": key_field, "value": _anchor_regex_batch(batch), "regex": True}],
            "start": 0,
            "length": len(batch),
            "window": 0,
            "fetchall": 1,
            "order": [],
            "order_cols": [],
        }
        resp = _post(client, f"{api_base}/api/{endpoint}", body)
        data = resp.json()
        for row in data.get("data", []):
            key_value = row.get(key_field)
            if key_value is not None:
                out[key_value] = row
        if sleep_seconds:
            time.sleep(sleep_seconds)
    return out


def fetch_genome_metadata(
    client: httpx.Client,
    genome_ids: list[str],
    batch_size: int = 100,
    api_base: str = OMDB_API_BASE,
    sleep_seconds: float = 0.3,
) -> dict[str, dict]:
    return _query_omdb_api(client, "genome-cols", "genome", genome_ids, batch_size, api_base, sleep_seconds)


def fetch_sample_metadata(
    client: httpx.Client,
    sample_ids: list[str],
    batch_size: int = 100,
    api_base: str = OMDB_API_BASE,
    sleep_seconds: float = 0.3,
) -> dict[str, dict]:
    return _query_omdb_api(client, "sample-cols", "sushisample", sample_ids, batch_size, api_base, sleep_seconds)


@dataclass
class EnrichSummary:
    family_id: str
    n_input_rows: int = 0
    n_target_ids: int = 0
    n_targets_matched_in_cluster_file: int = 0
    n_distinct_genomes_needed: int = 0
    n_genomes_resolved: int = 0
    n_distinct_samples_needed: int = 0
    n_samples_resolved: int = 0
    n_output_rows: int = 0


def enrich_unique_targets(
    unique_targets_path: Path,
    cluster_tsv_gz_path: Path,
    out_path: Path,
    max_genomes_per_target: int = 5,
    genome_batch_size: int = 100,
    sample_batch_size: int = 100,
    api_base: str = OMDB_API_BASE,
    sleep_seconds: float = 0.3,
) -> EnrichSummary:
    """Reads unique_targets_path (a <family>_unique_targets.tsv), resolves
    every target_id to its NR100 cluster's genome(s) (streaming
    cluster_tsv_gz_path once), then batch-queries OMDB's own API for
    per-genome GTDB taxonomy and per-sample location/ecosystem metadata,
    and writes out_path with one row per (target_id, genome)."""
    family_id = unique_targets_path.stem.removesuffix("_unique_targets")
    summary = EnrichSummary(family_id=family_id)

    with open(unique_targets_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        original_fieldnames = reader.fieldnames or []
        rows = list(reader)
    summary.n_input_rows = len(rows)

    wanted_target_ids = {r["target_id"] for r in rows}
    summary.n_target_ids = len(wanted_target_ids)

    target_to_cluster = stream_target_cluster_genomes(cluster_tsv_gz_path, wanted_target_ids, max_genomes_per_target)
    summary.n_targets_matched_in_cluster_file = len(target_to_cluster)

    all_genomes: set[str] = set()
    for info in target_to_cluster.values():
        all_genomes.update(info.genomes)
    summary.n_distinct_genomes_needed = len(all_genomes)

    with _client() as client:
        genome_metadata = fetch_genome_metadata(client, sorted(all_genomes), genome_batch_size, api_base, sleep_seconds)
        summary.n_genomes_resolved = len(genome_metadata)

        all_samples = {row["shushisample"] for row in genome_metadata.values() if row.get("shushisample")}
        summary.n_distinct_samples_needed = len(all_samples)

        sample_metadata = fetch_sample_metadata(client, sorted(all_samples), sample_batch_size, api_base, sleep_seconds)
        summary.n_samples_resolved = len(sample_metadata)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = original_fieldnames + CLUSTER_SUMMARY_COLUMNS + GENOME_ROW_COLUMNS + SAMPLE_ROW_COLUMNS
    n_output_rows = 0
    with open(out_path, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            info = target_to_cluster.get(row["target_id"])
            genomes = info.genomes if info else []
            cluster_common = {
                "n_members_in_cluster": info.n_members if info else "",
                "n_genomes_in_cluster_total": info.n_distinct_genomes if info else "",
                "n_genomes_in_cluster_shown": len(genomes),
            }
            if not genomes:
                out_row = {**row, **cluster_common}
                for col in GENOME_ROW_COLUMNS + SAMPLE_ROW_COLUMNS:
                    out_row[col] = False if col.endswith("_found") else ""
                writer.writerow(out_row)
                n_output_rows += 1
                continue

            for idx, genome in enumerate(genomes, start=1):
                g = genome_metadata.get(genome)
                genome_cols = {
                    "genome_index": idx,
                    "genome": genome,
                    "genome_metadata_found": g is not None,
                    "gtdb_domain": g.get("domain", "") if g else "",
                    "gtdb_phylum": g.get("phylum", "") if g else "",
                    "gtdb_class": g.get("class_", "") if g else "",
                    "gtdb_order": g.get("order", "") if g else "",
                    "gtdb_family": g.get("family", "") if g else "",
                    "gtdb_genus": g.get("genus", "") if g else "",
                    "gtdb_species": g.get("species", "") if g else "",
                    "gtdb_version": g.get("gtdb_version", "") if g else "",
                    "completeness": g.get("completeness", "") if g else "",
                    "contamination": g.get("contamination", "") if g else "",
                    "sample_id": g.get("shushisample", "") if g else "",
                    "study_id": g.get("sushistudy", "") if g else "",
                }
                sample_id = genome_cols["sample_id"]
                s = sample_metadata.get(sample_id) if sample_id else None
                sample_cols = {
                    "sample_metadata_found": s is not None,
                    "environment": s.get("environment", "") if s else "",
                    "ecosystem": s.get("ecosystem", "") if s else "",
                    "ecosystem_type": s.get("ecosystem_type", "") if s else "",
                    "ecosystem_name": s.get("ecosystem_name", "") if s else "",
                    "ecosystem_compartment": s.get("ecosystem_compartment", "") if s else "",
                    "sample_source": s.get("sample_source", "") if s else "",
                    "latitude_degN": s.get("latitude_degN", "") if s else "",
                    "longitude_degE": s.get("longitude_degE", "") if s else "",
                    "biosample": s.get("biosample", "") if s else "",
                    "bioproject": s.get("bioproject", "") if s else "",
                    "public_sample_link": s.get("public_sample_link", "") if s else "",
                }
                writer.writerow({**row, **cluster_common, **genome_cols, **sample_cols})
                n_output_rows += 1
    summary.n_output_rows = n_output_rows
    return summary
