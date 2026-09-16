"""Adds real numeric sample depth (meters) to an already-enriched
<family>_unique_targets_with_metadata.tsv, by fetching each sample's
underlying NCBI BioSample record and reading its harmonized "depth"
attribute -- confirmed live this exists for at least some OMDB samples
(e.g. SAMN10388028, a BATS time-series sample: harmonized_name="depth" ->
"0 m") even though OMDB's own genome/sample API (pipeline/omdb_metadata.py)
has no depth field at all.

Coverage is inherently PARTIAL and uneven across studies, for two
independent reasons, both kept visible rather than silently dropped:
1. Not every OMDB sample has a real NCBI BioSample accession in the first
   place -- omdb_metadata.py's "biosample" field is "Unknown_biosample"
   for those (e.g. JGI GOLD-only samples), and those never reach NCBI at
   all here.
2. Even a real BioSample record may simply not have reported depth --
   common for sediment/isolation-source-only submissions, much less
   common for water-column time-series studies (BATS, Tara-style) that
   follow MIxS environmental-package conventions more closely.

depth_raw (the exact NCBI attribute text) is always kept alongside
depth_m (this module's best-effort numeric parse of it), so a value this
parser got wrong or couldn't parse is visible and correctable rather than
an indistinguishable blank.
"""

from __future__ import annotations

import csv
import os
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import httpx
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_CONTACT = "REPLACE_WITH_CONTACT_EMAIL@example.org"
UNRESOLVABLE_BIOSAMPLE_MARKERS = {"", "unknown_biosample", "not applicable", "na", "none"}

# (label, lower_bound_inclusive, upper_bound_exclusive) -- upper=None means unbounded.
DEPTH_BINS: list[tuple[str, float, float | None]] = [
    ("0-50m", 0.0, 50.0),
    ("50-200m", 50.0, 200.0),
    ("200-1000m", 200.0, 1000.0),
    ("1000-4000m", 1000.0, 4000.0),
    (">4000m", 4000.0, None),
]

_MISSING_TOKENS = {"not collected", "not applicable", "missing", "n/a", "na", "not provided", "unknown", "restricted access"}
_SURFACE_RE = re.compile(r"surface", re.IGNORECASE)
_RANGE_RE = re.compile(r"(-?\d+(?:\.\d+)?)\s*m?\s*(?:-|to)\s*(-?\d+(?:\.\d+)?)\s*m?", re.IGNORECASE)
_SINGLE_RE = re.compile(r"(-?\d+(?:\.\d+)?)\s*m?\b", re.IGNORECASE)


def parse_depth_meters(raw: str | None) -> float | None:
    """Best-effort parse of an NCBI BioSample "depth" attribute's free-text
    value into meters. Handles "0 m", "200", "10-20 m" (midpoint), and
    "surface" (0.0); returns None for missing/unrecognized text rather
    than guessing."""
    if raw is None:
        return None
    text = raw.strip()
    if not text or text.lower() in _MISSING_TOKENS:
        return None
    if _SURFACE_RE.search(text):
        return 0.0
    m = _RANGE_RE.search(text)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        return (lo + hi) / 2.0
    m = _SINGLE_RE.search(text)
    if m:
        return float(m.group(1))
    return None


def depth_zone(depth_m: float | None) -> str:
    """Bins a parsed depth (meters) into DEPTH_BINS' labels. A negative
    depth (rare -- e.g. a splash-zone/above-sea-level artifact in some
    submissions) is clipped to 0 rather than dropped, since the practical
    difference from "at the surface" is negligible for this purpose.
    Returns "" (not a bin) when depth_m is None."""
    if depth_m is None:
        return ""
    d = max(depth_m, 0.0)
    for label, lo, hi in DEPTH_BINS:
        if d >= lo and (hi is None or d < hi):
            return label
    return ""  # unreachable given DEPTH_BINS' last bound is None (unbounded)


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
def _get(client: httpx.Client, url: str, params: dict) -> httpx.Response:
    resp = client.get(url, params=params)
    if resp.status_code == 429:
        wait_s = float(resp.headers.get("Retry-After", "5"))
        time.sleep(wait_s)
        resp.raise_for_status()
    resp.raise_for_status()
    return resp


def _chunks(items: list[str], size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _extract_depth_raw(biosample_xml_element: ET.Element) -> str | None:
    attrs = biosample_xml_element.find("Attributes")
    if attrs is None:
        return None
    for attr in attrs.findall("Attribute"):
        if attr.get("harmonized_name") == "depth":
            return (attr.text or "").strip()
    # fall back to a case-insensitive raw attribute_name match if harmonization missed it
    for attr in attrs.findall("Attribute"):
        if (attr.get("attribute_name") or "").strip().lower() == "depth":
            return (attr.text or "").strip()
    return None


def fetch_sample_depths(
    client: httpx.Client,
    biosample_accessions: list[str],
    batch_size: int = 100,
    sleep_seconds: float = 0.35,
    contact_email: str | None = None,
) -> dict[str, str | None]:
    """Batch-fetches BioSample XML records from NCBI's E-utilities
    (confirmed live: a comma-separated accession list in one efetch call
    correctly returns multiple full BioSample records, no need to
    pre-resolve accessions to internal UIDs via esearch first). Returns
    {accession: depth_raw or None} for every accession that had a
    resolvable BioSample record -- an accession simply absent from the
    result means NCBI had no record for it at all (distinct from having a
    record with no depth attribute, which maps to None)."""
    out: dict[str, str | None] = {}
    wanted = [
        a for a in sorted(set(biosample_accessions))
        if a and a.strip().lower() not in UNRESOLVABLE_BIOSAMPLE_MARKERS
    ]
    contact = contact_email or os.environ.get("PHA_REFERENCE_CONTACT_EMAIL", DEFAULT_CONTACT)
    for batch in _chunks(wanted, batch_size):
        params = {
            "db": "biosample",
            "id": ",".join(batch),
            "rettype": "full",
            "retmode": "xml",
            "tool": "PHA-Ocean-Atlas-pha-reference",
            "email": contact,
        }
        resp = _get(client, f"{NCBI_EUTILS_BASE}/efetch.fcgi", params)
        root = ET.fromstring(resp.text)
        for biosample_el in root.findall("BioSample"):
            accession = biosample_el.get("accession")
            if accession:
                out[accession] = _extract_depth_raw(biosample_el)
        if sleep_seconds:
            time.sleep(sleep_seconds)
    return out


def enrich_with_depth(
    metadata_path: Path,
    out_path: Path,
    batch_size: int = 100,
    sleep_seconds: float = 0.35,
) -> dict:
    """Reads a <family>_unique_targets_with_metadata.tsv (from
    pipeline/omdb_metadata.py), adds depth_raw/depth_m/depth_zone columns
    keyed off its existing "biosample" column, and writes out_path.
    Returns a small dict of coverage counts for reporting."""
    with open(metadata_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    biosample_accessions = [r.get("biosample", "") for r in rows]
    with _client() as client:
        depths_raw = fetch_sample_depths(client, biosample_accessions, batch_size=batch_size, sleep_seconds=sleep_seconds)

    out_fieldnames = fieldnames + ["depth_raw", "depth_m", "depth_zone"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_with_depth_value = 0
    with open(out_path, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=out_fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            accession = row.get("biosample", "")
            raw = depths_raw.get(accession)
            depth_m = parse_depth_meters(raw)
            if depth_m is not None:
                n_with_depth_value += 1
            writer.writerow({**row, "depth_raw": raw or "", "depth_m": depth_m if depth_m is not None else "", "depth_zone": depth_zone(depth_m)})

    return {
        "n_input_rows": len(rows),
        "n_distinct_biosamples": len({a for a in biosample_accessions if a and a.strip().lower() not in UNRESOLVABLE_BIOSAMPLE_MARKERS}),
        "n_biosamples_found_in_ncbi": len(depths_raw),
        "n_rows_with_parsed_depth": n_with_depth_value,
    }
