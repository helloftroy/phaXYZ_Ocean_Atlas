"""BRENDA full JSON download + per-EC-record streaming extraction.

Download: automates the exact same license-acceptance form submit a browser
does at https://www.brenda-enzymes.org/download.php (checkbox
"accept-license" + the "BRENDA JSON file" button) -- this is BRENDA's own
documented bulk-download route for the full dataset, not page-scraping of
individual enzyme entries. No login/account is required for this download.

Extraction: the decompressed JSON is a single ~700MB object
(`{"release": ..., "version": ..., "data": {<EC number>: {...}, ...}}`), too
large to comfortably `json.load()` into a plain dict repeatedly during
development. `iter_ec_records` streams through `data`'s key/value pairs with
ijson so memory stays bounded to one EC record at a time, and only the EC
numbers actually requested (a handful, from family_definitions.yaml) get
returned -- everything else is parsed-and-discarded in one pass, matching
what a real single-key lookup costs on a file this shape.
"""

from __future__ import annotations

import gzip
import io
import os
import shutil
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import httpx
import ijson

from phaatlas.config_loader import REPO_ROOT

BRENDA_DOWNLOAD_URL = "https://www.brenda-enzymes.org/download.php"
DEFAULT_CACHE_DIR = REPO_ROOT / "data" / "brenda"
DEFAULT_CONTACT = "REPLACE_WITH_CONTACT_EMAIL@example.org"


def _user_agent() -> str:
    contact = os.environ.get("PHA_REFERENCE_CONTACT_EMAIL", DEFAULT_CONTACT)
    return f"PHA-Ocean-Atlas-pha-reference/0.1 (research pipeline; contact: {contact})"


def download_json(cache_dir: Path = DEFAULT_CACHE_DIR, force: bool = False) -> Path:
    """Downloads and extracts the full BRENDA JSON dump, accepting the CC BY
    4.0 license exactly as the website's own download form does. Returns the
    path to the extracted brenda_<release>.json file. Cached: re-running
    without force=True reuses an already-extracted file.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    existing = sorted(cache_dir.glob("brenda_*.json"))
    if existing and not force:
        return existing[-1]

    with httpx.Client(headers={"User-Agent": _user_agent()}, timeout=180.0) as client:
        resp = client.post(
            BRENDA_DOWNLOAD_URL,
            data={"dlfile": "dl-json", "accept-license": "1"},
            follow_redirects=True,
        )
        resp.raise_for_status()
        content_disposition = resp.headers.get("content-disposition", "")
        if "json" not in content_disposition.lower() and resp.headers.get("content-type") != "application/x-gzip":
            raise RuntimeError(
                "BRENDA download did not return the expected JSON tarball "
                f"(content-disposition={content_disposition!r}, "
                f"content-type={resp.headers.get('content-type')!r}) -- "
                "the license-accept form on brenda-enzymes.org/download.php may have changed."
            )
        tar_bytes = resp.content

    with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".json")]
        if not members:
            raise RuntimeError("BRENDA download tarball contained no .json member")
        member = members[0]
        extracted_path = cache_dir / Path(member.name).name
        with tf.extractfile(member) as src, open(extracted_path, "wb") as dst:
            shutil.copyfileobj(src, dst)

    return extracted_path


@dataclass
class BrendaProteinRef:
    local_id: str
    organism: str | None
    source: str | None  # "swissprot" | "trembl" | "genbank" | None
    accessions: list[str]
    reference_ids: list[str]
    comment: str | None


@dataclass
class BrendaReference:
    local_id: str
    pmid: str | None
    doi: str | None
    title: str | None
    authors: list[str]
    journal: str | None
    year: int | None


@dataclass
class BrendaEntry:
    local_id: str
    references: list[str]
    comment: str | None
    value: str | None


@dataclass
class BrendaECRecord:
    ec_number: str
    recommended_name: str | None
    systematic_name: str | None
    synonyms: list[str]
    proteins: dict[str, BrendaProteinRef]
    references: dict[str, BrendaReference]
    natural_substrates_products: list[dict]
    substrates_products: list[dict]
    general_information: list[dict]
    protein_variants: list[dict]
    temperature_optimum: list[dict]
    temperature_range: list[dict]
    temperature_stability: list[dict]
    ph_optimum: list[dict]
    ph_range: list[dict]
    specific_activity: list[dict]
    km_value: list[dict]
    kcat_km_value: list[dict]
    turnover_number: list[dict]
    raw: dict = field(repr=False)


def _parse_proteins(raw_protein_obj: dict | None) -> dict[str, BrendaProteinRef]:
    out = {}
    for local_id, p in (raw_protein_obj or {}).items():
        out[local_id] = BrendaProteinRef(
            local_id=local_id,
            organism=p.get("organism"),
            source=p.get("source"),
            accessions=list(p.get("accessions") or []),
            reference_ids=list(p.get("references") or []),
            comment=p.get("comment") or None,
        )
    return out


def _parse_references(raw_reference_obj: dict | None) -> dict[str, BrendaReference]:
    out = {}
    for local_id, r in (raw_reference_obj or {}).items():
        out[local_id] = BrendaReference(
            local_id=local_id,
            pmid=str(r["pmid"]) if r.get("pmid") else None,
            doi=r.get("doi"),
            title=r.get("title"),
            authors=list(r.get("authors") or []),
            journal=r.get("journal"),
            year=r.get("year"),
        )
    return out


def parse_ec_record(ec_number: str, raw: dict) -> BrendaECRecord:
    return BrendaECRecord(
        ec_number=ec_number,
        recommended_name=raw.get("recommended_name"),
        systematic_name=raw.get("systematic_name"),
        synonyms=[s.get("value") for s in raw.get("synonyms") or [] if s.get("value")],
        proteins=_parse_proteins(raw.get("protein")),
        references=_parse_references(raw.get("reference")),
        natural_substrates_products=list(raw.get("natural_substrates_products") or []),
        substrates_products=list(raw.get("substrates_products") or []),
        general_information=list(raw.get("general_information") or []),
        protein_variants=list(raw.get("protein_variants") or []),
        temperature_optimum=list(raw.get("temperature_optimum") or []),
        temperature_range=list(raw.get("temperature_range") or []),
        temperature_stability=list(raw.get("temperature_stability") or []),
        ph_optimum=list(raw.get("ph_optimum") or []),
        ph_range=list(raw.get("ph_range") or []),
        specific_activity=list(raw.get("specific_activity") or []),
        km_value=list(raw.get("km_value") or []),
        kcat_km_value=list(raw.get("kcat_km_value") or []),
        turnover_number=list(raw.get("turnover_number") or []),
        raw=raw,
    )


def get_release(json_path: Path) -> str | None:
    """Reads just the top-level `release` field without loading the whole
    file -- ijson streams the first few tokens and stops."""
    with open(json_path, "rb") as f:
        for prefix, event, value in ijson.parse(f):
            if prefix == "release" and event == "string":
                return value
            if prefix.startswith("data"):
                break
    return None


def iter_ec_records(json_path: Path, ec_numbers: set[str]) -> Iterator[BrendaECRecord]:
    """Streams through the top-level `data` object, yielding a parsed
    BrendaECRecord for each requested EC number found. Single pass,
    memory-bounded to one EC record at a time regardless of how large the
    full file is.
    """
    remaining = set(ec_numbers)
    if not remaining:
        return
    with open(json_path, "rb") as f:
        for key, value in ijson.kvitems(f, "data"):
            if key in remaining:
                yield parse_ec_record(key, value)
                remaining.discard(key)
                if not remaining:
                    return
