"""UniProt REST client: exact-accession fetch (Route A) and family-specific
search (Route B, tiers 1/2). Returns full raw entry JSON plus a parsed
ParsedUniProtEntry -- both kept, per spec (`raw_payload_json` on
source_evidence retains everything; `protein`/`phenotype` only get the
subset the pipeline currently uses).

No API key needed for UniProt's REST API; a descriptive User-Agent
(PHA_REFERENCE_CONTACT_EMAIL) is sent as a politeness courtesy, same
convention as ../fair_ocean_agent's RateLimitedClient.
"""

from __future__ import annotations

import hashlib
import os
import time
from dataclasses import dataclass, field

import httpx
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
DEFAULT_CONTACT = "REPLACE_WITH_CONTACT_EMAIL@example.org"


def _user_agent() -> str:
    contact = os.environ.get("PHA_REFERENCE_CONTACT_EMAIL", DEFAULT_CONTACT)
    return f"PHA-Ocean-Atlas-pha-reference/0.1 (research pipeline; contact: {contact})"


def _client() -> httpx.Client:
    return httpx.Client(headers={"User-Agent": _user_agent()}, timeout=30.0)


@retry(
    reraise=True,
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((httpx.TransportError, httpx.HTTPStatusError)),
)
def _get(client: httpx.Client, url: str, params: dict | None = None) -> httpx.Response:
    resp = client.get(url, params=params)
    if resp.status_code == 429:
        wait_s = float(resp.headers.get("Retry-After", "5"))
        time.sleep(wait_s)
        resp.raise_for_status()
    resp.raise_for_status()
    return resp


@dataclass
class ParsedUniProtEntry:
    accession: str
    secondary_accessions: list[str]
    reviewed: bool
    annotation_score: str | None
    protein_existence: str | None
    protein_name: str | None
    alternative_protein_names: list[str]
    gene_name: str | None
    gene_synonyms: list[str]
    organism: str | None
    strain: str | None
    ncbi_taxid: str | None
    taxonomic_lineage: list[str]
    sequence: str | None
    sequence_length: int | None
    sequence_sha256: str | None
    ec_numbers: list[str]
    rhea_ids: list[str]
    function_text: str | None
    catalytic_activity_text: str | None
    pathway_text: str | None
    active_sites_text: str | None
    literature_pmids: list[str]
    literature_dois: list[str]
    raw: dict = field(repr=False)


def parse_uniprot_entry(entry: dict) -> ParsedUniProtEntry:
    accession = entry["primaryAccession"]
    secondary = list(entry.get("secondaryAccessions") or [])
    reviewed = str(entry.get("entryType", "")).startswith("UniProtKB reviewed")

    desc = entry.get("proteinDescription") or {}
    rec = desc.get("recommendedName") or {}
    protein_name = (rec.get("fullName") or {}).get("value")
    alt_names = []
    for alt in desc.get("alternativeNames") or []:
        v = (alt.get("fullName") or {}).get("value")
        if v:
            alt_names.append(v)
    ec_numbers = [e.get("value") for e in (rec.get("ecNumbers") or []) if e.get("value")]
    for alt in desc.get("alternativeNames") or []:
        ec_numbers.extend(e.get("value") for e in (alt.get("ecNumbers") or []) if e.get("value"))

    # Unreviewed (TrEMBL) entries almost never have a `recommendedName` --
    # UniProt puts their (uncurated, submitter-supplied) name under
    # `submissionNames` instead. Without this fallback, protein_name/
    # ec_numbers silently come back empty for nearly every ANNOTATED_UNREVIEWED
    # row (confirmed live: gene:phaR AND reviewed:false entries all had
    # proteinDescription == {"submissionNames": [...]} with no recommendedName
    # at all), which also starves the phaR synthase/regulator disambiguation
    # of the protein_name text it depends on.
    if not protein_name:
        for sub in desc.get("submissionNames") or []:
            v = (sub.get("fullName") or {}).get("value")
            if v:
                protein_name = v
                ec_numbers.extend(e.get("value") for e in (sub.get("ecNumbers") or []) if e.get("value"))
                break

    genes = entry.get("genes") or []
    gene_name = None
    gene_synonyms: list[str] = []
    if genes:
        gene_name = (genes[0].get("geneName") or {}).get("value")
        for g in genes:
            for syn in g.get("synonyms") or []:
                if syn.get("value"):
                    gene_synonyms.append(syn["value"])

    organism_obj = entry.get("organism") or {}
    organism = organism_obj.get("scientificName")
    ncbi_taxid = str(organism_obj.get("taxonId")) if organism_obj.get("taxonId") else None
    lineage = list(organism_obj.get("lineage") or [])
    strain = None
    for s in organism_obj.get("strains") or []:
        if s.get("name"):
            strain = s["name"]
            break

    seq_obj = entry.get("sequence") or {}
    sequence = seq_obj.get("value")
    sequence_length = seq_obj.get("length")
    sequence_sha256 = hashlib.sha256(sequence.encode("ascii")).hexdigest() if sequence else None

    function_text = None
    catalytic_activity_text = None
    pathway_text = None
    rhea_ids: list[str] = []
    for comment in entry.get("comments") or []:
        ctype = comment.get("commentType")
        if ctype == "FUNCTION":
            texts = [t.get("value") for t in comment.get("texts") or [] if t.get("value")]
            if texts:
                function_text = " ".join(texts) if function_text is None else function_text + " " + " ".join(texts)
        elif ctype == "CATALYTIC ACTIVITY":
            reaction = comment.get("reaction") or {}
            name = reaction.get("name")
            if name:
                catalytic_activity_text = name if catalytic_activity_text is None else catalytic_activity_text + "; " + name
            ec = reaction.get("ecNumber")
            if ec:
                ec_numbers.append(ec)
            for xref in reaction.get("reactionCrossReferences") or []:
                if xref.get("database") == "Rhea" and xref.get("id"):
                    rhea_ids.append(xref["id"])
        elif ctype == "PATHWAY":
            texts = [t.get("value") for t in comment.get("texts") or [] if t.get("value")]
            if texts:
                pathway_text = " ".join(texts) if pathway_text is None else pathway_text + " " + " ".join(texts)

    active_site_descs = []
    for feat in entry.get("features") or []:
        if feat.get("type") == "Active site":
            loc = feat.get("location") or {}
            pos = (loc.get("start") or {}).get("value")
            desc_text = feat.get("description") or ""
            active_site_descs.append(f"{pos}:{desc_text}".strip(":"))
    active_sites_text = "|".join(d for d in active_site_descs if d) or None

    pmids: list[str] = []
    dois: list[str] = []
    for ref in entry.get("references") or []:
        citation = ref.get("citation") or {}
        for xref in citation.get("citationCrossReferences") or []:
            if xref.get("database") == "PubMed" and xref.get("id"):
                pmids.append(str(xref["id"]))
            elif xref.get("database") == "DOI" and xref.get("id"):
                dois.append(str(xref["id"]))

    def _dedup(seq: list[str]) -> list[str]:
        seen = []
        for x in seq:
            if x not in seen:
                seen.append(x)
        return seen

    return ParsedUniProtEntry(
        accession=accession,
        secondary_accessions=secondary,
        reviewed=reviewed,
        annotation_score=str(entry.get("annotationScore")) if entry.get("annotationScore") is not None else None,
        protein_existence=entry.get("proteinExistence"),
        protein_name=protein_name,
        alternative_protein_names=_dedup(alt_names),
        gene_name=gene_name,
        gene_synonyms=_dedup(gene_synonyms),
        organism=organism,
        strain=strain,
        ncbi_taxid=ncbi_taxid,
        taxonomic_lineage=lineage,
        sequence=sequence,
        sequence_length=sequence_length,
        sequence_sha256=sequence_sha256,
        ec_numbers=_dedup(ec_numbers),
        rhea_ids=_dedup(rhea_ids),
        function_text=function_text,
        catalytic_activity_text=catalytic_activity_text,
        pathway_text=pathway_text,
        active_sites_text=active_sites_text,
        literature_pmids=_dedup(pmids),
        literature_dois=_dedup(dois),
        raw=entry,
    )


def fetch_by_accession(accession: str) -> ParsedUniProtEntry | None:
    with _client() as client:
        try:
            resp = _get(client, f"{UNIPROT_BASE}/{accession}.json")
        except httpx.HTTPStatusError as exc:
            if exc.response.status_code == 404:
                return None
            raise
        return parse_uniprot_entry(resp.json())


def search(
    query: str,
    max_results: int = 2000,
    page_size: int = 500,
) -> list[ParsedUniProtEntry]:
    """Runs a UniProt search query, paginating via the `next` cursor Link
    header, and returns full parsed entries (UniProt's search endpoint
    returns complete entry JSON per hit when no `fields=` filter is
    supplied, so no second per-accession fetch is needed).
    """
    results: list[ParsedUniProtEntry] = []
    url = f"{UNIPROT_BASE}/search"
    params: dict | None = {"query": query, "format": "json", "size": page_size}
    with _client() as client:
        while url and len(results) < max_results:
            resp = _get(client, url, params=params)
            data = resp.json()
            for entry in data.get("results", []):
                results.append(parse_uniprot_entry(entry))
                if len(results) >= max_results:
                    break
            next_url = resp.links.get("next", {}).get("url")
            url = next_url
            params = None  # cursor is already encoded in next_url
    return results


def release_info() -> str | None:
    """Best-effort UniProt release string, from the response headers of a
    trivial query -- UniProt doesn't expose this in the JSON body itself."""
    with _client() as client:
        try:
            resp = _get(client, f"{UNIPROT_BASE}/search", params={"query": "accession:P23608", "format": "json", "size": 1})
        except httpx.HTTPError:
            return None
        return resp.headers.get("x-uniprot-release") or resp.headers.get("X-UniProt-Release")
