"""Shared upsert helpers so ingest_brenda.py and ingest_uniprot.py never
duplicate a protein/family_assignment row -- both funnel through here.
"""

from __future__ import annotations

import datetime as dt

from sqlalchemy.orm import Session

from phaatlas.db.models import FamilyAssignment, Phenotype, Protein, SourceEvidence
from phaatlas.pipeline import policy
from phaatlas.sources.uniprot import ParsedUniProtEntry


def _pipe_union(existing: str | None, new_items: list[str]) -> str | None:
    items = [x for x in (existing or "").split("|") if x]
    for it in new_items:
        if it and it not in items:
            items.append(it)
    return "|".join(items) if items else None


def protein_id_for_accession(accession: str) -> str:
    return f"UNIPROT:{accession}"


def protein_id_for_brenda_local(ec_number: str, local_id: str) -> str:
    return f"BRENDA:{ec_number}:{local_id}"


def upsert_protein_from_uniprot(
    session: Session,
    entry: ParsedUniProtEntry,
    uniprot_release: str | None,
) -> Protein:
    protein_id = protein_id_for_accession(entry.accession)
    protein = session.get(Protein, protein_id)
    now = dt.datetime.now(dt.timezone.utc)
    if protein is None:
        protein = Protein(protein_id=protein_id, retrieved_at=now)
        session.add(protein)

    protein.uniprot_accession = entry.accession
    protein.secondary_accessions = "|".join(entry.secondary_accessions) or None
    protein.reviewed = entry.reviewed
    protein.annotation_score = entry.annotation_score
    protein.protein_existence = entry.protein_existence
    protein.protein_name = entry.protein_name
    protein.alternative_protein_names = "|".join(entry.alternative_protein_names) or None
    protein.gene_name = entry.gene_name
    protein.gene_synonyms = "|".join(entry.gene_synonyms) or None
    protein.organism = entry.organism
    protein.strain = entry.strain
    protein.ncbi_taxid = entry.ncbi_taxid
    protein.taxonomic_lineage = "|".join(entry.taxonomic_lineage) or None
    protein.sequence = entry.sequence
    protein.sequence_length = entry.sequence_length
    protein.sequence_sha256 = entry.sequence_sha256
    protein.ec_numbers = _pipe_union(protein.ec_numbers, entry.ec_numbers)
    protein.rhea_ids = _pipe_union(protein.rhea_ids, entry.rhea_ids)
    protein.function_text = entry.function_text or protein.function_text
    protein.catalytic_activity_text = entry.catalytic_activity_text or protein.catalytic_activity_text
    protein.pathway_text = entry.pathway_text or protein.pathway_text
    protein.active_sites_text = entry.active_sites_text or protein.active_sites_text
    protein.uniprot_release = uniprot_release or protein.uniprot_release
    protein.source_releases = _pipe_union(protein.source_releases, [f"UniProt:{uniprot_release}"] if uniprot_release else [])
    protein.updated_at = now

    session.flush()
    return protein


def get_or_create_brenda_local_protein(
    session: Session,
    ec_number: str,
    local_id: str,
    organism: str | None,
    brenda_release: str | None,
) -> Protein:
    """Used only when a BRENDA protein record has NO accession at all --
    still stored as its own protein row per spec ("do not deduplicate
    different accessions yet by deleting records"), just without a
    sequence.
    """
    protein_id = protein_id_for_brenda_local(ec_number, local_id)
    protein = session.get(Protein, protein_id)
    now = dt.datetime.now(dt.timezone.utc)
    if protein is None:
        protein = Protein(protein_id=protein_id, retrieved_at=now, organism=organism)
        session.add(protein)
    protein.organism = organism or protein.organism
    session.flush()
    attach_brenda_context_to_protein(session, protein.protein_id, ec_number, local_id, brenda_release, None)
    return protein


def attach_brenda_context_to_protein(
    session: Session,
    protein_id: str,
    ec_number: str,
    local_id: str,
    brenda_release: str | None,
    brenda_protein_source: str | None,
) -> None:
    """Stamps brenda_ec_numbers/brenda_record_ids/brenda_sources/source_releases
    onto a protein row -- used for BOTH the accession-linked (Route A) and
    no-accession BRENDA protein paths, so a protein reached via BRENDA
    always carries this provenance regardless of which path found it.
    """
    protein = session.get(Protein, protein_id)
    protein.brenda_ec_numbers = _pipe_union(protein.brenda_ec_numbers, [ec_number])
    protein.brenda_record_ids = _pipe_union(protein.brenda_record_ids, [f"{ec_number}:{local_id}"])
    if brenda_protein_source:
        protein.brenda_sources = _pipe_union(protein.brenda_sources, [brenda_protein_source])
    protein.source_releases = _pipe_union(protein.source_releases, [f"BRENDA:{brenda_release}"] if brenda_release else [])
    protein.updated_at = dt.datetime.now(dt.timezone.utc)
    session.flush()


def upsert_family_assignment(
    session: Session,
    protein_id: str,
    family_id: str,
    reason: str,
    tier: str,
    needs_manual_review: bool,
    pha_subfamily: str | None = None,
    confidence_notes: str | None = None,
) -> FamilyAssignment:
    existing = (
        session.query(FamilyAssignment)
        .filter_by(protein_id=protein_id, family_id=family_id)
        .one_or_none()
    )
    now = dt.datetime.now(dt.timezone.utc)
    if existing is None:
        existing = FamilyAssignment(
            protein_id=protein_id,
            family_id=family_id,
            pha_subfamily=pha_subfamily,
            retrieval_reasons=reason,
            evidence_tier=tier,
            confidence_notes=confidence_notes,
            needs_manual_review=needs_manual_review,
            assigned_at=now,
        )
        session.add(existing)
    else:
        existing.retrieval_reasons = _pipe_union(existing.retrieval_reasons, [reason])
        existing.evidence_tier = policy.better_tier(existing.evidence_tier, tier)
        # A protein confirmed GOLD via one route stays GOLD even if a later,
        # separate query also happens to match it through an ambiguous
        # route -- manual review is only for cases where NO route has yet
        # confirmed it.
        if existing.evidence_tier != policy.CANDIDATE_AMBIGUOUS:
            existing.needs_manual_review = False
        if pha_subfamily and not existing.pha_subfamily:
            existing.pha_subfamily = pha_subfamily
        if confidence_notes:
            existing.confidence_notes = _pipe_union(existing.confidence_notes, [confidence_notes])
        existing.updated_at = now
    session.flush()
    return existing


def add_source_evidence(
    session: Session,
    protein_id: str,
    source_name: str,
    source_release: str | None,
    source_record_id: str | None,
    ec_number: str | None,
    organism_text: str | None,
    literature_pmids: list[str],
    literature_dois: list[str],
    raw_payload: dict | None,
    notes: str | None = None,
) -> SourceEvidence:
    ev = SourceEvidence(
        protein_id=protein_id,
        source_name=source_name,
        source_release=source_release,
        source_record_id=source_record_id,
        ec_number=ec_number,
        organism_text=organism_text,
        literature_pmids="|".join(literature_pmids) or None,
        literature_dois="|".join(literature_dois) or None,
        raw_payload_json=_safe_json(raw_payload),
        notes=notes,
    )
    session.add(ev)
    session.flush()
    return ev


def add_phenotype(
    session: Session,
    protein_id: str,
    source_evidence_id: int | None,
    phenotype_type: str,
    source_field: str,
    normalized,  # NormalizedPhenotype
    literature_pmids: list[str] | None = None,
    literature_dois: list[str] | None = None,
    temperature_optimum: str | None = None,
    temperature_range: str | None = None,
    temperature_stability: str | None = None,
    ph_optimum: str | None = None,
    ph_range: str | None = None,
    specific_activity: str | None = None,
) -> Phenotype:
    ph = Phenotype(
        protein_id=protein_id,
        source_evidence_id=source_evidence_id,
        phenotype_type=phenotype_type,
        source_field=source_field,
        phenotype_raw=normalized.phenotype_raw,
        demonstrated_monomers=normalized.demonstrated_monomers,
        demonstrated_polymers=normalized.demonstrated_polymers,
        composition_text=normalized.composition_text,
        phenotype_evidence_level=normalized.phenotype_evidence_level,
        temperature_optimum=temperature_optimum,
        temperature_range=temperature_range,
        temperature_stability=temperature_stability,
        ph_optimum=ph_optimum,
        ph_range=ph_range,
        specific_activity=specific_activity,
        literature_pmids="|".join(literature_pmids) if literature_pmids else None,
        literature_dois="|".join(literature_dois) if literature_dois else None,
    )
    session.add(ph)
    session.flush()
    return ph


def _safe_json(payload: dict | None) -> str | None:
    if payload is None:
        return None
    import json

    return json.dumps(payload, default=str)
