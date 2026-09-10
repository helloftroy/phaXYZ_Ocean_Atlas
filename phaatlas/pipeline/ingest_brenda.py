"""BRENDA ingestion: EC-anchored starting point for phaA/phaB/phaC/phaJ
(the families with a clean or config-declared ec_brenda), family-config
driven and gated by require_pha_context exactly as specified -- e.g. an
EC 2.3.1.9 protein is only assigned to phaA if the record's own text
(protein comment, EC-level synonyms, or any substrate/product/general-info
entry that references this specific protein) actually mentions PHA/
polyhydroxyalkanoate context, never from the EC number alone.

Families with ec_brenda: null (phasins, regulators, synthase partners) are
skipped here entirely, per spec -- ingest_uniprot.py is their only route.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

from sqlalchemy.orm import Session

from phaatlas.config_loader import FamilyDefinition
from phaatlas.db.models import RetrievalRun
from phaatlas.phenotype.normalize import EXPERIMENTAL_LITERATURE, normalize_phenotype_text
from phaatlas.pipeline import policy
from phaatlas.pipeline.upsert import (
    add_phenotype,
    add_source_evidence,
    attach_brenda_context_to_protein,
    get_or_create_brenda_local_protein,
    upsert_family_assignment,
)
from phaatlas.pipeline.ingest_uniprot import record_uniprot_hit
from phaatlas.sources import brenda as brenda_source
from phaatlas.sources import uniprot as uniprot_source

SOFTWARE_VERSION = "phaatlas/0.1"


def _entries_for_protein(entries: list[dict], local_id: str) -> list[dict]:
    return [e for e in entries if local_id in (e.get("proteins") or [])]


def _entry_text(entry: dict) -> str:
    value = entry.get("value") or ""
    comment = entry.get("comment") or ""
    if value and value != "more = ?":
        return f"{value} ({comment})" if comment else value
    return comment or value


def _protein_has_pha_context(record: brenda_source.BrendaECRecord, local_id: str) -> tuple[bool, str | None]:
    protein = record.proteins[local_id]
    candidates = [protein.comment, record.recommended_name, " ".join(record.synonyms)]
    for field_entries in (
        record.natural_substrates_products,
        record.substrates_products,
        record.general_information,
        record.protein_variants,
    ):
        for e in _entries_for_protein(field_entries, local_id):
            candidates.append(_entry_text(e))
    blob_parts = [c for c in candidates if c]
    if policy.text_has_pha_context(*blob_parts):
        matched = next((c for c in blob_parts if policy.text_has_pha_context(c)), None)
        return True, matched
    return False, None


def _reference_pmids_dois(record: brenda_source.BrendaECRecord, reference_ids: list[str]) -> tuple[list[str], list[str]]:
    pmids, dois = [], []
    for rid in reference_ids:
        ref = record.references.get(rid)
        if ref is None:
            continue
        if ref.pmid:
            pmids.append(ref.pmid)
        if ref.doi:
            dois.append(ref.doi)
    return pmids, dois


def _add_phenotypes_for_protein(
    session: Session,
    protein_id: str,
    source_evidence_id: int | None,
    record: brenda_source.BrendaECRecord,
    local_id: str,
) -> None:
    field_map = [
        (record.natural_substrates_products, "natural_substrates_products", "substrate_specificity"),
        (record.substrates_products, "substrates_products", "substrate_specificity"),
        (record.general_information, "general_information", "general_information"),
        (record.protein_variants, "protein_variants", "variant_phenotype"),
    ]
    for entries, field_name, phenotype_type in field_map:
        for e in _entries_for_protein(entries, local_id):
            text = _entry_text(e)
            if not text:
                continue
            normalized = normalize_phenotype_text(text, EXPERIMENTAL_LITERATURE)
            pmids, dois = _reference_pmids_dois(record, e.get("references") or [])
            add_phenotype(
                session,
                protein_id=protein_id,
                source_evidence_id=source_evidence_id,
                phenotype_type=phenotype_type,
                source_field=field_name,
                normalized=normalized,
                literature_pmids=pmids,
                literature_dois=dois,
            )

    condition_field_map = [
        (record.temperature_optimum, "temperature_optimum"),
        (record.temperature_range, "temperature_range"),
        (record.temperature_stability, "temperature_stability"),
        (record.ph_optimum, "ph_optimum"),
        (record.ph_range, "ph_range"),
        (record.specific_activity, "specific_activity"),
    ]
    for entries, column_name in condition_field_map:
        for e in _entries_for_protein(entries, local_id):
            text = _entry_text(e)
            if not text:
                continue
            normalized = normalize_phenotype_text(text, EXPERIMENTAL_LITERATURE)
            pmids, dois = _reference_pmids_dois(record, e.get("references") or [])
            kwargs = {column_name: text}
            add_phenotype(
                session,
                protein_id=protein_id,
                source_evidence_id=source_evidence_id,
                phenotype_type="reaction_conditions",
                source_field=column_name,
                normalized=normalized,
                literature_pmids=pmids,
                literature_dois=dois,
                **kwargs,
            )


def ingest_brenda_for_family(
    session: Session,
    family: FamilyDefinition,
    json_path: Path,
    brenda_release: str | None,
) -> dict:
    summary = {"family_id": family.family_id, "ec_number": family.ec_brenda, "proteins_total": 0, "proteins_assigned": 0, "accessions_found": 0}
    if not family.ec_brenda:
        return summary

    run = RetrievalRun(
        database="BRENDA",
        database_release=brenda_release,
        family_id=family.family_id,
        query=family.ec_brenda,
        software_version=SOFTWARE_VERSION,
        status="running",
    )
    session.add(run)
    session.flush()

    try:
        records = list(brenda_source.iter_ec_records(json_path, {family.ec_brenda}))
    except Exception as exc:
        run.status = "failed"
        run.error_text = str(exc)
        run.finished_at = dt.datetime.now(dt.timezone.utc)
        session.commit()
        raise

    if not records:
        run.status = "completed"
        run.number_returned = 0
        run.finished_at = dt.datetime.now(dt.timezone.utc)
        session.commit()
        return summary

    record = records[0]
    summary["proteins_total"] = len(record.proteins)
    uniprot_release = uniprot_source.release_info()

    for local_id, brenda_protein in record.proteins.items():
        if family.require_pha_context:
            confirmed, _matched_text = _protein_has_pha_context(record, local_id)
            if not confirmed:
                continue
        summary["proteins_assigned"] += 1

        if brenda_protein.accessions:
            for accession in brenda_protein.accessions:
                entry = uniprot_source.fetch_by_accession(accession)
                if entry is None:
                    continue
                record_uniprot_hit(session, entry, family, "brenda_accession", uniprot_release)
                summary["accessions_found"] += 1
                protein_id = f"UNIPROT:{accession}"
                attach_brenda_context_to_protein(
                    session, protein_id, family.ec_brenda, local_id, brenda_release, brenda_protein.source
                )

                pmids, dois = _reference_pmids_dois(record, brenda_protein.reference_ids)
                evidence = add_source_evidence(
                    session,
                    protein_id=protein_id,
                    source_name="BRENDA",
                    source_release=brenda_release,
                    source_record_id=f"{family.ec_brenda}:{local_id}",
                    ec_number=family.ec_brenda,
                    organism_text=brenda_protein.organism,
                    literature_pmids=pmids,
                    literature_dois=dois,
                    raw_payload={"protein": vars(brenda_protein), "recommended_name": record.recommended_name},
                    notes=f"brenda_source={brenda_protein.source}",
                )
                _add_phenotypes_for_protein(session, protein_id, evidence.evidence_id, record, local_id)
        else:
            protein = get_or_create_brenda_local_protein(
                session, family.ec_brenda, local_id, brenda_protein.organism, brenda_release
            )
            tier, needs_review = policy.compute_tier_for_brenda(has_accession=False)
            upsert_family_assignment(
                session,
                protein_id=protein.protein_id,
                family_id=family.family_id,
                reason="brenda_no_accession",
                tier=tier,
                needs_manual_review=needs_review,
            )
            pmids, dois = _reference_pmids_dois(record, brenda_protein.reference_ids)
            evidence = add_source_evidence(
                session,
                protein_id=protein.protein_id,
                source_name="BRENDA",
                source_release=brenda_release,
                source_record_id=f"{family.ec_brenda}:{local_id}",
                ec_number=family.ec_brenda,
                organism_text=brenda_protein.organism,
                literature_pmids=pmids,
                literature_dois=dois,
                raw_payload={"protein": vars(brenda_protein), "recommended_name": record.recommended_name},
                notes=f"brenda_source={brenda_protein.source}",
            )
            _add_phenotypes_for_protein(session, protein.protein_id, evidence.evidence_id, record, local_id)

    run.status = "completed"
    run.number_returned = summary["proteins_assigned"]
    run.finished_at = dt.datetime.now(dt.timezone.utc)
    session.commit()
    return summary
