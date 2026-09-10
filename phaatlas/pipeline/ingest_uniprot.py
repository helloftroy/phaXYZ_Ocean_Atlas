"""UniProt ingestion: Route A (exact BRENDA-linked accessions) and Route B
(family-specific reviewed/unreviewed searches), family-config driven.

phaR is handled specially: phaR_synthase and phaR_regulator are the one
deliberate pair of families sharing an exact gene name (see
family_definitions.yaml), so their `gene:phaR` search results are resolved
jointly against both families' disambiguation_terms rather than run as two
independent per-family gene searches.
"""

from __future__ import annotations

import datetime as dt

from sqlalchemy.orm import Session

from phaatlas.config_loader import FamilyDefinition
from phaatlas.db.models import RetrievalRun
from phaatlas.pipeline import policy
from phaatlas.pipeline.upsert import (
    add_phenotype,
    add_source_evidence,
    protein_id_for_accession,
    upsert_family_assignment,
    upsert_protein_from_uniprot,
)
from phaatlas.phenotype.normalize import (
    GENERIC_REACTION_ANNOTATION,
    normalize_phenotype_text,
    uniprot_function_evidence_level,
)
from phaatlas.sources import uniprot as uniprot_source

SOFTWARE_VERSION = "phaatlas/0.1"


def record_uniprot_hit(
    session: Session,
    entry: uniprot_source.ParsedUniProtEntry,
    family: FamilyDefinition,
    reason: str,
    uniprot_release: str | None,
    pha_subfamily: str | None = None,
) -> None:
    if entry.gene_name and entry.gene_name in family.exclude_genes:
        return
    if any(g in family.exclude_genes for g in entry.gene_synonyms):
        return

    protein = upsert_protein_from_uniprot(session, entry, uniprot_release)

    func_level = uniprot_function_evidence_level(entry.function_text or "")
    specific_function_evidence = func_level != GENERIC_REACTION_ANNOTATION
    if reason == "brenda_accession":
        tier, needs_review = policy.compute_tier_for_brenda(has_accession=True)
    else:
        tier, needs_review = policy.compute_tier_for_uniprot(reason, entry.reviewed, specific_function_evidence)

    upsert_family_assignment(
        session,
        protein_id=protein.protein_id,
        family_id=family.family_id,
        reason=reason,
        tier=tier,
        needs_manual_review=needs_review,
        pha_subfamily=pha_subfamily,
    )

    evidence = add_source_evidence(
        session,
        protein_id=protein.protein_id,
        source_name="UniProt",
        source_release=uniprot_release,
        source_record_id=entry.accession,
        ec_number=(entry.ec_numbers[0] if entry.ec_numbers else None),
        organism_text=entry.organism,
        literature_pmids=entry.literature_pmids,
        literature_dois=entry.literature_dois,
        raw_payload=entry.raw,
        notes=f"retrieval_reason={reason}",
    )

    if entry.function_text:
        normalized = normalize_phenotype_text(entry.function_text, func_level)
        add_phenotype(
            session,
            protein_id=protein.protein_id,
            source_evidence_id=evidence.evidence_id,
            phenotype_type="function_annotation",
            source_field="uniprot_function",
            normalized=normalized,
            literature_pmids=entry.literature_pmids,
        )

    if entry.catalytic_activity_text:
        normalized = normalize_phenotype_text(entry.catalytic_activity_text, GENERIC_REACTION_ANNOTATION)
        add_phenotype(
            session,
            protein_id=protein.protein_id,
            source_evidence_id=evidence.evidence_id,
            phenotype_type="substrate_specificity",
            source_field="uniprot_catalytic_activity",
            normalized=normalized,
        )


def _run_query(
    session: Session,
    query: str,
    family_id: str,
    max_results: int,
) -> list[uniprot_source.ParsedUniProtEntry]:
    run = RetrievalRun(
        database="UniProt",
        family_id=family_id,
        query=query,
        software_version=SOFTWARE_VERSION,
        status="running",
    )
    session.add(run)
    session.flush()
    try:
        entries = uniprot_source.search(query, max_results=max_results)
        run.number_returned = len(entries)
        run.status = "completed"
    except Exception as exc:
        run.status = "failed"
        run.error_text = str(exc)
        run.finished_at = dt.datetime.now(dt.timezone.utc)
        session.flush()
        raise
    run.finished_at = dt.datetime.now(dt.timezone.utc)
    session.flush()
    return entries


def ingest_uniprot_for_family(
    session: Session,
    family: FamilyDefinition,
    max_results_per_query: int = 1000,
    skip_gene_queries: bool = False,
) -> dict:
    uniprot_release = uniprot_source.release_info()
    summary = {"family_id": family.family_id, "queries": 0, "hits": 0}

    queries: list[tuple[str, str]] = []
    if not skip_gene_queries:
        for gene in family.genes:
            queries.append((f'gene:{gene} AND reviewed:true', "uniprot_reviewed_gene_name"))
            queries.append((f'gene:{gene} AND reviewed:false', "uniprot_unreviewed_gene_name"))
    for term in family.protein_name_terms:
        queries.append((f'protein_name:"{term}" AND reviewed:true', "uniprot_reviewed_protein_name"))
        queries.append((f'protein_name:"{term}" AND reviewed:false', "uniprot_unreviewed_protein_name"))
    if family.broad_uniprot_search and family.ec_brenda:
        queries.append((f'ec:{family.ec_brenda} AND reviewed:true', "uniprot_reviewed_ec"))
        queries.append((f'ec:{family.ec_brenda} AND reviewed:false', "uniprot_unreviewed_ec"))

    for query, reason in queries:
        entries = _run_query(session, query, family.family_id, max_results_per_query)
        summary["queries"] += 1
        for entry in entries:
            record_uniprot_hit(session, entry, family, reason, uniprot_release)
            summary["hits"] += 1

    session.commit()
    return summary


def ingest_phaR_disambiguation(
    session: Session,
    phaR_synthase: FamilyDefinition,
    phaR_regulator: FamilyDefinition,
    max_results_per_query: int = 1000,
) -> dict:
    """The one deliberate shared-gene-name case: runs `gene:phaR` once (per
    reviewed/unreviewed tier) and classifies each hit against BOTH families'
    disambiguation_terms rather than two independent per-family gene
    searches -- see spec's "must be separate internal families" requirement.
    """
    uniprot_release = uniprot_source.release_info()
    summary = {"queries": 0, "hits_synthase": 0, "hits_regulator": 0, "hits_ambiguous": 0}

    for reviewed_flag, base_reason in ((True, "uniprot_reviewed_gene_name"), (False, "uniprot_unreviewed_gene_name")):
        query = f'gene:phaR AND reviewed:{"true" if reviewed_flag else "false"}'
        entries = _run_query(session, query, "phaR_synthase|phaR_regulator", max_results_per_query)
        summary["queries"] += 1

        for entry in entries:
            blob = " ".join(filter(None, [entry.protein_name, entry.function_text, entry.pathway_text])).lower()
            is_synthase = any(term.lower() in blob for term in phaR_synthase.disambiguation_terms)
            is_regulator = any(term.lower() in blob for term in phaR_regulator.disambiguation_terms)

            if is_synthase and not is_regulator:
                record_uniprot_hit(session, entry, phaR_synthase, base_reason, uniprot_release)
                summary["hits_synthase"] += 1
            elif is_regulator and not is_synthase:
                record_uniprot_hit(session, entry, phaR_regulator, base_reason, uniprot_release)
                summary["hits_regulator"] += 1
            else:
                # Neither or both disambiguation_terms matched -- record
                # under both families, flagged for manual review, rather
                # than guessing which role this phaR protein plays.
                note = "phaR gene hit could not be disambiguated between synthase/regulator roles from available text"
                for fam in (phaR_synthase, phaR_regulator):
                    protein = upsert_protein_from_uniprot(session, entry, uniprot_release)
                    upsert_family_assignment(
                        session,
                        protein_id=protein.protein_id,
                        family_id=fam.family_id,
                        reason=base_reason,
                        tier=policy.CANDIDATE_AMBIGUOUS,
                        needs_manual_review=True,
                        confidence_notes=note,
                    )
                    add_source_evidence(
                        session,
                        protein_id=protein.protein_id,
                        source_name="UniProt",
                        source_release=uniprot_release,
                        source_record_id=entry.accession,
                        ec_number=None,
                        organism_text=entry.organism,
                        literature_pmids=entry.literature_pmids,
                        literature_dois=entry.literature_dois,
                        raw_payload=entry.raw,
                        notes=note,
                    )
                summary["hits_ambiguous"] += 1

    session.commit()
    return summary


def ingest_uniprot_accessions(
    session: Session,
    accessions: list[str],
    family: FamilyDefinition,
    pha_subfamily: str | None = None,
) -> dict:
    """Route A: fetch the exact UniProt entry for every accession found via
    BRENDA (the highest-value linkage -- experimental BRENDA phenotype tied
    to an exact protein sequence)."""
    uniprot_release = uniprot_source.release_info()
    summary = {"family_id": family.family_id, "requested": len(accessions), "found": 0, "missing": 0}
    for accession in sorted(set(accessions)):
        entry = uniprot_source.fetch_by_accession(accession)
        if entry is None:
            summary["missing"] += 1
            continue
        record_uniprot_hit(session, entry, family, "brenda_accession", uniprot_release, pha_subfamily=pha_subfamily)
        summary["found"] += 1
    session.commit()
    return summary
