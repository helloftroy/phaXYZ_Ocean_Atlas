"""Shared evidence-tier / retrieval-reason policy used by both
ingest_brenda.py and ingest_uniprot.py, so the rules live in exactly one
place.

Evidence tiers (spec):
  GOLD_EXPERIMENTAL   specific PHA role/function backed by literature/BRENDA,
                      or a highly specific experimental UniProt annotation.
  CURATED_REVIEWED    reviewed Swiss-Prot PHA annotation, no protein-specific
                      experiment found.
  ANNOTATED_UNREVIEWED TrEMBL/unreviewed PHA annotation.
  CANDIDATE_AMBIGUOUS  matched only through a broad/shared-chemistry signal
                      (an EC class or gene name that also covers non-PHA
                      proteins) with no confirming PHA-context text found --
                      never force these into a higher tier.

Retrieval reasons in play (uniprot_reviewed_protein_name and
uniprot_unreviewed_ec extend the spec's own named examples by the same
gene-name/protein-name/ec pattern it already establishes):
  brenda_accession, brenda_no_accession,
  uniprot_reviewed_gene_name, uniprot_reviewed_protein_name, uniprot_reviewed_ec,
  uniprot_unreviewed_gene_name, uniprot_unreviewed_protein_name, uniprot_unreviewed_ec
"""

from __future__ import annotations

GOLD_EXPERIMENTAL = "GOLD_EXPERIMENTAL"
CURATED_REVIEWED = "CURATED_REVIEWED"
ANNOTATED_UNREVIEWED = "ANNOTATED_UNREVIEWED"
CANDIDATE_AMBIGUOUS = "CANDIDATE_AMBIGUOUS"

_TIER_RANK = {
    CANDIDATE_AMBIGUOUS: 1,
    ANNOTATED_UNREVIEWED: 2,
    CURATED_REVIEWED: 3,
    GOLD_EXPERIMENTAL: 4,
}

PHA_CONTEXT_KEYWORDS = (
    "polyhydroxyalkanoate",
    "polyhydroxybutyrate",
    "poly-3-hydroxybutyrate",
    "poly(3-hydroxy",
    "hydroxyalkanoate",
    "hydroxybutyrate",
    "hydroxyvalerate",
    " pha ",
    "(pha)",
    " phb ",
    "(phb)",
)


def text_has_pha_context(*texts: str | None) -> bool:
    blob = " ".join(t for t in texts if t).lower()
    if not blob:
        return False
    return any(kw in blob for kw in PHA_CONTEXT_KEYWORDS)


def better_tier(a: str, b: str) -> str:
    return a if _TIER_RANK[a] >= _TIER_RANK[b] else b


def compute_tier_for_brenda(has_accession: bool) -> tuple[str, bool]:
    # BRENDA evidence only ever reaches ingest_brenda.py after passing the
    # require_pha_context gate (see ingest_brenda.py) or belonging to a
    # family with a clean, PHA-specific EC (phaC) -- so by the time we're
    # here it's always literature-backed for this specific protein.
    tier = GOLD_EXPERIMENTAL
    return tier, False


def compute_tier_for_uniprot(
    reason: str,
    reviewed: bool,
    specific_function_evidence: bool,
) -> tuple[str, bool]:
    """reason is one of the uniprot_* retrieval reasons above.
    specific_function_evidence: True if this entry's UniProt FUNCTION
    comment names a specific monomer/polymer from the controlled vocabulary
    (see phenotype.normalize.uniprot_function_evidence_level).
    """
    if specific_function_evidence:
        return GOLD_EXPERIMENTAL, False
    if reason.endswith("_ec"):
        # EC-based UniProt hits are the broad/shared-chemistry case the
        # spec warns about directly -- always ambiguous pending review,
        # regardless of reviewed status, unless upgraded above.
        return CANDIDATE_AMBIGUOUS, True
    if reason.endswith("_gene_name") or reason.endswith("_protein_name"):
        if reviewed:
            return CURATED_REVIEWED, False
        return ANNOTATED_UNREVIEWED, False
    return CANDIDATE_AMBIGUOUS, True
