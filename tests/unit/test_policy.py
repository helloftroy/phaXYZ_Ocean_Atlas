from phaatlas.pipeline import policy


def test_ec_based_hits_are_always_candidate_ambiguous_unless_specific_function():
    tier, needs_review = policy.compute_tier_for_uniprot("uniprot_reviewed_ec", reviewed=True, specific_function_evidence=False)
    assert tier == policy.CANDIDATE_AMBIGUOUS
    assert needs_review is True


def test_ec_based_hit_upgraded_when_specific_function_evidence_present():
    tier, needs_review = policy.compute_tier_for_uniprot("uniprot_reviewed_ec", reviewed=True, specific_function_evidence=True)
    assert tier == policy.GOLD_EXPERIMENTAL
    assert needs_review is False


def test_reviewed_gene_name_hit_is_curated_reviewed():
    tier, needs_review = policy.compute_tier_for_uniprot("uniprot_reviewed_gene_name", reviewed=True, specific_function_evidence=False)
    assert tier == policy.CURATED_REVIEWED
    assert needs_review is False


def test_unreviewed_gene_name_hit_is_annotated_unreviewed():
    tier, needs_review = policy.compute_tier_for_uniprot("uniprot_unreviewed_gene_name", reviewed=False, specific_function_evidence=False)
    assert tier == policy.ANNOTATED_UNREVIEWED
    assert needs_review is False


def test_brenda_evidence_is_always_gold():
    tier, needs_review = policy.compute_tier_for_brenda(has_accession=True)
    assert tier == policy.GOLD_EXPERIMENTAL
    assert needs_review is False
    tier, needs_review = policy.compute_tier_for_brenda(has_accession=False)
    assert tier == policy.GOLD_EXPERIMENTAL


def test_better_tier_never_downgrades():
    assert policy.better_tier(policy.GOLD_EXPERIMENTAL, policy.CANDIDATE_AMBIGUOUS) == policy.GOLD_EXPERIMENTAL
    assert policy.better_tier(policy.CANDIDATE_AMBIGUOUS, policy.CURATED_REVIEWED) == policy.CURATED_REVIEWED


def test_pha_context_keyword_matching():
    assert policy.text_has_pha_context("accumulates polyhydroxybutyrate under nitrogen limitation")
    assert policy.text_has_pha_context("acetyl-CoA acetyltransferase", "involved in PHB biosynthesis")
    assert not policy.text_has_pha_context("generic acetyl-CoA acetyltransferase in fatty acid metabolism")
    assert not policy.text_has_pha_context(None, "")
