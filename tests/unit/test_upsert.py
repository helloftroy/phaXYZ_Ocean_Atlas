import pytest

from phaatlas.db.models import FamilyAssignment
from phaatlas.db.session import get_sessionmaker, init_db
from phaatlas.pipeline import policy
from phaatlas.pipeline.upsert import upsert_family_assignment


@pytest.fixture
def session(tmp_path):
    path = tmp_path / "test.sqlite"
    init_db(path)
    Session = get_sessionmaker(path)
    s = Session()
    from phaatlas.db.models import FamilyDefinitionRow, Protein

    s.add(Protein(protein_id="UNIPROT:P99999", uniprot_accession="P99999"))
    s.add(FamilyDefinitionRow(family_id="phaC", display_name="PhaC", role="synthase", config_json="{}"))
    s.add(FamilyDefinitionRow(family_id="phaA", display_name="PhaA", role="thiolase", config_json="{}"))
    s.commit()
    yield s
    s.close()


def test_family_assignment_widens_reasons_across_two_independent_hits(session):
    upsert_family_assignment(
        session,
        protein_id="UNIPROT:P99999",
        family_id="phaC",
        reason="uniprot_unreviewed_gene_name",
        tier=policy.ANNOTATED_UNREVIEWED,
        needs_manual_review=False,
    )
    upsert_family_assignment(
        session,
        protein_id="UNIPROT:P99999",
        family_id="phaC",
        reason="brenda_accession",
        tier=policy.GOLD_EXPERIMENTAL,
        needs_manual_review=False,
    )
    session.commit()

    row = session.query(FamilyAssignment).filter_by(protein_id="UNIPROT:P99999", family_id="phaC").one()
    assert row.evidence_tier == policy.GOLD_EXPERIMENTAL  # widened, never narrowed
    reasons = set(row.retrieval_reasons.split("|"))
    assert reasons == {"uniprot_unreviewed_gene_name", "brenda_accession"}


def test_ambiguous_hit_does_not_downgrade_a_confirmed_assignment(session):
    upsert_family_assignment(
        session,
        protein_id="UNIPROT:P99999",
        family_id="phaA",
        reason="uniprot_reviewed_gene_name",
        tier=policy.CURATED_REVIEWED,
        needs_manual_review=False,
    )
    # a later, separate EC-based query also happens to match this protein,
    # but ambiguously -- must not drag the confirmed tier back down
    upsert_family_assignment(
        session,
        protein_id="UNIPROT:P99999",
        family_id="phaA",
        reason="uniprot_reviewed_ec",
        tier=policy.CANDIDATE_AMBIGUOUS,
        needs_manual_review=True,
    )
    session.commit()

    row = session.query(FamilyAssignment).filter_by(protein_id="UNIPROT:P99999", family_id="phaA").one()
    assert row.evidence_tier == policy.CURATED_REVIEWED
    assert row.needs_manual_review is False


def test_only_one_assignment_row_per_protein_family_pair(session):
    for _ in range(3):
        upsert_family_assignment(
            session,
            protein_id="UNIPROT:P99999",
            family_id="phaC",
            reason="uniprot_reviewed_gene_name",
            tier=policy.CURATED_REVIEWED,
            needs_manual_review=False,
        )
    session.commit()
    count = session.query(FamilyAssignment).filter_by(protein_id="UNIPROT:P99999", family_id="phaC").count()
    assert count == 1
