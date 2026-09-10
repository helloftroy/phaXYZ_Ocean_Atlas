import shutil

import pytest

from phaatlas.db.models import FamilyAssignment, FamilyDefinitionRow, Protein
from phaatlas.db.session import get_sessionmaker, init_db
from phaatlas.pipeline import cluster95

MMSEQS_AVAILABLE = shutil.which("mmseqs") is not None


def _seed(session, protein_id, sequence, sha, tier, family="phaC"):
    session.add(Protein(protein_id=protein_id, uniprot_accession=protein_id.split(":")[-1], sequence=sequence, sequence_sha256=sha))
    session.add(
        FamilyAssignment(
            protein_id=protein_id,
            family_id=family,
            retrieval_reasons="uniprot_reviewed_gene_name",
            evidence_tier=tier,
            needs_manual_review=False,
        )
    )


@pytest.fixture
def session(tmp_path):
    path = tmp_path / "test.sqlite"
    init_db(path)
    Session = get_sessionmaker(path)
    s = Session()
    s.add(FamilyDefinitionRow(family_id="phaC", display_name="PhaC", role="synthase", config_json="{}"))
    s.commit()
    yield s
    s.close()


def test_pick_canonical_per_hash_groups_exact_duplicates_and_is_deterministic():
    members = [
        cluster95.FamilyMember("UNIPROT:B", "MKV", "hash1", "ANNOTATED_UNREVIEWED"),
        cluster95.FamilyMember("UNIPROT:A", "MKV", "hash1", "GOLD_EXPERIMENTAL"),
        cluster95.FamilyMember("UNIPROT:C", "MKL", "hash2", "CURATED_REVIEWED"),
    ]
    hash_to_ids, hash_to_canonical, hash_to_seq = cluster95._pick_canonical_per_hash(members)
    assert set(hash_to_ids["hash1"]) == {"UNIPROT:A", "UNIPROT:B"}
    # GOLD_EXPERIMENTAL beats ANNOTATED_UNREVIEWED as the canonical pick
    assert hash_to_canonical["hash1"] == "UNIPROT:A"
    assert hash_to_seq["hash1"] == "MKV"
    assert hash_to_ids["hash2"] == ["UNIPROT:C"]


@pytest.mark.skipif(not MMSEQS_AVAILABLE, reason="mmseqs binary not on PATH")
def test_cluster_family_end_to_end_writes_cluster95_columns(session):
    # Two exact-duplicate accessions of the same sequence (should collapse
    # to ONE mmseqs input record but both still get a cluster95 assignment
    # with cluster95_size counting both) plus one clearly distinct sequence
    # (should land in its own singleton cluster).
    seq_a = "MSTAVENLNKRFAAALAELQRSSVQPTAAQKAEATRFVQSFVQAAKADPAGAFAAAAQPLDSPAALAAYTAKLGSPMAQAWLTAQNYLENPQKVLEDQRFVVQENIRDNAKTLTGGAAARQ" * 2
    seq_b = "MTQTQTQTQTQAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAAVQAA" * 2
    import hashlib

    sha_a = hashlib.sha256(seq_a.encode()).hexdigest()
    _seed(session, "UNIPROT:AAA1", seq_a, sha_a, "GOLD_EXPERIMENTAL")
    _seed(session, "UNIPROT:AAA2", seq_a, sha_a, "ANNOTATED_UNREVIEWED")  # exact duplicate of AAA1
    sha_b = hashlib.sha256(seq_b.encode()).hexdigest()
    _seed(session, "UNIPROT:BBB1", seq_b, sha_b, "CURATED_REVIEWED")
    session.commit()

    summary = cluster95.cluster_family(session, "phaC", min_seq_id=0.95, min_cov=0.90)
    assert summary["members_with_sequence"] == 3
    assert summary["unique_sequences"] == 2
    assert summary["clusters"] == 2

    rows = {fa.protein_id: fa for fa in session.query(FamilyAssignment).filter_by(family_id="phaC").all()}
    # AAA1 (higher tier) is the canonical pick within its exact-duplicate group
    assert rows["UNIPROT:AAA1"].cluster95_representative == "UNIPROT:AAA1"
    assert rows["UNIPROT:AAA2"].cluster95_representative == "UNIPROT:AAA1"
    assert rows["UNIPROT:AAA1"].cluster95_id == rows["UNIPROT:AAA2"].cluster95_id
    assert rows["UNIPROT:AAA1"].cluster95_size == 2
    assert rows["UNIPROT:AAA2"].cluster95_size == 2
    assert rows["UNIPROT:BBB1"].cluster95_representative == "UNIPROT:BBB1"
    assert rows["UNIPROT:BBB1"].cluster95_size == 1


def test_cluster_family_no_members_is_a_noop(session):
    summary = cluster95.cluster_family(session, "phaC")
    assert summary["members_with_sequence"] == 0
