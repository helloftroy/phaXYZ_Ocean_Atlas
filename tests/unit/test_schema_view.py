import datetime as dt

import pytest

from phaatlas.db.models import FamilyAssignment, FamilyDefinitionRow, Phenotype, Protein
from phaatlas.db.session import get_sessionmaker, init_db


@pytest.fixture
def db_path(tmp_path):
    path = tmp_path / "test_pha_reference.sqlite"
    init_db(path)
    return path


def test_init_db_creates_all_core_tables_and_the_master_view(db_path):
    import sqlite3

    conn = sqlite3.connect(str(db_path))
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    views = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='view'")}
    conn.close()
    assert {"protein", "family_definition", "family_assignment", "source_evidence", "phenotype", "retrieval_run"} <= tables
    assert "protein_master_export" in views


def test_master_view_flattens_protein_and_family_assignment(db_path):
    Session = get_sessionmaker(db_path)
    session = Session()
    try:
        protein = Protein(
            protein_id="UNIPROT:P00000",
            uniprot_accession="P00000",
            protein_name="Test PHA synthase",
            gene_name="phaC",
            organism="Testus organismus",
            sequence="MKV",
            sequence_length=3,
            sequence_sha256="deadbeef",
        )
        session.add(protein)
        session.add(FamilyDefinitionRow(family_id="phaC", display_name="PhaC", role="synthase", config_json="{}"))
        session.flush()
        session.add(
            FamilyAssignment(
                protein_id=protein.protein_id,
                family_id="phaC",
                retrieval_reasons="uniprot_reviewed_gene_name",
                evidence_tier="CURATED_REVIEWED",
                needs_manual_review=False,
            )
        )
        session.commit()

        import sqlite3

        conn = sqlite3.connect(str(db_path))
        row = conn.execute(
            "SELECT pha_family, evidence_tier, gene_name, sequence FROM protein_master_export WHERE protein_id=?",
            ("UNIPROT:P00000",),
        ).fetchone()
        conn.close()
        assert row == ("phaC", "CURATED_REVIEWED", "phaC", "MKV")
    finally:
        session.close()


def test_master_view_aggregates_phenotype_fields_pipe_separated(db_path):
    Session = get_sessionmaker(db_path)
    session = Session()
    try:
        protein = Protein(protein_id="UNIPROT:P11111", uniprot_accession="P11111")
        session.add(protein)
        session.add(FamilyDefinitionRow(family_id="phaC", display_name="PhaC", role="synthase", config_json="{}"))
        session.flush()
        session.add(
            FamilyAssignment(
                protein_id=protein.protein_id,
                family_id="phaC",
                retrieval_reasons="brenda_accession",
                evidence_tier="GOLD_EXPERIMENTAL",
                needs_manual_review=False,
            )
        )
        session.add(
            Phenotype(
                protein_id=protein.protein_id,
                phenotype_type="substrate_specificity",
                source_field="natural_substrates_products",
                phenotype_raw="accumulates PHB",
                demonstrated_polymers="PHB",
                phenotype_evidence_level="experimental_literature",
            )
        )
        session.add(
            Phenotype(
                protein_id=protein.protein_id,
                phenotype_type="substrate_specificity",
                source_field="natural_substrates_products",
                phenotype_raw="also accumulates PHBV",
                demonstrated_polymers="PHBV",
                phenotype_evidence_level="experimental_literature",
            )
        )
        session.commit()

        import sqlite3

        conn = sqlite3.connect(str(db_path))
        row = conn.execute(
            "SELECT demonstrated_polymers FROM protein_master_export WHERE protein_id=?",
            ("UNIPROT:P11111",),
        ).fetchone()
        conn.close()
        values = set(row[0].split("|"))
        assert values == {"PHB", "PHBV"}
    finally:
        session.close()
