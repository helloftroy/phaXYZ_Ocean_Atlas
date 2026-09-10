from phaatlas.pipeline.export import MASTER_COLUMNS, _dedup_pipe
from phaatlas.db.models import FamilyAssignment, FamilyDefinitionRow, Phenotype, Protein
from phaatlas.db.session import get_sessionmaker, init_db
from phaatlas.pipeline import export as export_pipeline


def test_dedup_pipe_preserves_order_and_drops_duplicates():
    assert _dedup_pipe("3HB|4HB|3HB|PHB") == "3HB|4HB|PHB"


def test_dedup_pipe_handles_none_and_empty():
    assert _dedup_pipe(None) is None
    assert _dedup_pipe("") == ""


def test_master_columns_match_spec_exactly():
    expected = [
        "protein_id", "pha_family", "pha_subfamily", "uniprot_accession", "secondary_accessions",
        "reviewed", "protein_name", "alternative_protein_names", "gene_name", "gene_synonyms",
        "organism", "strain", "ncbi_taxid", "taxonomic_lineage", "sequence_length", "sequence_sha256",
        "sequence", "protein_existence", "annotation_score", "ec_numbers", "rhea_ids", "function_text",
        "catalytic_activity_text", "pathway_text", "brenda_ec_numbers", "brenda_record_ids", "brenda_sources",
        "retrieval_reason", "evidence_tier", "literature_pmids", "literature_dois", "demonstrated_monomers",
        "demonstrated_polymers", "polymer_composition_text", "substrate_specificity_text",
        "phenotype_evidence_level", "phenotype_references", "temperature_optimum", "temperature_range",
        "temperature_stability", "ph_optimum", "ph_range", "specific_activity", "source_releases", "retrieved_at",
    ]
    assert MASTER_COLUMNS == expected


def test_export_csv_and_fasta_round_trip(tmp_path):
    db_path = tmp_path / "t.sqlite"
    init_db(db_path)
    Session = get_sessionmaker(db_path)
    session = Session()
    session.add(
        Protein(
            protein_id="UNIPROT:P22222",
            uniprot_accession="P22222",
            gene_name="phaC",
            organism="Testus organismus",
            sequence="MKVLLLLLL" * 10,  # >60 chars, exercises FASTA line-wrapping
            sequence_length=90,
            sequence_sha256="abc123",
        )
    )
    session.add(FamilyDefinitionRow(family_id="phaC", display_name="PhaC", role="synthase", config_json="{}"))
    session.flush()
    session.add(
        FamilyAssignment(
            protein_id="UNIPROT:P22222",
            family_id="phaC",
            retrieval_reasons="uniprot_reviewed_gene_name",
            evidence_tier="CURATED_REVIEWED",
            needs_manual_review=False,
        )
    )
    session.commit()
    session.close()

    csv_path = tmp_path / "master.csv"
    faa_path = tmp_path / "ref.faa"
    n_rows = export_pipeline.export_master_csv(db_path, csv_path)
    n_seqs = export_pipeline.export_fasta(db_path, faa_path)

    assert n_rows == 1
    assert n_seqs == 1
    faa_text = faa_path.read_text()
    assert faa_text.startswith(">UNIPROT:P22222")
    # every sequence line except the header/last should be exactly 60 chars
    lines = [l for l in faa_text.splitlines() if not l.startswith(">")]
    assert all(len(l) <= 60 for l in lines)
    assert "".join(lines) == "MKVLLLLLL" * 10
