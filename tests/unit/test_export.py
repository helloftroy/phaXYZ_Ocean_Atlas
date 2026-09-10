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
        "retrieval_reason", "evidence_tier", "cluster95_id", "cluster95_representative", "cluster95_size",
        "literature_pmids", "literature_dois", "demonstrated_monomers",
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


def test_export_fasta_all_and_nr95_and_cluster_tsv(tmp_path):
    db_path = tmp_path / "t2.sqlite"
    init_db(db_path)
    Session = get_sessionmaker(db_path)
    session = Session()
    session.add(FamilyDefinitionRow(family_id="phaC", display_name="PhaC", role="synthase", config_json="{}"))
    session.add(Protein(protein_id="UNIPROT:R1", uniprot_accession="R1", sequence="MKV", sequence_sha256="h1"))
    session.add(Protein(protein_id="UNIPROT:M1", uniprot_accession="M1", sequence="MKL", sequence_sha256="h2"))
    session.flush()
    # R1 is its own cluster representative (size 2); M1 is also in that
    # SAME cluster despite a different sequence -- exercises that nr95 only
    # emits the row whose protein_id equals its own cluster95_representative.
    session.add(
        FamilyAssignment(
            protein_id="UNIPROT:R1", family_id="phaC", retrieval_reasons="uniprot_reviewed_gene_name",
            evidence_tier="GOLD_EXPERIMENTAL", needs_manual_review=False,
            cluster95_id="phaC:UNIPROT:R1", cluster95_representative="UNIPROT:R1", cluster95_size=2,
        )
    )
    session.add(
        FamilyAssignment(
            protein_id="UNIPROT:M1", family_id="phaC", retrieval_reasons="uniprot_reviewed_gene_name",
            evidence_tier="CURATED_REVIEWED", needs_manual_review=False,
            cluster95_id="phaC:UNIPROT:R1", cluster95_representative="UNIPROT:R1", cluster95_size=2,
        )
    )
    session.commit()
    session.close()

    all_faa = tmp_path / "all.faa"
    nr95_faa = tmp_path / "nr95.faa"
    cluster_tsv = tmp_path / "cluster95.tsv"
    n_all = export_pipeline.export_fasta_all(db_path, all_faa)
    n_nr95 = export_pipeline.export_fasta_nr95(db_path, nr95_faa)
    n_tsv = export_pipeline.export_cluster95_tsv(db_path, cluster_tsv)

    assert n_all == 2
    assert n_nr95 == 1  # only R1, since it's the only row where protein_id == cluster95_representative
    assert ">UNIPROT:R1" in nr95_faa.read_text()
    assert ">UNIPROT:M1" not in nr95_faa.read_text()
    assert n_tsv == 2
    tsv_text = cluster_tsv.read_text()
    assert "UNIPROT:R1\tphaC\tphaC:UNIPROT:R1\tUNIPROT:R1\t2" in tsv_text
    assert "UNIPROT:M1\tphaC\tphaC:UNIPROT:R1\tUNIPROT:R1\t2" in tsv_text
