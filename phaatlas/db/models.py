"""SQLAlchemy ORM models for PHA_reference/pha_reference.sqlite.

Normalized core tables (family_definition, protein, family_assignment,
source_evidence, phenotype, retrieval_run) plus the flattened
protein_master_export VIEW (created separately in db/views.py, not an ORM
model since it's read-only and SQLite-native).

Design notes:
  - protein.protein_id is a stable internal surrogate key, NOT the UniProt
    accession, because some BRENDA protein records have no UniProt/GenBank
    accession at all and must still get their own row (spec: "do not
    deduplicate different accessions yet by deleting records"). When a
    UniProt accession IS known, protein_id == "UNIPROT:<accession>" so the
    same protein encountered from both BRENDA and UniProt naturally
    resolves to one row. When only a BRENDA protein-local id is known,
    protein_id == "BRENDA:<ec>:<local_id>" -- and gets promoted (a second
    row is NOT created) if a later import discovers its UniProt accession
    via BRENDA's own `accessions` field.
  - Every table that stores free text preserves the field's origin
    (source_field) so provenance survives all the way to the master CSV.
  - Multiplicity is preserved deliberately: a protein can have many
    source_evidence rows (one per BRENDA reference / UniProt fetch) and
    many phenotype rows (one per literature statement) -- only
    protein_master_export flattens these with pipe-separated joins.
"""

from __future__ import annotations

import datetime as dt

from sqlalchemy import (
    Boolean,
    DateTime,
    ForeignKey,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


def _now() -> dt.datetime:
    return dt.datetime.now(dt.timezone.utc)


class FamilyDefinitionRow(Base):
    """One row per internal PHA family, mirrored from config/family_definitions.yaml.

    Re-synced (upserted) every time the CLI runs, so the database always
    carries a record of the exact config that produced it, even if the
    YAML file is edited later.
    """

    __tablename__ = "family_definition"

    family_id: Mapped[str] = mapped_column(String, primary_key=True)
    display_name: Mapped[str] = mapped_column(String, nullable=False)
    role: Mapped[str] = mapped_column(String, nullable=False)
    ec_brenda: Mapped[str | None] = mapped_column(String, nullable=True)
    require_pha_context: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    config_json: Mapped[str] = mapped_column(Text, nullable=False)
    updated_at: Mapped[dt.datetime] = mapped_column(DateTime, nullable=False, default=_now)


class Protein(Base):
    __tablename__ = "protein"

    protein_id: Mapped[str] = mapped_column(String, primary_key=True)

    uniprot_accession: Mapped[str | None] = mapped_column(String, nullable=True, unique=True)
    secondary_accessions: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated
    reviewed: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    annotation_score: Mapped[str | None] = mapped_column(String, nullable=True)
    protein_existence: Mapped[str | None] = mapped_column(String, nullable=True)

    protein_name: Mapped[str | None] = mapped_column(Text, nullable=True)
    alternative_protein_names: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated
    gene_name: Mapped[str | None] = mapped_column(String, nullable=True)
    gene_synonyms: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated

    organism: Mapped[str | None] = mapped_column(String, nullable=True)
    strain: Mapped[str | None] = mapped_column(String, nullable=True)
    ncbi_taxid: Mapped[str | None] = mapped_column(String, nullable=True)
    taxonomic_lineage: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated

    sequence: Mapped[str | None] = mapped_column(Text, nullable=True)
    sequence_length: Mapped[int | None] = mapped_column(Integer, nullable=True)
    sequence_sha256: Mapped[str | None] = mapped_column(String, nullable=True)

    ec_numbers: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated
    rhea_ids: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated
    function_text: Mapped[str | None] = mapped_column(Text, nullable=True)
    catalytic_activity_text: Mapped[str | None] = mapped_column(Text, nullable=True)
    pathway_text: Mapped[str | None] = mapped_column(Text, nullable=True)
    active_sites_text: Mapped[str | None] = mapped_column(Text, nullable=True)

    brenda_ec_numbers: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated
    brenda_record_ids: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated "ec:local_id"
    brenda_sources: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated (swissprot/trembl/genbank/...)

    source_releases: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated, e.g. "UniProt:2026_04|BRENDA:2026.1"
    uniprot_release: Mapped[str | None] = mapped_column(String, nullable=True)
    retrieved_at: Mapped[dt.datetime] = mapped_column(DateTime, nullable=False, default=_now)
    updated_at: Mapped[dt.datetime] = mapped_column(DateTime, nullable=False, default=_now, onupdate=_now)

    family_assignments: Mapped[list["FamilyAssignment"]] = relationship(back_populates="protein")
    source_evidence: Mapped[list["SourceEvidence"]] = relationship(back_populates="protein")
    phenotypes: Mapped[list["Phenotype"]] = relationship(back_populates="protein")


class FamilyAssignment(Base):
    """Assignment of one protein to one internal PHA family.

    One row per (protein_id, family_id). retrieval_reasons/evidence_tier
    are recomputed (widened, never narrowed) every time an independent
    search re-discovers the same protein for the same family -- see
    pipeline/ingest_uniprot.py's upsert_assignment.
    """

    __tablename__ = "family_assignment"
    __table_args__ = (UniqueConstraint("protein_id", "family_id", name="uq_family_assignment_protein_family"),)

    assignment_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    protein_id: Mapped[str] = mapped_column(ForeignKey("protein.protein_id"), nullable=False)
    family_id: Mapped[str] = mapped_column(ForeignKey("family_definition.family_id"), nullable=False)
    pha_subfamily: Mapped[str | None] = mapped_column(String, nullable=True)

    retrieval_reasons: Mapped[str] = mapped_column(Text, nullable=False)  # pipe-separated
    evidence_tier: Mapped[str] = mapped_column(String, nullable=False)
    confidence_notes: Mapped[str | None] = mapped_column(Text, nullable=True)
    needs_manual_review: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    assigned_at: Mapped[dt.datetime] = mapped_column(DateTime, nullable=False, default=_now)
    updated_at: Mapped[dt.datetime] = mapped_column(DateTime, nullable=False, default=_now, onupdate=_now)

    # Populated by `pha-reference cluster95` (pipeline/cluster95.py), scoped
    # to this (protein, family) row rather than to `protein` directly --
    # the one protein that legitimately belongs to two families (an
    # unresolved phaR hit, see ingest_phaR_disambiguation) can land in two
    # different clusters, one per family's own MMseqs2 run. NULL for any
    # row whose protein has no sequence (a BRENDA-only, no-accession
    # record) or that predates the first clustering run -- never backfilled
    # by deleting/recreating rows, only ever updated in place.
    cluster95_id: Mapped[str | None] = mapped_column(String, nullable=True)
    cluster95_representative: Mapped[str | None] = mapped_column(String, nullable=True)  # protein_id of the cluster representative
    cluster95_size: Mapped[int | None] = mapped_column(Integer, nullable=True)

    protein: Mapped["Protein"] = relationship(back_populates="family_assignments")


class SourceEvidence(Base):
    """One BRENDA or UniProt source record backing a protein.

    raw_payload_json retains the full untouched source record (the BRENDA
    EC-level protein/reference sub-object, or the full UniProt JSON entry)
    so nothing is discarded even though only specific fields are lifted
    into `protein`/`phenotype` -- see spec's "do not discard fields merely
    because we are not using them immediately".
    """

    __tablename__ = "source_evidence"

    evidence_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    protein_id: Mapped[str] = mapped_column(ForeignKey("protein.protein_id"), nullable=False)

    source_name: Mapped[str] = mapped_column(String, nullable=False)  # "BRENDA" | "UniProt"
    source_release: Mapped[str | None] = mapped_column(String, nullable=True)
    source_record_id: Mapped[str | None] = mapped_column(String, nullable=True)
    ec_number: Mapped[str | None] = mapped_column(String, nullable=True)
    organism_text: Mapped[str | None] = mapped_column(String, nullable=True)

    literature_pmids: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated
    literature_dois: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated

    raw_payload_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    notes: Mapped[str | None] = mapped_column(Text, nullable=True)

    retrieved_at: Mapped[dt.datetime] = mapped_column(DateTime, nullable=False, default=_now)

    protein: Mapped["Protein"] = relationship(back_populates="source_evidence")


class Phenotype(Base):
    """One literature/annotation statement about polymer/monomer phenotype,
    substrate specificity, or reaction conditions for a protein.

    Crucial evidence rule (see spec): `demonstrated_monomers` /
    `demonstrated_polymers` are populated ONLY when phenotype_evidence_level
    is 'experimental_literature' or 'curated_specific_function'. A generic
    Rhea/catalytic-activity annotation ('generic_reaction_annotation') or a
    family-membership-only inference ('inferred_annotation') must leave
    those two fields NULL even if the raw text superficially mentions a
    monomer name -- see phenotype/normalize.py.
    """

    __tablename__ = "phenotype"

    phenotype_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    protein_id: Mapped[str] = mapped_column(ForeignKey("protein.protein_id"), nullable=False)
    source_evidence_id: Mapped[int | None] = mapped_column(ForeignKey("source_evidence.evidence_id"), nullable=True)

    phenotype_type: Mapped[str] = mapped_column(String, nullable=False)  # monomer_composition | substrate_specificity | temperature | ph | specific_activity | ...
    source_field: Mapped[str] = mapped_column(String, nullable=False)  # e.g. "natural_substrates_products", "uniprot_function"

    phenotype_raw: Mapped[str] = mapped_column(Text, nullable=False)
    demonstrated_monomers: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated codes
    demonstrated_polymers: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated codes
    composition_text: Mapped[str | None] = mapped_column(Text, nullable=True)

    phenotype_evidence_level: Mapped[str] = mapped_column(String, nullable=False)

    temperature_optimum: Mapped[str | None] = mapped_column(String, nullable=True)
    temperature_range: Mapped[str | None] = mapped_column(String, nullable=True)
    temperature_stability: Mapped[str | None] = mapped_column(String, nullable=True)
    ph_optimum: Mapped[str | None] = mapped_column(String, nullable=True)
    ph_range: Mapped[str | None] = mapped_column(String, nullable=True)
    specific_activity: Mapped[str | None] = mapped_column(String, nullable=True)

    literature_pmids: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated
    literature_dois: Mapped[str | None] = mapped_column(Text, nullable=True)  # pipe-separated

    created_at: Mapped[dt.datetime] = mapped_column(DateTime, nullable=False, default=_now)

    protein: Mapped["Protein"] = relationship(back_populates="phenotypes")


class RetrievalRun(Base):
    """Provenance record for every retrieval/import run. Never overwritten
    or deleted -- append-only history of "what source/query/version produced
    what, when" (spec: "Do not silently overwrite older provenance.").
    """

    __tablename__ = "retrieval_run"

    run_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    started_at: Mapped[dt.datetime] = mapped_column(DateTime, nullable=False, default=_now)
    finished_at: Mapped[dt.datetime | None] = mapped_column(DateTime, nullable=True)

    database: Mapped[str] = mapped_column(String, nullable=False)  # "BRENDA" | "UniProt"
    database_release: Mapped[str | None] = mapped_column(String, nullable=True)
    family_id: Mapped[str | None] = mapped_column(String, nullable=True)
    query: Mapped[str | None] = mapped_column(Text, nullable=True)
    number_returned: Mapped[int | None] = mapped_column(Integer, nullable=True)
    software_version: Mapped[str | None] = mapped_column(String, nullable=True)
    status: Mapped[str] = mapped_column(String, nullable=False, default="running")  # running | completed | failed
    error_text: Mapped[str | None] = mapped_column(Text, nullable=True)
