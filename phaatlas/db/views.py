"""SQL for the protein_master_export VIEW.

One row per (protein, pha_family) pair -- normally one row per protein
accession, per spec, except for the rare genuine dual-family case (e.g. a
phaR hit resolved to both phaR_synthase and phaR_regulator pending manual
review), which intentionally produces two rows rather than silently picking
one. Child-table fields (phenotype, literature refs) are flattened with
GROUP_CONCAT('|'); duplicate exact-string values are collapsed via a
DISTINCT sub-select first (SQLite's GROUP_CONCAT doesn't support DISTINCT
together with a custom separator), though dedup at the individual
pipe-token level -- e.g. the same monomer code appearing in two different
phenotype statements -- is intentionally left to exports/export.py's CSV
writer, which also does that token-level dedup for the columns copied
straight off `protein` (ec_numbers, gene_synonyms, etc).
"""

PROTEIN_MASTER_EXPORT_VIEW_SQL = """
CREATE VIEW protein_master_export AS
WITH evidence_pmids AS (
    SELECT protein_id, GROUP_CONCAT(pmids, '|') AS literature_pmids
    FROM (SELECT DISTINCT protein_id, literature_pmids AS pmids FROM source_evidence WHERE literature_pmids IS NOT NULL AND literature_pmids != '')
    GROUP BY protein_id
),
evidence_dois AS (
    SELECT protein_id, GROUP_CONCAT(dois, '|') AS literature_dois
    FROM (SELECT DISTINCT protein_id, literature_dois AS dois FROM source_evidence WHERE literature_dois IS NOT NULL AND literature_dois != '')
    GROUP BY protein_id
),
pheno_monomers AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS demonstrated_monomers
    FROM (SELECT DISTINCT protein_id, demonstrated_monomers AS v FROM phenotype WHERE demonstrated_monomers IS NOT NULL AND demonstrated_monomers != '')
    GROUP BY protein_id
),
pheno_polymers AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS demonstrated_polymers
    FROM (SELECT DISTINCT protein_id, demonstrated_polymers AS v FROM phenotype WHERE demonstrated_polymers IS NOT NULL AND demonstrated_polymers != '')
    GROUP BY protein_id
),
pheno_composition AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS polymer_composition_text
    FROM (SELECT DISTINCT protein_id, composition_text AS v FROM phenotype WHERE composition_text IS NOT NULL AND composition_text != '')
    GROUP BY protein_id
),
pheno_substrate AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS substrate_specificity_text
    FROM (SELECT DISTINCT protein_id, phenotype_raw AS v FROM phenotype WHERE phenotype_type = 'substrate_specificity')
    GROUP BY protein_id
),
pheno_evidence_level AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS phenotype_evidence_level
    FROM (SELECT DISTINCT protein_id, phenotype_evidence_level AS v FROM phenotype)
    GROUP BY protein_id
),
pheno_refs AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS phenotype_references
    FROM (
        SELECT DISTINCT protein_id, literature_pmids AS v FROM phenotype WHERE literature_pmids IS NOT NULL AND literature_pmids != ''
        UNION
        SELECT DISTINCT protein_id, literature_dois AS v FROM phenotype WHERE literature_dois IS NOT NULL AND literature_dois != ''
    )
    GROUP BY protein_id
),
pheno_temp_opt AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS temperature_optimum
    FROM (SELECT DISTINCT protein_id, temperature_optimum AS v FROM phenotype WHERE temperature_optimum IS NOT NULL AND temperature_optimum != '')
    GROUP BY protein_id
),
pheno_temp_range AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS temperature_range
    FROM (SELECT DISTINCT protein_id, temperature_range AS v FROM phenotype WHERE temperature_range IS NOT NULL AND temperature_range != '')
    GROUP BY protein_id
),
pheno_temp_stab AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS temperature_stability
    FROM (SELECT DISTINCT protein_id, temperature_stability AS v FROM phenotype WHERE temperature_stability IS NOT NULL AND temperature_stability != '')
    GROUP BY protein_id
),
pheno_ph_opt AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS ph_optimum
    FROM (SELECT DISTINCT protein_id, ph_optimum AS v FROM phenotype WHERE ph_optimum IS NOT NULL AND ph_optimum != '')
    GROUP BY protein_id
),
pheno_ph_range AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS ph_range
    FROM (SELECT DISTINCT protein_id, ph_range AS v FROM phenotype WHERE ph_range IS NOT NULL AND ph_range != '')
    GROUP BY protein_id
),
pheno_spec_act AS (
    SELECT protein_id, GROUP_CONCAT(v, '|') AS specific_activity
    FROM (SELECT DISTINCT protein_id, specific_activity AS v FROM phenotype WHERE specific_activity IS NOT NULL AND specific_activity != '')
    GROUP BY protein_id
)
SELECT
    p.protein_id AS protein_id,
    fa.family_id AS pha_family,
    fa.pha_subfamily AS pha_subfamily,
    p.uniprot_accession AS uniprot_accession,
    p.secondary_accessions AS secondary_accessions,
    p.reviewed AS reviewed,
    p.protein_name AS protein_name,
    p.alternative_protein_names AS alternative_protein_names,
    p.gene_name AS gene_name,
    p.gene_synonyms AS gene_synonyms,
    p.organism AS organism,
    p.strain AS strain,
    p.ncbi_taxid AS ncbi_taxid,
    p.taxonomic_lineage AS taxonomic_lineage,
    p.sequence_length AS sequence_length,
    p.sequence_sha256 AS sequence_sha256,
    p.sequence AS sequence,
    p.protein_existence AS protein_existence,
    p.annotation_score AS annotation_score,
    p.ec_numbers AS ec_numbers,
    p.rhea_ids AS rhea_ids,
    p.function_text AS function_text,
    p.catalytic_activity_text AS catalytic_activity_text,
    p.pathway_text AS pathway_text,
    p.brenda_ec_numbers AS brenda_ec_numbers,
    p.brenda_record_ids AS brenda_record_ids,
    p.brenda_sources AS brenda_sources,
    fa.retrieval_reasons AS retrieval_reason,
    fa.evidence_tier AS evidence_tier,
    ep.literature_pmids AS literature_pmids,
    ed.literature_dois AS literature_dois,
    pm.demonstrated_monomers AS demonstrated_monomers,
    pp.demonstrated_polymers AS demonstrated_polymers,
    pc.polymer_composition_text AS polymer_composition_text,
    ps.substrate_specificity_text AS substrate_specificity_text,
    pel.phenotype_evidence_level AS phenotype_evidence_level,
    pr.phenotype_references AS phenotype_references,
    pto.temperature_optimum AS temperature_optimum,
    ptr.temperature_range AS temperature_range,
    pts.temperature_stability AS temperature_stability,
    pho.ph_optimum AS ph_optimum,
    phr.ph_range AS ph_range,
    psa.specific_activity AS specific_activity,
    p.source_releases AS source_releases,
    p.retrieved_at AS retrieved_at
FROM protein p
LEFT JOIN family_assignment fa ON fa.protein_id = p.protein_id
LEFT JOIN evidence_pmids ep ON ep.protein_id = p.protein_id
LEFT JOIN evidence_dois ed ON ed.protein_id = p.protein_id
LEFT JOIN pheno_monomers pm ON pm.protein_id = p.protein_id
LEFT JOIN pheno_polymers pp ON pp.protein_id = p.protein_id
LEFT JOIN pheno_composition pc ON pc.protein_id = p.protein_id
LEFT JOIN pheno_substrate ps ON ps.protein_id = p.protein_id
LEFT JOIN pheno_evidence_level pel ON pel.protein_id = p.protein_id
LEFT JOIN pheno_refs pr ON pr.protein_id = p.protein_id
LEFT JOIN pheno_temp_opt pto ON pto.protein_id = p.protein_id
LEFT JOIN pheno_temp_range ptr ON ptr.protein_id = p.protein_id
LEFT JOIN pheno_temp_stab pts ON pts.protein_id = p.protein_id
LEFT JOIN pheno_ph_opt pho ON pho.protein_id = p.protein_id
LEFT JOIN pheno_ph_range phr ON phr.protein_id = p.protein_id
LEFT JOIN pheno_spec_act psa ON psa.protein_id = p.protein_id
"""
