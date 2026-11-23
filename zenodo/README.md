# PHA Ocean Atlas — data deposit

Derived data behind a survey of polyhydroxyalkanoate (PHA) biosynthesis genes
across **274,282 ocean metagenome-assembled and isolate genomes** from the Ocean
Microbiomics Database (OMDB v2.0). The deposit covers the full PHA gene
complement, with the PHA synthase (**phaC**) annotated in the most depth:
sequence evidence, catalytic-site scoring, synthase class, clustering at four
identity thresholds, predicted structure, and predicted thermostability.

**39 files, 293 MB compressed.** Large tables are gzipped; `MANIFEST.tsv` lists
every file with its byte size, data-row count, SHA-256 checksum and a one-line
description.

An interactive browser for these data is at
**https://helloftroy.github.io/phaXYZ_Ocean_Atlas/** — it loads a condensed
extract of this deposit and lets you filter genomes, phaC proteins and phaC
clusters by taxonomy, habitat, depth, gene content and annotation quality, map
occurrences, and export any filtered subset as CSV.

Code that produced everything here: https://github.com/helloftroy/phaXYZ_Ocean_Atlas

## Headline numbers

| | |
|---|---:|
| Genomes surveyed | 274,282 |
| Genomes with at least one PHA-family gene | 241,864 |
| phaC-positive genomes | 38,744 |
| Distinct phaC protein sequences (QC-passing) | 33,770 |
| phaC clusters at 70% identity | 4,345 |
| Novel phaC recovered from apparently phaC-negative genomes | 255 |
| ESMFold models | 2,108 |

## Contents

**`01_gene_annotations/`** — one table per PHA family (phaA, phaB, phaC, phaD,
phaE, phaF, phaG, phaI, phaJ, phaP, phaR, phaY, phaZ), each row a
protein–genome pair with GTDB taxonomy, study and sample metadata, sampling
coordinates, and the best-matching reference with identity and coverage. Plus
`genome_family_matrix.tsv.gz` (genome × family copy number, with a derived
pathway-architecture label) and `omdb_all_genome_taxonomy.tsv.gz` (taxonomy for
all 274,282 genomes, including PHA-negative ones, for use as a denominator).

**`02_clustering/`** — MMseqs2 cluster assignments for all phaC proteins at 50,
60, 70 and 90% identity; cluster representative sequences; and
`phaC_all_genomes_from_nr100_clusters.tsv.gz`, the uncapped protein-to-genome
membership list. **Use this file for any genome count** (see Caveats).

**`03_catalytic_site/`** — per-protein catalytic-triad scoring against a
625-column phaC profile HMM (the residues aligned to the Cys/Asp/His columns and
a completeness flag), class I/II/III HMM scores with trusted-cutoff calls, and a
sequence-evidence table with lengths, coverages and query-to-reference length
ratios.

**`04_structure/`** — Foldseek best hits of each predicted model against folded
references, with TM-scores and coverages; fold verdicts for the 255 novel
candidates; substrate-channel geometry for a multi-copy *Desulfoluna* genome
alongside the same measurement on PhaC crystal structures of known chain-length
preference; and `esmfold_models.tar.gz`, the 2,108 ESMFold PDB models.

**`05_thermostability/`** — TemStaPro predictions per phaC protein at eight
temperature thresholds from 40 to 80 °C.

**`06_novel_candidates/`** — the 255 phaC genes recovered from genomes the
primary annotation called phaC-negative despite carrying more than five other
PHA genes, with sequences, the full candidate list with verdicts, triad scoring
for every candidate examined, and the 4,170 target genomes searched.

**`07_environment/`** — sample provenance and coordinates for every genome, and
phaC-positive genomes joined to World Ocean Atlas 2023 annual climatological
means for temperature, salinity and dissolved oxygen at their sampling
coordinate and depth.

**`08_figure_data/`** — the companion table behind each published figure, so any
plotted value can be read directly without rerunning the pipeline.

## Caveats

Please read these before using the tables; each one has bitten us.

1. **Genome counts: use the uncapped membership list.** The metadata join that
   produced the per-family annotation tables and `genome_family_matrix.tsv`
   records at most 5 genomes per identical protein sequence. For phaC this
   recorded 7,298 genuine carriers as phaC-negative — the matrix gives 31,464
   phaC-positive genomes, while the correct figure is **38,744**, from
   `02_clustering/phaC_all_genomes_from_nr100_clusters.tsv.gz`. The same cap
   applies to every other family, so **all per-family genome counts derived from
   the matrix are lower bounds**. The size of that undercount was measured only
   for phaC.

2. **Environmental joins use the capped subset.** The WOA23 tables in
   `07_environment/` were built from the capped genome set (31,443 rows) and so
   do not cover the cap-hidden carriers. Depth is also sparse: roughly a third of
   phaC-positive genomes have usable coordinates and depth together.

3. **Catalytic-triad completeness is a floor, not a measurement.** It is scored
   from HMM alignment columns, and glocal alignment collapses divergent regions
   into insert states, so divergent proteins produce false negatives. Where a
   structure exists, geometric triad detection finds intact active sites in
   proteins this flag calls incomplete.

4. **Structures and thermostability are predictions.** ESMFold models are not
   experimental structures; side chains are modelled and absolute pocket volumes
   should be read comparatively, not absolutely. TemStaPro output is a predicted
   label, not a melting temperature. Nothing here has been validated
   biochemically.

5. **Structural scores are not available for every protein.** qTM-based fold
   confidence exists for the 15,867 representatives after 95% dereplication, not
   all 33,770 proteins. A missing value means not measured, not low.

6. **PhaC and PhaZ share a fold.** PHA synthases and PHA depolymerases are both
   α/β-hydrolases with a nucleophile–His–acid triad. An intact triad and a
   confident fold match show a folded, catalytically plausible enzyme; they do
   not by themselves establish synthase function. Proteins where the depolymerase
   evidence is competitive are flagged in the sequence-evidence table.

## Not included

Raw OMDB genomes and protein catalogues (obtain from the Ocean Microbiomics
Database), World Ocean Atlas 2023 gridded fields (NOAA NCEI), and the curated
reference database, which is rebuildable from BRENDA and UniProt with the
pipeline in the code repository.

## Provenance

- Genomes and proteins: Ocean Microbiomics Database v2.0 (`OMDBv2.0_AA_G_NR100`).
- Taxonomy: GTDB, as distributed with OMDB.
- Oceanographic context: World Ocean Atlas 2023 1° annual objectively analysed
  climatological means (NOAA NCEI).
- Coastlines used in figures and the web map: Natural Earth 110 m land (public
  domain).

## Citation

> *Manuscript in preparation — citation to be added on publication.*

Please cite both this deposit and the paper. Individual tables carry the caveats
above; if you reuse one in isolation, carry the relevant caveat with it.

## Licence

Data files: CC BY 4.0. Analysis code in the linked repository: see that
repository's licence. Third-party inputs retain their own terms.
