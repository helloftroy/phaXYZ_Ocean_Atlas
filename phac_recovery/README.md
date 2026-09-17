# phaC recovery: searching for the missed synthase

Goal: for the 13,522 genomes identified as carrying >5 other PHA pathway
genes but no phaC hit (see `figures/scripts/plot_phac_negative_taxonomy.py`
and `plot_phac_negative_treemap.py`), search the genomic neighborhood
around their known other-pha-gene hits with a permissive phaC HMM, since
these are likely divergent phaC variants our original mmseqs2 seed search
missed by sequence similarity, not genomes truly lacking a synthase.

Why a neighborhood search, not whole-genome: PHA genes are typically
operonic, so restricting to each genome's own flanking ORFs (±10 genes,
naturally the whole scaffold if it has fewer) both keeps 13k+ genomes
tractable and raises confidence any hit is really part of this pathway.

## Order of operations

**1. Local (already done, artifacts in `fair_ocean_agent/phac_recovery/`)**

```
.venv/bin/python phac_recovery/build_local_manifests.py
```

Produces (from data already sitting in `fair_ocean_agent/`, no cluster
needed): `qualifying_genomes.txt` (13,522), `other_pha_hits.tsv` (251,786
genome/family/target_id anchor rows), `target_ids.txt` (170,086 distinct
target_ids needing sequence extraction), `genome_download_manifest.tsv`
(genome → its OMDB `.genes.faa.gz` URL, from the small ~7.5MB
`OMDBv2.0_data.tsv` catalog, all 13,522 confirmed present).

Also already done: `phac_recovery/hmm/PF07167.hmm` (Pfam PhaC N-terminal
domain, fetched from EBI) and `phac_recovery/hmm/phaC_custom.hmm`
(full-length, built from this repo's own 822-sequence curated phaC
diversity via `phac_recovery/hmm/build_custom_phac_hmm.sh`) — both
committed to the repo, both validated live against known true phaC vs.
random ORFs (clean separation, including a correctly-lower-scoring
divergent archaeal phaC — exactly the kind of case a relaxed threshold
is meant to still catch).

**2. scp the manifests up to the cluster**

```
scp fair_ocean_agent/phac_recovery/target_ids.txt \
    fair_ocean_agent/phac_recovery/other_pha_hits.tsv \
    fair_ocean_agent/phac_recovery/genome_download_manifest.tsv \
    <cluster>:<repo>/phac_recovery/
```

**3. Three sbatch steps, in order, on the cluster**

```
sbatch --account=191001-364393 cluster/run_phac_recovery_extract_targets.sbatch
sbatch --account=191001-364393 cluster/run_phac_recovery_extract_neighborhoods.sbatch
sbatch --account=191001-364393 cluster/run_phac_recovery_hmmsearch.sbatch
```

- **extract_targets**: pulls real sequences for the 170,086 target_ids from
  the mmseqs `target_db` already built for the original search (`mmseqs
  createsubdb` + `convert2fasta` — no need to touch the 49GB NR100.faa.gz
  or the 4.2GB cluster.tsv). CPU-only, `service` partition.
- **extract_neighborhoods**: downloads each of the 13,522 genomes' own gene
  calls directly from OMDB (predictable per-genome URL), locates each
  genome's known anchors by **exact sequence match** (not the cluster.tsv
  member list — NR100 is 100%-identity dereplication, so a genome that
  truly carries target_id T has a byte-identical copy of T's sequence
  somewhere in its own gene calls; a straight dict lookup resolves it
  without ever touching the 4.2GB member file), and writes the combined
  ±10-gene neighborhoods to one FASTA. **Needs internet** (13,522
  individual HTTPS fetches) — `service` partition, long walltime.
- **hmmsearch**: searches the neighborhoods with both HMMs at a relaxed
  `--domE 10` threshold. One-time prerequisite: `conda install -n
  pha-reference -c bioconda -c conda-forge hmmer` (not in the env by
  default).

**4. scp the results back**

```
scp <cluster>:<repo>/phac_recovery/candidate_phaC_hits.tsv fair_ocean_agent/phac_recovery/
scp <cluster>:<repo>/phac_recovery/extraction_report.tsv fair_ocean_agent/phac_recovery/
```

`candidate_phaC_hits.tsv`: `gene_id, genome, profile, score, evalue` — one
row per domain hit per HMM profile. A gene hit by *both* PF07167 and
phaC_custom is a stronger candidate than one hit by only one. `genome` is
recovered directly from `gene_id` (`<genome>-scaffold_...`), so this joins
straight back onto `genome_family_matrix.tsv` / the taxonomy work already
done for further curation (e.g. does a candidate hit concentrate in the
same Bacteroidota lineage flagged earlier?).

## What's NOT handled yet

- No attempt yet to cross-reference hits back against the taxonomic
  clustering from `plot_phac_negative_taxonomy.py` (e.g., "did the
  22-genera Bacteroidota/ABJGZF group actually turn up a consistent
  candidate phaC, confirming the single-lineage hypothesis?"). That's the
  natural next step once `candidate_phaC_hits.tsv` exists.
- No manual curation/validation of candidate hits (e.g. checking they
  have plausible catalytic residues, aren't a truncated fragment, etc.)
  — a permissive `--domE 10` threshold will produce false positives by
  design; treat this as a candidate list to triage, not a final answer.
