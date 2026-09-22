# PHA Gene Family Atlas: Verified Results

This document covers the ocean-metagenome PHA (polyhydroxyalkanoate) pathway gene search — 15 gene families, searched against a global ocean gene catalog, with every recruiting reference sequence individually verified against UniProt and InterPro before being counted. It documents the current, verified state of the dataset going forward.

## 1. Setup

### 1.1 Reference sequences: how each family's search bait was assembled

Each of the 15 PHA pathway gene families (phaA, phaB, phaC, phaD, phaE, phaF, phaG, phaI, phaJ, phaP, phaQ, phaR_regulator, phaR_synthase, phaY, phaZ) is defined in `config/family_definitions.yaml` by: its known UniProt gene-name aliases (e.g. phaA/phbA), a set of free-text protein-name search phrases (e.g. "acetyl-CoA acetyltransferase," "acetoacetyl-CoA reductase"), and — where the family has one — a BRENDA EC-number seed (e.g. phaA=EC 2.3.1.9, phaB=EC 1.1.1.36). Two families that share the literal gene name "phaR" for genuinely different proteins — the class III/IV PHA synthase small subunit (phaR_synthase) and the PHA-responsive transcriptional autorepressor (phaR_regulator) — are kept as deliberately separate families, disambiguated by required co-occurring text (phaR_synthase requires "synthase"/"class IV"/"PhaRC" context; phaR_regulator requires "regulator"/"repressor"/"transcriptional" context).

Candidate sequences were pulled from BRENDA (by EC number, where applicable) and UniProt (by gene name and protein-name phrase), then reduced to non-redundant representatives at 95% identity / 90% coverage (`mmseqs easy-cluster --min-seq-id 0.95 -c 0.90 --cov-mode 0`) — one query FASTA per family, headers carrying only the protein accession.

### 1.2 Target database

The search target is **OMDBv2.0_AA_G_NR100** — the Ocean Microbiomics Database's amino-acid gene catalog, dereplicated to 100% identity. Target IDs throughout this project (e.g. `OMDBv2.0_AA_G_NR100_000103673917`) trace directly to entries in this catalog, each ultimately resolvable back to a specific metagenome-assembled genome (MAG) or isolate genome, its GTDB taxonomy, and its sample's collection metadata (location, depth, ecosystem type).

### 1.3 Search

Each family's query set was searched against the target database with MMseqs2 (GPU-accelerated build):

```
mmseqs search <query_db> <omdb_target_db> <result_db> <tmp> \
    -s 7.0 -e 1e-5 -c 0.0 --max-seqs 10000 --alignment-mode 3 --gpu 1

mmseqs convertalis <query_db> <omdb_target_db> <result_db> <hits.tsv> \
    --sub-mat aa:blosum62.out --format-mode 4 --translation-table 1 \
    --gap-open aa:11 --gap-extend aa:1 \
    --format-output query,target,evalue,bits,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov
```

Sensitivity 7.0 (MMseqs2's most sensitive setting) and coverage 0.0 at search time were deliberate: coverage is filtered downstream per-analysis rather than at the alignment step, so fragmentary hits (partial gene calls, contig edges — common in MAG data) aren't silently discarded before they can even be evaluated. Every family's target hits were then collapsed to one row per target (`<family>_unique_targets_with_metadata.tsv`): the best-scoring reference match (`best_query`), its identity/coverage/E-value, and the target's resolved genome/taxonomy/sample metadata.

### 1.4 Cross-referencing every reference for accuracy

Before any hit was counted, the reference sequence that recruited it (`best_query`) was independently verified — not assumed correct because it originated from the curated BRENDA/UniProt retrieval above. For every one of the roughly 9,600 distinct reference sequences that ever recruited a hit across all 15 families:

1. **UniProt annotation retrieved directly** (protein name, organism, sequence length) — not taken from the original retrieval's own labeling, which turned out to be unreliable for a meaningful share of references (a real, recurring bacterial gene-naming collision: a multi-subunit ion-transport operon, unrelated to PHA metabolism, happens to also be historically abbreviated "Pha" with subunits lettered A–G, colliding one-for-one with this project's own phaA–phaG gene names).
2. **InterPro/Pfam domain composition checked for every reference, not just a sample** (not name text alone) — this is the step that catches a reference whose *name* looks right but whose actual domain architecture doesn't match. Each of the ~9,600 references was resolved to its full InterPro match set and checked against a confirmed-good marker-domain signature for its family (e.g. IPR011283 specifically for phaB's acetoacetyl-CoA reductase, IPR002155 for phaA's thiolase fold, IPR010134 for the real PhaR autorepressor) and against a confirmed-bad list of unrelated domain families (AraC-type regulators, the antiporter-subunit domains, phenylacetate-pathway enzymes, polyketide/fatty-acid synthase domains, and others). This caught references with generic, non-specific names (e.g. "PhaR protein," "PhaD protein," "PhaF2 protein," even one literally named "Pesticidal protein Cry1Ba") that were, on inspection, a structurally unrelated protein family wearing a misleading label — and, in the other direction, recovered several genuine references that generic wording alone hadn't confidently matched.
3. **Cross-checked against this project's own `config/family_definitions.yaml`** identity for that family, so a real PHA-pathway protein recruited under the *wrong* family (e.g. a genuine depolymerase reference filed as phaD) was also caught, not just outright non-PHA contamination.

References that failed this check were added to an exclusion list (`FAMILY_BAD_QUERIES` in `phaatlas/pipeline/pathway_architecture.py`) and every target they recruited was removed from that family's counts. The genome × gene-family matrix and every count in this document reflect that corrected, verified reference set.

## 2. Verified family scope

Genomes with at least one verified hit, per family, out of the full ocean-genome dataset searched:

| Family | Role | Genomes (verified) |
|---|---|---|
| phaZ | PHA depolymerase | 156,495 |
| phaB | Acetoacetyl-CoA reductase | 174,573 |
| phaA | β-ketothiolase | 184,489 |
| phaJ | (R)-specific enoyl-CoA hydratase | 86,056 |
| phaC | PHA synthase | 68,424 |
| phaR_regulator | PHA-responsive transcriptional regulator | 20,292 |
| phaF | Granule-associated phasin/regulator | 19,902 |
| phaY | Intracellular PHA oligomer hydrolase | 18,375 |
| phaP | Phasin (granule surface protein) | 16,097 |
| phaG | Fatty-acid-synthesis-to-PHA precursor transacylase | 4,995 |
| phaE | Class III synthase partner subunit | 3,167 |
| phaI | Phasin | 5,153 |
| phaQ | Regulatory protein | 6,323 |
| phaR_synthase | Class III/IV synthase partner subunit | 464 |
| phaD | Regulatory protein | 823* |

\* phaD's verified count is small because essentially none of the reference sequences originally retrieved for it turned out to be genuine matches to its defined identity ("PHA regulatory protein PhaD") — every high-volume candidate resolved, on inspection, to a different protein. The 823 genomes retained here come from the long tail of lower-volume references that did check out; this family should be treated as sparse/exploratory going forward rather than well-characterized, pending a fresh, more targeted reference search specifically for phaD.

phaM (granule/nucleoid protein) returned zero hits in the original search and is not yet meaningfully represented in this dataset.

## 3. phaC evidence tiers: catalytic triad, HMM, and verified pathway context

For all 68,424 verified phaC genomes, five evidence tiers were built, strongest evidence first and mutually exclusive:

1. **Catalytic triad complete** — the Cys-Asp-His catalytic triad (confirmed via P23608/Cupriavidus necator PhaC1 numbering: Cys319 in a G-x-C-x-G lipase-box-like motif, Asp480, His508) is present at all three expected positions, located by projecting every target onto the project's own 625-column PhaC profile HMM via `hmmalign`.
2. **HMM-supported (no triad)** — hits at least one of 6 independently-built PhaC profile HMMs (Pfam PF07167, this project's own 822-sequence model, NCBIFam TIGR01838/01839/01836 for Class I/II/III, PANTHER PTHR36837) but the triad wasn't resolvable (usually because the alignment doesn't reach that far — see below).
3. **No HMM/triad, ≥5 other PHA genes** — no direct protein-level confirmation, but the genome carries 5 or more *other, individually verified* PHA pathway genes (§1.4's reference audit applied to every family, not just phaC).
4. **No HMM/triad, 1–4 other PHA genes** — same, with fewer supporting genes.
5. **phaC only** — no HMM hit, no triad, and no other verified PHA pathway gene in the genome at all.

![evidence tiers](phac_verified_triad_hmm_pathway_groups.png)

| Tier | Genomes | % |
|---|---|---|
| Catalytic triad complete | 28,737 | 42.0% |
| HMM-supported (no triad) | 3,790 | 5.5% |
| No HMM/triad, ≥5 other PHA genes | 4,620 | 6.8% |
| No HMM/triad, 1–4 other PHA genes | 30,050 | 43.9% |
| phaC only | 1,227 | 1.8% |

Two things stand out immediately against the pre-audit version of this same breakdown. First, the highest-confidence tier is much larger than it looked before verification (triad-complete alone is 42.0%, versus the old HMM-supported figure of 23.5%) — a direct effect of the reference cleanup concentrating genuine signal rather than diluting it across contaminated hits. Second, the "≥5 other genes" tier collapsed from 37.9% of the dataset (pre-audit) to 6.8% now — because most of the genomes that used to qualify were only reaching 5+ genes by counting phaD/phaE/phaG hits that turned out to be a different protein entirely (§1.4, and see `PHA_ALL_FAMILIES_REFERENCE_AUDIT.md`). The "1-4 other genes" tier is now the single largest bucket (43.9%) — a large population of phaC calls that remain genuinely uncertain at the direct-evidence level even after full verification, not resolved one way or the other by this analysis. Only 1.8% of genomes have phaC standing completely alone with no other supporting evidence of any kind.

### Taxonomy

| Tier | Distinct phyla | Distinct genera | % Pseudomonadota |
|---|---|---|---|
| Catalytic triad complete | 42 | 2,005 | 86.1% |
| HMM-supported (no triad) | 44 | 747 | 63.1% |
| ≥5 other PHA genes | 38 | 925 | 71.8% |
| 1–4 other PHA genes | 82 | 2,415 | 60.6% |
| phaC only | 28 | 248 | 66.9% |

![phylum composition](phac_verified_group_phylum_composition.png)

The taxonomic pattern found before the audit survives correction, essentially unchanged: the tier with the strongest direct protein evidence (catalytic triad) is the most Pseudomonadota-dominated (86.1% — the phylum nearly every curated PhaC reference comes from), while the "1-4 other genes" tier — now the largest single group, and genuinely unresolved at the protein level — is both the most taxonomically diverse by far (82 distinct phyla, more than double any other tier) and the least Pseudomonadota-heavy (60.6%), picking up real share in Bacteroidota (16.3% of that tier), Marinisomatota, SAR324, and Verrucomicrobiota. This is consistent with the same interpretation as before: genomes lacking direct confirmation skew toward lineages current PhaC references simply cover less well, not toward a biologically different (less real) population.

### Location / habitat — the sponge signal reverses

This is the one pattern that did **not** survive verification. Before the audit, marine-sponge-tissue genomes looked enriched in the high-gene-count ("≥5 other genes") tier. With verified gene counts, that's no longer true — it's reversed:

| Tier | Marine sponge tissue | vs. dataset-wide baseline (6.5%) |
|---|---|---|
| Catalytic triad complete | 8.7% | **1.32x enriched** |
| HMM-supported (no triad) | 7.1% | 1.09x |
| ≥5 other PHA genes | 4.0% | **0.62x depleted** |
| 1–4 other PHA genes | 5.1% | 0.78x |
| phaC only | 0.6% | 0.09x depleted |

Sponge-associated genomes are now enriched in the *highest*-confidence tier (catalytic triad complete, 1.32x baseline) and depleted in the tier that used to look sponge-rich. The most likely explanation: the old "≥5 other genes ⇒ sponge-enriched" signal was an artifact of exactly the contamination this audit removed — phaD, phaE, and phaG lost 67-99% of their genomes to the antiporter-locus naming collision (§1.4), and whatever taxonomic/habitat skew existed among the specific mislabeled reference sequences responsible would have shown up as a spurious pattern in the uncorrected pathway-richness counts. This is a concrete example of why the reference audit mattered beyond just shrinking numbers: at least one previously-reported ecological pattern in this dataset was not real.

Depth shows no comparably dramatic pattern: median resolved depth is 10-20 m across every tier, and every tier reaches down to hadal-trench depths (max 9,700-10,900 m) at similar rates (~5-7% of depth-resolved genomes below 1,000 m in every tier, phaC-only included). Depth/water-column position does not track evidence tier the way sponge association does.

**Files:** `figures/scripts/plot_phac_triad_hmm_pathway_groups.py`, `figures/scripts/plot_phac_verified_group_phylum.py`; `/tmp/phac_verified_triad_hmm_group.pkl` (genome → tier label, 68,424 entries) and `/tmp/verified_phac_genome_targets.pkl` (genome → verified target_ids) if continuing this analysis.

## 4. Sequence evidence against the full audited reference set

Every verified phaC target was searched against all 1,875 audited references (§1.4's ON_TARGET/FUSION set — not just the 14-sequence "trusted" benchmark used in earlier passes) with MMseqs2 (`-s 7.0 -e 1e-5 -c 0.0`), and reduced to its single best hit by bit score. For each of the 82,846 verified targets this gives: protein length, the matched reference's own accession and length, %identity, alignment length, E-value, query/target coverage, and the query-to-reference length ratio — the full table is `catalytic_domain/phac_sequence_evidence.tsv`.

![identity vs coverage](phac_pident_vs_coverage.png)

The scatter above shows every point at once but overplots into a solid mass at this n (81,514 points) — the violin view below separates identity from coverage per tier and makes the actual shape of each distribution, not just its density, legible:

![evidence violins](phac_evidence_violins.png)

The two metrics tell different stories. **%identity drops hard** across tiers — median 53% (triad-complete) → 48% (HMM-supported) → 28% → 26% → 26% (the three tiers with no direct protein-level support). **Query coverage barely moves** — 95% → 94% → 82% → 83% → 84%. If the unresolved tiers were spurious short-fragment matches, coverage would have collapsed alongside identity; it didn't.

Six metrics, median per tier, from the same 81,514-target evidence table:

| Tier | n (with hit) | Median %identity | Median qcov | Median tcov | Median alignment length | Median candidate length | Median reference length | Median length ratio (candidate/reference) |
|---|---|---|---|---|---|---|---|---|
| Catalytic triad complete | 47,251 | 53.2% | 95.4% | 89.5% | 389 aa | 426 aa | 561 aa | 1.00 |
| HMM-supported (no triad) | 5,293 | 48.0% | 94.2% | 52.3% | 292 aa | 355 aa | 559 aa | 0.61 |
| ≥5 other PHA genes | 4,873 | 27.7% | 81.8% | 64.7% | 351 aa | 399 aa | 428 aa | 0.92 |
| 1–4 other PHA genes | 23,207 | 26.2% | 82.9% | 78.6% | 329 aa | 388 aa | 383 aa | 0.99 |
| phaC only | 890 | 25.7% | 84.3% | 77.2% | 326 aa | 386 aa | 397 aa | 0.98 |

(No-hit-at-all rate per tier, for completeness: 0.3% / 0.6% / 1.9% / 4.2% / 3.7%, same order.)

Two things stand out beyond the identity/coverage split itself. First, the **length ratios for the three unresolved tiers are all close to 1.0** (0.92–0.99) — candidate and matched reference are essentially the same length, which is exactly what you'd expect from a real, if divergent, full-length homolog and not from a short domain or motif match riding on a much longer reference. Second, the **HMM-supported-but-no-triad tier is the one genuine outlier**, with a length ratio of only 0.61 (candidate median 355 aa vs. reference median 559 aa) and a correspondingly low target coverage (52.3%, well below its own 94.2% query coverage) — this tier's candidates are themselves usually *shorter* than the references they match, which is the direct, mechanical explanation for why the catalytic triad wasn't resolvable for them (§3): the alignment often doesn't reach far enough into the C-terminal region to test the Asp/His positions, not because the protein is unrelated.

The broader, 1,875-sequence reference set changes the picture from earlier, narrower comparisons: even the pathway-only tiers (no HMM, no triad) reach **~82-84% median query coverage**, not the ~7-10% seen against the old 14-sequence trusted set. With enough reference diversity, most candidates *do* find a substantially full-length match somewhere in the confirmed-genuine phaC population — the gap between tiers is almost entirely in **identity** (~48-53% for HMM/triad-supported vs. ~26-28% for pathway-only), not coverage. That's a materially different, more specific signature than "fragment-only match": it looks like broad but distant homology, consistent with the working interpretation elsewhere in this document — divergent-but-real phaC, not spurious short alignments. Coverage below 20% (the "50% identity, 8% coverage" weak-evidence case named directly on the plot) is rare in every tier (0.1-0.9%), confirming that pattern is the exception, not the norm, across the whole verified dataset.

**Files:** `catalytic_domain/phac_sequence_evidence.tsv` (82,846 rows), `figures/scripts/plot_phac_pident_coverage.py`, `figures/scripts/plot_phac_evidence_violins.py`, `catalytic_domain/audited_search/` (mmseqs2 databases/results), `/tmp/audited_phac_refs.faa` (1,875-sequence audited reference FASTA, fetched fresh from UniProt).
## 5. Ecology, geography, and taxonomy: the full figure set, rebuilt

Every figure below was originally built earlier in this project's history, before the reference-query audit (§1.4). All are rebuilt here against the verified 68,424-genome phaC set — either automatically (most of them import `figures/scripts/_phac_qc.py`, now updated to the full 67-accession exclusion list, so rerunning the same script picks up the correction with no other changes) or, for a few, by regenerating an intermediate file first. Two findings below are not just "the same result on cleaner data" — they changed materially, and are flagged as such.

### 5.1 Taxonomic host-specificity

**UpSet plot** (`figures/phaC_upset.png`) — gene-family co-occurrence pattern across all 15 PHA families, from the corrected `architecture_summary.tsv` (§2). Top phaC-negative combination is still `ABJZFRReg` (595 genomes) — a genome with the core biosynthesis/mobilization/regulation genes but no synthase, the "1-4 other genes" pattern from §3 recurring at the whole-pathway level.

**Host-phylum specificity, radial plots** (`figures/phac_taxonomy_radial.png`, `_detailed.png`) — for each phaC cluster (70%-identity mmseqs grouping), is it found in genomes from only one bacterial/archaeal phylum, or does it span multiple? This is the figure that changed the most dramatically of the whole batch: **96.8% of clusters (11,396 of 11,775) are now single-phylum, up from the pre-audit 83%.** Extreme cross-phylum spread turns out to have been disproportionately a contamination signature — a cluster built partly or wholly from mislabeled antiporter/isochorismate-synthase/etc. hits will spuriously look "promiscuous" because those contaminating gene families are broadly distributed across the tree of life, unlike real phaC lineages which tend to track host phylogeny more closely. 92 distinct phyla appear in the detailed version, all shown individually (no "other" bucket).

**Promiscuous-cluster deep dives** (`figures/promiscuous_cluster_*_radial.png`) — the two clusters singled out pre-audit as the most taxonomically promiscuous (`...131018176`: 244 genera/10 phyla; `...095227834`: 198 genera/15 phyla) turned out to be almost entirely the contamination just described: `...131018176` collapsed to 5 genomes/4 genera/1 phylum once wrong-gene targets were removed, and `...095227834` **disappeared completely** (zero genomes left). Re-derived the current top-2 by distinct-phyla count fresh from the corrected data: `...246448549` (110 genomes, 43 genera, 9 phyla, Pseudomonadota-dominated, includes both REDSEA-S09-B13 and UBA868 among its host genera) and `...022439766` (230 genomes, 24 genera, 8 phyla, Thermoproteota/archaea-dominated, Nitrosopelagicus-heavy). Real promiscuous clusters do exist — just far less extreme (max 9 phyla, not 15) and different ones than before.

### 5.2 Geography

**Global prevalence heatmap** (`figures/phac_pct_global_heatmap.png`) — % of genomes carrying phaC, gridded at 2°×2° over every location OMDB screened (not just phaC-positive samples, so this is a true rate, not a density map). 254,872 genomes with usable coordinates, 23.7% phaC-positive overall; 815 grid cells shown after the ≥15-genomes/≥1-sample completeness filter.

**Regional specialists** (`figures/phaC_regional_specialists.png`) — clusters geographically concentrated in one place (`geo_mean_resultant_length > 0.85`, a circular-statistics tightness measure) rather than cosmopolitan. Of the 16 originally hand-picked clusters, **5 no longer exist at all post-audit** — their entire membership was wrong-gene contamination: the second Sydney-AU (SHLQ01) cluster, the S. Atlantic hydrothermal-vent Thermococcus/Pyrococcus cluster, the tropical N. Pacific Crocosphaera cluster, the Baltic Sea BACL27 cluster, and the S. Pacific vent Hydrogenivirga cluster. The remaining 11 all still pass the same tightness threshold and are plotted — spanning Sydney AU, NE Pacific abyssal (Moritella), sub-Arctic and true-high-Arctic Davis Strait/near-pole (UBA4427, two separate clusters), Norwegian shelf (UBA4582), Philippine Sea (JACZQZ01), mid N. Atlantic (Robiginitomaculum_A), Hawaii abyssal (Rs1), Southern Ocean (Algiphilus), Amazon-influenced W. Atlantic (Nanopelagicaceae/UBA7398), and mid-depth N. Pacific (UBA8229). No new candidates were curated in to replace the 5 lost ones — a fresh candidate pass would be needed to restore full geographic coverage if that's wanted later.

**Global clusters, re-picked** (`figures/phaC_marinobacter_global.png`, `phaC_sulfitobacter_global.png`, `figures/phaC_pseudomonas_e_global.png`) — UBA868 (`figures/phaC_uba868_global.png`, `_overlay.png`, kept for reference) was the original focus organism for "which genus has multiple genuinely-global phaC clusters," picked earlier because it dominated the two clusters flagged as most taxonomically promiscuous (§5.1) — but since both of those turned out to be almost entirely contamination, UBA868 itself is a weak choice post-audit: only 6 qualifying global clusters (≥30 genomes, geographic concentration R<0.4, max pairwise spread >10,000km), 711 total genomes.

Re-scanned every genus in the corrected `phaC_cluster0.7_cluster_ecology.tsv` against the same criteria and found three real, well-characterized marine genera with substantially more/bigger genuinely-global clusters:

| Genus | Qualifying global clusters | Total genomes | Character |
|---|---|---|---|
| Marinobacter | 9 | 1,587 | Most uniform spread of the three — R as low as 0.08–0.13 in several clusters, not just far-flung hotspots |
| Sulfitobacter | 6 | 2,931 | Fewest clusters but by far the largest scale — top cluster alone has 1,034 genomes |
| Pseudomonas_E | 9 | 2,034 | Ties Marinobacter for cluster count, largest total genome count |

All three small-multiples figures were built with `figures/scripts/plot_genus_global_clusters.py <genus>` (generalized from the UBA868-specific script). **Sulfitobacter was selected as the headline pick** — single-overlay version built with `figures/scripts/plot_genus_global_overlay.py Sulfitobacter` (`figures/phaC_sulfitobacter_overlay.png`, generalized the same way from the UBA868 overlay script): the 6 Sulfitobacter-dominant clusters with >200 genomes (of 39 total) plotted together on one map, one color/marker per cluster, largest-cluster-first so smaller clusters stay visible on top. All 6 clusters visibly span every major ocean basin with real overlap between them, not just each independently far-flung — the actual "shares a range" story the figure is named for.

### 5.3 Habitat

**Prevalence by ocean habitat** (`figures/phac_pct_by_ocean_habitat.png`) — % phaC-positive per habitat category (≥200 genomes each), true-denominator approach (same genome universe as the global map). 61,130/257,321 genomes phaC-positive overall (23.8%). Clear host-association enrichment at the top: whale-fall bone biofilm (78.5%), estuarine sediment (62.1%), hydrozoa tissue (58.2%), sea ice (57.4%), algae tissue (55.0%), coral tissue (51.3%) — all well above the open-water baseline (21.5%). Sponge tissue sits at 38.2% here (dataset-wide, all evidence tiers combined) — for the tier-specific sponge-enrichment finding (which *did* change substantially post-audit), see §3's habitat discussion.

### 5.4 Depth

**REDSEA-S09-B13 case study** (`figures/redsea_depth_clusters.png`) — does this specific genus's phaC cluster membership shift with depth? Its largest cluster (`...123205028`, 89 genomes) spans the full range this genus was sampled at, from 0m to 5,100m, with no obvious depth-driven cluster turnover — the same cluster persists from surface to deep.

**Hadal trench depths** (`figures/hadal_cluster_depth_span.png`, `hadal_genera_depth_clusters.png`) — is phaC found at hadal-trench depths (>6000m) genuinely distinct, or just a depth-extreme sample of shallower lineages? **29 verified genomes** reach hadal depth (up from 27 pre-audit — this side of the analysis barely moved, consistent with phaC's own reference set having checked out clean in the audit). These 29 fall into **15 distinct 70%-identity clusters** across 10 genera (Nitrosopelagicus, Nitrosopumilus, DTSX01, NORP23, DUCF01, UBA9611, Casp-alpha2, JAYZAP01, Phenylobacterium, UBA9659, GCA-002726655). The genera-vs-cluster figure makes the key point directly: where a genus's hadal and shallow-water members share the same color (same cluster), hadal phaC is not a novel lineage — it's the same cluster reaching down to extreme depth. This rebuild required reconstructing the underlying all-vs-all self-search from scratch (§1.4/§4's audited-reference search doesn't answer this — it's phaC-vs-phaC, not phaC-vs-reference; see `figures/scripts/build_selfsearch_derived_tables.py`, run against `catalytic_domain/selfsearch/` , an mmseqs2 all-vs-all among the 82,846 verified target sequences).

**Genomes whose synthase is highly divergent from the rest of the database** (`figures/phac_divergent35_vs_depth.png`) — same all-vs-all self-search, restricted to genomes whose best non-self hit is <35% identity anywhere in the verified dataset. **284 genome-level candidates** (vs. 280 pre-audit), **222 distinct target_ids** (vs. 221) that form **222 singleton clusters** under the same phaC_cluster0.7 settings — i.e. still independently divergent, not one hidden lineage, spanning **203 genera across 31 phyla** (vs. 202/35 pre-audit). This result was essentially unchanged by the audit, for the same reason as the hadal figures: phaC's own reference set didn't need correcting.

### 5.5 Temperature

**WOA23 annual temperature** (`figures/phaC_woa_temperature_depth.png`) — genomes joined to NOAA WOA23 1°-grid annual climatological temperature at their sample depth. 26,131 of 26,206 attempted joins matched (75 unmatched — no nearby WOA grid cell with data at that depth). 2,608 of 25,831 plotted rows (10.1%) come from sub-zero-annual-temperature water.

**TemStaPro stability vs. ocean temperature** (`figures/phaC_temstapro_vs_ocean_temperature.png`) — per-sequence thermostability predictions (from TemStaPro, a protein language-model-based stability classifier — predictions themselves don't change with which genomes are now excluded, since they're a function of sequence alone; only the metadata join and which points get plotted changed here) against each genome's local ocean temperature. 12,150 rows plotted after QC (45,154 removed — the majority of the raw prediction set, since TemStaPro was originally run on a broader sequence pool than just the verified-phaC / depth-resolved subset this figure plots).

## 6. Structure prediction: resolving the tiers with no direct sequence evidence

§3 established five evidence tiers for phaC. Two are already confident (catalytic triad complete, 42.0%; HMM-supported, 5.5%) — direct protein-level evidence, no further work needed. The other three (≥5 other genes 6.8%, 1-4 other genes 43.9%, phaC only 1.8% — 52.5% of the whole verified dataset) have no direct evidence either way: real pathway context in the first two, nothing at all in the third. Sequence-only methods (HMM, mmseqs homology, catalytic-triad projection) have been run against these already (§3-4) and can't resolve them further — the next available independent check is structure: does the *folded* protein look like a real PhaC catalytic domain, regardless of what its sequence identity says.

### 6.1 Scale and approach

Folding all 68,424 verified genomes individually would mean predicting ~48,000 structures for the three unresolved tiers alone — a real undertaking even on a GPU cluster. Deduping to one representative sequence per phaC_cluster0.7 cluster (70%-identity groups) cuts this to a much more tractable **4,992 sequences** (one per cluster touched by an unresolved-tier genome, length-filtered to 150-700aa to exclude obvious fragments) — clusters, not genomes, are the natural unit here anyway, since near-identical sequences within a cluster would fold near-identically. Combined with a **428-sequence positive-control sample** (stratified by phylum, drawn from clusters touched only by triad-complete/HMM-supported genomes — confirms the whole ESMFold→Foldseek approach actually recovers the right answer on cases already known to be real, before trusting it on the uncertain ones) and the **1,875-sequence audited reference set** (§1.4/§4 — these get folded too, to serve as the Foldseek comparison targets, built from the same verified pool as everything else in this project rather than an arbitrary external structure database), the total scope is **~7,300 structures**. That's a reasonable scale for ESMFold on a single GPU node — not an AlphaFold2/ColabFold-scale undertaking, which is exactly why ESMFold (a single forward pass, no MSA search, seconds-to-low-minutes per structure) is the right first tool rather than starting with the slower, MSA-dependent alternatives.

### 6.2 Pipeline

- `structure_prediction/build_fold_inputs.py` — builds the three FASTA sets described above from the same verified-tier data as §3 (`/tmp/phac_verified_triad_hmm_group.pkl`, `/tmp/verified_phac_genome_targets.pkl`) plus the cluster assignment table. Deterministic (fixed random seed) so the positive-control sample is reproducible.
- `structure_prediction/run_esmfold.py` — the folding driver. Uses HuggingFace `transformers`' `EsmForProteinFolding` rather than Meta's original `fair-esm[esmfold]` package deliberately: the original needs `openfold`'s custom CUDA attention kernels compiled against the exact local CUDA/PyTorch build, a common source of cluster-install pain; the `transformers` port reimplements the same model in plain PyTorch with no custom kernel step. Writes one PDB file per sequence plus a summary TSV (length, mean pLDDT, pTM, status). Resumable — skips any target_id whose PDB already exists, so a retried/re-submitted array task only fills gaps.
- `cluster/install_esmfold.sh` — environment setup. Adds ESMFold's packages (torch, transformers, accelerate, einops) to the *existing* `pha-reference` env (the one `source cluster/env_activate.sh` activates) rather than creating a separate one — only 4 new packages, not enough to justify a second env and a second activation step in the normal workflow. Model-weight cache (`$HF_HOME`, ~2.7GB) and the fallback dedicated-env path (only if the shared env ever hits a real conflict) both default to `/scratch/morrill/users/hmp278`, not `$HOME` — avoids blowing a small home-directory quota.
- `cluster/run_esmfold.sbatch` — SLURM array job, one task per batch slice of one input FASTA; run once per input file (uncertain / positive-control / reference). Sources `cluster/env_activate.sh` directly (same env as everything else in this project) and sets `$HF_HOME` to the same scratch path.

**Honest caveat:** none of this has been run end-to-end on an actual GPU node yet — written from documented `transformers`/ESMFold usage patterns, not confirmed live (this session has no direct cluster access, per the established pattern all session of writing/pushing code for the user to run and report back errors on). Treat the first real submission as a shakedown run, not a sure thing — watch for install friction in `install_esmfold.sh` and tune `N_BATCHES`/chunk-size in `run_esmfold.sbatch` once real per-sequence timing is in hand.

### 6.3 Next steps (not yet started)

1. `source cluster/env_activate.sh` then run `cluster/install_esmfold.sh` once, then submit `cluster/run_esmfold.sbatch` three times (once per input FASTA).
2. Once structures land: build a Foldseek database from the folded reference set, then search the uncertain-tier and positive-control structures against it (structural alignment, not sequence) — this is the actual "does it look like a real PhaC" test.
3. Positive controls should score well against the reference Foldseek DB if the pipeline is working; if the uncertain tiers score comparably, that's real structural evidence for divergent-but-genuine phaC (a third, independent line of evidence beyond §3's HMM/triad and §4's sequence homology). If they score poorly despite the positive controls succeeding, that's real evidence some of the "1-4 other genes"/"≥5 other genes" tier is not phaC after all.
4. For any candidates ESMFold handles poorly (very low pLDDT, likely for the most divergent/unusual-length sequences) or where the Foldseek verdict is ambiguous, escalate to ColabFold/AlphaFold2 (MSA-based, slower but often more accurate for divergent sequences with few close homologs) for that specific subset only — not the whole 7,300, which is why doing the cheap ESMFold pass first matters.

**Files:** `structure_prediction/build_fold_inputs.py`, `run_esmfold.py`, `uncertain_cluster_representatives.faa` (4,992), `positive_control_representatives.faa` (428), `audited_reference_set.faa` (1,875), `fold_manifest.tsv`; `cluster/install_esmfold.sh`, `cluster/run_esmfold.sbatch`.
