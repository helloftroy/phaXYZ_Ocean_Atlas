"""Phylum-controlled habitat enrichment test for phaC prevalence.

plot_phac_pct_by_ocean_habitat.py showed several habitats (whale-fall bone
biofilm, estuarine sediment, hydrozoa tissue, sea ice, algae tissue, coral
tissue, ...) sitting well above the dataset-wide phaC-positive rate -- but
that is a raw percentage with no statistical test and, more importantly, no
control for taxonomic composition. If e.g. "Animal bone biofilm" happens to
be dominated by a phylum that is independently phaC-rich everywhere, the
habitat's high raw rate could be entirely a composition artifact rather than
a real habitat effect. This script asks: does each habitat's enrichment
survive once you control for which phyla live there?

Three complementary tests are run per habitat (each answers a slightly
different question, deliberately not collapsed into one number):

  1. Raw, unstratified Fisher's exact test (habitat vs rest of dataset) --
     reproduces the existing figure's percentages as a Fisher OR/p-value, for
     reference against the phylum-controlled numbers below.
  2. Cochran-Mantel-Haenszel (CMH) test, stratified by phylum -- the actual
     phylum-controlled significance test. Combines habitat-vs-rest evidence
     WITHIN each phylum into one pooled odds ratio and p-value. If the CMH
     odds ratio/p-value stays strong while phylum composition varies a lot
     between the habitat and the rest of the ocean, that is real evidence
     the habitat effect is not just a composition artifact.
     Formula: standard Mantel-Haenszel combined OR and chi-square (Mantel &
     Haenszel 1959); implemented in `_stats_utils.mantel_haenszel` (shared
     with plot_phac_mcl_precursor_strategy.py's phaG/phaJ-vs-habitat test)
     rather than via statsmodels (not installed everywhere this project's
     scripts run) -- see that module for the formula, short enough to
     audit by eye. Here, a_i/b_i/c_i/d_i are the phaC-positive/negative x
     in-habitat/outside-habitat counts within phylum stratum i.
  3. Phylum-standardized expected rate (direct standardization): what phaC
     rate would this habitat show if each of its phyla had the SAME phaC
     rate they show everywhere else in the ocean, weighted by this
     habitat's own actual phylum mix? Observed-minus-expected is the most
     intuitive single number for "how much of the raw enrichment survives
     composition-adjustment" -- reported alongside the CMH test, not
     instead of it.

A separate per-phylum breakdown TSV is also written (habitat x phylum) as a
robustness check: does the direction/significance of enrichment hold up
within each individual major phylum, not just in the pooled CMH number.

Coverage caveat (real, not hypothetical -- see plot_phac_pct_by_ocean_habitat.py's
own docstring): genome_family_matrix.tsv (the only file with both phylum and
genome ID) is missing ~31,000 genomes present in the fresh genome universe
(omdb_all_genomes_with_locations.tsv), entirely missing the "Cold seep
sediment" category, and has an older/mismatched ecosystem_compartment column
that MUST NOT be used for habitat labels (that mismatch already burned one
earlier figure -- see the habitat script's docstring). This script only
trusts genome_family_matrix.tsv for gtdb_phylum, joined strictly by genome
ID; habitat and the phaC-positive genome set both come from the same fresh
files the habitat figure uses. Genomes with no phylum match are excluded
from tests 2 and 3 (impossible to stratify by an unknown phylum) but are
still counted in test 1 and reported explicitly as a coverage number per
habitat, not silently dropped.

Usage:
    python figures/scripts/phac_habitat_phylum_controlled_test.py

Outputs:
    figures/phac_habitat_phylum_enrichment_test.tsv       -- one row per habitat
    figures/phac_habitat_phylum_enrichment_breakdown.tsv  -- one row per habitat x phylum
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc
from _stats_utils import mantel_haenszel

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

# same denylist as plot_phac_pct_by_ocean_habitat.py -- keep in sync
EXCLUDE = {
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}
MIN_N = 200          # same threshold as the habitat figure: habitats smaller than this aren't tested
MIN_STRATUM_N = 1    # strata with n_i<=1 are skipped (undefined CMH variance)
MIN_PHYLUM_BREAKDOWN_N = 20  # per-phylum robustness rows below this are still written but flagged

# ---- phaC-positive genome IDs (same QC filter as the habitat figure) ----
phac_genome_ids = set()
with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if _phac_qc.is_bad(row['best_query']):
            continue
        phac_genome_ids.add(row['genome'])

# ---- phylum per genome (only source with genome+phylum together) ----
genome_phylum = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        genome_phylum[row['genome']] = row['gtdb_phylum']

# ---- genome universe + habitat label, from the fresh canonical file ----
# genome_id -> (habitat, phaC_positive: bool, phylum: str or None)
genomes = []
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        c = row['ecosystem_compartment']
        if c in EXCLUDE:
            continue
        g = row['genome']
        genomes.append((g, c, g in phac_genome_ids, genome_phylum.get(g)))

n_total = len(genomes)
n_phac_total = sum(1 for _, _, p, _ in genomes if p)
n_with_phylum = sum(1 for _, _, _, ph in genomes if ph is not None)
overall_pct = 100 * n_phac_total / n_total
print(f'{n_phac_total:,} / {n_total:,} genomes phaC-positive overall ({overall_pct:.1f}%)')
print(f'{n_with_phylum:,} / {n_total:,} genomes ({100*n_with_phylum/n_total:.1f}%) have a phylum match '
      f'via genome_family_matrix.tsv -- only these contribute to the phylum-controlled tests below.')

by_habitat = defaultdict(list)
for g, c, p, ph in genomes:
    by_habitat[c].append((p, ph))

habitats = [c for c, rows in by_habitat.items() if len(rows) >= MIN_N]
print(f'{len(habitats)} habitats with >= {MIN_N} genomes\n')


habitat_rows = []
breakdown_rows = []

for hab in habitats:
    rows = by_habitat[hab]
    n_hab = len(rows)
    phac_hab = sum(1 for p, _ in rows if p)
    pct_hab = 100 * phac_hab / n_hab
    n_hab_with_phylum = sum(1 for _, ph in rows if ph is not None)

    # ---- test 1: raw, unstratified Fisher (habitat vs rest of whole dataset) ----
    a = phac_hab
    b = n_hab - phac_hab
    c = n_phac_total - phac_hab
    d = (n_total - n_hab) - c
    raw_or, raw_p = fisher_exact([[a, b], [c, d]], alternative='two-sided')

    # ---- build per-phylum strata, restricted to genomes with a phylum match ----
    # phylum -> [a,b,c,d] within that phylum: in-habitat vs outside-habitat, both restricted to phylum-known genomes
    phylum_counts_hab = defaultdict(lambda: [0, 0])   # phylum -> [phac+, phac-] within habitat
    phylum_counts_out = defaultdict(lambda: [0, 0])   # phylum -> [phac+, phac-] outside habitat (same dataset)
    for g, c_, p, ph in genomes:
        if ph is None:
            continue
        bucket = phylum_counts_hab if c_ == hab else phylum_counts_out
        bucket[ph][0 if p else 1] += 1

    all_phyla = set(phylum_counts_hab) | set(phylum_counts_out)
    strata = []
    for ph in all_phyla:
        ha, hb = phylum_counts_hab.get(ph, [0, 0])
        oa, ob = phylum_counts_out.get(ph, [0, 0])
        strata.append((ha, hb, oa, ob))
        n_ph_hab = ha + hb
        n_ph_out = oa + ob
        if n_ph_hab == 0:
            continue  # phylum not present in this habitat at all -- nothing to report for it here
        pct_ph_hab = 100 * ha / n_ph_hab if n_ph_hab else float('nan')
        pct_ph_out = 100 * oa / n_ph_out if n_ph_out else float('nan')
        if n_ph_hab >= 2 and n_ph_out >= 2:
            ph_or, ph_p = fisher_exact([[ha, hb], [oa, ob]], alternative='two-sided')
        else:
            ph_or, ph_p = float('nan'), float('nan')
        breakdown_rows.append({
            'habitat': hab, 'phylum': ph,
            'n_in_habitat': n_ph_hab, 'phac_in_habitat': ha, 'pct_in_habitat': f'{pct_ph_hab:.1f}',
            'n_outside_habitat': n_ph_out, 'phac_outside_habitat': oa, 'pct_outside_habitat': f'{pct_ph_out:.1f}',
            'fisher_or': f'{ph_or:.3f}' if ph_or == ph_or else '', 'fisher_p': f'{ph_p:.3g}' if ph_p == ph_p else '',
            'low_n_flag': 'yes' if n_ph_hab < MIN_PHYLUM_BREAKDOWN_N else '',
        })

    or_mh, cmh_stat, cmh_p, n_strata = mantel_haenszel(strata)

    # ---- test 3: phylum-standardized expected rate ----
    # expected_pct = what this habitat's phaC rate would be if each of its own
    # phyla had the SAME phaC rate as that phylum shows OUTSIDE this habitat,
    # weighted by this habitat's own actual phylum composition.
    if n_hab_with_phylum > 0:
        expected_hits = 0.0
        for ph, (ha, hb) in phylum_counts_hab.items():
            n_ph_hab = ha + hb
            oa, ob = phylum_counts_out.get(ph, [0, 0])
            n_ph_out = oa + ob
            rate_out = (oa / n_ph_out) if n_ph_out > 0 else (n_phac_total / n_total)  # fall back to dataset baseline if phylum unseen outside habitat
            expected_hits += n_ph_hab * rate_out
        expected_pct = 100 * expected_hits / n_hab_with_phylum
        observed_pct_phylum_subset = 100 * sum(ha for ha, _ in phylum_counts_hab.values()) / n_hab_with_phylum
    else:
        expected_pct = observed_pct_phylum_subset = float('nan')

    habitat_rows.append({
        'habitat': hab, 'n_total': n_hab, 'n_phac': phac_hab, 'pct_raw': f'{pct_hab:.2f}',
        'overall_pct': f'{overall_pct:.2f}',
        'raw_fisher_or': f'{raw_or:.3f}', 'raw_fisher_p': f'{raw_p:.3g}',
        'n_with_phylum': n_hab_with_phylum, 'pct_phylum_coverage': f'{100*n_hab_with_phylum/n_hab:.1f}',
        'cmh_or': f'{or_mh:.3f}' if or_mh is not None else '',
        'cmh_chi2': f'{cmh_stat:.2f}' if cmh_stat is not None else '',
        'cmh_p': f'{cmh_p:.3g}' if cmh_p is not None else '',
        'n_strata_used': n_strata,
        'observed_pct_phylum_subset': f'{observed_pct_phylum_subset:.2f}' if observed_pct_phylum_subset == observed_pct_phylum_subset else '',
        'expected_pct_if_composition_only': f'{expected_pct:.2f}' if expected_pct == expected_pct else '',
    })

habitat_rows.sort(key=lambda r: -float(r['pct_raw']))

OUT.mkdir(parents=True, exist_ok=True)
habitat_out = OUT / 'phac_habitat_phylum_enrichment_test.tsv'
with open(habitat_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(habitat_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(habitat_rows)
print(f'wrote {habitat_out} ({len(habitat_rows)} habitats)')

breakdown_out = OUT / 'phac_habitat_phylum_enrichment_breakdown.tsv'
with open(breakdown_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(breakdown_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(breakdown_rows)
print(f'wrote {breakdown_out} ({len(breakdown_rows)} habitat x phylum rows)')

print()
print(f'{"habitat":32s} {"raw %":>7s} {"raw p":>10s} {"cmh OR":>8s} {"cmh p":>10s} {"expected %":>11s} {"strata":>7s}')
for r in habitat_rows:
    print(f'{r["habitat"]:32s} {r["pct_raw"]:>6s}% {r["raw_fisher_p"]:>10s} {r["cmh_or"]:>8s} '
          f'{r["cmh_p"]:>10s} {r["expected_pct_if_composition_only"]:>10s}% {r["n_strata_used"]:>7}')
