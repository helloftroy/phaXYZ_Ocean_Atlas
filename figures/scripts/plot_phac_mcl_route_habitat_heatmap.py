"""mcl-PHA precursor-supply route, focused view: log2 enrichment of each
route by habitat, as a heatmap -- a zoomed-in companion to Panel C of
plot_phac_mcl_precursor_strategy.py (figures/phac_mcl_precursor_strategy.png),
which shows the same routes as stacked-percentage bars. The stacked-bar view
makes within-habitat composition easy to read but enrichment/depletion
relative to baseline hard to see at a glance across many habitats at once;
this heatmap trades that for the opposite: exact per-habitat composition is
gone, but "which habitats favor which route, and by how much" is immediate.

Four columns, each a per-gene presence check (not mutually exclusive, unlike
section 8's Neither/J-only/G-only/Both quadrant split):
  phaAB   -- canonical scl-PHA route present (phaA AND phaB)
  phaJ    -- beta-oxidation-linked mcl route present (phaJ), regardless of G
  phaG    -- FAS-linked mcl route present (phaG), regardless of J
  Both J&G -- both mcl routes present together

Baseline is the dataset-wide rate across all phaC-positive genomes (not all
screened marine genomes) -- these four columns are gene-content questions
that are only meaningful for a genome that already carries phaC (see
plot_phac_mcl_precursor_strategy.py's own docstring on this), so "enriched"
here means "more common in this habitat than in the average phaC-positive
genome," not "more common than in seawater generally."

Cell value: log2((habitat rate + 0.5) / (habitat n + 1)) minus the same
Haldane-Anscombe-corrected log2 baseline rate -- the +0.5/+1 pseudocount
avoids -inf for the (rare) habitat/route combination with zero positives,
matching a standard small-count correction rather than leaving those cells
undefined.

Rows are ordered by average-linkage hierarchical clustering on each
habitat's 4-column enrichment vector, so habitats with similar route
profiles land next to each other instead of in an arbitrary or purely
alphabetical order.

Significance: same phylum-controlled Cochran-Mantel-Haenszel test as
plot_phac_mcl_precursor_strategy.py's habitat table (via _stats_utils.py),
run independently for each of the four route columns -- marked on the
heatmap as */** for p<0.05/p<0.001.

Usage:
    python figures/scripts/plot_phac_mcl_route_habitat_heatmap.py

Outputs:
    figures/phac_mcl_route_habitat_heatmap.png / .pdf
    figures/phac_mcl_route_habitat_heatmap.tsv
"""
import csv
import math
import sys
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stats_utils import mantel_haenszel

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

EXCLUDE_HABITATS = {  # keep in sync with plot_phac_mcl_precursor_strategy.py
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}
MIN_HABITAT_N = 200   # min genomes screened in a habitat overall, for it to be considered at all
MIN_PHAC_N = 20        # min phaC-positive genomes in a habitat, for a stable route mix to be reported

ROUTES = ['phaAB', 'phaJ', 'phaG', 'Both J&G']


def route_flags(rec):
    return {
        'phaAB': rec['ab'],
        'phaJ': rec['j'],
        'phaG': rec['g'],
        'Both J&G': rec['j'] and rec['g'],
    }


# ---------------------------------------------------------------------
# 1. per-genome A/B/G/J presence, phylum -- phaC-positive genomes only
#    (identical loading to plot_phac_mcl_precursor_strategy.py)
# ---------------------------------------------------------------------
genomes = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if int(row['n_phaC']) <= 0:
            continue
        genomes[row['genome']] = {
            'ab': int(row['n_phaA']) > 0 and int(row['n_phaB']) > 0,
            'g': int(row['n_phaG']) > 0,
            'j': int(row['n_phaJ']) > 0,
            'phylum': row['gtdb_phylum'],
        }

n_total = len(genomes)
print(f'{n_total:,} phaC-positive genomes')

baseline_rate = {}
for route in ROUTES:
    n_pos = sum(1 for r in genomes.values() if route_flags(r)[route])
    baseline_rate[route] = (n_pos + 0.5) / (n_total + 1)
    print(f'  baseline {route:9s} {n_pos:6,d}/{n_total:,}  ({100*n_pos/n_total:.2f}%)')

# ---------------------------------------------------------------------
# 2. habitat labels (same source/gating as plot_phac_mcl_precursor_strategy.py)
# ---------------------------------------------------------------------
genome_habitat = {}
habitat_totals_all = defaultdict(int)
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        c = row['ecosystem_compartment']
        if c in EXCLUDE_HABITATS:
            continue
        habitat_totals_all[c] += 1
        if row['genome'] in genomes:
            genome_habitat[row['genome']] = c

valid_habitats = [c for c, n in habitat_totals_all.items() if n >= MIN_HABITAT_N]

# ---------------------------------------------------------------------
# 3. per habitat x route: log2 enrichment + phylum-controlled CMH test
# ---------------------------------------------------------------------
rows_out = []
enrichment = {}   # habitat -> {route: log2 enrichment}
sig = {}           # habitat -> {route: p-value}

for hab in valid_habitats:
    hab_genomes = [g for g, h in genome_habitat.items() if h == hab]
    if len(hab_genomes) < MIN_PHAC_N:
        continue
    recs = [genomes[g] for g in hab_genomes]
    n = len(recs)

    enrichment[hab] = {}
    sig[hab] = {}
    row = {'habitat': hab, 'n_phac_positive': n}

    for route in ROUTES:
        n_pos = sum(1 for r in recs if route_flags(r)[route])
        hab_rate = (n_pos + 0.5) / (n + 1)
        log2_enr = math.log2(hab_rate / baseline_rate[route])
        enrichment[hab][route] = log2_enr

        # phylum-controlled CMH test: is this route's presence associated
        # with being in this habitat (vs. the rest of the phaC-positive,
        # habitat-labeled dataset) once phylum is controlled for -- same
        # machinery as plot_phac_mcl_precursor_strategy.py's phaG test,
        # just run for all four route columns instead of only phaG.
        phylum_in = defaultdict(lambda: [0, 0])
        phylum_out = defaultdict(lambda: [0, 0])
        for g, h in genome_habitat.items():
            rec = genomes[g]
            bucket = phylum_in if h == hab else phylum_out
            bucket[rec['phylum']][0 if route_flags(rec)[route] else 1] += 1
        strata = []
        for p in set(phylum_in) | set(phylum_out):
            ia, ib = phylum_in.get(p, [0, 0])
            oa, ob = phylum_out.get(p, [0, 0])
            strata.append((ia, ib, oa, ob))
        or_mh, cmh_stat, cmh_p, n_strata = mantel_haenszel(strata)
        sig[hab][route] = cmh_p

        row[f'pct_{route}'] = f'{100*n_pos/n:.2f}'
        row[f'log2_enrichment_{route}'] = f'{log2_enr:.3f}'
        row[f'cmh_or_{route}'] = f'{or_mh:.3f}' if or_mh is not None else ''
        row[f'cmh_p_{route}'] = f'{cmh_p:.3g}' if cmh_p is not None else ''

    rows_out.append(row)

rows_out.sort(key=lambda r: -r['n_phac_positive'])
out_tsv = OUT / 'phac_mcl_route_habitat_heatmap.tsv'
with open(out_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(rows_out)
print(f'\nwrote {out_tsv} ({len(rows_out)} habitats)')

# ---------------------------------------------------------------------
# 4. row order via hierarchical clustering on each habitat's 4-route vector
# ---------------------------------------------------------------------
habitats = [r['habitat'] for r in rows_out]
mat = np.array([[enrichment[h][route] for route in ROUTES] for h in habitats])
if len(habitats) >= 3:
    order = leaves_list(linkage(mat, method='average', metric='euclidean'))
else:
    order = np.arange(len(habitats))
habitats_ordered = [habitats[i] for i in order]
mat_ordered = mat[order]

n_by_hab = {r['habitat']: r['n_phac_positive'] for r in rows_out}

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig_h = 2.0 + 0.42 * len(habitats_ordered)
fig, ax = plt.subplots(figsize=(9.6, fig_h), dpi=300)

# color normalized PER COLUMN, not globally: phaAB/phaJ have high baseline
# prevalence (85.6%/65.2%) so their log2 enrichment is mathematically capped
# near +-0.6, while phaG/Both J&G have low baseline prevalence (5.5%/4.4%)
# and can swing several log2 units on a handful of genomes -- a shared color
# scale would make the first two columns wash out to near-white even where
# they carry real, significant signal (e.g. phaJ depletion in estuarine/
# brackish water). Cell text always shows the true log2 value regardless.
col_vmax = np.abs(mat_ordered).max(axis=0)
col_vmax[col_vmax == 0] = 1.0
mat_norm = mat_ordered / col_vmax
im = ax.imshow(mat_norm, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

ax.set_xticks(range(len(ROUTES)))
ax.set_xticklabels(ROUTES, fontsize=11, fontweight='bold')
ax.set_xticks(np.arange(-0.5, len(ROUTES), 1), minor=True)
ax.set_yticks(range(len(habitats_ordered)))
ax.set_yticklabels([f'{h} (n={n_by_hab[h]:,})' for h in habitats_ordered], fontsize=9.3)
ax.set_yticks(np.arange(-0.5, len(habitats_ordered), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=1.6)
ax.tick_params(which='minor', length=0)
ax.tick_params(which='major', length=0)
for spine in ax.spines.values():
    spine.set_visible(False)

for i, hab in enumerate(habitats_ordered):
    for j, route in enumerate(ROUTES):
        val = mat_ordered[i, j]
        p = sig[hab][route]
        stars = '**' if (p is not None and p < 0.001) else ('*' if (p is not None and p < 0.05) else '')
        txt_color = 'white' if abs(mat_norm[i, j]) > 0.55 else '#20302C'
        ax.text(j, i, f'{val:+.2f}{stars}', ha='center', va='center', fontsize=8.6,
                color=txt_color, fontweight='bold' if stars else 'normal')

cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03, ticks=[-1, 0, 1])
cbar.ax.set_yticklabels(['most\ndepleted', '0', 'most\nenriched'], fontsize=8)
cbar.set_label('color: rank within column only\n(read log$_2$ values in cells)', fontsize=8.6)

ax.set_title('mcl-PHA precursor route enrichment by habitat', fontsize=14, fontweight='bold', pad=14, loc='left')
fig.text(0.02, 0.99,
          f'Cell text = log2(habitat rate / dataset-wide phaC-positive rate), Haldane-Anscombe corrected. '
          f'Color is normalized per column (not globally) since phaAB/phaJ have high baseline prevalence '
          f'(85.6%/65.2%) and so are mathematically capped near +-0.6, while phaG/Both J&G (baseline '
          f'5.5%/4.4%) can swing several log2 units -- a shared scale would wash out the first two columns.\n'
          f'{len(habitats_ordered)} habitats with ≥{MIN_PHAC_N} phaC-positive genomes shown. '
          f'* p<0.05, ** p<0.001 (phylum-controlled CMH test). Rows ordered by hierarchical clustering of route profile.',
          transform=fig.transFigure, ha='left', va='top', fontsize=7.6, color='#5B6E70', wrap=True)

fig.subplots_adjust(left=0.30, right=0.87, top=0.87, bottom=0.03)

out_path = OUT / 'phac_mcl_route_habitat_heatmap.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_mcl_route_habitat_heatmap.pdf', facecolor='white')
print('saved pdf too')
