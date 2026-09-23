"""% of genomes carrying phaC, broken down by ocean habitat (ecosystem_compartment).

Earlier depth/location panels (plot_overview.py Panel C) only had COUNTS of
phaC-positive genomes -- no true denominator, since the phaC hit file alone
doesn't say how many genomes were screened per habitat. omdb_all_genomes_with_locations.tsv
(274,282 genomes, freshly pulled) is the full OMDB genome universe with
per-genome ecosystem_compartment + lat/lon, so it supplies that denominator
and lets us compute a real percentage, not just a raw count.

Restricted to genuine marine habitats -- OMDB also includes freshwater,
terrestrial, and lab-control samples (its "ecosystem" field is uninformative,
always "ocean" as a project label even for those), so those are excluded via
an explicit denylist rather than trusting any single column to mean "is this
the ocean."

IMPORTANT: numerator and denominator are joined by GENOME ID, not by habitat
label. The phaC hits file (phaC_unique_targets_with_metadata.tsv) carries its
own, older ecosystem_compartment column from an earlier metadata pull, and it
disagrees with the freshly-pulled omdb_all_genomes_with_locations.tsv for
several categories ("Marine sponge tissue" vs "Marine Porifera tissue",
"Hydrothermal vent plume" vs "Hydrothermal vent fluid/plume", "Hydrozoa
tissue" vs "Marine Hydrozoa tissue", "Cold seep sediment" missing entirely
from the older file). Matching on that stale label produced real-looking but
wrong 0% bars for those habitats. Using only the fresh file's labels (for
both which genomes exist AND what habitat they're in) and just checking phaC
genome-ID membership against it sidesteps the mismatch entirely.

QC: excludes phaC hits whose best_query is UNIPROT:C7BNH2, a mislabeled
reference (actually isochorismate synthase, not a PHA synthase) -- see
_phac_qc.py.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = Path('/Users/hellpark/multimodal_seusmbol/PHA_Ocean_Atlas/PHA_bioprospecting/omdb_search/results')
ACCENT = '#1B7A6E'
LOW_COLOR = '#C9622D'

# non-marine / non-habitat categories present in OMDB's ecosystem_compartment
# (freshwater, lab controls/synthetic spike-ins, terrestrial, microcosm
# experiments) -- excluded so "parts of the ocean" means actual ocean habitats
EXCLUDE = {
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}

# shorter display labels for compartments whose raw names are long/technical
DISPLAY = {
    'Seawater': 'Open water (seawater)',
    'Marine sediment': 'Seafloor sediment',
    'Marine Porifera tissue': 'Sponge tissue',
    'Brackish sea water': 'Brackish water',
    'Marine biofilm': 'Biofilm',
    'Hydrothermal vent fluid/plume': 'Hydrothermal vent (fluid/plume)',
    'Cold seep sediment': 'Cold seep sediment',
    'Hydrothermal vent sediment': 'Hydrothermal vent sediment',
    'Marine Hydrozoa tissue': 'Hydrozoa tissue',
    'Estuarine water': 'Estuarine water',
    'Coral tissue': 'Coral tissue',
    'Sea ice': 'Sea ice',
    'Marine algae thallus': 'Algae tissue',
    'Estuarine sediment': 'Estuarine sediment',
    'Cold seep water': 'Cold seep water',
    'Marine seep water': 'Seep water',
    'Animal bone biofilm': 'Whale-fall bone biofilm',
    'Subseafloor aquifer water': 'Subseafloor aquifer',
}

# ---- phaC-positive genome IDs (label-free -- just membership) ----
phac_genome_ids = set()
with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if _phac_qc.is_bad(row['best_query']):
            continue
        phac_genome_ids.add(row['genome'])

# ---- denominator AND habitat label both come from the fresh canonical
# file only -- numerator is just "is this genome ID in phac_genome_ids" ----
total_by_compartment = defaultdict(int)
phac_by_compartment = defaultdict(int)
n_genomes_total = 0
n_phac_total = 0
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        c = row['ecosystem_compartment']
        if c in EXCLUDE:
            continue
        total_by_compartment[c] += 1
        n_genomes_total += 1
        if row['genome'] in phac_genome_ids:
            phac_by_compartment[c] += 1
            n_phac_total += 1

overall_pct = 100 * n_phac_total / n_genomes_total

MIN_N = 200  # drop habitats too small to give a stable percentage
rows = []
for c, total in total_by_compartment.items():
    if total < MIN_N:
        continue
    hits = phac_by_compartment.get(c, 0)
    rows.append((c, DISPLAY.get(c, c), hits, total, 100 * hits / total))

rows.sort(key=lambda r: r[4])
print(f'{n_phac_total:,} / {n_genomes_total:,} genomes are phaC-positive overall ({overall_pct:.1f}%)')
print(f'{len(rows)} marine habitats shown (>= {MIN_N} genomes each)')
for _, label, hits, total, pct in rows[::-1]:
    print(f'  {label:32s} {pct:5.1f}%  ({hits:,}/{total:,})')

# ---- phylum-composition-expected rate, from the CMH-controlled test
# (phac_habitat_phylum_controlled_test.py) -- 'habitat' column there uses the
# same raw ecosystem_compartment strings as this script's own compartment
# key, so the join is exact, not label-matched. Run that script first if this
# file doesn't exist yet.
expected_pct_by_c = {}
test_path = ROOT / 'figures/phac_habitat_phylum_enrichment_test.tsv'
if test_path.exists():
    with open(test_path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['expected_pct_if_composition_only']:
                expected_pct_by_c[row['habitat']] = float(row['expected_pct_if_composition_only'])
missing_expected = [label for c, label, *_ in rows if c not in expected_pct_by_c]
if missing_expected:
    print(f'NOTE: no phylum-controlled expected%% for: {missing_expected} -- run '
          f'figures/scripts/phac_habitat_phylum_controlled_test.py first (or re-run it '
          f'if habitat thresholds/inputs changed since it was last run).')

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(11.5, 0.62 * len(rows) + 2), dpi=300)

labels = [r[1] for r in rows]
pcts = [r[4] for r in rows]
expected = [expected_pct_by_c.get(r[0]) for r in rows]
obs_colors = [ACCENT if p >= overall_pct else LOW_COLOR for p in pcts]
EXPECTED_COLOR = '#9AA6A0'

y = list(range(len(rows)))
h = 0.36
obs_bars = ax.barh([yi + h / 2 + 0.02 for yi in y], pcts, height=h, color=obs_colors, zorder=3, label='Observed % phaC-positive')
exp_bars = ax.barh([yi - h / 2 - 0.02 for yi in y], [e if e is not None else 0 for e in expected], height=h,
                    color=EXPECTED_COLOR, hatch='////', edgecolor='white', linewidth=0.6, zorder=3,
                    label='Expected % from phylum composition alone (CMH-controlled)')

ax.axvline(overall_pct, color='#5B6E70', linewidth=1.3, linestyle='--', zorder=2)
ax.text(overall_pct, len(rows) - 0.15, f'  dataset overall: {overall_pct:.1f}%', color='#5B6E70',
        fontsize=9.5, va='center', ha='left', style='italic')

for yi, (c, label, hits, total, pct) in zip(y, rows):
    ax.text(pct + 0.6, yi + h / 2 + 0.02, f'{pct:.1f}%', va='center', fontsize=8.7, color='#20302C', fontweight='bold')
    e = expected_pct_by_c.get(c)
    if e is not None:
        ax.text(e + 0.6, yi - h / 2 - 0.02, f'{e:.1f}%', va='center', fontsize=8.2, color='#5B6E70')

ax.set_yticks(y)
ax.set_yticklabels([f'{label}\n(n={total:,})' for _, label, _, total, _ in rows], fontsize=9.3)

ax.set_xlabel('% of genomes with phaC', fontsize=11.5)
ax.set_xlim(0, max(pcts) * 1.28)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
ax.grid(axis='x', color='#E4E8E5', linewidth=0.8, zorder=0)
ax.set_axisbelow(True)
ax.legend(loc='lower right', fontsize=9, frameon=False)

ax.set_title('phaC prevalence across ocean habitats', fontsize=15.5, fontweight='bold', loc='left', pad=14)
fig.text(0.01, 0.005,
          f'{n_phac_total:,} of {n_genomes_total:,} screened OMDB genomes ({overall_pct:.1f}%) carry a phaC hit. '
          f'Habitats with <{MIN_N} genomes and non-marine/control samples excluded. "Expected %" = phylum-composition-standardized rate '
          f'(each habitat’s own phylum mix, applied to that phylum’s phaC rate elsewhere in the ocean) from the Cochran-Mantel-Haenszel '
          f'phylum-controlled test -- see figures/PHA_CLEAN_RESULTS.md §5.3.1 and figures/phac_habitat_phylum_enrichment_test.tsv.',
          ha='left', fontsize=8.0, color='#5B6E70', wrap=True)
fig.tight_layout(rect=[0, 0.045, 1, 1])

out_path = OUT / 'phac_pct_by_ocean_habitat.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_pct_by_ocean_habitat.pdf', facecolor='white')
