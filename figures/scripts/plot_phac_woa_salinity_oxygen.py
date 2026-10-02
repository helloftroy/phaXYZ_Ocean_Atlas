"""Dissolved O2 vs. salinity for phaC-positive genomes, both from WOA23
(map_phac_genomes_to_woa_salinity_oxygen.py) -- answers a direct question
about the physiochemical envelope phaC-positive genomes are sampled from,
not from any single in-situ measurement but from the WOA23 objectively-
analyzed climatological mean at each genome's own sample lat/lon/depth.

Population: 10,734/31,443 QC-passing phaC-positive genomes (34.1%) --
WOA matching needs lat/lon/depth, and only 34.2% of the QC-passing
population has depth populated at all (see the mapping script's own
printed stats). This is not a sampling choice made here, it is a real
ceiling on what this join can ever cover -- stated directly on the
figure rather than left implicit.

Color = dominant habitat per genome (ecosystem_compartment, joined from
genome_family_matrix.tsv -- not present in the WOA output file itself),
same categorical convention/palette as this project's other recent
habitat-colored figures, so colors are consistent across the atlas.

Usage:
    python figures/scripts/plot_phac_woa_salinity_oxygen.py
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

EXCLUDE_HABITATS = {
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}
CATEGORICAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
OTHER_COLOR = '#B9B7AE'
HYPOXIA_THRESHOLD_UMOL_KG = 63   # ~2 mg/L, the standard marine-hypoxia cutoff

# ---------------------------------------------------------------------
# 1. load WOA-matched genomes + habitat
# ---------------------------------------------------------------------
genome_habitat = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        h = row['ecosystem_compartment']
        if h not in EXCLUDE_HABITATS:
            genome_habitat[row['genome']] = h

rows = []
with open(FA / 'phaC_genomes_woa23_salinity_oxygen.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['woa23_salinity_psu'] and row['woa23_oxygen_umol_kg']:
            rows.append({
                'genome': row['genome'],
                'salinity': float(row['woa23_salinity_psu']),
                'oxygen': float(row['woa23_oxygen_umol_kg']),
                'habitat': genome_habitat.get(row['genome'], 'Mixed/Unknown'),
            })
print(f'{len(rows):,} genomes with both WOA23 salinity and oxygen matched')

habitat_counts = Counter(r['habitat'] for r in rows)
top_habitats = [h for h, _ in habitat_counts.most_common(len(CATEGORICAL))]
habitat_color = {h: CATEGORICAL[i] for i, h in enumerate(top_habitats)}
print('habitat counts:', habitat_counts.most_common())

n_hypoxic = sum(1 for r in rows if r['oxygen'] < HYPOXIA_THRESHOLD_UMOL_KG)
print(f'{n_hypoxic:,}/{len(rows):,} ({100*n_hypoxic/len(rows):.1f}%) in hypoxic water (<{HYPOXIA_THRESHOLD_UMOL_KG} umol/kg O2)')

# ---------------------------------------------------------------------
# 2. plot
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11.5})
TEXT_MUTED = '#5B6E70'
EDGE_WHITE = '#FCFCFB'

fig, ax = plt.subplots(figsize=(12.5, 9), dpi=300)

ax.axhspan(0, HYPOXIA_THRESHOLD_UMOL_KG, color='#F3EFE8', zorder=0)
ax.text(40.3, HYPOXIA_THRESHOLD_UMOL_KG / 2, f'hypoxic (<{HYPOXIA_THRESHOLD_UMOL_KG} μmol/kg)', fontsize=9,
        color='#8B7E63', ha='right', va='center', style='italic')

for h in top_habitats:
    pts = [r for r in rows if r['habitat'] == h]
    ax.scatter([r['salinity'] for r in pts], [r['oxygen'] for r in pts],
               s=22, color=habitat_color[h], alpha=0.75, edgecolors=EDGE_WHITE, linewidths=0.3, zorder=3)
other_pts = [r for r in rows if r['habitat'] not in habitat_color]
if other_pts:
    ax.scatter([r['salinity'] for r in other_pts], [r['oxygen'] for r in other_pts],
               s=22, color=OTHER_COLOR, alpha=0.75, edgecolors=EDGE_WHITE, linewidths=0.3, zorder=3)

ax.set_xlabel('Salinity (PSU, WOA23 annual climatological mean)', fontsize=12.5)
ax.set_ylabel('Dissolved O₂ (μmol/kg, WOA23 annual climatological mean)', fontsize=12.5)
ax.grid(True, color='#EBEEEC', linewidth=0.7, zorder=0)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('#C3C2B7')
ax.spines['bottom'].set_color('#C3C2B7')

legend_handles = [Patch(facecolor=habitat_color[h], edgecolor='none', label=h) for h in top_habitats]
if other_pts:
    legend_handles.append(Patch(facecolor=OTHER_COLOR, edgecolor='none', label='Other'))
ax.legend(handles=legend_handles, loc='lower left', fontsize=9.3, frameon=False,
          title='Dominant habitat', title_fontsize=9.8, handlelength=1.1)


# 31,443 QC-passing phaC-positive genomes / 34.2% with lat-lon-depth populated:
# printed directly by map_phac_genomes_to_woa_salinity_oxygen.py's own run,
# not recomputed here (re-deriving the full QC-passing genome set a second
# time would just duplicate that script's work for the same answer).
N_QC_PASSING_GENOMES = 31443
PCT_WITH_DEPTH = 34.2
stats_text = (
    f'{len(rows):,} phaC-positive genomes with a WOA23 salinity + O₂ match\n'
    f'({100*len(rows)/N_QC_PASSING_GENOMES:.0f}% of the {N_QC_PASSING_GENOMES:,} QC-passing phaC-positive genomes --\n'
    f'only {PCT_WITH_DEPTH:.0f}% of those have lat/lon/depth populated at all)\n'
    f'{n_hypoxic:,} ({100*n_hypoxic/len(rows):.0f}%) fall in hypoxic water'
)
ax.text(0.015, 0.985, stats_text, transform=ax.transAxes, ha='left', va='top', fontsize=9.3,
        color='#3A4A46', linespacing=1.5,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#F7F8F6', edgecolor='#D8DCD9', linewidth=0.8))

fig.suptitle('Dissolved O₂ vs. salinity, phaC-positive genomes', fontsize=18, fontweight='bold',
              x=0.065, ha='left', y=0.985)
fig.text(0.065, 0.945,
         'WOA23 1°-grid annual climatological mean at each genome’s own sample lat/lon/depth — a climatology, not a raw in-situ reading.',
         ha='left', va='top', fontsize=10.3, color=TEXT_MUTED)

fig.tight_layout(rect=[0.01, 0.01, 0.99, 0.90])

out_path = OUT / 'phac_woa_salinity_oxygen.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_woa_salinity_oxygen.pdf', facecolor='white')
