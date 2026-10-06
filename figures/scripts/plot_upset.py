"""PhaC-anchored gene-architecture UpSet plot, publication PNG/PDF.

Two column groups:
  - left group:  top 30 architectures that DO contain phaC, ranked by genome count
  - right group: top 10 architectures that are missing phaC but still carry
    more than 5 other PHA genes -- i.e. a substantial accessory-gene suite
    with no detected synthase.

Every genome now known to carry a phaC is excluded from that right group,
including the ones the capped metadata table hid and the ones the
missed-synthase search turned up (see the block below). What is left on the
right is the population the figure is actually asserting: PHA-gene-rich
genomes where a dedicated whole-proteome search for a synthase still came back
empty.

Row order and the left per-gene total bars are both scoped to the
phaC-positive population only (the primary group) so their meaning doesn't
shift depending on how many phaC-negative columns are appended on the right.

Source data: PHA_bioprospecting/omdb_search/results/architecture_summary.tsv.
Regenerate that file via
`pha-reference pathway-architecture-summary` (see phaatlas/cli.py) if the
underlying genome_family_matrix.tsv changes.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phaatlas.pipeline.pathway_architecture import short_code, default_family_order

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FA = Path('/Users/hellpark/multimodal_seusmbol/PHA_Ocean_Atlas/PHA_bioprospecting/omdb_search/results')
OUT = Path(__file__).resolve().parent.parent

# Architectures are rebuilt here from genome_family_matrix.tsv, but phaC presence is
# NOT taken from the matrix's own n_phaC. This figure's right-hand columns claim to
# show genomes that carry an accessory-gene suite and genuinely have no synthase, so
# every genome known to carry one has to be moved out of that group first. Two
# separate reasons a real phaC carrier is missing from the matrix:
#   1. the matrix comes from a metadata table that lists at most 5 genomes per
#      identical protein, hiding 7,298 carriers (PHA_CLEAN_RESULTS.md section 2);
#   2. the missed-synthase search (section 13) found a triad-complete phaC in 623
#      further genomes -- some QC false negatives, some proteins the atlas never
#      held at all.
# Both are added here. The other 14 families still come from the capped matrix.
ROOT = FA.parent.parent.parent
import _phac_qc  # noqa: E402
_bad = _phac_qc.load_bad_targets()
phac_genomes = set()
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    _r = csv.reader(f, delimiter='\t')
    next(_r)
    for _t, _g in _r:
        if _t not in _bad:
            phac_genomes.add(_g)
_n_capped = len(phac_genomes)

_recovered = set()
with open(ROOT / 'phac_recovery/novel_phac_candidate_verdicts.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        _recovered.add(row['genome'])
_n_new = len(_recovered - phac_genomes)
phac_genomes |= _recovered
print(f'{_n_capped:,} phaC-positive genomes from the full membership list; the missed-synthase search found a '
      f'triad-complete phaC in {len(_recovered):,} genomes, {_n_new:,} of them not already in that list '
      f'-> {len(phac_genomes):,} total')
_arch = Counter()
_seen = set()
_n_rescued = 0
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        fams = {k[2:] for k, v in row.items() if k.startswith('n_pha') and int(v) > 0}
        _seen.add(row['genome'])
        if row['genome'] in phac_genomes:
            _n_rescued += 'phaC' not in fams
            fams.add('phaC')
        else:
            fams.discard('phaC')
        if fams:
            _arch[frozenset(fams)] += 1
for _g in phac_genomes - _seen:   # phaC-positive genomes with no row in the matrix at all
    _arch[frozenset({'phaC'})] += 1
    _n_rescued += 1
print(f'{len(phac_genomes):,} phaC-positive genomes in total; {_n_rescued:,} of them were phaC-negative or '
      f'absent in the capped matrix this figure was previously drawn from')
rows = [{'families_present': ','.join(sorted(k)), 'n_genomes': str(v), 'architecture': '+'.join(sorted(k))}
        for k, v in _arch.items()]

phac_rows = [r for r in rows if 'phaC' in r['families_present'].split(',')]
phac_rows.sort(key=lambda r: -int(r['n_genomes']))

missing_c_rows = [r for r in rows if 'phaC' not in r['families_present'].split(',') and len(r['families_present'].split(',')) > 5]
missing_c_rows.sort(key=lambda r: -int(r['n_genomes']))

# How much of the right-hand population was actually put through the whole-proteome
# search, so the figure can say "no synthase found" without overclaiming.
_searched = set(open(ROOT / 'phac_recovery/qualifying_genomes.txt').read().split())
_right_genomes = set()
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        _fams = {k[2:] for k, v in row.items() if k.startswith('n_pha') and int(v) > 0} - {'phaC'}
        if row['genome'] not in phac_genomes and len(_fams) > 5:
            _right_genomes.add(row['genome'])
n_right = len(_right_genomes)
n_right_searched = len(_right_genomes & _searched)
print(f'{n_right:,} genomes in the phaC-negative, >5-other-genes group; {n_right_searched:,} of them were '
      f'searched for a synthase across their whole proteome')

TOP_N_PHAC = 30
TOP_N_MISSING = 10
top_phac = phac_rows[:TOP_N_PHAC]
top_missing = missing_c_rows[:TOP_N_MISSING]

all_shown = top_phac + top_missing
n_genomes_shown = [int(r['n_genomes']) for r in all_shown]
combo_families = [set(r['families_present'].split(',')) for r in all_shown]
is_missing_group = [False] * len(top_phac) + [True] * len(top_missing)

# per-family totals + row order scoped to the phaC-positive population only
family_totals = Counter()
for r in phac_rows:
    for fam in r['families_present'].split(','):
        family_totals[fam] += int(r['n_genomes'])

family_order = default_family_order()
present_families = [f for f in family_order if family_totals.get(f, 0) > 0]
present_families.sort(key=lambda f: -family_totals[f])

def label_for(fam):
    code = short_code(fam)
    if code == 'RReg':
        return 'phaR-reg'
    if code == 'RSyn':
        return 'phaR-syn'
    return fam

n_rows = len(present_families)
n_cols = len(all_shown)

# one distinct color per gene, ordered to match present_families -- used for
# the left total-bars, the row-label text, and that row's filled dots, so a
# reader can trace "which gene is this row" by color alone, not just position.
ROW_COLORS = [
    '#C9622D',  # phaC -- anchor gene, warm terracotta
    '#1E6E7A',  # phaB -- teal
    '#8B5FBF',  # phaA -- violet
    '#3E8914',  # phaZ -- moss green
    '#C2A83E',  # phaJ -- ochre
    '#4A7FB5',  # phaD -- slate blue
    '#B33951',  # phaF -- rose/berry
    '#D98E04',  # phaG -- amber
    '#3E9E8C',  # phaE -- teal-green (distinct from phaB's teal)
    '#9E3B3B',  # phaR(reg) -- brick red
    '#6B4226',  # phaP -- umber brown
    '#5B7C99',  # phaY -- steel blue
    '#A8763E',  # phaI -- tan
    '#7A7A7A',  # phaR(syn) -- neutral grey
    '#556B4F',  # phaQ -- deep pine
]
row_color = {fam: ROW_COLORS[i % len(ROW_COLORS)] for i, fam in enumerate(present_families)}

PHAC_POS_COLOR = '#C9622D'
PHAC_NEG_COLOR = '#9E3B3B'
DOT_OFF = '#1B2A2E'
DOT_OFF_ALPHA = 0.15
LINE_POS = '#D9A98C'
LINE_NEG = '#D6A3A3'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.edgecolor': '#3A4442',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
})

fig = plt.figure(figsize=(19, 9.6), dpi=300)
gs = fig.add_gridspec(
    2, 2,
    width_ratios=[1.55, n_cols * 0.30],
    height_ratios=[1.35, n_rows * 0.165],
    wspace=0.16, hspace=0.04,
    top=0.87, bottom=0.115, left=0.10, right=0.99,
)

ax_corner = fig.add_subplot(gs[0, 0]); ax_corner.axis('off')
ax_top = fig.add_subplot(gs[0, 1])
ax_left = fig.add_subplot(gs[1, 0])
ax_matrix = fig.add_subplot(gs[1, 1], sharex=ax_top)

# ---- top bar chart: genome counts per architecture combo ----
xs = range(n_cols)
bar_colors = [PHAC_NEG_COLOR if m else PHAC_POS_COLOR for m in is_missing_group]
ax_top.bar(xs, n_genomes_shown, width=0.62, color=bar_colors, zorder=3)
for i, v in enumerate(n_genomes_shown):
    ax_top.text(i, v + max(n_genomes_shown) * 0.015, f'{v:,}', ha='center', va='bottom', fontsize=7.2,
                rotation=90 if v < max(n_genomes_shown) * 0.25 else 0)
ax_top.set_ylabel('Genomes', fontsize=10.5)
for s in ('top', 'right', 'left'):
    ax_top.spines[s].set_visible(False)
ax_top.set_xticks([])
ax_top.set_ylim(0, max(n_genomes_shown) * 1.18)
ax_top.grid(axis='y', color='#E4E7E2', linewidth=0.6, zorder=0)
ax_top.set_axisbelow(True)

# divider between the two column groups + small group labels
divider_x = TOP_N_PHAC - 0.5
ax_top.axvline(divider_x, color='#3A4442', linewidth=1.1, linestyle=(0, (3, 2)), zorder=4)
ax_top.text(TOP_N_PHAC / 2 - 0.5, max(n_genomes_shown) * 1.12, 'phaC-positive architectures',
            ha='center', fontsize=10, fontweight='bold', color=PHAC_POS_COLOR)
ax_top.text(TOP_N_PHAC + TOP_N_MISSING / 2 - 0.5, max(n_genomes_shown) * 1.12, 'no synthase found, ≥6 other genes',
            ha='center', fontsize=10, fontweight='bold', color=PHAC_NEG_COLOR)

# ---- left panel: per-gene totals (phaC+ population), color-coded per gene ----
ys = range(n_rows)
totals = [family_totals[f] for f in present_families]
colors_left = [row_color[f] for f in present_families]
ax_left.barh(list(ys)[::-1], totals, height=0.60, color=colors_left, zorder=3)
ax_left.set_xscale('log')
ax_left.set_xlim(max(totals) * 2.2, min(totals) * 0.7)  # padded + inverted: longest bar stops well short of the matrix edge
ax_left.set_yticks([])
for s in ('top', 'left', 'bottom', 'right'):
    ax_left.spines[s].set_visible(False)
ax_left.set_xlabel('Genomes among phaC+\n(log scale)', fontsize=9)
ax_left.tick_params(labelsize=8)
ax_left.set_ylim(-0.6, n_rows - 0.4)

# ---- dot matrix ----
for i, fam in enumerate(present_families):
    row_y = n_rows - 1 - i
    ax_matrix.axhline(row_y, color='#EEF0EC', linewidth=8, zorder=0)

ax_matrix.axvline(divider_x, color='#3A4442', linewidth=1.1, linestyle=(0, (3, 2)), zorder=1)

for col, fams in enumerate(combo_families):
    line_color = LINE_NEG if is_missing_group[col] else LINE_POS
    present_idx = [n_rows - 1 - present_families.index(f) for f in present_families if f in fams]
    if len(present_idx) > 1:
        ax_matrix.plot([col, col], [min(present_idx), max(present_idx)], color=line_color, linewidth=2.4, zorder=2)
    for i, fam in enumerate(present_families):
        row_y = n_rows - 1 - i
        on = fam in fams
        ax_matrix.scatter(
            col, row_y,
            s=100 if on else 60,
            c=row_color[fam] if on else DOT_OFF,
            alpha=1.0 if on else DOT_OFF_ALPHA,
            zorder=3, linewidths=0,
        )

ax_matrix.set_xlim(-0.6, n_cols - 0.4)
ax_matrix.set_ylim(-0.6, n_rows - 0.4)
ax_matrix.set_xticks([])
ax_matrix.set_yticks(range(n_rows))
ax_matrix.set_yticklabels([label_for(f) for f in present_families[::-1]], style='italic', fontsize=11, fontweight='bold')
for tick_label, fam in zip(ax_matrix.get_yticklabels(), present_families[::-1]):
    tick_label.set_color(row_color[fam])
ax_matrix.tick_params(left=False, pad=8)
for spine in ax_matrix.spines.values():
    spine.set_visible(False)

fig.suptitle('PHA Pathway Gene Architectures: Presence, Co-occurrence, and PhaC-Negative Outliers',
             fontsize=16, fontweight='bold', x=0.5, y=0.975)
fig.text(0.5, 0.925,
          f'Left: top {TOP_N_PHAC} architectures among {sum(int(r["n_genomes"]) for r in phac_rows):,} phaC-positive genomes '
          f'({len(phac_rows):,} distinct architectures total)  |  '
          f'Right: top {TOP_N_MISSING} of {len(missing_c_rows):,} architectures covering {n_right:,} genomes with >5 other PHA genes and no synthase found by any method',
          ha='center', fontsize=10, color='#5B6E70')

footnote = (
    "Architecture = the set of PHA pathway genes co-occurring in one genome (pipeline/pathway_architecture.py). Row order and left-panel "
    "totals are scoped to phaC-positive genomes only. Dot color = gene identity (matches row label and left bar); hollow/faint dot = absent.\n"
    f"Right-hand columns (red) carry a substantial accessory-gene suite and no synthase that any method has found: {n_right_searched:,} of these {n_right:,} genomes had their whole proteome searched with validated synthase queries.\n"
    "phaC presence combines the full protein-to-genome membership list with the triad-complete hits from that search; the other 14 genes use a table capped at 5 genomes per identical protein, so are lower bounds."
)
fig.text(0.5, 0.058, footnote, ha='center', va='top', fontsize=8.1, color='#5B6E70')

out_path = OUT / 'phaC_upset.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phaC_upset.pdf', facecolor='white')
print('saved pdf too')

print()
print('row order:', [label_for(f) for f in present_families])
print('phaC-negative outlier top combo:', top_missing[0]['architecture'], top_missing[0]['n_genomes'], 'genomes' if top_missing else 'none')
