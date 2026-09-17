"""Treemap of phaC-negative, PHA-gene-rich genomes (>5 other PHA genes,
no phaC -- the same qualifying criterion as plot_phac_negative_taxonomy.py,
used here only as a filter, not as the grouping variable). Two-level
hierarchy: phylum blocks (colored, sized by total qualifying genomes),
each subdivided into individual genus tiles (sized by that genus's
qualifying genome count) -- captures both the single-huge-genus case
(Luminiphilus, 1,255 genomes) and the many-smaller-genera-adding-up case
(22 Bacteroidota genera, ~800 genomes total, none individually huge) in
one figure.

Genera with fewer than MIN_GENUS_SIZE qualifying genomes are folded into
one "other genera" tile per phylum (still real genomes, just too many
distinct tiny genera to label individually); phyla with fewer than
MIN_PHYLUM_SIZE total qualifying genomes are folded into one "Other phyla"
block.

Squarified treemap layout implemented directly (Bruls, Huizing & van Wijk
2000) -- short enough that a plain implementation is more transparent
than a new dependency for one figure.
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = Path(__file__).resolve().parent.parent

MIN_GENUS_SIZE = 15
MIN_PHYLUM_SIZE = 60

PHYLUM_COLORS = {
    'Pseudomonadota': '#C9622D', 'Bacteroidota': '#1E6E7A', 'Actinomycetota': '#8B5FBF',
    'Chloroflexota': '#3E8914', 'Planctomycetota': '#C2A83E', 'Acidobacteriota': '#B33951',
    'Desulfobacterota': '#6B4226', 'Bdellovibrionota': '#6FA88A', 'SAR324': '#4A7FB5',
    'Asgardarchaeota': '#5B7C99', 'Myxococcota': '#A8763E', 'Latescibacterota': '#9E3B3B',
    'Bacillota': '#556B4F', 'Gemmatimonadota': '#D98E04', 'Verrucomicrobiota': '#3E9E8C',
}
OTHER_PHYLUM_COLOR = '#7A7A7A'
OTHER_GENUS_ALPHA = 0.35


# ---- squarified treemap layout (Bruls, Huizing, van Wijk 2000) ----
def _layoutrow(sizes, x, y, dx, dy):
    width = sum(sizes) / dy
    rects, cy = [], y
    for s in sizes:
        h = s / width
        rects.append({'x': x, 'y': cy, 'dx': width, 'dy': h})
        cy += h
    return rects


def _layoutcol(sizes, x, y, dx, dy):
    height = sum(sizes) / dx
    rects, cx = [], x
    for s in sizes:
        w = s / height
        rects.append({'x': cx, 'y': y, 'dx': w, 'dy': height})
        cx += w
    return rects


def _worst_ratio(sizes, x, y, dx, dy):
    rects = _layoutrow(sizes, x, y, dx, dy) if dx >= dy else _layoutcol(sizes, x, y, dx, dy)
    return max(max(r['dx'] / r['dy'], r['dy'] / r['dx']) for r in rects)


def squarify(sizes, x, y, dx, dy):
    sizes = list(sizes)
    if not sizes:
        return []
    if len(sizes) == 1:
        return _layoutrow(sizes, x, y, dx, dy) if dx >= dy else _layoutcol(sizes, x, y, dx, dy)
    i = 1
    while i < len(sizes) and _worst_ratio(sizes[:i], x, y, dx, dy) >= _worst_ratio(sizes[:i + 1], x, y, dx, dy):
        i += 1
    current, remaining = sizes[:i], sizes[i:]
    if dx >= dy:
        row = _layoutrow(current, x, y, dx, dy)
        width = sum(current) / dy
        leftover = (x + width, y, dx - width, dy)
    else:
        row = _layoutcol(current, x, y, dx, dy)
        height = sum(current) / dx
        leftover = (x, y + height, dx, dy - height)
    return row + squarify(remaining, *leftover)


def normalize(sizes, area):
    total = sum(sizes)
    return [s * area / total for s in sizes]


# ---- data ----
rows = list(csv.DictReader(open(FA / 'genome_family_matrix.tsv', newline=''), delimiter='\t'))
qualifying = [r for r in rows if int(r['n_phaC']) == 0 and int(r['n_families_present']) > 5]

phylum_total = Counter(r['gtdb_phylum'] for r in qualifying if r['gtdb_phylum'])
genus_total = defaultdict(Counter)
for r in qualifying:
    p, g = r['gtdb_phylum'], r['gtdb_genus']
    if p and g:
        genus_total[p][g] += 1

CANVAS_W, CANVAS_H = 100.0, 62.0

phyla_shown = [(p, n) for p, n in phylum_total.most_common() if n >= MIN_PHYLUM_SIZE]
other_phyla_total = sum(n for p, n in phylum_total.items() if n < MIN_PHYLUM_SIZE)
phylum_entries = [(p, n) for p, n in phyla_shown]
if other_phyla_total > 0:
    phylum_entries.append(('Other phyla', other_phyla_total))
phylum_entries.sort(key=lambda x: -x[1])

phylum_sizes = normalize([n for _, n in phylum_entries], CANVAS_W * CANVAS_H)
phylum_rects = squarify(phylum_sizes, 0, 0, CANVAS_W, CANVAS_H)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(17, 10.2), dpi=300)

all_tiles = []  # for genus-level tile info, used for reporting/labels

for (phylum, total_n), prect in zip(phylum_entries, phylum_rects):
    color = OTHER_PHYLUM_COLOR if phylum == 'Other phyla' else PHYLUM_COLORS.get(phylum, OTHER_PHYLUM_COLOR)

    if phylum == 'Other phyla':
        ax.add_patch(mpatches.Rectangle((prect['x'], prect['y']), prect['dx'], prect['dy'],
                                          facecolor=color, alpha=0.5, edgecolor='white', linewidth=1.6, zorder=1))
        if prect['dx'] > 4 and prect['dy'] > 3:
            ax.text(prect['x'] + prect['dx'] / 2, prect['y'] + prect['dy'] / 2, f'Other phyla\n(n={total_n})',
                     ha='center', va='center', fontsize=8, color='white', fontweight='bold', zorder=2)
        continue

    # Reserve a header band at the top of this phylum's block for its name --
    # otherwise the header text is drawn directly on top of whichever genus
    # tile happens to land in that corner, colliding with that tile's own label.
    header_h = max(1.4, min(2.8, prect['dy'] * 0.13))
    header_y = prect['y'] + prect['dy'] - header_h
    ax.add_patch(mpatches.Rectangle((prect['x'], header_y), prect['dx'], header_h,
                                      facecolor=color, alpha=0.95, edgecolor='white', linewidth=1.1, zorder=1))
    header_fontsize = 11.5 if prect['dx'] > 14 else (9 if prect['dx'] > 7 else 7)
    ax.text(prect['x'] + 0.6, header_y + header_h / 2, phylum, ha='left', va='center',
             fontsize=header_fontsize, color='white', fontweight='bold', zorder=3)

    genus_counts = genus_total[phylum]
    big = [(g, n) for g, n in genus_counts.items() if n >= MIN_GENUS_SIZE]
    big.sort(key=lambda x: -x[1])
    other_n = sum(n for g, n in genus_counts.items() if n < MIN_GENUS_SIZE)
    entries = list(big)
    if other_n > 0:
        entries.append((f'{len(genus_counts) - len(big)} smaller genera', other_n))

    genus_area_dy = prect['dy'] - header_h
    sizes = normalize([n for _, n in entries], prect['dx'] * genus_area_dy)
    rects = squarify(sizes, prect['x'], prect['y'], prect['dx'], genus_area_dy)

    for (label, n), r in zip(entries, rects):
        is_other = label.endswith('smaller genera')
        ax.add_patch(mpatches.Rectangle((r['x'], r['y']), r['dx'], r['dy'],
                                          facecolor=color, alpha=(OTHER_GENUS_ALPHA if is_other else 0.85),
                                          edgecolor='white', linewidth=1.1, zorder=1))
        area = r['dx'] * r['dy']
        if area > 2.2 and r['dx'] > 2.5 and r['dy'] > 1.8:
            fontsize = 6.3 if area < 6 else (7.6 if area < 20 else 9.5)
            txt = label if is_other else f'{label}\n{n}'
            ax.text(r['x'] + r['dx'] / 2, r['y'] + r['dy'] / 2, txt, ha='center', va='center',
                     fontsize=fontsize, color='white', fontweight='bold', zorder=2, linespacing=1.3)
        all_tiles.append((phylum, label, n))

ax.set_xlim(0, CANVAS_W)
ax.set_ylim(0, CANVAS_H)
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

fig.suptitle('Taxonomy of phaC-Negative, PHA-Gene-Rich Genomes', fontsize=18, fontweight='bold', y=0.975)
fig.text(0.5, 0.925,
          f'{len(qualifying):,} genomes with >5 other PHA genes but no phaC hit, grouped by phylum then genus — tile area = genome count',
          ha='center', fontsize=11, color='#5B6E70')

footnote = (
    f"Genera with <{MIN_GENUS_SIZE} qualifying genomes are folded into one faint \"N smaller genera\" tile per phylum (still real\n"
    f"genomes, just too many distinct tiny genera to label individually); phyla with <{MIN_PHYLUM_SIZE} total are folded into \"Other phyla.\""
)
fig.text(0.5, 0.02, footnote, ha='center', va='bottom', fontsize=8.6, color='#5B6E70')

fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.07)

out_path = OUT / 'phac_negative_treemap.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_negative_treemap.pdf', facecolor='white')
print('saved pdf too')

print()
print(f'{len(phylum_entries)} phylum blocks, {len(all_tiles)} genus tiles drawn')
print('largest genus tiles:')
for phylum, label, n in sorted(all_tiles, key=lambda t: -t[2])[:15]:
    print(f'  {phylum:<18} {label:<28} n={n}')
