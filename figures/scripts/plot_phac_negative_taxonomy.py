"""Do phaC-negative-but-PHA-gene-rich genomes (>5 other PHA genes, no
phaC) cluster taxonomically -- pointing at one/few divergent synthase
lineages our search missed -- or are they scattered across many
unrelated taxa (more likely independent, noisy losses)? Both cases exist
in the data; this figure is built specifically to tell them apart, one
genus per point:

  x = penetrance = % of that genus's TOTAL PHA-gene genomes that are
      phaC-negative. Low = a minority exception within an otherwise
      normal phaC+ genus (scattered/noisy). High = this looks like the
      genus's normal biology, not an exception.
  y = architecture consistency = % of that genus's qualifying genomes
      that share its single most common architecture. Low = many
      different gene combinations, no one fixed replacement pathway
      (scattered). High = one specific accessory-gene "kit" recurs
      consistently (consistent with one fixed lineage-level event, not
      independent noise).

Upper-right = strongest single-novel-lineage candidates (worth a
targeted divergent-phaC search first). Lower-left = scattered/likely
independent low-frequency occurrences.
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = Path(__file__).resolve().parent.parent

MIN_QUALIFYING = 20  # per-genus floor to keep the plot to genera with enough N to trust a %

rows = list(csv.DictReader(open(FA / 'genome_family_matrix.tsv', newline=''), delimiter='\t'))
qualifying = [r for r in rows if int(r['n_phaC']) == 0 and int(r['n_families_present']) > 5]

genus_total = Counter(r['gtdb_genus'] for r in rows if r['gtdb_genus'])
genus_qual = Counter(r['gtdb_genus'] for r in qualifying if r['gtdb_genus'])
genus_qual_arch = defaultdict(Counter)
genus_phylum = {}
for r in qualifying:
    g = r['gtdb_genus']
    if not g:
        continue
    genus_qual_arch[g][r['architecture']] += 1
    genus_phylum[g] = r['gtdb_phylum']

genera = [g for g, n in genus_qual.items() if n >= MIN_QUALIFYING]

PHYLUM_COLORS = {
    'Pseudomonadota': '#C9622D', 'Bacteroidota': '#1E6E7A', 'Actinomycetota': '#8B5FBF',
    'Chloroflexota': '#3E8914', 'Planctomycetota': '#C2A83E', 'Acidobacteriota': '#B33951',
    'SAR324': '#4A7FB5', 'Desulfobacterota': '#6B4226', 'Bdellovibrionota': '#6FA88A',
}
OTHER_COLOR = '#7A7A7A'

points = []
for g in genera:
    n_qual = genus_qual[g]
    n_total = genus_total[g]
    penetrance = 100 * n_qual / n_total
    top_arch, top_n = genus_qual_arch[g].most_common(1)[0]
    consistency = 100 * top_n / n_qual
    phylum = genus_phylum[g]
    points.append({
        'genus': g, 'n_qual': n_qual, 'n_total': n_total, 'penetrance': penetrance,
        'consistency': consistency, 'top_arch': top_arch, 'phylum': phylum,
        'n_distinct_arch': len(genus_qual_arch[g]),
    })

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(13, 9.5), dpi=300)

ax.axvspan(70, 101, ymin=0.68, ymax=1.0, color='#E8B25A', alpha=0.12, zorder=0)
ax.text(99, 97, 'strongest single-lineage\ncandidates', ha='right', va='top', fontsize=9.5,
        style='italic', color='#8A5A1E', fontweight='bold')

for p in points:
    color = PHYLUM_COLORS.get(p['phylum'], OTHER_COLOR)
    size = 30 + 5.5 * p['n_qual'] ** 0.5
    ax.scatter(p['penetrance'], p['consistency'], s=size, color=color, alpha=0.75,
               edgecolor='white', linewidth=0.7, zorder=3)

# label the most interesting extremes: top-right (clustered) and largest bubbles overall
label_candidates = sorted(points, key=lambda p: -(p['penetrance'] * p['consistency']))[:8]
label_candidates += sorted(points, key=lambda p: -p['n_qual'])[:4]
seen_labels = set()
for p in label_candidates:
    if p['genus'] in seen_labels:
        continue
    seen_labels.add(p['genus'])
    ax.annotate(p['genus'], (p['penetrance'], p['consistency']), fontsize=8.3, fontweight='bold',
                color='#20302C', xytext=(5, 4), textcoords='offset points', zorder=4)

ax.set_xlim(0, 102)
ax.set_ylim(0, 102)
ax.set_xlabel('Penetrance: % of this genus\'s PHA-gene genomes that lack phaC', fontsize=11.5)
ax.set_ylabel('Architecture consistency: % of its phaC-negative genomes\nsharing one dominant gene combination', fontsize=11.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(color='#E4E7E2', linewidth=0.6, zorder=0)
ax.set_axisbelow(True)

legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=9, label=p)
                   for p, c in PHYLUM_COLORS.items()]
legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=OTHER_COLOR, markersize=9, label='Other phylum'))
leg = ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=9.5, title='Phylum', title_fontsize=10)
leg.get_title().set_fontweight('bold')

fig.suptitle('phaC-Negative, Gene-Rich Genomes: Which Ones Point to a Single Missed Lineage?',
             fontsize=15.5, fontweight='bold', y=0.985)
fig.text(0.5, 0.925,
          f'{len(genera)} genera with ≥{MIN_QUALIFYING} phaC-negative/>5-other-gene genomes (of {len(genus_qual)} total, {len(qualifying):,} qualifying genomes) — bubble size ~ √(genome count)',
          ha='center', fontsize=10, color='#5B6E70')

footnote = (
    "Upper-right = high penetrance + one consistent replacement architecture -- the strongest candidates for a single divergent/unrecognized\n"
    "phaC lineage worth a targeted search. Lower-left = low penetrance (a minority exception in an otherwise phaC+ genus) and/or many\n"
    "different architectures (no one fixed pathway) -- more consistent with scattered, independent occurrences than one missed lineage."
)
fig.text(0.5, 0.015, footnote, ha='center', va='bottom', fontsize=8.3, color='#5B6E70')

fig.subplots_adjust(left=0.09, right=0.78, top=0.87, bottom=0.13)

out_path = OUT / 'phac_negative_taxonomy.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_negative_taxonomy.pdf', facecolor='white')
print('saved pdf too')

print()
print('=== full table, sorted by penetrance x consistency ===')
for p in sorted(points, key=lambda p: -(p['penetrance'] * p['consistency'])):
    print(f"{p['genus']:<20} phylum={p['phylum']:<16} n_qual={p['n_qual']:>5}  penetrance={p['penetrance']:>5.1f}%  "
          f"consistency={p['consistency']:>5.1f}%  n_distinct_arch={p['n_distinct_arch']:>3}  top_arch={p['top_arch']}")
