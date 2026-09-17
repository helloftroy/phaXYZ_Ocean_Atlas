import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# FA = where HPC results get scp'd to locally (raw TSV/FASTA/tree inputs live here, not in the repo)
FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
# OUT = this repo's figures directory -- generated outputs, not FA, so `git status` in PHA_Ocean_Atlas shows them
OUT = Path(__file__).resolve().parent.parent

rows = list(csv.DictReader(open(FA / 'phaC_site_ordination.tsv', newline=''), delimiter='\t'))
mantel = {r['metric']: r for r in csv.DictReader(open(FA / 'phaC_mantel_test.tsv', newline=''), delimiter='\t')}
jr, jp = float(mantel['jaccard']['r']), float(mantel['jaccard']['p_value'])

PCTVAR_PC1, PCTVAR_PC2 = 2.6, 1.9  # from live PCoA rerun, confirmed bit-identical to saved coords

biome_counts = Counter(r['ecosystem_type'] for r in rows)
# collapse rare categories (<10 sites) into "Other" for a clean legend
RARE_CUTOFF = 10
def biome_of(r):
    b = r['ecosystem_type'] or 'Unknown'
    return b if biome_counts[b] >= RARE_CUTOFF else 'Other'

biome_order_raw = [b for b, n in biome_counts.most_common() if n >= RARE_CUTOFF]
biomes = biome_order_raw + (['Other'] if any(biome_counts[b] < RARE_CUTOFF for b in biome_counts) else [])

PALETTE = {
    'Marine': '#9AA5A0',       # dominant background class -- muted, recedes
    'Porifera': '#C9622D',
    'Cnidaria': '#8B5FBF',
    'Lacustrine': '#3E8914',
    'Algae': '#C2A83E',
    'Riverine': '#4A7FB5',
    'Subsurface aquatic': '#B33951',
    'Marine cryosphere': '#1E6E7A',
    'Other': '#B9C3BB',
}
LABELS = {
    'Marine': 'Marine (open water)',
    'Porifera': 'Sponge-associated',
    'Cnidaria': 'Coral/cnidarian-associated',
    'Lacustrine': 'Lacustrine (lake)',
    'Algae': 'Algae-associated',
    'Riverine': 'Riverine',
    'Subsurface aquatic': 'Subsurface aquatic',
    'Marine cryosphere': 'Marine cryosphere (sea ice)',
    'Other': 'Other / mixed (n<10 each)',
}

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.edgecolor': '#3A4442',
    'axes.linewidth': 0.9,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'svg.fonttype': 'none',
})

fig, ax = plt.subplots(figsize=(8.5, 7.2), dpi=300)

x = {b: [] for b in biomes}
y = {b: [] for b in biomes}
for r in rows:
    b = biome_of(r)
    x[b].append(float(r['jaccard_pcoa_x']))
    y[b].append(float(r['jaccard_pcoa_y']))

# plot dominant "Marine" background class first, small & muted; overlay rarer biomes on top, larger
draw_order = sorted(biomes, key=lambda b: -biome_counts.get(b, 0) if b != 'Other' else -1e9)
draw_order = [b for b in draw_order if b == 'Marine'] + [b for b in draw_order if b != 'Marine']

for b in draw_order:
    if not x[b]:
        continue
    is_bg = (b == 'Marine')
    ax.scatter(
        x[b], y[b],
        s=10 if is_bg else 26,
        c=PALETTE.get(b, '#999'),
        alpha=0.35 if is_bg else 0.85,
        linewidths=0.3 if not is_bg else 0,
        edgecolors='#20302C' if not is_bg else 'none',
        label=f'{LABELS.get(b,b)} (n={len(x[b])})',
        zorder=2 if is_bg else 3,
        rasterized=True,
    )

ax.axhline(0, color='#C7D0CB', linewidth=0.7, zorder=1)
ax.axvline(0, color='#C7D0CB', linewidth=0.7, zorder=1)

ax.set_xlabel(f'PCoA axis 1 ({PCTVAR_PC1:.1f}% of variance)', fontsize=12)
ax.set_ylabel(f'PCoA axis 2 ({PCTVAR_PC2:.1f}% of variance)', fontsize=12)
ax.set_title('Site ordination by PhaC gene-cluster composition\n(Jaccard dissimilarity, PCoA)', fontsize=14.5, fontweight='bold', pad=14)

handles, labels_ = ax.get_legend_handles_labels()
# keep legend order = draw_order but Marine first as already is
order_map = {f'{LABELS.get(b,b)} (n={len(x[b])})': i for i, b in enumerate(draw_order)}
pairs = sorted(zip(handles, labels_), key=lambda hl: order_map.get(hl[1], 999))
handles, labels_ = zip(*pairs)
leg = ax.legend(handles, labels_, loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=9.5, markerscale=1.3, title='Biome (ecosystem_type)', title_fontsize=10)
leg.get_title().set_fontweight('bold')

footnote = (
    f"n = {len(rows)} sites (richest 2,000 of qualifying sites, min. 3 distinct clusters/site)  |  "
    f"total variance explained (2 axes) = {PCTVAR_PC1+PCTVAR_PC2:.1f}%\n"
    f"Mantel test (Jaccard dissimilarity vs. geographic distance): r = {jr:.3f}, p = {jp:.3f} "
    f"(999 permutations) — significant but weak isolation-by-distance signal"
)
fig.text(0.5, -0.02, footnote, ha='center', va='top', fontsize=8.7, color='#5B6E70', wrap=True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=10)

fig.tight_layout(rect=[0, 0.03, 0.78, 1])
out_path = OUT / 'jaccard_pcoa_biome.png'
fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
print('saved', out_path)

# also a vector PDF for true publication use (journals generally want vector, not raster, for scatter figures)
fig.savefig(OUT / 'jaccard_pcoa_biome.pdf', bbox_inches='tight', facecolor='white')
print('saved pdf too')
