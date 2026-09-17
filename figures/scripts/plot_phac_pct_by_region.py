"""Global choropleth: % of genomes carrying phaC, by 10x10-degree region.

Denominator note (important, not hidden): the only genome universe
available locally is genome_family_matrix.tsv, which is genomes with
>=1 hit against ANY of the 16 searched PHA-pathway families -- not the
full raw OMDB/GOPC genome catalog (which isn't downloaded locally). So
this is "% of PHA-gene-carrying genomes that specifically carry phaC",
not "% of all ocean bacteria" -- a real, meaningful ratio, just not the
broadest possible one. 251,858 genomes total, 90,846 (36.1%) carry phaC.

Cells with <MIN_GENOMES_PER_CELL genomes are masked out (light hatched
grey) rather than colored, since a handful of genomes can swing a %
wildly (e.g. 1/1 = 100%) and would misrepresent regions as extreme.
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = Path(__file__).resolve().parent.parent

GRID_DEG = 10
MIN_GENOMES_PER_CELL = 20

lat_bins = list(range(-90, 90, GRID_DEG))
lon_bins = list(range(-180, 180, GRID_DEG))
lat_idx = {b: i for i, b in enumerate(lat_bins)}
lon_idx = {b: i for i, b in enumerate(lon_bins)}

cell_total = defaultdict(int)
cell_phac = defaultdict(int)

with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        lat_raw, lon_raw = row['latitude_degN'], row['longitude_degE']
        if not lat_raw or not lon_raw:
            continue
        lat, lon = float(lat_raw), float(lon_raw)
        cell = (math.floor(lat / GRID_DEG) * GRID_DEG, math.floor(lon / GRID_DEG) * GRID_DEG)
        cell_total[cell] += 1
        if int(row['n_phaC']) > 0:
            cell_phac[cell] += 1

n_rows, n_cols = len(lat_bins), len(lon_bins)
pct_grid = np.full((n_rows, n_cols), np.nan)
n_grid = np.zeros((n_rows, n_cols))
for (clat, clon), total in cell_total.items():
    if clat not in lat_idx or clon not in lon_idx:
        continue
    i, j = lat_idx[clat], lon_idx[clon]
    n_grid[i, j] = total
    if total >= MIN_GENOMES_PER_CELL:
        pct_grid[i, j] = 100 * cell_phac[(clat, clon)] / total

max_pct = np.nanmax(pct_grid)
min_pct = np.nanmin(pct_grid)

CMAP = mcolors.LinearSegmentedColormap.from_list('phac_pct', ['#F7F0DD', '#E8B25A', '#C9622D', '#8A2E1E'])

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig = plt.figure(figsize=(17, 9), dpi=300)
ax = fig.add_axes([0.03, 0.10, 0.94, 0.76], projection=ccrs.PlateCarree())
ax.set_global()
ax.add_feature(cfeature.LAND, facecolor='#EDE9DD', edgecolor='#B9AF98', linewidth=0.4, zorder=2)
ax.gridlines(draw_labels=False, linewidth=0.3, color='#C7D0CB', linestyle=(0, (1, 3)), zorder=3)
for spine in ax.spines.values():
    spine.set_edgecolor('#3A4442')
    spine.set_linewidth(0.9)

# grey/hatched base layer for insufficient-data cells, drawn first so the ocean feature (below) still shows through gaps
insufficient = (n_grid > 0) & np.isnan(pct_grid)
lon_edges = lon_bins + [180]
lat_edges = lat_bins + [90]
masked_pct = np.ma.masked_invalid(pct_grid)

ax.add_feature(cfeature.OCEAN, facecolor='#F2F4F1', zorder=0)
mesh = ax.pcolormesh(lon_edges, lat_edges, masked_pct, transform=ccrs.PlateCarree(),
                      cmap=CMAP, vmin=min_pct, vmax=max_pct, zorder=1, edgecolors='white', linewidth=0.3)

# hatch insufficient-data-but-nonzero cells so they're visibly "no reliable estimate", not "0%"
for i in range(n_rows):
    for j in range(n_cols):
        if insufficient[i, j]:
            ax.add_patch(plt.Rectangle((lon_bins[j], lat_bins[i]), GRID_DEG, GRID_DEG, transform=ccrs.PlateCarree(),
                                        facecolor='#D8DED7', hatch='////', edgecolor='#B9C3BB', linewidth=0.2, zorder=1))

cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', shrink=0.75, pad=0.02, aspect=22)
cbar.set_label('% of PHA-gene genomes carrying phaC', fontsize=10.5)
cbar.ax.tick_params(labelsize=9)

fig.suptitle('Does phaC Prevalence Vary by Region?', fontsize=18, fontweight='bold', y=0.975)
fig.text(0.5, 0.905,
          f'% of PHA-pathway-gene-carrying genomes ({251858:,} total) that specifically carry phaC, by {GRID_DEG}°×{GRID_DEG}° region — '
          f'colored 0 to {max_pct:.0f}% (the observed max)',
          ha='center', fontsize=10.5, color='#5B6E70')

footnote = (
    f"Denominator is genomes with >=1 hit against any of the 16 searched PHA-pathway families (not the full raw genome catalog, which\n"
    f"isn't available locally) — this is \"% of PHA-gene genomes that specifically carry phaC,\" not \"% of all ocean bacteria.\" Hatched grey\n"
    f"cells have fewer than {MIN_GENOMES_PER_CELL} genomes and are masked rather than colored, since a handful of genomes can swing a % wildly."
)
fig.text(0.5, 0.02, footnote, ha='center', va='bottom', fontsize=8.4, color='#5B6E70')

out_path = OUT / 'phaC_pct_by_region.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phaC_pct_by_region.pdf', facecolor='white')
print('saved pdf too')

print(f'n cells with >= {MIN_GENOMES_PER_CELL} genomes: {np.sum(~np.isnan(pct_grid))}')
print(f'pct range shown: {min_pct:.1f}% - {max_pct:.1f}%')
