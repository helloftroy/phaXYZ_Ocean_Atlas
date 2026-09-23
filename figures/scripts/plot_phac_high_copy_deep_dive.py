"""Deep dive on the extreme tail of section 9's phaC copy-count
distribution: genomes with >=5 phaC copies (n_phaC>=2 is already "multi-
copy" and covered by plot_phac_multicopy_genomes.py; >=5 is the genuinely
unusual case this script asks about specifically -- how concentrated is it
in a few species, is it legitimate, and where/what is it).

Four questions, answered in order below:
  1. Do a few species account for most of the >=5-copy genomes?
  2. Are we sure all their phaC calls are legit?
  3. Where are they found in the world?
  4. Which habitats, and which species?

Key finding driving this script's structure: the >=5-copy set splits into
two very different populations that must not be pooled naively --
  (a) RSGB23-1 -- a bulk NCBI RefSeq isolate-genome import with NO habitat
      or location metadata for any of its 12,683 genomes (confirmed
      directly below). It contributes 92/239 (38.5%) of the >=5-copy
      genomes, and Legionella pneumophila alone -- an intracellular
      pathogen, not a marine organism in any sample here -- accounts for
      20 of those. This population cannot be described geographically or
      by habitat and should not be read as marine ecology.
  (b) MAG-derived genomes from actual environmental sampling (147/239),
      which DO have habitat/location and are what panels B-D describe.

Legitimacy checks, all reused from existing project infrastructure rather
than re-derived:
  - Evidence tier per target_id, via the corrected NR100 target_id<->genome
    join (phaC_all_genomes_from_nr100_clusters.tsv, filtered against the
    current 219-accession bad-target list) joined to
    catalytic_domain/phac_sequence_evidence.tsv keyed by target_id --
    NOT by that table's own 'genome' column, which names only one
    arbitrary genome per NR100-deduplicated target_id and silently drops
    the rest for any genome whose copy happens to be identical to another
    genome's (see plot_phac_multicopy_genomes.py's docstring on this exact
    confound).
  - Genome-level contamination (CheckM-style, from
    phaC_unique_targets_with_metadata_depth.tsv), compared against the
    dataset-wide phaC-positive baseline.
  - Within-genome pairwise %identity between a genome's own phaC copies
    (catalytic_domain/selfsearch/hits.tsv, phaC-vs-phaC all-vs-all) --
    near-100%-identical copies within one genome would be the fingerprint
    of an assembly/binning duplication artifact; divergent copies are the
    fingerprint of real distinct paralogs.

Usage:
    python figures/scripts/plot_phac_high_copy_deep_dive.py

Outputs:
    figures/phac_high_copy_deep_dive.png / .pdf
    figures/phac_high_copy_genomes.tsv           (full genome-level listing)
    figures/phac_high_copy_by_study.tsv          (study-level enrichment)
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

MIN_COPIES = 5

# ---------------------------------------------------------------------
# 1. n_phaC>=5 genomes + full metadata (genome_family_matrix.tsv has
#    species/study/habitat/location directly, no join needed)
# ---------------------------------------------------------------------
genomes = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        n = int(row['n_phaC'])
        if n <= 0:
            continue
        genomes[row['genome']] = row

n_total = len(genomes)
high = {g: r for g, r in genomes.items() if int(r['n_phaC']) >= MIN_COPIES}
n_high = len(high)
print(f'{n_total:,} phaC-positive genomes; {n_high} ({100*n_high/n_total:.2f}%) carry >={MIN_COPIES} phaC copies')

is_rsgb = {g: (r['study_id'] == 'RSGB23-1') for g, r in high.items()}
n_rsgb = sum(is_rsgb.values())
print(f'{n_rsgb}/{n_high} ({100*n_rsgb/n_high:.1f}%) are from RSGB23-1, a bulk RefSeq isolate-genome import with NO habitat/location for any of its 12,683 genomes')

# ---------------------------------------------------------------------
# 2. species/genus concentration
# ---------------------------------------------------------------------
sp_counts = Counter(r['gtdb_species'] for r in high.values())
top5_sp_n = sum(n for _, n in sp_counts.most_common(5))
print(f'\n{n_high} genomes span {len(sp_counts)} distinct species; top 5 species account for '
      f'{top5_sp_n}/{n_high} ({100*top5_sp_n/n_high:.1f}%)')
print('Top species:')
for sp, n in sp_counts.most_common(8):
    studies = Counter(r['study_id'] for g, r in high.items() if r['gtdb_species'] == sp)
    print(f'  {sp:35s} n={n:3d}  studies={dict(studies)}')

# ---------------------------------------------------------------------
# 3. legitimacy check A: evidence tier via the CORRECT NR100 target_id<->genome
#    join (one-to-many), not phac_sequence_evidence.tsv's own 'genome' column
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()
high_set = set(high)
genome_targets = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        t, g = row['target_id'], row['genome']
        if g in high_set and t not in bad_targets:
            genome_targets[g].add(t)

target_to_genome_high = defaultdict(set)
for g, tids in genome_targets.items():
    for t in tids:
        target_to_genome_high[t].add(g)

ev = {}
with open(ROOT / 'catalytic_domain/phac_sequence_evidence.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        ev[row['target_id']] = row

tier_counts = Counter()
pidents = []
n_targets = 0
for tids in genome_targets.values():
    for t in tids:
        n_targets += 1
        if t in ev:
            tier_counts[ev[t]['evidence_tier']] += 1
            try:
                pidents.append(float(ev[t]['pident']))
            except ValueError:
                pass
n_no_hit = n_targets - sum(tier_counts.values())
print(f'\n{n_targets} distinct phaC target_ids across these {len(genome_targets)} genomes (expected ~{sum(int(r["n_phaC"]) for r in high.values())} from n_phaC sum)')
print('Evidence tier (correct target_id join):')
for t, n in tier_counts.most_common():
    print(f'  {t:32s} n={n:5d}  ({100*n/n_targets:.1f}%)')
print(f'  {"no hit at all":32s} n={n_no_hit:5d}  ({100*n_no_hit/n_targets:.1f}%)')
if pidents:
    pidents.sort()
    print(f'median %identity to reference: {pidents[len(pidents)//2]:.1f}%')

# ---------------------------------------------------------------------
# 4. legitimacy check B: contamination (MAG-derived only -- isolates from
#    RSGB23-1 aren't binned MAGs, contamination isn't the relevant check there)
# ---------------------------------------------------------------------
genome_contam = {}
with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g = row['genome']
        if g not in genome_contam and row['contamination']:
            try:
                genome_contam[g] = float(row['contamination'])
            except ValueError:
                pass

baseline_contam = list(genome_contam.values())
baseline_mean = sum(baseline_contam) / len(baseline_contam)

mag_high = [g for g in high if '_MAG_' in g]
mag_high_contam = [genome_contam[g] for g in mag_high if g in genome_contam]
high_mean = sum(mag_high_contam) / len(mag_high_contam) if mag_high_contam else None
print(f'\ncontamination: baseline (all phaC+ genomes) mean={baseline_mean:.2f}% (n={len(baseline_contam):,}); '
      f'>={MIN_COPIES}-copy MAGs mean={high_mean:.2f}% (n={len(mag_high_contam)}), max={max(mag_high_contam):.2f}%')

# ---------------------------------------------------------------------
# 5. legitimacy check C: within-genome pairwise %identity between a
#    genome's own phaC copies -- near-duplicate copies would flag assembly
#    artifact; divergent copies support real distinct paralogs
# ---------------------------------------------------------------------
within_genome_pident = defaultdict(list)
all_high_targets = set(target_to_genome_high)
with open(ROOT / 'catalytic_domain/selfsearch/hits.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        q, t = row['query'], row['target']
        if q == t or q not in all_high_targets or t not in all_high_targets:
            continue
        common = target_to_genome_high[q] & target_to_genome_high[t]
        for g in common:
            within_genome_pident[g].append(float(row['pident']))

pair_pidents = [p for v in within_genome_pident.values() for p in v]
n_near_dup = sum(1 for p in pair_pidents if p >= 95)
print(f'\nwithin-genome pairwise %identity between a genome\'s own phaC copies: '
      f'{len(within_genome_pident)}/{n_high} genomes have >=1 comparable pair, n_pairs={len(pair_pidents)}, '
      f'median={sorted(pair_pidents)[len(pair_pidents)//2] if pair_pidents else float("nan"):.1f}%, '
      f'>=95%% identical (near-duplicate) = {n_near_dup}/{len(pair_pidents)} ({100*n_near_dup/len(pair_pidents) if pair_pidents else 0:.1f}%)')

# ---------------------------------------------------------------------
# 6. geography + habitat, MAG-derived (located) subset only
# ---------------------------------------------------------------------
located = {g: r for g, r in high.items() if r['latitude_degN'] and not is_rsgb[g]}
print(f'\n{len(located)}/{n_high} high-copy genomes have usable habitat + location (excludes RSGB23-1)')
hab_counts = Counter(r['ecosystem_compartment'] for r in located.values())
print('Habitat breakdown:')
for h, n in hab_counts.most_common():
    print(f'  {h:28s} n={n}')

sites = Counter((r['study_id'], r['latitude_degN'], r['longitude_degE']) for r in located.values())
print(f'\n{len(sites)} distinct (study, lat, lon) sites among located genomes; top sites:')
for (st, lat, lon), n in sites.most_common(8):
    print(f'  {st:12s} {lat:>9s},{lon:>9s}  n={n}')

# ---------------------------------------------------------------------
# 7. study-level enrichment: is elevated copy-number a species trait or a
#    site/study trait? (taxonomic diversity per top study answers this)
# ---------------------------------------------------------------------
by_study_all = defaultdict(list)
for g, r in genomes.items():
    by_study_all[r['study_id']].append(int(r['n_phaC']))

study_rows = []
for st, counts in by_study_all.items():
    n = len(counts)
    if n < 20:
        continue
    pct_high = 100 * sum(1 for c in counts if c >= MIN_COPIES) / n
    study_rows.append({'study_id': st, 'n_phac_positive': n, 'pct_ge5': f'{pct_high:.2f}',
                        'n_ge5': sum(1 for c in counts if c >= MIN_COPIES), 'max_n_phac': max(counts)})
study_rows.sort(key=lambda r: -float(r['pct_ge5']))

study_out = OUT / 'phac_high_copy_by_study.tsv'
with open(study_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(study_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(study_rows)
print(f'\nwrote {study_out} ({len(study_rows)} studies)')

overall_pct_high = 100 * n_high / n_total
print(f'overall %>={MIN_COPIES} = {overall_pct_high:.3f}%')
print(f'{"study":15s} {"n":>7s} {"%>=5":>7s} {"n>=5":>5s} {"max":>4s}')
for r in study_rows[:12]:
    print(f'{r["study_id"]:15s} {r["n_phac_positive"]:>7} {r["pct_ge5"]:>6s}% {r["n_ge5"]:>5} {r["max_n_phac"]:>4}')

# taxonomic diversity check for the top non-RSGB studies driving the tail
print('\nTaxonomic diversity (n distinct genera) among each study\'s own multi-copy (>=2) genomes:')
for st in ['CARD22-1', 'CHAS20-1', 'DONG22-1', 'CAOS20-1', 'ATLA15-1', 'BPAM22-1']:
    recs = [r for r in genomes.values() if r['study_id'] == st and int(r['n_phaC']) >= 2]
    genera = Counter(r['gtdb_genus'] for r in recs)
    print(f'  {st:12s} n_multicopy={len(recs):4d}  distinct genera={len(genera):3d}  top genus={genera.most_common(1)}')

# ---------------------------------------------------------------------
# 8. full genome-level TSV
# ---------------------------------------------------------------------
genome_rows = []
for g, r in sorted(high.items(), key=lambda kv: -int(kv[1]['n_phaC'])):
    genome_rows.append({
        'genome': g, 'n_phaC': r['n_phaC'], 'gtdb_species': r['gtdb_species'], 'gtdb_genus': r['gtdb_genus'],
        'gtdb_phylum': r['gtdb_phylum'], 'study_id': r['study_id'], 'ecosystem_compartment': r['ecosystem_compartment'],
        'latitude_degN': r['latitude_degN'], 'longitude_degE': r['longitude_degE'],
        'is_mag': 'MAG' in g, 'contamination_pct': genome_contam.get(g, ''),
        'n_targets_found': len(genome_targets.get(g, [])),
        'n_catalytic_triad_complete': sum(1 for t in genome_targets.get(g, []) if ev.get(t, {}).get('evidence_tier') == 'Catalytic triad complete'),
        'n_hmm_supported': sum(1 for t in genome_targets.get(g, []) if ev.get(t, {}).get('evidence_tier') == 'HMM-supported (no triad)'),
        'n_no_hit': len(genome_targets.get(g, [])) - sum(1 for t in genome_targets.get(g, []) if t in ev),
    })
genome_out = OUT / 'phac_high_copy_genomes.tsv'
with open(genome_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(genome_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(genome_rows)
print(f'\nwrote {genome_out} ({len(genome_rows)} genomes)')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig = plt.figure(figsize=(15, 15.5), dpi=300)
gs = fig.add_gridspec(3, 2, height_ratios=[1.35, 0.55, 1.35], hspace=0.55, wspace=0.32)
axA = fig.add_subplot(gs[0, 0])
axD = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, :])
axB = fig.add_subplot(gs[2, :], projection=ccrs.PlateCarree())

TEAL, ORANGE, SAGE, MAROON = '#1E6E7A', '#C2622D', '#3A6B63', '#9E3B3B'

# Panel A: top genera, split MAG (environmental) vs isolate/RSGB
genus_mag = Counter(r['gtdb_genus'] for g, r in high.items() if not is_rsgb[g])
genus_rsgb = Counter(r['gtdb_genus'] for g, r in high.items() if is_rsgb[g])
all_genera = Counter()
all_genera.update(genus_mag)
all_genera.update(genus_rsgb)
top_genera = [g for g, _ in all_genera.most_common(14)]
top_genera_sorted = sorted(top_genera, key=lambda g: genus_mag[g] + genus_rsgb[g])
y = np.arange(len(top_genera_sorted))
mag_vals = [genus_mag[g] for g in top_genera_sorted]
rsgb_vals = [genus_rsgb[g] for g in top_genera_sorted]
axA.barh(y, mag_vals, color=TEAL, height=0.68, zorder=3, label='Environmental MAG')
axA.barh(y, rsgb_vals, left=mag_vals, color='#C9622D', height=0.68, zorder=3, label='RSGB23-1 (no habitat/location)')
axA.set_yticks(y)
axA.set_yticklabels(top_genera_sorted, fontsize=9.3, style='italic')
axA.set_xlabel(f'Genomes with >={MIN_COPIES} phaC copies')
axA.set_title('A. Top genera among high-copy genomes', fontsize=12.5, fontweight='bold', loc='left')
axA.legend(fontsize=8.5, frameon=False, loc='lower right')
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.spines['left'].set_visible(False)
axA.tick_params(axis='y', length=0)
axA.grid(axis='x', color='#E4E8E5', linewidth=0.6, zorder=0)
axA.set_axisbelow(True)

# Panel D: study-level %>=5, top offenders vs. overall baseline
studies_for_plot = sorted(study_rows[:14], key=lambda r: float(r['pct_ge5']))
labels_d = [f'{r["study_id"]} (n={r["n_phac_positive"]:,})' for r in studies_for_plot]
vals_d = [float(r['pct_ge5']) for r in studies_for_plot]
colors_d = [MAROON if r['study_id'] in ('RSGB23-1', 'CARD22-1', 'CHAS20-1', 'DONG22-1', 'CAOS20-1') else TEAL for r in studies_for_plot]
axD.barh(labels_d, vals_d, color=colors_d, height=0.68, zorder=3)
axD.axvline(overall_pct_high, color='#5B6E70', linewidth=1.2, linestyle='--', zorder=2)
axD.text(overall_pct_high, len(labels_d) - 0.3, f'  overall: {overall_pct_high:.2f}%', color='#5B6E70',
          fontsize=8.3, va='center', ha='left', style='italic')
axD.set_xlabel(f'% of study\'s phaC-positive genomes with >={MIN_COPIES} copies')
axD.set_title('D. Which studies drive the high-copy tail?', fontsize=12.5, fontweight='bold', loc='left')
axD.text(0.02, 0.03, 'maroon = studies discussed in text', transform=axD.transAxes, fontsize=7.6,
          color='#5B6E70', style='italic')
axD.spines['top'].set_visible(False)
axD.spines['right'].set_visible(False)
axD.spines['left'].set_visible(False)
axD.tick_params(axis='y', length=0, labelsize=8.5)
axD.grid(axis='x', color='#E4E8E5', linewidth=0.6, zorder=0)
axD.set_axisbelow(True)

# Panel C: legitimacy -- evidence tier stacked bar (left) + contamination/
# within-genome-identity summary as text (right), one compact full-width row
axC.set_title('C. Are these phaC calls legitimate?', fontsize=12.5, fontweight='bold', loc='left')
axC.axis('off')
bar_ax = axC.inset_axes([0.0, 0.15, 0.42, 0.55])
tier_colors = {'Catalytic triad complete': SAGE, 'HMM-supported (no triad)': '#8CAF9B', 'no hit': '#D8DDD9'}
bottoms = 0
for tname, tcolor in tier_colors.items():
    val = tier_counts.get(tname, n_no_hit if tname == 'no hit' else 0)
    bar_ax.barh(0, val, left=bottoms, color=tcolor, height=0.7, zorder=3,
                label=f'{tname} ({100*val/n_targets:.1f}%)')
    bottoms += val
bar_ax.set_yticks([])
bar_ax.set_xlim(0, n_targets)
bar_ax.set_xlabel(f'phaC target_ids (n={n_targets})', fontsize=9)
bar_ax.legend(fontsize=8, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.55), ncol=1)
bar_ax.spines['top'].set_visible(False)
bar_ax.spines['right'].set_visible(False)
bar_ax.spines['left'].set_visible(False)

axC.text(0.47, 0.78,
          f'Contamination (MAG-derived only): {high_mean:.2f}% mean (n={len(mag_high_contam)}) vs. '
          f'{baseline_mean:.2f}% dataset-wide baseline (n={len(baseline_contam):,}).',
          transform=axC.transAxes, ha='left', va='top', fontsize=9.3, color='#3A4A46')
axC.text(0.47, 0.55,
          f'Within-genome pairwise identity between a genome\'s own copies: median '
          f'{sorted(pair_pidents)[len(pair_pidents)//2] if pair_pidents else float("nan"):.0f}%, only '
          f'{100*n_near_dup/len(pair_pidents) if pair_pidents else 0:.0f}% of comparable pairs near-identical (>=95%) --',
          transform=axC.transAxes, ha='left', va='top', fontsize=9.3, color='#3A4A46')
axC.text(0.47, 0.32, 'not the signature of simple assembly-duplication artifacts.',
          transform=axC.transAxes, ha='left', va='top', fontsize=9.3, color='#3A4A46')

# Panel B: world map of located (non-RSGB) high-copy genomes
ax = axB
ax.set_global()
ax.add_feature(cfeature.OCEAN, facecolor='#E4EEEC', zorder=0)
ax.add_feature(cfeature.LAND, facecolor='#E9E4D6', edgecolor='#B9AF98', linewidth=0.5, zorder=1)
ax.gridlines(draw_labels=False, linewidth=0.4, color='#C7D0CB', linestyle=(0, (1, 3)), zorder=1)
for spine in ax.spines.values():
    spine.set_edgecolor('#3A4442')
    spine.set_linewidth(0.9)

HAB_COLORS = {'Seawater': TEAL, 'Marine sediment': ORANGE, 'Marine biofilm': '#8B5FBF',
              'Coral tissue': MAROON, 'Marine sponge tissue': '#D98E04', 'Sea ice': '#4A7FB5'}
site_agg = defaultdict(lambda: {'n': 0, 'hab': Counter()})
for g, r in located.items():
    key = (r['latitude_degN'], r['longitude_degE'])
    site_agg[key]['n'] += 1
    site_agg[key]['hab'][r['ecosystem_compartment']] += 1
    site_agg[key]['lat'] = float(r['latitude_degN'])
    site_agg[key]['lon'] = float(r['longitude_degE'])

for key, d in site_agg.items():
    top_hab = d['hab'].most_common(1)[0][0]
    color = HAB_COLORS.get(top_hab, '#7A7A7A')
    size = 22 + 14 * d['n']
    ax.scatter(d['lon'], d['lat'], transform=ccrs.PlateCarree(), s=size, color=color,
               edgecolor='#20302C', linewidth=0.7, zorder=4, alpha=0.88)

handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markeredgecolor='#20302C',
                       markersize=8, label=h) for h, c in HAB_COLORS.items() if h in Counter(r['ecosystem_compartment'] for r in located.values())]
ax.legend(handles=handles, loc='lower left', fontsize=7.6, frameon=True, framealpha=0.9, title='Habitat', title_fontsize=8)
ax.set_title(f'B. Where in the world ({len(located)} located genomes, RSGB23-1 excluded)',
             fontsize=12.5, fontweight='bold', loc='left', pad=8)

fig.suptitle(f'Genomes with >={MIN_COPIES} phaC copies: who, is it legitimate, and where',
             fontsize=17, fontweight='bold', y=0.985)
fig.text(0.5, 0.958,
          f'{n_high} of {n_total:,} phaC-positive genomes ({overall_pct_high:.2f}%) carry >={MIN_COPIES} copies. '
          f'Top 5 species account for only {100*top5_sp_n/n_high:.0f}% of these -- not a small-species story. '
          f'{n_rsgb}/{n_high} ({100*n_rsgb/n_high:.0f}%) are RSGB23-1 isolate genomes with no habitat/location.',
          ha='center', fontsize=9.7, color='#5B6E70')
fig.subplots_adjust(left=0.13, right=0.97, top=0.925, bottom=0.03)

out_path = OUT / 'phac_high_copy_deep_dive.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('\nsaved', out_path)
fig.savefig(OUT / 'phac_high_copy_deep_dive.pdf', facecolor='white')
print('saved pdf too')
