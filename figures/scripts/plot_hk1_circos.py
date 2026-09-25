"""Chord diagram for "Unknown HK1" (section 9.5's runner-up, unnamed-
species candidate with 167 phaC-positive genomes and 2 median copies per
genome) -- the second half of piece 5 of the multi-part request that
CARD22-1's diagram (plot_card22_desulfoluna_circos.py) already covers.

HK1 is far too large to show all 167 genomes legibly in one chart (unlike
CARD22-1's 11), so this is scoped to the genomes actually worth looking
at: those with >=3 structurally-confirmed (qtmscore>=0.5, missing data
kept by default -- see _circos_chart.py's struct_qtm handling) phaC
copies, after dropping this project's own QC exclusions. That is 21 of
167 phaC-positive HK1 genomes (155 have >=1 confirmed copy, 87 have >=2).

Those 21 split cleanly by ecosystem, not by study: 14 are sponge/coral
symbionts (9 KELL22-1 + 3 PANK22-1 + 2 LUOR22-1, Porifera tissue; 1 TPAC,
coral tissue), matching section 9.5's note that this lineage's sponge-
symbiont prevalence is what makes it worth resolving structurally; the
other 7 are free-living sediment MAGs (2 DONG22-1 + 1 SOGI22-1 + 1
CHAS20-1, marine sediment; 2 SILV23-1, freshwater lake sediment) -- so
sector bands are colored by that host/free-living split (with sub-shades
per specific habitat) rather than by one-color-per-study as in the
CARD22-1 diagram, since 8 studies of 1-9 genomes each would be a much
less informative legend than this 2-way ecological split.

Usage:
    python figures/scripts/plot_hk1_circos.py

Outputs:
    figures/hk1_circos.png / .pdf
    figures/hk1_circos_nodes.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc
from _circos_chart import build_circos_chart, FA

ROOT = Path(__file__).resolve().parent.parent.parent

bad_targets = _phac_qc.load_bad_targets()

struct_qtm = {}
with open(ROOT / 'structural_evidence_best_hit.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        struct_qtm[row['candidate_id']] = float(row['qtmscore'])

genomes_meta = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['gtdb_species'] == 'Unknown HK1' and int(row['n_phaC']) >= 1:
            genomes_meta[row['genome']] = row

genome_of_target = {}
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        genome_of_target[row['target_id']] = row['genome']

confirmed = defaultdict(set)
for t, g in genome_of_target.items():
    if g not in genomes_meta or t in bad_targets:
        continue
    if struct_qtm.get(t, 0.5) >= 0.5:
        confirmed[g].add(t)

GENOME_ORDER = sorted((g for g, s in confirmed.items() if len(s) >= 3),
                       key=lambda g: (genomes_meta[g]['sample_source'], g))

HABITAT_COLOR = {
    'Marine sponge tissue': '#C77AB0', 'Geodia barretti tissue': '#9A5AA6', 'Agelas tissue': '#7B68C4',
    'Porites panamensis tissue': '#E2954F',
    'Marine sediment': '#4FA8A0', 'Freshwater lake sediment': '#5D8FD1',
}
BAND_COLOR = {g: HABITAT_COLOR[genomes_meta[g]['sample_source']] for g in GENOME_ORDER}
BAND_LEGEND = [
    ('#C77AB0', 'Marine sponge tissue (KELL22-1)'), ('#9A5AA6', 'Geodia barretti tissue (LUOR22-1)'),
    ('#7B68C4', 'Agelas tissue (PANK22-1)'), ('#E2954F', 'Coral tissue, Porites panamensis (TPAC)'),
    ('#4FA8A0', 'Marine sediment (DONG22-1, SOGI22-1, CHAS20-1)'),
    ('#5D8FD1', 'Freshwater lake sediment (SILV23-1)'),
]


def short_id(g):
    parts = g.split('_')
    return f'{parts[1][-4:]}-{parts[-1][-3:]}'


GENOME_LABEL = {g: f'{genomes_meta[g]["study_id"]} {short_id(g)}' for g in GENOME_ORDER}


def cluster_label(c):
    return f'paralog {c[-6:]}'


build_circos_chart(
    genome_order=GENOME_ORDER, genome_label=GENOME_LABEL, genome_band_color=BAND_COLOR, band_legend=BAND_LEGEND,
    cluster_label_fn=cluster_label, out_stem='hk1_circos', bad_targets=bad_targets,
    struct_qtm=struct_qtm, qtm_min=0.5, figsize=15, require_triad_complete=True,
    title='phaC paralogs across "Unknown HK1": the 21 genomes with >=3 structurally-confirmed copies',
    subtitle='Each sector = one genome (of 167 phaC-positive HK1 genomes, scoped to those with >=3 confirmed copies); each dot = one\n'
             'triad-complete phaC copy (qtmscore>=0.5). Chords connect copies sharing the same 70%-identity paralog cluster;\n'
             'opacity/width = exact pairwise %identity. Sector color = host tissue vs. free-living sediment habitat.',
)
