"""Recreates the phaC pairwise-%identity heatmap (same method/style as
figures/scripts/plot_hk1_phac_deep_dive.py, section 9.5.1) for the two
original single-genome deep dives -- CARD22-1_SAMN24292811_MAG_00000010
(section 9.4, Desulfoluna) and the Modicisalibacter zincidurans
representative genome ZHEN20-1_SAMN07748058_MAG_00000118 (section 9.5) --
now filtered to triad-complete copies only, matching the circos rebuild
(section 9.10 addendum) rather than the full raw copy count those
sections originally reported.

Per-genome target lists come straight from phaC_all_genomes_from_nr100_clusters.tsv
(current QC applied), then narrowed with the same hybrid triad-complete
check circos now uses (_triad_filter.py): CARD22-1's 7 raw copies drop to
6 (...041187467 fails structurally); Modicisalibacter's 5 stay at 5 (all
pass) -- not a uniform "fewer everywhere" story, just whatever the
structural check actually finds for each genome.

Usage:
    python figures/scripts/plot_deep_dive_identity_heatmaps.py

Outputs:
    figures/card22_desulfoluna_identity_heatmap.png / .pdf
    figures/modicisalibacter_identity_heatmap.png / .pdf
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc
import _triad_filter

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

GENOMES = [
    ('CARD22-1_SAMN24292811_MAG_00000010', 'card22_desulfoluna_identity_heatmap',
     'CARD22-1_SAMN24292811_MAG_00000010\n(Desulfoluna, section 9.4)'),
    ('ZHEN20-1_SAMN07748058_MAG_00000118', 'modicisalibacter_identity_heatmap',
     'ZHEN20-1_SAMN07748058_MAG_00000118\n(Modicisalibacter zincidurans, section 9.5)'),
]

bad_targets = _phac_qc.load_bad_targets()

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.mode = 'global'


def pident(a, b):
    aln = aligner.align(a, b)[0]
    s1, s2 = str(aln[0]), str(aln[1])
    m = sum(1 for x, y in zip(s1, s2) if x == y and x != '-')
    al = sum(1 for x, y in zip(s1, s2) if x != '-' and y != '-')
    return 100 * m / al if al else 0.0


plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5})

for genome, out_stem, title in GENOMES:
    raw_targets = []
    with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['genome'] == genome and row['target_id'] not in bad_targets:
                raw_targets.append(row['target_id'])

    triad_complete_ids = _triad_filter.load_triad_complete_set(raw_targets)
    targets = sorted(t for t in raw_targets if t in triad_complete_ids)
    print(f'{genome}: {len(raw_targets)} QC-passing targets -> {len(targets)} triad-complete')

    seqs = {}
    wanted = set(targets)
    cur_id, cur_seq = None, []
    with open(FA / 'phaC_cluster_sequences.faa') as f:
        for line in f:
            if line.startswith('>'):
                if cur_id in wanted:
                    seqs[cur_id] = ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id in wanted:
            seqs[cur_id] = ''.join(cur_seq)

    n = len(targets)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            mat[i, j] = 100.0 if i == j else pident(seqs[targets[i]], seqs[targets[j]])
    print(f'  {n*(n-1)//2} pairwise identities computed')

    short_labels = [f'{t[-6:]}\n({len(seqs[t])}aa)' for t in targets]

    fig, ax = plt.subplots(figsize=(8.5, 7.5), dpi=300)
    im = ax.imshow(mat, cmap='YlGnBu', vmin=0, vmax=100)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(short_labels, fontsize=8.5, rotation=45, ha='right')
    ax.set_yticklabels(short_labels, fontsize=8.5)
    for i in range(n):
        for j in range(n):
            val = mat[i, j]
            color = 'white' if val > 60 else '#20302C'
            ax.text(j, i, f'{val:.0f}', ha='center', va='center', fontsize=8.5, color=color)
    ax.set_title(f'phaC pairwise %identity: {title}\n({n} triad-complete copies)', fontsize=11.5, fontweight='bold')
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('%identity (BLOSUM62 global alignment)')
    fig.text(0.5, 0.01, f'Triad-complete only (structural geometry, section 9.12) -- {len(raw_targets) - n} of {len(raw_targets)} raw QC-passing copies dropped.',
              ha='center', fontsize=8, color='#5B6E70')
    fig.tight_layout(rect=[0, 0.04, 1, 0.88])
    fig.savefig(OUT / f'{out_stem}.png', dpi=300, facecolor='white')
    fig.savefig(OUT / f'{out_stem}.pdf', facecolor='white')
    print('saved', OUT / f'{out_stem}.png')
    plt.close(fig)
