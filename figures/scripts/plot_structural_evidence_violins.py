"""Violin-plot view of the candidate-vs-reference structural evidence from
Foldseek (PHA_CLEAN_RESULTS.md section 6): positive-control clusters
(triad-complete/HMM-supported genomes, already trusted) vs. the uncertain-
tier cluster representatives (no direct sequence evidence), each folded
with ESMFold and searched against the folded audited reference set.

Deliberately NOT drawing a cutoff line on these plots. The whole point of
this figure is to look at where the two distributions actually sit before
picking one -- see structure_prediction/build_structural_evidence_table.py's
own docstring and PHA_CLEAN_RESULTS.md section 6 for why a cutoff isn't
chosen yet. If positive controls cluster high and uncertain candidates
split into a high-scoring sub-population plus a low-scoring one, that
split point (not some a priori threshold) is what should set the cutoff
later.

Candidates with status='no_hit' (nothing above even the loose -e 10 used
in run_foldseek_search.sh) are excluded from the violins themselves (there
is no numeric score to plot) but their count is reported in the subtitle
per group -- silently dropping them would understate how much of the
uncertain tier failed to structurally match anything at all.
"""
import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
IN_PATH = ROOT / 'structure_prediction/foldseek_out/structural_evidence_best_hit.tsv'
OUT = ROOT / 'figures'

GROUPS = ['positive_control', 'uncertain']
GROUP_LABELS = ['Positive controls\n(triad/HMM-confident)', 'Uncertain tier\n(no direct sequence evidence)']
COLORS = ['#1E6E7A', '#C2A83E']

METRICS = [
    ('qtmscore', 'Query-normalized TM-score', (0, 1)),
    ('ttmscore', 'Target-normalized TM-score', (0, 1)),
    ('alntmscore', 'Alignment-normalized TM-score', (0, 1)),
    ('structural_pident', 'Structural-alignment %identity', None),
    ('qcov_struct', 'Query structural coverage', None),
    ('tcov_struct', 'Target structural coverage', None),
]


def load():
    if not IN_PATH.exists():
        raise SystemExit(
            f'{IN_PATH} not found -- run structure_prediction/run_foldseek_search.sh then '
            f'structure_prediction/build_structural_evidence_table.py first.'
        )
    by_group = {g: defaultdict(list) for g in GROUPS}
    n_hit = {g: 0 for g in GROUPS}
    n_no_hit = {g: 0 for g in GROUPS}
    with open(IN_PATH) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            g = row['group']
            if g not in by_group:
                continue
            if row['status'] != 'hit':
                n_no_hit[g] += 1
                continue
            n_hit[g] += 1
            for metric, _, _ in METRICS:
                val = row[metric]
                if val == '':
                    continue
                by_group[g][metric].append(float(val))
    return by_group, n_hit, n_no_hit


def draw_violin(ax, data_by_group, metric, title, ylim):
    data = [data_by_group[g].get(metric, []) for g in GROUPS]
    positions = range(len(GROUPS))
    non_empty = [i for i, d in enumerate(data) if d]
    if non_empty:
        parts = ax.violinplot([data[i] for i in non_empty], positions=[positions[i] for i in non_empty],
                               showmedians=False, showextrema=False, widths=0.7)
        for idx, pc in zip(non_empty, parts['bodies']):
            pc.set_facecolor(COLORS[idx])
            pc.set_alpha(0.75)
            pc.set_edgecolor('#2A2A2A')
            pc.set_linewidth(0.6)
        bp = ax.boxplot([data[i] for i in non_empty], positions=[positions[i] for i in non_empty], widths=0.1,
                         patch_artist=True, showfliers=False,
                         medianprops=dict(color='white', linewidth=2),
                         boxprops=dict(facecolor='#2A2A2A', alpha=0.85, edgecolor='none'),
                         whiskerprops=dict(color='#2A2A2A', linewidth=1.2),
                         capprops=dict(color='#2A2A2A', linewidth=1.2))
        for i in non_empty:
            med = np.median(data[i])
            ax.text(i, med, f' {med:.2f}', va='center', ha='left', fontsize=8.5, fontweight='bold', color='#111')
    ax.set_xticks(list(positions))
    ax.set_xticklabels(GROUP_LABELS, fontsize=8.5)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=8)
    if ylim:
        ax.set_ylim(*ylim)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#E4E7E2', linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def main():
    by_group, n_hit, n_no_hit = load()

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), dpi=300)

    for ax, (metric, title, ylim) in zip(axes.flat, METRICS):
        draw_violin(ax, by_group, metric, title, ylim)

    fig.suptitle('ESMFold + Foldseek: Candidate Structures vs. Audited Reference Set', fontsize=15, fontweight='bold', y=0.99)
    subtitle = '  |  '.join(
        f'{label.split(chr(10))[0]}: n={n_hit[g]} with a best hit, {n_no_hit[g]} no hit'
        for g, label in zip(GROUPS, GROUP_LABELS)
    )
    fig.text(0.5, 0.945, subtitle, ha='center', fontsize=9.3, color='#5B6E70')

    footnote = (
        'One representative structure per phaC_cluster0.7 cluster, searched against the 1,875-structure folded reference set (foldseek easy-search, TM-align mode, -e 10).\n'
        'Best hit per candidate chosen by highest alntmscore. No cutoff has been chosen yet -- see PHA_CLEAN_RESULTS.md section 6 for how these distributions inform that choice.'
    )
    fig.text(0.5, 0.01, footnote, ha='center', va='bottom', fontsize=7.8, color='#5B6E70')

    fig.subplots_adjust(left=0.06, right=0.98, top=0.87, bottom=0.1, wspace=0.28, hspace=0.35)

    out_path = OUT / 'structural_evidence_violins.png'
    fig.savefig(out_path, dpi=300, facecolor='white')
    print('saved', out_path)
    fig.savefig(OUT / 'structural_evidence_violins.pdf', facecolor='white')


if __name__ == '__main__':
    main()
