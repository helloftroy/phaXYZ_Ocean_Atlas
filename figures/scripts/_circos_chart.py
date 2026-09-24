"""Shared chord-diagram builder, factored out of
plot_modicisalibacter_circos.py once its layout (label rotation, legend
placement, sector geometry) was validated -- see that script and
PHA_CLEAN_RESULTS.md section 9.9. One sector per genome, one node per
phaC copy, chords connect same-70%-cluster copies across different
genomes, chord opacity/width = exact pairwise %identity, node color =
paralog cluster.
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path(__file__).resolve().parent.parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
OUT = Path(__file__).resolve().parent.parent

PALETTE = ['#1E6E7A', '#C2622D', '#3A6B63', '#9E3B3B', '#7A5FA0', '#C9A227', '#4A7FB5',
           '#B33951', '#6FA88A', '#8B4513', '#2E86AB', '#D98E04', '#5C8A8A', '#A65A2A']

_aligner = PairwiseAligner()
_aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
_aligner.open_gap_score = -10
_aligner.extend_gap_score = -0.5
_aligner.mode = 'global'


def _pident(a, b):
    aln = _aligner.align(a, b)[0]
    s1, s2 = str(aln[0]), str(aln[1])
    m = sum(1 for x, y in zip(s1, s2) if x == y and x != '-')
    al = sum(1 for x, y in zip(s1, s2) if x != '-' and y != '-')
    return 100 * m / al if al else 0.0


def _polar(r, deg):
    rad = np.radians(deg)
    return r * np.cos(rad), r * np.sin(rad)


def build_circos_chart(genome_order, genome_label, genome_band_color, band_legend,
                        cluster_label_fn, out_stem, title, subtitle,
                        bad_targets, struct_qtm=None, qtm_min=0.5, figsize=13):
    """
    genome_order: list of genome ids, sector order around the circle
    genome_label: dict genome -> display label (must be pre-disambiguated,
        e.g. include a unique suffix -- duplicate labels will overlap silently)
    genome_band_color: dict genome -> hex color for its sector background band
    band_legend: list of (color, label) for the sector-color legend
    cluster_label_fn: function(cluster_id) -> display string for the paralog legend
    out_stem: output filename stem (writes figures/<out_stem>.png/.pdf/_nodes.tsv)
    bad_targets: the current _phac_qc.load_bad_targets() set
    struct_qtm: optional dict target_id -> qtmscore; if given, targets WITH
        an entry below qtm_min are dropped -- targets with NO entry (never
        folded/tested, e.g. not part of whatever ESMFold run struct_qtm
        came from) are always kept, since missing structural evidence is
        not the same as failing the structural check and should not be
        silently treated as a rejection (a genome left with zero targets
        after this is dropped from genome_order too)
    """
    genome_targets = defaultdict(set)
    with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['genome'] in genome_order and row['target_id'] not in bad_targets:
                genome_targets[row['genome']].add(row['target_id'])

    if struct_qtm is not None:
        n_dropped = 0
        n_no_data = 0
        for g in list(genome_targets):
            kept = {t for t in genome_targets[g] if struct_qtm.get(t, qtm_min) >= qtm_min}
            n_no_data += sum(1 for t in genome_targets[g] if t not in struct_qtm)
            n_dropped += len(genome_targets[g]) - len(kept)
            genome_targets[g] = kept
        print(f'structural filter (qtm>={qtm_min}): dropped {n_dropped} confirmed-below-threshold targets; '
              f'{n_no_data} targets had no structural data at all and were kept by default (not treated as a failure)')

    genome_order = [g for g in genome_order if genome_targets.get(g)]
    all_tids = set()
    for s in genome_targets.values():
        all_tids |= s
    print(f'{len(genome_order)} genomes (with >=1 kept target), '
          f'{sum(len(v) for v in genome_targets.values())} protein-genome instances, {len(all_tids)} distinct target_ids')

    member_to_cluster = {}
    with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
        for row in csv.reader(f, delimiter='\t'):
            if row[1] in all_tids:
                member_to_cluster[row[1]] = row[0]

    seqs = {}
    cur_id, cur_seq = None, []
    with open(FA / 'phaC_cluster_sequences.faa') as f:
        for line in f:
            if line.startswith('>'):
                if cur_id in all_tids:
                    seqs[cur_id] = ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id in all_tids:
            seqs[cur_id] = ''.join(cur_seq)

    nodes = []
    for g in genome_order:
        for t in sorted(genome_targets[g]):
            nodes.append({'genome': g, 'target_id': t, 'cluster': member_to_cluster.get(t, t)})

    cluster_ids = sorted(set(n['cluster'] for n in nodes), key=lambda c: -sum(1 for n in nodes if n['cluster'] == c))
    cluster_color = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(cluster_ids)}

    distinct_tids = sorted(seqs)
    pident_cache = {}
    for i, t1 in enumerate(distinct_tids):
        for t2 in distinct_tids[i + 1:]:
            pident_cache[frozenset((t1, t2))] = _pident(seqs[t1], seqs[t2])
    print(f'{len(pident_cache)} pairwise identities computed')

    nodes_out = OUT / f'{out_stem}_nodes.tsv'
    with open(nodes_out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['genome', 'genome_label', 'target_id', 'cluster', 'cluster_label', 'protein_length_aa'])
        for n in nodes:
            w.writerow([n['genome'], genome_label[n['genome']], n['target_id'], n['cluster'],
                        cluster_label_fn(n['cluster']), len(seqs.get(n['target_id'], ''))])
    print(f'wrote {nodes_out}')

    # ---- geometry ----
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5})
    fig, ax = plt.subplots(figsize=(figsize, figsize), dpi=300, subplot_kw={'aspect': 'equal'})

    n_genomes = len(genome_order)
    GAP_DEG = min(3.0, 300.0 / max(n_genomes, 1))  # shrink the gap automatically for many sectors
    R_OUTER, R_NODE, R_LABEL = 1.0, 0.93, 1.08

    total_gap = GAP_DEG * n_genomes
    counts = {g: max(len(genome_targets[g]), 1) for g in genome_order}
    total_count = sum(counts.values())
    avail_deg = 360.0 - total_gap

    sector_span, sector_start = {}, {}
    start = 90.0
    for g in genome_order:
        span = avail_deg * counts[g] / total_count
        sector_start[g] = start
        sector_span[g] = span
        start -= (span + GAP_DEG)

    node_angle = {}
    for g in genome_order:
        tids_g = sorted(genome_targets[g])
        span = sector_span[g]
        n = len(tids_g)
        for i, t in enumerate(tids_g):
            frac = (i + 0.5) / n
            node_angle[(g, t)] = sector_start[g] - frac * span

    label_fontsize = 8.3 if n_genomes <= 20 else max(4.5, 8.3 - 0.05 * (n_genomes - 20))
    node_size = 70 if n_genomes <= 20 else max(18, 70 - 2 * (n_genomes - 20))

    for g in genome_order:
        wedge = mpatches.Wedge((0, 0), R_OUTER + 0.05, sector_start[g] - sector_span[g], sector_start[g],
                                 width=0.09, facecolor=genome_band_color.get(g, '#E5E5E5'), edgecolor='#8B8F8C',
                                 linewidth=0.5, zorder=2)
        ax.add_patch(wedge)
        mid = sector_start[g] - sector_span[g] / 2
        x, y = _polar(R_LABEL, mid)
        mid_norm = (mid + 180) % 360 - 180
        rot = mid_norm if -90 < mid_norm <= 90 else mid_norm + 180
        ha = 'left' if -90 < mid_norm <= 90 else 'right'
        ax.text(x, y, genome_label[g], rotation=rot, ha=ha, va='center', fontsize=label_fontsize,
                 rotation_mode='anchor', color='#20302C')

    drawn = set()
    n_chords = 0
    for i, n1 in enumerate(nodes):
        for n2 in nodes[i + 1:]:
            if n1['genome'] == n2['genome'] or n1['cluster'] != n2['cluster']:
                continue
            key = frozenset((n1['target_id'], n2['target_id'], n1['genome'], n2['genome']))
            if key in drawn:
                continue
            drawn.add(key)
            pid = 100.0 if n1['target_id'] == n2['target_id'] else pident_cache.get(frozenset((n1['target_id'], n2['target_id'])), 0)
            x1, y1 = _polar(R_NODE, node_angle[(n1['genome'], n1['target_id'])])
            x2, y2 = _polar(R_NODE, node_angle[(n2['genome'], n2['target_id'])])
            path = MplPath([(x1, y1), (0, 0), (x2, y2)], [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
            alpha = 0.12 + 0.7 * (pid / 100)
            lw = (0.3 + 2.2 * (pid / 100)) * (1.0 if n_genomes <= 20 else 0.6)
            patch = mpatches.PathPatch(path, facecolor='none', edgecolor=cluster_color[n1['cluster']], linewidth=lw, alpha=alpha, zorder=1)
            ax.add_patch(patch)
            n_chords += 1
    print(f'{n_chords} chords drawn')

    for n in nodes:
        x, y = _polar(R_NODE, node_angle[(n['genome'], n['target_id'])])
        ax.scatter([x], [y], s=node_size, color=cluster_color[n['cluster']], edgecolor='#20302C', linewidth=0.6, zorder=3)

    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.axis('off')

    legend_handles = [mpatches.Patch(color=cluster_color[c], label=cluster_label_fn(c)) for c in cluster_ids
                       if sum(1 for n in nodes if n['cluster'] == c) > 1]
    legend_handles.append(mpatches.Patch(color='#B0B6B2', label='singleton (no shared cluster in this set)'))
    legend1 = ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.28, -0.04), fontsize=8,
                          frameon=False, ncol=1, title='Paralog group (node color)', title_fontsize=8.5)
    ax.add_artist(legend1)

    band_handles = [mpatches.Patch(color=c, label=l) for c, l in band_legend]
    legend2 = ax.legend(handles=band_handles, loc='upper center', bbox_to_anchor=(0.72, -0.04), fontsize=8, frameon=False,
               title='Genome group', title_fontsize=8.5)

    # subtitle/title y positions must leave room for however many lines the
    # subtitle has -- fixed y=0.945/0.98 (validated on the 2-line Modicisalibacter/
    # CARD22-1 subtitles) collided with the title for HK1's 3-line subtitle, so
    # both are pushed apart proportionally to the extra line count.
    n_subtitle_lines = subtitle.count('\n') + 1
    suptitle_text = fig.suptitle(title, fontsize=15, fontweight='bold', y=0.98 + 0.006 * (n_subtitle_lines - 2))
    subtitle_text = fig.text(0.5, 0.945 - 0.02 * (n_subtitle_lines - 2), subtitle, ha='center', fontsize=9, color='#5B6E70')

    # bbox_inches='tight' does not reliably auto-discover legends placed via
    # bbox_to_anchor outside the axes' own data limits (confirmed live: both
    # legends below the circle were silently cropped, worse as more paralog-
    # cluster rows were added) -- bbox_extra_artists explicitly tells savefig
    # to include them (and the subtitle text, same issue) in the tight-bbox
    # calculation, which pad_inches alone cannot fix since the crop itself
    # was wrong, not just under-padded.
    extra_artists = (legend1, legend2, subtitle_text, suptitle_text)
    out_path = OUT / f'{out_stem}.png'
    fig.savefig(out_path, dpi=300, facecolor='white', bbox_inches='tight', pad_inches=0.3, bbox_extra_artists=extra_artists)
    print('\nsaved', out_path)
    fig.savefig(OUT / f'{out_stem}.pdf', facecolor='white', bbox_inches='tight', pad_inches=0.3, bbox_extra_artists=extra_artists)
    print('saved pdf too')
    return nodes, n_chords
