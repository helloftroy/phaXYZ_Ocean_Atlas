"""Superpose one novel phaC candidate on the reference it structurally matches.

The example is the clearest case in the 255 (PHA_CLEAN_RESULTS.md section 13.4):
a Haliangium (Myxococcota) protein at 26.7% sequence identity to any validated
phaC, matching a Mesobacillus persicus (Bacillota) PHA synthase at TM 0.907 over
97.9% of that reference, with only 15.6% identity inside the structural
alignment itself. Two different phyla, no usable sequence signal, same fold.

Structural alignment is computed here rather than taken from Foldseek, which
returns scores but not the residue correspondence or the transform. The method
is the standard iterative one:

  1. seed the correspondence with a BLOSUM62 global sequence alignment;
  2. superpose the aligned CA pairs by Kabsch;
  3. rebuild the correspondence by Needleman-Wunsch over the TM-score distance
     term 1 / (1 + (d_ij/d0)^2) of the superposed coordinates;
  4. repeat to convergence, keeping the best-scoring iteration.

Step 1 is only a seed. At 15.6% identity the sequence alignment is close to
meaningless on its own, which is the point of iterating -- the reported TM-score
is a check on whether this reproduces Foldseek's number from a sequence seed it
should not be able to rely on.

Usage:
    python structure_prediction/superpose_novel_phac_example.py

This script is the independent measurement, not the published figure. The
figure is ray-traced by render_novel_phac_superposition.pml and laid out by
compose_superposition_figure.py; the Ca-trace rendering here exists so the
TM-score, the residue correspondence and the RMSD come from an implementation
that does not share any code with Foldseek or with PyMOL's cealign.

Outputs:
    figures/novel_phac_superposition_trace.png / .pdf   (diagnostic)
    structure_prediction/novel_phac_superposition_alignment.tsv
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'catalytic_domain'))

CANDIDATE_ID = 'RSGB23-1_GCF-000024805-V1_GENO_10000001-scaffold_1_5425'
CANDIDATE_PDB = ROOT / f'{CANDIDATE_ID}.pdb'
REFERENCE_PDB = ROOT / 'references.pdb'
REFERENCE_ACC = 'A0A1H8CR08'

CANDIDATE_COLOR, REFERENCE_COLOR = '#1c5cab', '#C2622D'
TRIAD_COLOR = '#0D9488'
TEXT_DARK, TEXT_MUTED = '#1E2630', '#5B6670'

# Triads from catalytic_domain/find_structural_triad.py. These are PDB residue
# numbers, not 0-based array indices -- find_triad keys off parse_pdb_residues'
# resseq. Verified against the residue names: 219/403/432 are CYS/ASP/HIS in the
# candidate and 134/290/319 are CYS/ASP/HIS in the reference. An earlier version
# indexed the coordinate array with these directly and so marked the residue
# after each one.
CANDIDATE_TRIAD = {'Cys': 219 - 1, 'Asp': 403 - 1, 'His': 432 - 1}
REFERENCE_TRIAD = {'Cys': 134 - 1, 'Asp': 290 - 1, 'His': 319 - 1}

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def read_ca(path):
    """CA coordinates and one-letter sequence, in residue order."""
    coords, seq = [], []
    seen = set()
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM') or line[12:16].strip() != 'CA':
                continue
            resseq = int(line[22:26])
            if resseq in seen:
                continue
            seen.add(resseq)
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            seq.append(THREE_TO_ONE.get(line[17:20].strip(), 'X'))
    return np.asarray(coords), ''.join(seq)


def kabsch(mobile, target):
    """Rotation and translation putting `mobile` onto `target` (both N x 3)."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    covariance = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(covariance)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rotation, tc - rotation @ mc


def d0_for(length):
    return max(1.24 * (length - 15) ** (1 / 3) - 1.8, 0.5)


def tm_score(distances, length_norm):
    d0 = d0_for(length_norm)
    return float(np.sum(1.0 / (1.0 + (distances / d0) ** 2)) / length_norm)


def needleman_wunsch(score, gap=-0.6):
    """Global alignment maximising `score`, returning index pairs."""
    n, m = score.shape
    best = np.zeros((n + 1, m + 1))
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)
    best[1:, 0] = np.arange(1, n + 1) * gap
    best[0, 1:] = np.arange(1, m + 1) * gap
    trace[1:, 0], trace[0, 1:] = 1, 2
    for i in range(1, n + 1):
        diag = best[i - 1, :-1] + score[i - 1]
        for j in range(1, m + 1):
            options = (diag[j - 1], best[i - 1, j] + gap, best[i, j - 1] + gap)
            k = int(np.argmax(options))
            best[i, j] = options[k]
            trace[i, j] = k
    pairs, i, j = [], n, m
    while i > 0 or j > 0:
        k = trace[i, j]
        if k == 0:
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif k == 1:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def seed_alignment(seq_a, seq_b):
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score, aligner.extend_gap_score, aligner.mode = -11, -1, 'global'
    alignment = aligner.align(seq_a, seq_b)[0]
    pairs = []
    for (a_start, a_end), (b_start, b_end) in zip(*alignment.aligned):
        pairs.extend(zip(range(a_start, a_end), range(b_start, b_end)))
    return pairs


def structural_align(coords_a, seq_a, coords_b, seq_b, length_norm, rounds=20):
    pairs = seed_alignment(seq_a, seq_b)
    best = (-1.0, pairs, np.eye(3), np.zeros(3))
    d0 = d0_for(length_norm)
    for _ in range(rounds):
        if len(pairs) < 3:
            break
        rotation, translation = kabsch(coords_a[[i for i, _ in pairs]], coords_b[[j for _, j in pairs]])
        moved = coords_a @ rotation.T + translation
        distances = np.linalg.norm(moved[[i for i, _ in pairs]] - coords_b[[j for _, j in pairs]], axis=1)
        score = tm_score(distances, length_norm)
        if score > best[0]:
            best = (score, list(pairs), rotation, translation)
        full = np.linalg.norm(moved[:, None, :] - coords_b[None, :, :], axis=2)
        new_pairs = needleman_wunsch(1.0 / (1.0 + (full / d0) ** 2))
        if new_pairs == pairs:
            break
        pairs = new_pairs
    return best


candidate_xyz, candidate_seq = read_ca(CANDIDATE_PDB)
reference_xyz, reference_seq = read_ca(REFERENCE_PDB)
print(f'candidate {len(candidate_seq)} residues, reference {len(reference_seq)} residues')

score, pairs, rotation, translation = structural_align(
    candidate_xyz, candidate_seq, reference_xyz, reference_seq, length_norm=len(reference_seq))
moved = candidate_xyz @ rotation.T + translation
distances = np.linalg.norm(moved[[i for i, _ in pairs]] - reference_xyz[[j for _, j in pairs]], axis=1)
close = distances <= 5.0
identical = sum(1 for (i, j) in pairs if candidate_seq[i] == reference_seq[j])
rmsd = float(np.sqrt(np.mean(distances[close] ** 2))) if close.any() else float('nan')
print(f'TM-score (normalised by the {len(reference_seq)}-residue reference): {score:.3f}')
print(f'{len(pairs)} aligned pairs, {int(close.sum())} within 5 A, RMSD over those {rmsd:.2f} A')
print(f'sequence identity within the structural alignment: {100 * identical / len(pairs):.1f}%')

with open(ROOT / 'structure_prediction/novel_phac_superposition_alignment.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['candidate_index', 'candidate_residue', 'reference_index', 'reference_residue', 'distance_A'])
    for (i, j), d in zip(pairs, distances):
        w.writerow([i + 1, candidate_seq[i], j + 1, reference_seq[j], f'{d:.2f}'])

# ------------------------------------------------------------------ figure
plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig = plt.figure(figsize=(10.5, 4.5), dpi=300)

for panel, (elev, azim) in enumerate(((18, -62), (74, -62))):
    ax = fig.add_subplot(1, 2, panel + 1, projection='3d')
    ax.plot(*moved.T, color=CANDIDATE_COLOR, linewidth=1.9, alpha=0.95, zorder=3)
    ax.plot(*reference_xyz.T, color=REFERENCE_COLOR, linewidth=1.9, alpha=0.85, zorder=2)
    for name, idx in CANDIDATE_TRIAD.items():
        ax.scatter(*moved[idx], s=70, color=TRIAD_COLOR, edgecolors='white', linewidths=0.8, zorder=6)
    for name, idx in REFERENCE_TRIAD.items():
        ax.scatter(*reference_xyz[idx], s=70, color=TRIAD_COLOR, edgecolors='white',
                   linewidths=0.8, marker='D', zorder=6)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    # Equal cube around both chains so the two panels are the same scale, then
    # zoom: mplot3d frames a cube with a wide default margin, which left the
    # protein occupying about a third of each panel.
    span = np.vstack([moved, reference_xyz])
    centre, reach = span.mean(axis=0), (span.max(axis=0) - span.min(axis=0)).max() / 2
    for setter, c in ((ax.set_xlim, centre[0]), (ax.set_ylim, centre[1]), (ax.set_zlim, centre[2])):
        setter(c - reach, c + reach)
    ax.set_box_aspect((1, 1, 1), zoom=1.55)

handles = [
    plt.Line2D([0], [0], color=CANDIDATE_COLOR, linewidth=2.4,
               label=f'Haliangium candidate, {len(candidate_seq)} aa (Myxococcota)'),
    plt.Line2D([0], [0], color=REFERENCE_COLOR, linewidth=2.4,
               label=f'{REFERENCE_ACC} PHA synthase, {len(reference_seq)} aa (Bacillota)'),
    plt.Line2D([0], [0], color=TRIAD_COLOR, linestyle='', marker='o', markersize=8,
               markeredgecolor='white', label='catalytic triad (Cys, Asp, His)'),
]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=9.6,
           bbox_to_anchor=(0.5, 0.004), handlelength=2.0, columnspacing=2.4)
fig.text(0.5, 0.098,
         f'TM-score {score:.3f} over {len(pairs)} aligned residues  ·  '
         f'{100 * identical / len(pairs):.1f}% identity within the alignment  ·  '
         f'backbone RMSD {rmsd:.2f} Å over the {int(close.sum())} pairs within 5 Å',
         ha='center', fontsize=9.6, color=TEXT_MUTED)

fig.subplots_adjust(left=0.0, right=1.0, top=1.06, bottom=0.20, wspace=-0.06)
fig.savefig(ROOT / 'figures/novel_phac_superposition_trace.png', dpi=300, facecolor='white')
fig.savefig(ROOT / 'figures/novel_phac_superposition_trace.pdf', facecolor='white')
print('saved', ROOT / 'figures/novel_phac_superposition_trace.png')
