"""Compare the substrate-binding pocket of each phaC copy in the Desulfoluna
multi-copy MAG (CARD22-1_SAMN24292811_MAG_00000010).

PHA synthase class is only a crude proxy for monomer preference -- classes I,
III and IV all take short-chain-length (C3-C5) substrates and only class II
takes medium-chain-length. What actually sets the preference is the hydrophobic
channel that the acyl chain of (R)-3-hydroxyacyl-CoA occupies, next to the
catalytic cysteine. This measures that channel directly in the ESMFold models,
so the six copies can be compared with each other before any docking.

Method is LIGSITE-style, implemented here rather than taken from a pocket
server so the geometry is reproducible from the repo:

  1. grid the space around the catalytic Cys sulfur;
  2. mark grid points that clash with a protein heavy atom (van der Waals
     radius plus a 1.4 A water probe);
  3. for each free point, cast rays along 3 axes and 4 cube diagonals and count
     how many of the 7 find protein on BOTH sides -- a protein-solvent-protein
     event. Points with at least 5 of 7 are enclosed rather than merely
     outside the protein, which is what separates a pocket from bulk solvent;
  4. measure the enclosed free volume in distance SHELLS around the catalytic
     sulfur, rather than flood-filling one connected component.

Shells, not a connected component, because the component approach is brittle:
whether the cavity next to the sulfur is one region or several depends on a
single voxel bridging two lobes, which swung one copy's answer between 6 and
291 A^3. Shells also map onto the question directly. An extended acyl chain
reaches roughly 1.25 A per carbon from the thioester sulfur, so 0-6 A is the
space a C4 substrate (3-hydroxybutyryl-CoA) occupies, 6-11 A is what C6-C8
needs beyond it, and 11-16 A is C10-C12 territory. A copy with buried
hydrophobic volume only in the first shell is built for short-chain-length
monomers; one with volume out to the third could take medium-chain-length.

Caveat carried into the output: these are predicted apo models. Side chains in
the cavity are modelled, not observed, so absolute volumes are soft. The
comparison ACROSS the six copies is the result; the absolute numbers are not.

Usage:
    python structure_prediction/analyze_phac_substrate_pocket.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'catalytic_domain'))
from find_structural_triad import find_triad, load_sequences, parse_pdb_residues  # noqa: E402

OUT = ROOT / 'structure_prediction/phac_substrate_pocket.tsv'

SPACING = 0.7          # A; fine enough to resolve a channel, coarse enough to stay quick
PROBE = 1.4            # A, water
BOX = 30.0             # A cube around the catalytic sulfur
MIN_PSP = 5            # of 7 ray directions, to call a free point enclosed
LINING_CUTOFF = 5.0    # A from a cavity point to count a residue as lining

VDW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80, 'H': 1.20}
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
      'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
      'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

COPIES = [
    ('184878', 708, 'I', 'weak class call; phaC/phaZ ambiguous'),
    ('185226', 578, 'I', ''),
    ('185910', 577, 'III', 'pLDDT 64 -- low-confidence model'),
    ('187467', 541, 'I', 'no closed triad; not a synthase'),
    ('187468', 380, 'III', 'confident class III; no partner subunit'),
    ('187880', 372, 'III', 'phaE immediately upstream'),
    ('189118', 736, 'I', 'weak class call; phaC/phaZ ambiguous'),
]


def read_atoms(path):
    """Heavy atoms: coordinates, element, residue number, residue name."""
    xyz, elements, resnum, resname = [], [], [], []
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if element == 'H':
                continue
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            elements.append(element)
            resnum.append(int(line[22:26]))
            resname.append(line[17:20].strip())
    return np.asarray(xyz), elements, np.asarray(resnum), resname


def catalytic_sulfur(path, sequence):
    """Coordinates of the catalytic nucleophile's side-chain atom.

    Canonical cysteine first, and taken as soon as it forms a tight triad;
    Ser/Thr are only an alternative when no Cys does. This is the same
    precedence figures/scripts/_triad_filter.py uses, and it matters: ranking
    all three nucleophiles by distance instead picked a surface Ser424 over the
    real Cys175 in copy 184878, and the pocket measured around it was 6 A^3 of
    nothing 23 A from the protein centre.
    """
    residues = parse_pdb_residues(path)
    best, nucleophile = None, None
    for candidate in ('C', 'S', 'T'):
        triad = find_triad(sequence, residues, nucleophile=candidate)
        if triad is None:
            continue
        tight = (triad['nuc_his_dist_A'] <= 4.5 and triad['his_asp_dist_A'] <= 4.5)
        if best is None:
            best, nucleophile = triad, candidate
        if tight:
            best, nucleophile = triad, candidate
            break
    if best is None:
        return None, None, None
    wanted = {'C': 'SG', 'S': 'OG', 'T': 'OG1'}[nucleophile]
    with open(path) as f:
        for line in f:
            if (line.startswith('ATOM') and int(line[22:26]) == best['nuc_pos']
                    and line[12:16].strip() == wanted):
                return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), best, nucleophile
    return None, best, nucleophile


SHELLS = [(0.0, 6.0, 'C3-C5'), (6.0, 11.0, 'C6-C8'), (11.0, 16.0, 'C9-C12')]


def enclosed_grid(xyz, elements, centre):
    """Boolean grid of buried, solvent-free points around `centre`, plus origin."""
    n = int(BOX / SPACING)
    origin = centre - BOX / 2
    axis = np.arange(n) * SPACING
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing='ij'), axis=-1).reshape(-1, 3) + origin

    tree = cKDTree(grid)
    occupied = np.zeros(len(grid), dtype=bool)
    for i in cKDTree(xyz).query_ball_point(centre, BOX * 0.87 + 3.0):
        occupied[tree.query_ball_point(xyz[i], VDW.get(elements[i], 1.7) + PROBE)] = True
    occupied = occupied.reshape(n, n, n)

    directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)]
    psp = np.zeros((n, n, n), dtype=np.int8)
    for d in directions:
        before = np.zeros((n, n, n), dtype=bool)
        after = np.zeros((n, n, n), dtype=bool)
        for step in range(1, n):
            before |= shift(occupied, tuple(s * step for s in d))
            after |= shift(occupied, tuple(-s * step for s in d))
        psp += (before & after).astype(np.int8)
    return (~occupied) & (psp >= MIN_PSP), origin


def shift(array, offset):
    """Boolean array shifted by `offset`, padding with False."""
    out = np.zeros_like(array)
    src = tuple(slice(max(0, -o), array.shape[i] - max(0, o)) for i, o in enumerate(offset))
    dst = tuple(slice(max(0, o), array.shape[i] - max(0, -o)) for i, o in enumerate(offset))
    if any(s.stop <= s.start for s in src):
        return out
    out[dst] = array[src]
    return out


rows = []
for short, length, model_class, note in COPIES:
    target = f'OMDBv2.0_AA_G_NR100_000041{short}'
    path = ROOT / f'{target}.pdb'
    if not path.exists():
        print(f'{short}: {path.name} not found -- skipped')
        continue
    sequence = load_sequences([target]).get(target, '')
    xyz, elements, resnum, resname = read_atoms(path)
    sulfur, triad, nucleophile = catalytic_sulfur(path, sequence)
    if sulfur is None:
        print(f'{short}: no triad found -- skipped')
        continue
    tight = triad['nuc_his_dist_A'] <= 4.5 and triad['his_asp_dist_A'] <= 4.5

    mask, origin = enclosed_grid(xyz, elements, sulfur)
    # Keep only buried space that is CONNECTED back to the catalytic nucleophile.
    # Without this the outer shells count any void that happens to lie 11-16 A
    # away, including ones walled off from the active site, which would read as
    # room for a long acyl chain that the substrate cannot actually reach.
    labels, _ = ndimage.label(mask)
    n = mask.shape[0]
    idx = np.argwhere(mask)
    if len(idx):
        d_to_nuc = np.linalg.norm(idx * SPACING + origin - sulfur, axis=1)
        reachable = {int(labels[tuple(c)]) for c in idx[d_to_nuc <= 5.0]} - {0}
        mask = np.isin(labels, list(reachable)) if reachable else np.zeros_like(mask)
    points = np.argwhere(mask) * SPACING + origin
    distance = np.linalg.norm(points - sulfur, axis=1) if len(points) else np.zeros(0)
    atom_tree = cKDTree(xyz)

    row = dict(copy=short, length_aa=length, model_class=model_class,
               triad=f"{nucleophile}{triad['nuc_pos']}/D{triad['asp_pos']}/H{triad['his_pos']}",
               triad_closed=str(bool(tight)))
    for lo, hi, label in SHELLS:
        shell = points[(distance >= lo) & (distance < hi)]
        volume = len(shell) * SPACING ** 3
        if len(shell):
            lining = set()
            hits = cKDTree(shell).query_ball_point(xyz, LINING_CUTOFF)
            for i, h in enumerate(hits):
                if h:
                    lining.add((int(resnum[i]), resname[i]))
            letters = [THREE_TO_ONE.get(name, 'X') for _, name in sorted(lining)]
            hydropathy = float(np.mean([KD.get(c, 0.0) for c in letters])) if letters else 0.0
        else:
            letters, hydropathy = [], 0.0
        row[f'volume_{label}_A3'] = round(volume, 1)
        row[f'kd_{label}'] = round(hydropathy, 2)
        row[f'n_lining_{label}'] = len(letters)
    row['note'] = note
    rows.append(row)
    print(f"{short} ({model_class}, {length}aa)  triad {row['triad']:>16} closed={str(bool(tight)):5s}  " +
          '  '.join(f"{lab} {row[f'volume_{lab}_A3']:6.1f} A^3 KD{row[f'kd_{lab}']:+5.2f}"
                    for _, _, lab in SHELLS))

with open(OUT, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
print(f'\nwrote {OUT}')
