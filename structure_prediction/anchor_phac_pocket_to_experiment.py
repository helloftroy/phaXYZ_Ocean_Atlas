"""Calibrate the Desulfoluna phaC channel measurements against solved PHA
synthase structures of KNOWN substrate preference.

analyze_phac_substrate_pocket.py measures a channel but cannot say what a given
volume means. These references supply the scale, because their monomer
preference is established biochemically rather than inferred:

  5HZ2         Cupriavidus necator (Ralstonia eutropha) PhaC1, class I, wild
               type -- and UniProt P23608, the same protein whose numbering
               this project's catalytic-triad columns are defined against. The
               canonical short-chain-length synthase: 3-hydroxybutyryl-CoA and
               C4-C5, essentially no medium-chain-length activity.
  5T6O         The same enzyme as the C319A mutant, catalytic cysteine removed.
               A negative control: it should and does return no channel.
  9KNK         Aeromonas caviae PhaC, class I, full length. The textbook
               broad-specificity exception -- it natively makes
               P(3HB-co-3HHx), so it accepts a C6 monomer.
  5XAV         Chromobacterium sp. USM2 PhaC, class I. Highest resolution of
               the set (1.48 A) and a useful third point.

So 5HZ2 and 9KNK bracket the question: a Desulfoluna copy whose channel
resembles the Aeromonas one more than the Cupriavidus one is the copy worth
docking a C6-C8 substrate into.

Sequences come from the structures' own ATOM records rather than a FASTA,
since these are not in this project's databases. Only chain A is used; these
are oligomers and a second copy would double every cavity that sits on the
interface.

Caveat kept in the output: the references are crystal structures and the
Desulfoluna copies are ESMFold models. Predicted models are generally less
tightly packed, so if anything this comparison flatters the models' volumes.
Treat the reference numbers as a scale, not as a threshold.

Usage:
    python structure_prediction/anchor_phac_pocket_to_experiment.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'catalytic_domain'))
sys.path.insert(0, str(ROOT / 'structure_prediction'))
from analyze_phac_substrate_pocket import (KD, LINING_CUTOFF, SHELLS, SPACING,  # noqa: E402
                                           THREE_TO_ONE, catalytic_sulfur, enclosed_grid, read_atoms)
from find_structural_triad import load_sequences  # noqa: E402

from scipy import ndimage  # noqa: E402

REF_DIR = ROOT / 'structure_prediction/reference_pdb'
OUT = ROOT / 'structure_prediction/phac_pocket_vs_experiment.tsv'

REFERENCES = [
    ('5HZ2', 'Cupriavidus necator PhaC1', 'scl (C4-C5)', 'canonical short-chain synthase'),
    # C319A: the catalytic cysteine is engineered away, so no triad and no
    # channel is the right answer for this entry. Kept as a second negative
    # control alongside copy 187467.
    ('5T6O', 'Cupriavidus necator PhaC1 C319A', 'inactive mutant', 'catalytic Cys mutated to Ala'),
    ('5XAV', 'Chromobacterium sp. USM2 PhaC', 'scl (C4-C5)', 'highest resolution, 1.48 A'),
    ('9KNK', 'Aeromonas caviae PhaC', 'scl + C6 (3HHx)', 'broad specificity; makes P(3HB-co-3HHx)'),
]
COPIES = [('184878', 'I'), ('185226', 'I'), ('185910', 'III'),
          ('187467', 'I'), ('187468', 'III'), ('187880', 'III'), ('189118', 'I')]


def single_chain(path, chain='A'):
    """Write chain `chain` to a scratch file and return its path plus a sequence
    INDEXED BY RESIDUE NUMBER.

    find_triad assumes sequence position i+1 is PDB residue i+1, which holds for
    ESMFold models numbered 1..N but not for crystal structures, which carry
    their own numbering and have gaps at disordered loops. Returning the
    observed residues concatenated made every reference report "no triad found".
    Padding absent positions with X keeps find_triad's assumption true, and X
    matches no nucleophile, Asp or His so the padding cannot create a hit.
    """
    out = path.with_suffix(f'.chain{chain}.pdb')
    by_number = {}
    with open(path) as fin, open(out, 'w') as fout:
        for line in fin:
            if not line.startswith('ATOM') or line[21] != chain:
                continue
            if line[16] not in (' ', 'A'):          # altloc A only
                continue
            if line[26] != ' ':                     # skip insertion codes
                continue
            fout.write(line)
            by_number.setdefault(int(line[22:26]), THREE_TO_ONE.get(line[17:20].strip(), 'X'))
    if not by_number:
        return out, ''
    highest = max(by_number)
    sequence = ''.join(by_number.get(i, 'X') for i in range(1, highest + 1))
    return out, sequence


def profile(path, sequence, label):
    xyz, elements, resnum, resname = read_atoms(path)
    sulfur, triad, nucleophile = catalytic_sulfur(path, sequence)
    if sulfur is None:
        print(f'{label}: no triad found -- skipped')
        return None
    tight = triad['nuc_his_dist_A'] <= 4.5 and triad['his_asp_dist_A'] <= 4.5

    mask, origin = enclosed_grid(xyz, elements, sulfur)
    labels, _ = ndimage.label(mask)
    idx = np.argwhere(mask)
    if len(idx):
        d = np.linalg.norm(idx * SPACING + origin - sulfur, axis=1)
        reachable = {int(labels[tuple(c)]) for c in idx[d <= 5.0]} - {0}
        mask = np.isin(labels, list(reachable)) if reachable else np.zeros_like(mask)
    points = np.argwhere(mask) * SPACING + origin
    distance = np.linalg.norm(points - sulfur, axis=1) if len(points) else np.zeros(0)

    row = {'triad': f"{nucleophile}{triad['nuc_pos']}/D{triad['asp_pos']}/H{triad['his_pos']}",
           'triad_closed': str(bool(tight))}
    for lo, hi, name in SHELLS:
        shell = points[(distance >= lo) & (distance < hi)]
        lining = set()
        if len(shell):
            for i, hits in enumerate(cKDTree(shell).query_ball_point(xyz, LINING_CUTOFF)):
                if hits:
                    lining.add((int(resnum[i]), resname[i]))
        letters = [THREE_TO_ONE.get(n, 'X') for _, n in sorted(lining)]
        row[f'volume_{name}_A3'] = round(len(shell) * SPACING ** 3, 1)
        row[f'kd_{name}'] = round(float(np.mean([KD.get(c, 0.0) for c in letters])) if letters else 0.0, 2)
    return row


rows = []
print('reference structures (known substrate preference)')
for pdb_id, organism, preference, note in REFERENCES:
    path = REF_DIR / f'{pdb_id}.pdb'
    if not path.exists():
        print(f'  {pdb_id}: not downloaded -- skipped')
        continue
    chain_path, sequence = single_chain(path)
    row = profile(chain_path, sequence, pdb_id)
    if row is None:
        continue
    row.update(entry=pdb_id, kind='crystal', organism=organism, substrate=preference, note=note)
    rows.append(row)
    print(f"  {pdb_id}  {organism:34s} {preference:16s} {row['triad']:>18s}  " +
          '  '.join(f"{n} {row[f'volume_{n}_A3']:6.1f}/{row[f'kd_{n}']:+.2f}" for _, _, n in SHELLS))

print('\nDesulfoluna copies (ESMFold models)')
for short, model_class in COPIES:
    target = f'OMDBv2.0_AA_G_NR100_000041{short}'
    path = ROOT / f'{target}.pdb'
    if not path.exists():
        continue
    sequence = load_sequences([target]).get(target, '')
    row = profile(path, sequence, short)
    if row is None:
        continue
    row.update(entry=short, kind='ESMFold model', organism='Desulfoluna sp. (CARD22-1 MAG 010)',
               substrate='unknown', note=f'class {model_class} by HMM')
    rows.append(row)
    print(f"  {short}  class {model_class:3s} {' ':30s} {row['triad']:>18s}  " +
          '  '.join(f"{n} {row[f'volume_{n}_A3']:6.1f}/{row[f'kd_{n}']:+.2f}" for _, _, n in SHELLS))

fields = ['entry', 'kind', 'organism', 'substrate', 'triad', 'triad_closed']
fields += [f'{p}_{n}' for _, _, n in SHELLS for p in ('volume', 'kd')]
fields = ['entry', 'kind', 'organism', 'substrate', 'triad', 'triad_closed'] + \
         [f'volume_{n}_A3' for _, _, n in SHELLS] + [f'kd_{n}' for _, _, n in SHELLS] + ['note']
with open(OUT, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fields, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
print(f'\nwrote {OUT}')
