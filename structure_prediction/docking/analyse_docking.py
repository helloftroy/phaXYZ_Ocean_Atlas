"""Turn the raw Vina affinities into a chain-length preference per protein.

A raw docking score cannot be compared across chain lengths. Vina's affinity
grows more negative with ligand size almost regardless of fit -- the C12
ligand has 21 heavy atoms against the C4 ligand's 13, so it makes more contacts
wherever it lands. Reading a per-protein series of raw scores as a preference
curve would mostly be reading ligand size.

Two corrections, reported side by side:

  ligand efficiency   affinity divided by heavy atom count. The usual first
                      move, and still imperfect -- efficiency tends to fall
                      with size for reasons unrelated to the pocket.
  reference-relative  for EACH ligand separately, each receptor's affinity
                      minus the mean affinity the two strictly short-chain
                      reference enzymes get for that same ligand. Because the
                      comparison is within one ligand, size cancels exactly.
                      This is the number to read.

A negative reference-relative value means the receptor binds that chain length
better than Cupriavidus necator and Chromobacterium USM2 PhaC do.

The panel validates itself. Aeromonas caviae PhaC accepts a C6 monomer and the
two references do not, so its reference-relative curve should improve from C4
towards C6-C8. If it does not, the readout is not measuring what it claims and
the Desulfoluna numbers should not be interpreted -- the script says so rather
than leaving it to the reader.

Usage:
    python structure_prediction/docking/analyse_docking.py
"""
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'docking_scores.tsv'
OUT = HERE / 'docking_preference.tsv'

SCL_REFERENCES = ['5HZ2', '5XAV']
MCL_REFERENCE = '9KNK'

rows = [r for r in csv.DictReader(open(SOURCE, newline=''), delimiter='\t') if r['affinity_kcal_mol']]
for r in rows:
    r['affinity'] = float(r['affinity_kcal_mol'])
    r['carbons'] = int(r['carbons'])

ligands = sorted({r['ligand'] for r in rows}, key=lambda n: int(n.split('_')[0][1:]))
receptors = {}
for r in rows:
    receptors.setdefault(r['receptor'], r)

baseline = {}
for ligand in ligands:
    vals = [r['affinity'] for r in rows if r['ligand'] == ligand and r['receptor'] in SCL_REFERENCES]
    baseline[ligand] = st.mean(vals) if vals else None

by_receptor = defaultdict(dict)
for r in rows:
    by_receptor[r['receptor']][r['ligand']] = r

out_rows = []
for name, info in receptors.items():
    for ligand in ligands:
        r = by_receptor[name].get(ligand)
        if r is None:
            continue
        heavy = int(r['heavy_atoms']) if r['heavy_atoms'] else None
        relative = r['affinity'] - baseline[ligand] if baseline[ligand] is not None else None
        out_rows.append(dict(receptor=name, kind=r['kind'], label=r['label'], substrate=r['substrate'],
                             ligand=ligand, carbons=r['carbons'], heavy_atoms=heavy or '',
                             affinity_kcal_mol=round(r['affinity'], 2),
                             ligand_efficiency=round(r['affinity'] / heavy, 3) if heavy else '',
                             vs_scl_reference=round(relative, 2) if relative is not None else ''))

with open(OUT, 'w', newline='') as f:
    w = csv.DictWriter(f, delimiter='\t', fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)

header = f"{'receptor':10s} {'substrate':12s} " + ' '.join(f'{l.split("_")[0]:>7s}' for l in ligands)
print('Reference-relative affinity (kcal/mol; negative = binds better than the scl references)')
print(header)
print('-' * len(header))
for name in receptors:
    vals = [next((o['vs_scl_reference'] for o in out_rows if o['receptor'] == name and o['ligand'] == l), '')
            for l in ligands]
    print(f"{name:10s} {receptors[name]['substrate'][:12]:12s} " +
          ' '.join(f'{v:>7.2f}' if v != '' else f'{"-":>7s}' for v in vals))

# How far apart are two enzymes that belong to the SAME class? That gap is the
# floor below which a difference between receptors means nothing. The two
# short-chain references are both strictly scl, so any separation between them
# is method noise, not biology.
print()
floor = []
for ligand in ligands:
    a, b = (by_receptor[r][ligand]['affinity'] for r in SCL_REFERENCES)
    floor.append(abs(a - b))
print(f'Noise floor: the two strictly short-chain references disagree by '
      f'{st.mean(floor):.2f} kcal/mol on average (range {min(floor):.2f}-{max(floor):.2f}).')
print('Differences between receptors smaller than that should not be interpreted.')

print('\nChain length at which each receptor scores best, and the C4 to C8 shift:')
for name in receptors:
    vals = {int(o['carbons']): o['vs_scl_reference'] for o in out_rows
            if o['receptor'] == name and o['vs_scl_reference'] != ''}
    if not vals:
        continue
    best = min(vals, key=vals.get)
    shift = vals.get(8, float('nan')) - vals.get(4, float('nan'))
    print(f'  {name:10s} best at C{best:<3d} ({vals[best]:+.2f})   C4 to C8 shift {shift:+.2f} '
          f"{'(towards longer chains)' if shift < -0.3 else '(towards shorter chains)' if shift > 0.3 else '(flat)'}")

# Panel check, done against the references INDIVIDUALLY rather than their mean.
# Checking against the mean is too lenient and hides the thing that matters: if
# two enzymes of the same known class disagree with each other by as much as the
# different-class enzyme disagrees with them, the readout is not separating
# substrate preference. The shift is measured on raw affinities because the
# reference-relative values are mirror images of each other by construction
# (they are defined against their own mean) and cannot show that spread.
raw = defaultdict(dict)
for r in rows:
    raw[r['receptor']][r['carbons']] = r['affinity']
shift = {name: v[8] - v[4] for name, v in raw.items() if 4 in v and 8 in v}
scl_shifts = [shift[n] for n in SCL_REFERENCES if n in shift]
mcl_shift = shift.get(MCL_REFERENCE)

print('\nPanel check -- can the readout separate enzymes whose preference is known?')
for name in SCL_REFERENCES:
    print(f'  {name} (strictly short-chain): C4 to C8 shift {shift[name]:+.2f} kcal/mol')
print(f'  {MCL_REFERENCE} (accepts C6):        C4 to C8 shift {mcl_shift:+.2f} kcal/mol')
margin = min(scl_shifts) - mcl_shift
noise = st.mean(floor)
if mcl_shift < min(scl_shifts) and margin > noise:
    print(f'  PASSES: the C6-accepting enzyme is separated from both same-class references '
          f'by {margin:.2f}, above the {noise:.2f} noise floor.')
else:
    print(f'  FAILS: the C6-accepting enzyme leads by only {margin:.2f} kcal/mol over the nearest '
          f'short-chain reference, inside the {noise:.2f} noise floor.')
    print(f'  The two strictly short-chain references differ from each other by '
          f'{abs(scl_shifts[0] - scl_shifts[1]):.2f}, more than that margin, so this docking cannot')
    print('  reliably rank chain-length preference. Read the per-copy numbers as a hypothesis to')
    print('  test, not as a prediction, and weight the pocket-geometry result above them.')

print(f'\nwrote {OUT}')
