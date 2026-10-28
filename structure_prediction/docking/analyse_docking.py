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

mcl = {o['ligand']: o['vs_scl_reference'] for o in out_rows if o['receptor'] == MCL_REFERENCE}
if mcl:
    short = mcl.get('C4_3HB')
    mid = min(v for k, v in mcl.items() if k in ('C6_3HHx', 'C8_3HO'))
    verdict = 'PASSES' if mid < short else 'FAILS'
    print(f"\nPanel check: {MCL_REFERENCE} (accepts C6) scores {short:+.2f} at C4 and {mid:+.2f} at its best "
          f"of C6/C8 relative to the short-chain references -- {verdict}.")
    if verdict == 'FAILS':
        print('The readout does not separate the enzyme with known medium-chain activity from the ones\n'
              'without it, so the Desulfoluna numbers above should not be read as chain-length preference.')
print(f'\nwrote {OUT}')
