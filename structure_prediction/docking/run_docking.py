"""Dock the (R)-3-hydroxyacyl-NAC series into each Desulfoluna phaC copy and
into the reference synthases, and record the affinity for every pair.

Run build_ligands.py first. Needs the phac_docking conda environment, which
carries AutoDock Vina and Open Babel:

    /Users/hellpark/miniforge3/envs/phac_docking/bin/python \
        structure_prediction/docking/run_docking.py

The search box is centred on each protein's own catalytic nucleophile, located
geometrically by catalytic_domain/find_structural_triad.py rather than by
sequence position, so the box sits on the real active site in every structure
regardless of numbering. 24 A is wide enough to let a C12 chain extend without
being large enough to let the ligand wander into unrelated surface grooves.

Two things this deliberately does NOT claim. The docking is non-covalent while
the real mechanism forms a thioester with the catalytic cysteine, so a pose is
a proxy for substrate admission, not for catalysis. And affinity scales with
ligand size, so a raw score cannot be compared across chain lengths -- that is
what the reference enzymes are in the panel for, and what analyse_docking.py
corrects for.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
HERE = Path(__file__).resolve().parent
LIGANDS = HERE / 'ligands'
WORK = HERE / 'work'
WORK.mkdir(exist_ok=True)
OUT = HERE / 'docking_scores.tsv'

CONDA_BIN = Path('/Users/hellpark/miniforge3/envs/phac_docking/bin')
OBABEL = CONDA_BIN / 'obabel'
VINA = CONDA_BIN / 'vina'

BOX = 24.0
EXHAUSTIVENESS = 24
SEED = 42
N_MODES = 9

SITES = HERE / 'active_sites.json'


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def receptor_pdbqt(name, pdb_path):
    out = WORK / f'{name}.pdbqt'
    if out.exists():
        return out
    # -xr makes a rigid receptor; -p 7.4 assigns protonation at physiological pH
    r = sh([str(OBABEL), str(pdb_path), '-O', str(out), '-xr', '-p', '7.4'])
    if not out.exists() or out.stat().st_size == 0:
        raise SystemExit(f'{name}: receptor conversion failed\n{r.stderr[-800:]}')
    return out


def ligand_pdbqt(name, sdf_path):
    out = WORK / f'lig_{name}.pdbqt'
    if out.exists():
        return out
    r = sh([str(OBABEL), str(sdf_path), '-O', str(out)])
    if not out.exists() or out.stat().st_size == 0:
        raise SystemExit(f'{name}: ligand conversion failed\n{r.stderr[-800:]}')
    return out


def dock(receptor, ligand, centre, tag):
    pose = WORK / f'{tag}.pdbqt'
    cmd = [str(VINA), '--receptor', str(receptor), '--ligand', str(ligand),
           '--center_x', f'{centre[0]:.3f}', '--center_y', f'{centre[1]:.3f}',
           '--center_z', f'{centre[2]:.3f}',
           '--size_x', str(BOX), '--size_y', str(BOX), '--size_z', str(BOX),
           '--exhaustiveness', str(EXHAUSTIVENESS), '--seed', str(SEED),
           '--num_modes', str(N_MODES), '--out', str(pose)]
    r = sh(cmd)
    best = None
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == '1':
            try:
                best = float(parts[1])
            except ValueError:
                pass
            break
    if best is None:
        print(f'  {tag}: no score parsed\n{r.stdout[-400:]}{r.stderr[-400:]}')
    return best


def heavy_atoms(sdf_path):
    """Heavy-atom count from the SDF counts line, used to flag the size confound
    that analyse_docking.py corrects for."""
    lines = sdf_path.read_text().splitlines()
    n_atoms = int(lines[3][:3])
    elements = [lines[4 + i][31:34].strip() for i in range(n_atoms)]
    return sum(1 for e in elements if e != 'H')


sites = json.loads(SITES.read_text())
ligands = sorted(LIGANDS.glob('*.sdf'), key=lambda p: int(p.stem.split('_')[0][1:]))
HEAVY = {p.stem: heavy_atoms(p) for p in ligands}
rows = []
print(f'{len(sites)} receptors x {len(ligands)} ligands')
for name, info in sites.items():
    receptor = receptor_pdbqt(name, ROOT / info['pdb'])
    centre = info['centre']
    for sdf in ligands:
        lig = ligand_pdbqt(sdf.stem, sdf)
        affinity = dock(receptor, lig, centre, f'{name}__{sdf.stem}')
        rows.append(dict(receptor=name, kind=info['kind'], label=info['label'],
                         substrate=info.get('substrate', ''), ligand=sdf.stem,
                         carbons=int(sdf.stem.split('_')[0][1:]),
                         heavy_atoms=HEAVY[sdf.stem],
                         affinity_kcal_mol=affinity))
        print(f'  {name:10s} {sdf.stem:10s} {affinity if affinity is not None else "FAIL"}')

with open(OUT, 'w', newline='') as f:
    w = csv.DictWriter(f, delimiter='\t', fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f'\nwrote {OUT}')
