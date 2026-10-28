"""Build the (R)-3-hydroxyacyl-N-acetylcysteamine thioester series to dock.

The physiological substrate is (R)-3-hydroxyacyl-CoA, which is far too large
and flexible to dock reliably -- coenzyme A alone carries about 25 rotatable
bonds, well past what Vina samples. The standard truncation replaces CoA with
N-acetylcysteamine, and that is not a convenience: 3-hydroxybutyryl-NAC and its
longer homologues are the actual substrate analogues used in PHA synthase
activity assays, so the fragment being docked is a real experimental compound
rather than an arbitrary cut.

Chain length is varied on the acyl tail, which is the part the hydrophobic
channel has to accommodate:

    C4  3HB   3-hydroxybutyryl     the standard short-chain-length monomer
    C5  3HV   3-hydroxyvaleryl     also short-chain-length
    C6  3HHx  3-hydroxyhexanoyl    what Aeromonas caviae PhaC accepts
    C8  3HO   3-hydroxyoctanoyl    medium-chain-length
    C10 3HD   3-hydroxydecanoyl    medium-chain-length
    C12 3HDD  3-hydroxydodecanoyl  long end of medium-chain-length

Stereochemistry is checked rather than assumed: RDKit assigns the CIP code
after embedding and the script asserts R at the 3-hydroxy carbon, because the
enzymes are strictly (R)-specific and a silently inverted centre would dock a
substrate they cannot use.

Usage:
    python structure_prediction/docking/build_ligands.py
"""
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
OUT = HERE / 'ligands'
OUT.mkdir(exist_ok=True)

# tail carbons beyond the 3-hydroxy centre -> (name, n_carbons_total)
SERIES = [
    ('C4_3HB', 'C', 4), ('C5_3HV', 'CC', 5), ('C6_3HHx', 'CCC', 6),
    ('C8_3HO', 'CCCCC', 8), ('C10_3HD', 'CCCCCCC', 10), ('C12_3HDD', 'CCCCCCCCC', 12),
]
NAC = 'SCCNC(C)=O'      # N-acetylcysteamine thioester leaving group


def build(name, tail, carbons):
    smiles = f'{tail}[C@@H](O)CC(=O){NAC}'
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f'{name}: SMILES failed to parse: {smiles}')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise SystemExit(f'{name}: 3D embedding failed')
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

    Chem.AssignStereochemistryFrom3D(mol)
    centres = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
    codes = [code for _, code in centres]
    heavy = mol.GetNumHeavyAtoms()
    rotatable = AllChem.CalcNumRotatableBonds(Chem.RemoveHs(mol))
    path = OUT / f'{name}.sdf'
    Chem.SDWriter(str(path)).write(mol)
    return codes, heavy, rotatable, path


print(f"{'ligand':10s} {'C':>3} {'heavy':>6} {'rot':>4}  {'CIP':>6}")
problems = []
for name, tail, carbons in SERIES:
    codes, heavy, rotatable, path = build(name, tail, carbons)
    ok = codes == ['R']
    if not ok:
        problems.append((name, codes))
    print(f"{name:10s} {carbons:3d} {heavy:6d} {rotatable:4d}  {','.join(codes) or '-':>6}  {'OK' if ok else 'WRONG'}")

if problems:
    print('\nStereochemistry is not (R) at the 3-hydroxy carbon for: '
          + ', '.join(f'{n} ({c})' for n, c in problems))
    print('These enzymes are strictly (R)-specific, so flip [C@@H] to [C@H] and re-run.')
    sys.exit(1)
print(f'\nwrote {len(SERIES)} ligands to {OUT}')
