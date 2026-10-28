"""Locate each receptor's catalytic nucleophile and write the docking box centres.

The box is centred on the nucleophile's own side-chain atom, found geometrically
rather than by sequence position, so it lands on the real active site in crystal
structures and predicted models alike.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / 'catalytic_domain'))
sys.path.insert(0, str(ROOT / 'structure_prediction'))
from analyze_phac_substrate_pocket import catalytic_sulfur  # noqa: E402
from anchor_phac_pocket_to_experiment import single_chain  # noqa: E402
from find_structural_triad import load_sequences  # noqa: E402

OUT = Path(__file__).resolve().parent / 'active_sites.json'

COPIES = [('184878', 'I', 708), ('185226', 'I', 578), ('185910', 'III', 577),
          ('187467', 'I', 541), ('187468', 'III', 380), ('187880', 'III', 372),
          ('189118', 'I', 736)]
REFS = [('5HZ2', 'Cupriavidus necator PhaC1', 'scl only'),
        ('5XAV', 'Chromobacterium sp. USM2 PhaC', 'scl only'),
        ('9KNK', 'Aeromonas caviae PhaC', 'scl + C6')]

sites = {}
for short, model_class, length in COPIES:
    target = f'OMDBv2.0_AA_G_NR100_000041{short}'
    pdb = ROOT / f'{target}.pdb'
    if not pdb.exists():
        print(f'{short}: missing'); continue
    seq = load_sequences([target]).get(target, '')
    centre, triad, nuc = catalytic_sulfur(pdb, seq)
    tight = triad['nuc_his_dist_A'] <= 4.5 and triad['his_asp_dist_A'] <= 4.5
    sites[short] = dict(pdb=pdb.name, kind='model', label=f'{short} ({length} aa, class {model_class})',
                        substrate='unknown', centre=[round(float(v), 3) for v in centre],
                        triad=f"{nuc}{triad['nuc_pos']}/D{triad['asp_pos']}/H{triad['his_pos']}",
                        triad_closed=bool(tight))
    print(f"{short:8s} {nuc}{triad['nuc_pos']:<5d} closed={tight}  centre {centre.round(2)}")

for pdb_id, organism, substrate in REFS:
    src = ROOT / f'structure_prediction/reference_pdb/{pdb_id}.pdb'
    chain, seq = single_chain(src)
    centre, triad, nuc = catalytic_sulfur(chain, seq)
    tight = triad['nuc_his_dist_A'] <= 4.5 and triad['his_asp_dist_A'] <= 4.5
    sites[pdb_id] = dict(pdb=str(chain.relative_to(ROOT)), kind='crystal',
                         label=f'{pdb_id} — {organism}', substrate=substrate,
                         centre=[round(float(v), 3) for v in centre],
                         triad=f"{nuc}{triad['nuc_pos']}/D{triad['asp_pos']}/H{triad['his_pos']}",
                         triad_closed=bool(tight))
    print(f"{pdb_id:8s} {nuc}{triad['nuc_pos']:<5d} closed={tight}  centre {centre.round(2)}")

OUT.write_text(json.dumps(sites, indent=2))
print(f'\nwrote {OUT} ({len(sites)} receptors)')
