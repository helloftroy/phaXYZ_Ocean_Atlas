"""Geometric catalytic-triad finder, working directly from an ESMFold PDB
structure + its own raw sequence -- no HMM alignment involved. Written to
fill a real gap: the two existing structural audit tables (figures/
structural_no_hmm_triad_active_site_audit.tsv and
_no_cys_alternative_triad_audit.tsv) already did exactly this for the
19 no_hmm_triad_support candidates, but the script that produced them is
not present in this repo -- this reimplements the same idea (found by
reverse-engineering their columns) as a small, reusable, generic
function, so it can also be run on sequences that DO have an HMM-called
triad (the deep-dive genomes' own copies, the divergent35 set) to get
each triad's actual 1-based sequence position for structure rendering,
which no existing table records.

Method: scan the raw sequence for the G-x-Cys-x-G lipase-box motif
(same regex as catalytic_domain/scan_lipase_box_raw.py); if none match,
fall back to every Cys in the sequence. For every (Cys, Asp, His)
combination, measure the two catalytically meaningful hydrogen-bond
distances (Cys SG<->His NE2/ND1, His NE2/ND1<->Asp OD1/OD2, falling back
to CA-CA if a side-chain atom is missing -- ESMFold output has full side
chains, so this is rare) and keep the combination that minimizes their
sum. This is a plausibility search, not a guarantee -- a real catalytic
triad is a short, specific hydrogen-bonded chain (His-Asp usually
~2.5-3.5A in a genuine active site), so a "best" combination that is
still tens of Angstroms apart means no real triad exists in that
structure, which is itself the answer for the no-triad population.
"""
import re
from pathlib import Path

MOTIF = re.compile(r'G.C.G')

ROOT = Path(__file__).resolve().parent.parent
PDB_DIR = ROOT / 'PDB_files'
SEQ_FAA = ROOT / 'PHA_bioprospecting/omdb_search/results/phaC_cluster_sequences.faa'


def load_sequences(target_ids):
    wanted = set(target_ids)
    seqs = {}
    cur_id, cur_seq = None, []
    with open(SEQ_FAA) as f:
        for line in f:
            if line.startswith('>'):
                if cur_id in wanted:
                    seqs[cur_id] = ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id in wanted:
            seqs[cur_id] = ''.join(cur_seq)
    return seqs


def parse_pdb_residues(path):
    """resi (1-based int) -> {'resn': str, 'plddt': float, atoms: {name: (x,y,z)}}"""
    residues = {}
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            resn = line[17:20].strip()
            resi = int(line[22:26])
            atom = line[12:16].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            b = float(line[60:66])
            d = residues.setdefault(resi, {'resn': resn, 'atoms': {}, 'plddt': None})
            d['atoms'][atom] = (x, y, z)
            if atom == 'CA':
                d['plddt'] = b
    return residues


def _dist(a, b):
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5


def _atom_dist(res_i, res_j, names_i, names_j):
    """min distance over the given candidate atom names on each side, CA fallback"""
    best = None
    for ni in names_i:
        if ni not in res_i['atoms']:
            continue
        for nj in names_j:
            if nj not in res_j['atoms']:
                continue
            d = _dist(res_i['atoms'][ni], res_j['atoms'][nj])
            if best is None or d < best:
                best = d
    if best is None:
        best = _dist(res_i['atoms']['CA'], res_j['atoms']['CA'])
    return best


NUCLEOPHILE_ATOMS = {'C': ['SG'], 'S': ['OG'], 'T': ['OG1']}


def find_triad(seq, residues, nucleophile='C'):
    """Returns a dict describing the best-scoring (nucleophile, Asp, His)
    combination, or None if the sequence has no residue of the requested
    nucleophile type at all. nucleophile='C' is the canonical PhaC search
    (Cys); 'S'/'T' repeats the same geometric search for a Ser/Thr
    nucleophile instead, the same alternative-mechanism check
    figures/structural_no_cys_alternative_triad_audit.tsv already ran for
    the 15 no-cysteine no_hmm_triad_support candidates.
    """
    nuc_atoms = NUCLEOPHILE_ATOMS[nucleophile]
    motif_nuc = {m.start() + 3 for m in re.finditer('G.[' + nucleophile + '].G', seq)}  # 1-based motif position
    nuc_pos = [i + 1 for i, c in enumerate(seq) if c == nucleophile]
    if not nuc_pos:
        return None
    asp_pos = [i + 1 for i, c in enumerate(seq) if c == 'D']
    his_pos = [i + 1 for i, c in enumerate(seq) if c == 'H']
    if not asp_pos or not his_pos:
        return None

    best = None
    for c in nuc_pos:
        if c not in residues:
            continue
        for h in his_pos:
            if h not in residues:
                continue
            nuc_his = _atom_dist(residues[c], residues[h], nuc_atoms, ['NE2', 'ND1'])
            for d in asp_pos:
                if d not in residues:
                    continue
                his_asp = _atom_dist(residues[h], residues[d], ['NE2', 'ND1'], ['OD1', 'OD2'])
                score = nuc_his + his_asp
                if best is None or score < best['score']:
                    nuc_asp = _atom_dist(residues[c], residues[d], nuc_atoms, ['OD1', 'OD2'])
                    best = {
                        'score': score, 'nuc_pos': c, 'asp_pos': d, 'his_pos': h,
                        'nuc_in_motif': c in motif_nuc,
                        'nuc_his_dist_A': nuc_his, 'his_asp_dist_A': his_asp, 'nuc_asp_dist_A': nuc_asp,
                        'order_nuc_asp_his': c < d < h,
                        'nuc_plddt': residues[c]['plddt'], 'asp_plddt': residues[d]['plddt'], 'his_plddt': residues[h]['plddt'],
                    }
    return best


def find_triad_for_target(target_id, nucleophile='C'):
    pdb_path = PDB_DIR / f'{target_id}.pdb'
    seqs = load_sequences([target_id])
    if target_id not in seqs or not pdb_path.exists():
        return None
    residues = parse_pdb_residues(pdb_path)
    return find_triad(seqs[target_id], residues, nucleophile=nucleophile)


if __name__ == '__main__':
    import sys
    for tid in sys.argv[1:]:
        r = find_triad_for_target(tid)
        print(tid, r)
