"""Renders ESMFold PDB structures with the catalytic triad highlighted,
for a hand-picked set of targets: a few copies each from the three
"deep dive" multi-copy genomes (CARD22-1_SAMN24292811_MAG_00000010,
section 9.4; ZHEN20-1_SAMN07748058_MAG_00000118, the Modicisalibacter
zincidurans representative, section 9.5), a few of the most sequence-
divergent validated phaC candidates (figures/phac_divergent35_validated.tsv,
via figures/divergent35_structural_triad_audit.tsv), and two
"no plausible triad found" examples from the no_hmm_triad_support
population (figures/structural_no_hmm_triad_active_site_audit.tsv) for
visual contrast. Triad positions come from find_structural_triad.py.

Must be run under the pymol_env conda environment (pymol-open-source's
pip wheel is broken on this machine -- see the docstring note this
script's own commit message/PHA_CLEAN_RESULTS.md section leaves for the
libpng dylib issue; `conda create -n pymol_env -c conda-forge
pymol-open-source` is what actually works):

    /Users/hellpark/miniforge3/envs/pymol_env/bin/python3 \
        catalytic_domain/render_triad_structures.py

Cartoon color = AlphaFold-style per-residue pLDDT confidence (ESMFold
stores pLDDT, 0-1 scale, directly in the PDB B-factor column): dark blue
>0.9 (very high), light blue 0.7-0.9 (confident), yellow 0.5-0.7 (low),
orange <=0.5 (very low) -- the same four-band scheme AlphaFold's own
viewers use, so these are readable as "trustworthy fold" vs not on
sight. Triad residues are drawn as sticks: nucleophile (Cys/Ser/Thr)
gold, Asp red, His marine blue, each labeled with residue+number.

Outputs: figures/structures/triad_highlights/<target_id>_<tag>_full.png
(whole structure) and _active_site.png (zoomed on the triad), for every
job in JOBS below.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from find_structural_triad import find_triad_for_target, PDB_DIR

import pymol
from pymol import cmd, util

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'figures/structures/triad_highlights'
OUT.mkdir(parents=True, exist_ok=True)

# (target_id, output tag, short title for reference, nucleophile to search for)
JOBS = [
    # CARD22-1_SAMN24292811_MAG_00000010 (section 9.4)
    ('OMDBv2.0_AA_G_NR100_000041187880', 'card22_classIII_operon', 'CARD22-1: Class III, phaE-phaC-phaJ operon member', 'C'),
    ('OMDBv2.0_AA_G_NR100_000041185910', 'card22_classIII_91pct', 'CARD22-1: Class III, 91.1% identity to reference', 'C'),
    ('OMDBv2.0_AA_G_NR100_000041184878', 'card22_phaCZ_ambiguous', 'CARD22-1: 708aa, phaC/phaZ-ambiguous by sequence', 'C'),
    # ZHEN20-1_SAMN07748058_MAG_00000118 (Modicisalibacter zincidurans, section 9.5)
    ('OMDBv2.0_AA_G_NR100_000054867173', 'modici_classI_confirmed', 'Modicisalibacter: Class I, trusted-cutoff confirmed', 'C'),
    ('OMDBv2.0_AA_G_NR100_000054866076', 'modici_universal_ambiguous', 'Modicisalibacter: pan-genome-universal, phaC/phaZ-ambiguous', 'C'),
    ('OMDBv2.0_AA_G_NR100_000246008371', 'modici_738aa_atypical', 'Modicisalibacter: 738aa, structurally atypical paralog', 'C'),
    # Most sequence-divergent validated candidates (section 9.6-era; see figures/phac_divergent35_validated.tsv)
    ('OMDBv2.0_AA_G_NR100_000218417001', 'divergent_24pct_cys', '24.5% identity to nearest reference, canonical Cys triad', 'C'),
    ('OMDBv2.0_AA_G_NR100_000148739328', 'divergent_26pct_ser', '26.3% identity, no Cys -- real Ser-nucleophile triad instead', 'S'),
    ('OMDBv2.0_AA_G_NR100_000033698400', 'divergent_26pct_ser2', '26.4% identity, no Cys -- second Ser-nucleophile example', 'S'),
    # no_hmm_triad_support: no plausible triad found by any check (contrast cases)
    ('OMDBv2.0_AA_G_NR100_000000013407', 'notriad_no_cysteine', 'no_hmm_triad_support: zero cysteines, canonical triad impossible', None),
    ('OMDBv2.0_AA_G_NR100_000098889079', 'notriad_cys_far_apart', 'no_hmm_triad_support: has a Cys, but nearest His/Asp are 10-13A away', 'C'),
]

NUC_COLOR = {'C': 'yellow', 'S': 'orange', 'T': 'brightorange'}
NUC_NAME = {'C': 'Cys', 'S': 'Ser', 'T': 'Thr'}


def color_by_plddt(sel):
    # PyMOL's selection algebra has no <= / >= operator (raises "Selector-Error:
    # Invalid Number" -- confirmed live). A first attempt used compound "b > x
    # and b < y" ranges with strict inequalities on both sides, which left gaps
    # at residues sitting exactly on a bin boundary (ESMFold's B-factor pLDDT is
    # rounded to 2dp, so exact 0.50/0.70/0.90 hits are common) -- those residues
    # matched none of the four rules and fell back to PyMOL's default green,
    # confirmed live as stray green patches on an otherwise all-blue cartoon.
    # Fixed by coloring the whole selection dark blue first, then repeatedly
    # narrowing with single-sided "<" cuts so every residue is covered by
    # exactly the last (most restrictive) rule that applies -- no gaps possible.
    cmd.color('0x0053D6', sel)
    cmd.color('0x65CBF3', f'{sel} and b < 0.9')
    cmd.color('0xFFDB13', f'{sel} and b < 0.7')
    cmd.color('0xFF7D45', f'{sel} and b < 0.5')


def render_one(target_id, tag, title, nucleophile):
    pdb_path = PDB_DIR / f'{target_id}.pdb'
    if not pdb_path.exists():
        print(f'SKIP {target_id}: no PDB file')
        return

    obj = 'mol'
    cmd.reinitialize()
    cmd.load(str(pdb_path), obj)
    cmd.hide('everything')
    cmd.bg_color('white')
    cmd.show('cartoon')
    cmd.set('cartoon_transparency', 0.0)
    color_by_plddt(obj)
    cmd.set('ray_opaque_background', 1)
    cmd.set('antialias', 2)
    cmd.orient()

    triad = None
    if nucleophile is not None:
        triad = find_triad_for_target(target_id, nucleophile)

    if triad is not None:
        nuc_resi, asp_resi, his_resi = triad['nuc_pos'], triad['asp_pos'], triad['his_pos']
        nuc_col = NUC_COLOR[nucleophile]
        sel_nuc, sel_asp, sel_his = f'resi {nuc_resi}', f'resi {asp_resi}', f'resi {his_resi}'
        for sel, col in [(sel_nuc, nuc_col), (sel_asp, 'red'), (sel_his, 'marine')]:
            cmd.show('sticks', sel)
            cmd.color(col, sel)
            util.cnc(sel)
        cmd.set('stick_radius', 0.25)
        cmd.set('label_size', 22)
        cmd.set('label_color', 'black')
        cmd.set('label_outline_color', 'white')

    # no labels in the full-structure view -- at this zoom level the three
    # labels sit close enough in screen space to overlap into an unreadable
    # smear (confirmed live on the first render); labels are added only for
    # the zoomed active-site shot below, where there is room for them.
    cmd.png(str(OUT / f'{target_id}_{tag}_full.png'), width=1400, height=1400, dpi=300, ray=1)
    print(f'wrote {tag}_full.png ({title})')

    if triad is not None:
        cmd.orient(f'resi {triad["nuc_pos"]}+{triad["asp_pos"]}+{triad["his_pos"]}')
        cmd.zoom(f'resi {triad["nuc_pos"]}+{triad["asp_pos"]}+{triad["his_pos"]}', buffer=9)
        for resi, label in [(nuc_resi, f'{NUC_NAME[nucleophile]}{nuc_resi}'), (asp_resi, f'Asp{asp_resi}'), (his_resi, f'His{his_resi}')]:
            cmd.label(f'resi {resi} and name CA', f'"{label}"')
        cmd.png(str(OUT / f'{target_id}_{tag}_active_site.png'), width=1400, height=1400, dpi=300, ray=1)
        print(f'wrote {tag}_active_site.png')
        dist_txt = (f'  nuc-his={triad["nuc_his_dist_A"]:.2f}A  his-asp={triad["his_asp_dist_A"]:.2f}A  '
                    f'motif={triad["nuc_in_motif"]}  plddt(nuc/asp/his)={triad["nuc_plddt"]:.2f}/{triad["asp_plddt"]:.2f}/{triad["his_plddt"]:.2f}')
        print(dist_txt)
    else:
        print('  no triad drawn (none found or nucleophile=None requested)')


if __name__ == '__main__':
    sys.stdout = sys.stderr  # pymol -qc silently swallows real stdout (confirmed live)
    pymol.finish_launching(['pymol', '-qc'])
    for target_id, tag, title, nucleophile in JOBS:
        render_one(target_id, tag, title, nucleophile)
    print('\ndone --', len(JOBS), 'jobs')
