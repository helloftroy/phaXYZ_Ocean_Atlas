# Ray-traced superposition of the Haliangium novel phaC candidate on the
# Mesobacillus persicus PHA synthase it structurally matches
# (PHA_CLEAN_RESULTS.md section 13.4).
#
# Run headless, from the repo root:
#   /Users/hellpark/miniforge3/envs/pymol_env/bin/pymol -cq \
#       structure_prediction/render_novel_phac_superposition.pml
#
# Writes two transparent-background renders into figures/_superposition_raw/,
# which structure_prediction/compose_superposition_figure.py then labels.
#
# cealign, not align or super: at 15.9% identity within the structural
# alignment there is no sequence signal to align on, and Combinatorial
# Extension is the method meant for exactly that case. It also gives the
# transform independently of the Kabsch/Needleman-Wunsch one computed in
# superpose_novel_phac_example.py, so the two serve as a check on each other.

load RSGB23-1_GCF-000024805-V1_GENO_10000001-scaffold_1_5425.pdb, candidate
load references.pdb, reference

cealign reference, candidate

hide everything
show cartoon, candidate or reference
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_transparency, 0.30, reference

set_color cand_blue, [0.110, 0.361, 0.670]
set_color ref_orange, [0.761, 0.384, 0.176]
set_color triad_teal, [0.051, 0.580, 0.533]
color cand_blue, candidate
color ref_orange, reference

# Catalytic triads. Residue numbers are the PDB numbering, verified against
# the residue names: candidate Cys219/Asp403/His432, reference
# Cys134/Asp290/His319.
select cand_triad, candidate and resi 219+403+432
select ref_triad, reference and resi 134+290+319
show sticks, (cand_triad or ref_triad) and not (name C+N+O)
color triad_teal, cand_triad or ref_triad
util.cnc("(cand_triad or ref_triad) and not elem C")
set stick_radius, 0.22
set valence, 0

bg_color white
set ray_opaque_background, 0
set antialias, 2
set orthoscopic, on
set depth_cue, 0
set ray_shadows, 0
set ambient, 0.22
set direct, 0.50
set reflect, 0.22
set specular, 0.18
set spec_power, 60
set cartoon_highlight_color, grey70

# cealign's own RMSD and alignment length, for the record and as a check on
# the independent Kabsch superposition in superpose_novel_phac_example.py.
python
from pymol import cmd
result = cmd.cealign('reference', 'candidate')
print('CEALIGN RMSD %.3f over %d aligned residues' % (result['RMSD'], result['alignment_length']))
for a, b, lab in ((219, 134, 'Cys'), (403, 290, 'Asp'), (432, 319, 'His')):
    d = cmd.get_distance('candidate and resi %d and name CA' % a, 'reference and resi %d and name CA' % b)
    print('  %s CA-CA candidate %d to reference %d: %.2f A' % (lab, a, b, d))
print('  catalytic Cys SG-SG: %.2f A' % cmd.get_distance('candidate and resi 219 and name SG',
                                                         'reference and resi 134 and name SG'))
python end

orient reference
zoom reference, 2

png figures/_superposition_raw/view1.png, width=2000, height=1900, dpi=300, ray=1

turn y, 90
png figures/_superposition_raw/view2.png, width=2000, height=1900, dpi=300, ray=1

# Third view: the active site, to show the two triads superposing. The
# surrounding cartoon is pushed well back so the sticks read; the clipping
# planes are widened because the default slab cut through the side chains.
orient cand_triad
zoom cand_triad or ref_triad, 4.5
clip slab, 55
set cartoon_transparency, 0.90, reference
set cartoon_transparency, 0.84, candidate
set stick_radius, 0.20
set label_size, 26
set label_color, grey20
set label_font_id, 7
set label_position, (1.6, 1.4, 2.0)
label candidate and resi 219+403+432 and name CA, "%s%s" % (one_letter[resn], resi)
png figures/_superposition_raw/view3.png, width=2000, height=1900, dpi=300, ray=1
