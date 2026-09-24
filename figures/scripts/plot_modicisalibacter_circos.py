"""Chord diagram (circos-style, not a genome-circular plot): one sector
per genome, that genome's own phaC copies placed as nodes around its
sector, and chords connecting phaC copies from DIFFERENT genomes that
belong to the same 70%-identity paralog cluster -- i.e. "this genome's
copy of paralog X connects to that genome's copy of the same paralog X".
Chord opacity/width scales with the exact pairwise %identity between the
two specific proteins (computed directly, not read off the cluster
threshold); node color encodes which paralog cluster a protein belongs to
(the same 70%-cluster framework used throughout section 9). Prototype
case, picked deliberately small and already well-understood:
Modicisalibacter zincidurans (section 9.5's 6 genomes) plus 4 "outside"
genomes -- one representative each from Cobetia, Vreelandella, Halomonas,
and Marinobacter, the four genera section 9.5's addendum found sharing
Modicisalibacter's paralog clusters. The actual chord-diagram-building
logic lives in _circos_chart.py (factored out once this prototype's
layout issues -- label rotation, legend placement, duplicate genome
labels -- were found and fixed here); this script just defines the
genome set. Same method is applied to HK1/CARD22-1 in
plot_hk1_circos.py/plot_card22_desulfoluna_circos.py.

Outside-genome selection: for each of the two Modicisalibacter clusters
with real cross-genus membership (...602748, Halomonadaceae-specific;
...783925, the broader cross-family lineage), the highest-completeness/
lowest-contamination genome of each contributing genus was picked (ties
broken toward completeness) -- KOPF15-1_SAMEA7392432_MAG_00000114
(Cobetia), OBRI24-1_SAMEA115094377_MAG_00000066 (Vreelandella),
RSGB23-1_GCF-003547075-V1_GENO_10000001 (Halomonas, an NCBI isolate),
SCHR19-1_SAMN09405896_MAG_00000020 (Marinobacter). Each contributes its
FULL phaC complement to the diagram, not just the one target that placed
it in a shared cluster -- so their own additional, Modicisalibacter-
unconnected paralogs show up too (5 more distinct clusters appear this
way, unconnected singletons in the diagram, which is itself informative:
Modicisalibacter's own paralog repertoire is not simply a subset of what
these other genera carry).

Usage:
    python figures/scripts/plot_modicisalibacter_circos.py

Outputs:
    figures/modicisalibacter_circos.png / .pdf
    figures/modicisalibacter_circos_nodes.tsv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc
from _circos_chart import build_circos_chart

MODZ_GENOMES = [
    'DUAR15-1_SAMN06266142_MAG_00000002', 'RSGB23-1_GCF-000731955-V1_GENO_10000001',
    'SANC23-1_SAMEA110646468_MAG_00000065', 'TARA_SAMEA2623562_MAG_00000025',
    'TARA_SAMEA2623564_MAG_00000048', 'ZHEN20-1_SAMN07748058_MAG_00000118',
]
OUTSIDE_GENOMES = {
    'Cobetia': 'KOPF15-1_SAMEA7392432_MAG_00000114',
    'Vreelandella': 'OBRI24-1_SAMEA115094377_MAG_00000066',
    'Halomonas': 'RSGB23-1_GCF-003547075-V1_GENO_10000001',
    'Marinobacter': 'SCHR19-1_SAMN09405896_MAG_00000020',
}
GENOME_ORDER = MODZ_GENOMES + list(OUTSIDE_GENOMES.values())
GENOME_LABEL = {g: f'M.z. {g.split("_")[0]}-{g.split("_")[1][-4:]}' for g in MODZ_GENOMES}
GENOME_LABEL.update({g: f'{genus} ({g.split("_")[0]})' for genus, g in OUTSIDE_GENOMES.items()})

BAND_COLOR = {**{g: '#DCEDEA' for g in MODZ_GENOMES},
              OUTSIDE_GENOMES['Cobetia']: '#F3E3D3', OUTSIDE_GENOMES['Vreelandella']: '#EFD9E8',
              OUTSIDE_GENOMES['Halomonas']: '#E3E9F7', OUTSIDE_GENOMES['Marinobacter']: '#EAF0DC'}
BAND_LEGEND = [('#DCEDEA', 'Modicisalibacter zincidurans (6 genomes)'), ('#F3E3D3', 'Cobetia (outside genus)'),
               ('#EFD9E8', 'Vreelandella (outside genus)'), ('#E3E9F7', 'Halomonas (outside genus)'),
               ('#EAF0DC', 'Marinobacter (outside genus)')]

KNOWN_CLUSTER_NAME = {
    'OMDBv2.0_AA_G_NR100_000148602748': 'paralog 602748 (Class I)',
    'OMDBv2.0_AA_G_NR100_000139783925': 'paralog 783925 (phaC/phaZ-ambig.)',
    'OMDBv2.0_AA_G_NR100_000061282509': 'paralog 282509 (98.8% struct. id.)',
    'OMDBv2.0_AA_G_NR100_000139753017': 'paralog 753017 (738aa, atypical)',
}


def cluster_label(c):
    return KNOWN_CLUSTER_NAME.get(c, f'other ({c[-6:]})')


bad_targets = _phac_qc.load_bad_targets()
build_circos_chart(
    genome_order=GENOME_ORDER, genome_label=GENOME_LABEL, genome_band_color=BAND_COLOR, band_legend=BAND_LEGEND,
    cluster_label_fn=cluster_label, out_stem='modicisalibacter_circos', bad_targets=bad_targets,
    title='phaC paralogs across Modicisalibacter zincidurans and its closest outside relatives',
    subtitle='Each sector = one genome; each dot = one phaC copy. Chords connect copies from different genomes '
             'sharing the same 70%-identity paralog cluster;\nchord opacity/width = exact pairwise %identity between those two proteins.',
)
