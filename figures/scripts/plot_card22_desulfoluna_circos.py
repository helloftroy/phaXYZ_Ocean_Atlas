"""Chord diagram for CARD22-1's two focal Desulfoluna genomes (section
9.4) plus every other Desulfoluna genome sharing their paralog clusters --
piece 5 of the multi-part request following section 9.7's structural
results, applying the method validated in plot_modicisalibacter_circos.py
to the actual target case.

Unlike Modicisalibacter, no outside-genus genomes are included here: a
direct check found all 7 of these two genomes' 70%-clusters are
Desulfoluna-only (no Cobetia/Vreelandella/Halomonas/Marinobacter-style
cross-family sharing) -- so the natural "more genomes to show" for this
case is more Desulfoluna, not another genus. That search found 9 more:
3 more CARD22-1 genomes (same Red Sea coral-tissue site, a different
unnamed species "Desulfoluna sp022360815"), 4 PELI21-1 genomes (marine
sediment, "Desulfoluna sp013619155", a third site entirely), and 2 NCBI
RefSeq isolate references (D. butyratoxydans, D. spongiiphila) -- 11
Desulfoluna genomes total, 3 studies, 3 named/unnamed species.

Structural filter applied first, per this request's own instruction:
targets with qtmscore<0.5 in structural_evidence_best_hit.tsv are
dropped before building the diagram (this dropped section 9.4's two
originally-flagged phaC/phaZ-ambiguous copies' WEAKER sibling target in
each of the two focal genomes -- 185910 and 220052, both qtm~0.16 -- but
NOT the ambiguous copies themselves, 184878/189118, which scored 0.76/0.91
and are kept; the phaC/phaZ ambiguity from section 9.4 is a separate,
still-open question from structural confidence, and these two do fold
confidently as *some* real structure even though which family they
belong to is unresolved).

Usage:
    python figures/scripts/plot_card22_desulfoluna_circos.py

Outputs:
    figures/card22_desulfoluna_circos.png / .pdf
    figures/card22_desulfoluna_circos_nodes.tsv
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc
from _circos_chart import build_circos_chart, FA

FOCAL_GENOMES = ['CARD22-1_SAMN24292811_MAG_00000010', 'CARD22-1_SAMN24292812_MAG_00000014']
OTHER_GENOMES = [
    'CARD22-1_SAMN24292799_MAG_00000026', 'CARD22-1_SAMN24292833_MAG_00000013', 'CARD22-1_SAMN24292837_MAG_00000005',
    'PELI21-1_SAMN14421543_MAG_00000002', 'PELI21-1_SAMN14421545_MAG_00000005',
    'PELI21-1_SAMN14421546_MAG_00000001', 'PELI21-1_SAMN14421547_MAG_00000003',
    'RSGB23-1_GCF-900699765-V1_GENO_10000001', 'RSGB23-1_GCF-902498735-V1_GENO_10000001',
]
GENOME_ORDER = FOCAL_GENOMES + OTHER_GENOMES

genomes_meta = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['genome'] in GENOME_ORDER:
            genomes_meta[row['genome']] = row


def short_species(g):
    sp = genomes_meta[g]['gtdb_species']
    if sp.startswith('Unknown'):
        return 'sp.?'
    return sp.replace('Desulfoluna ', '').replace(' ', '')[:14]


GENOME_LABEL = {g: f'{g.split("_")[0]} {short_species(g)} ({g.split("_")[-1][-3:]})' for g in GENOME_ORDER}

STUDY_COLOR = {
    'CARD22-1': '#DCEDEA', 'PELI21-1': '#F3E3D3', 'RSGB23-1': '#E3E9F7',
}
BAND_COLOR = {g: STUDY_COLOR[genomes_meta[g]['study_id']] for g in GENOME_ORDER}
BAND_LEGEND = [
    ('#DCEDEA', 'CARD22-1 (Red Sea coral tissue, 5 genomes, 2 species)'),
    ('#F3E3D3', 'PELI21-1 (marine sediment, 4 genomes, 1 species)'),
    ('#E3E9F7', 'RSGB23-1 (NCBI isolate references, 2 genomes, 2 species)'),
]

CLUSTER_SHORT = {}  # filled in after loading, just use last-6-digits labeling


def cluster_label(c):
    return f'paralog {c[-6:]}'


bad_targets = _phac_qc.load_bad_targets()

struct_qtm = {}
with open('structural_evidence_best_hit.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        struct_qtm[row['candidate_id']] = float(row['qtmscore'])

build_circos_chart(
    genome_order=GENOME_ORDER, genome_label=GENOME_LABEL, genome_band_color=BAND_COLOR, band_legend=BAND_LEGEND,
    cluster_label_fn=cluster_label, out_stem='card22_desulfoluna_circos', bad_targets=bad_targets,
    struct_qtm=struct_qtm, qtm_min=0.5,
    title='phaC paralogs across Desulfoluna: CARD22-1 and every other genome sharing its clusters',
    subtitle='Each sector = one genome; each dot = one structurally-confirmed (qtmscore≥0.5) phaC copy. Chords connect copies from\n'
             'different genomes sharing the same 70%-identity paralog cluster; opacity/width = exact pairwise %identity.',
)
