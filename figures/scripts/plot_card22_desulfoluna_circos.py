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

Copies shown are the triad-complete ones, the same rule the identity
heatmap and the contig map use, so all three panels of the figure report
the same number of phaC per genome (see _deep_dive_copies). An earlier
version also required qtmscore>=0.5, which removed 185910 and 220052
(both qtm~0.16) from the two CARD22-1 genomes. That filter is not applied
any more: 185910 is triad-complete, class III, and 91.1% identical over
99% coverage to a characterised phaC, so a qTM that low indicts the
structural comparison rather than the protein -- and it was being dropped
while 184878 and 189118 were kept, which are the copies whose phaC-versus-
phaZ assignment is actually unresolved.

Only the deep-dive genome itself is outlined. GENOME_ORDER still leads
with both CARD22-1 genomes from section 9.4, but outlining both made the
panel look like it had two subjects.

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

# The one genome this figure is about -- the only sector drawn with an outline.
DEEP_DIVE_GENOME = 'CARD22-1_SAMN24292811_MAG_00000010'
# Both section 9.4 genomes still lead the ring; the second is context, not a
# second subject.
FOCAL_GENOMES = [DEEP_DIVE_GENOME, 'CARD22-1_SAMN24292812_MAG_00000014']
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
    """Species alone, no genus: every genome here is Desulfoluna, so repeating it
    11 times around the ring only costs label width. Unnamed GTDB placeholders
    (Desulfoluna sp022360815) keep their 6-digit stem, which is what tells the
    two unnamed species apart."""
    sp = genomes_meta[g]['gtdb_species']
    if sp.startswith('Unknown'):
        return 'sp.'
    sp = sp.replace('Desulfoluna ', '')
    if sp.startswith('sp'):
        return 'sp' + sp[2:8]
    return sp[:9] + '.' if len(sp) > 10 else sp


# Study prefix, species, and the MAG's own 3-digit suffix. The full labels ran
# off both sides of the figure at the larger font the panel needs.
GENOME_LABEL = {g: f'{g.split("_")[0][:6]} {short_species(g)} {g.split("_")[-1][-3:]}' for g in GENOME_ORDER}

STUDY_COLOR = {
    'CARD22-1': '#4FA8A0', 'PELI21-1': '#E2954F', 'RSGB23-1': '#5D8FD1',
}
BAND_COLOR = {g: STUDY_COLOR[genomes_meta[g]['study_id']] for g in GENOME_ORDER}
BAND_LEGEND = [
    ('#4FA8A0', 'CARD22-1 (Red Sea coral tissue, 5 genomes, 2 species)'),
    ('#E2954F', 'PELI21-1 (marine sediment, 4 genomes, 1 species)'),
    ('#5D8FD1', 'RSGB23-1 (NCBI isolate references, 2 genomes, 2 species)'),
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
    struct_qtm=struct_qtm, qtm_min=0.0, require_triad_complete=True, focal_genomes={DEEP_DIVE_GENOME},
    title='phaC paralogs across Desulfoluna: CARD22-1 and every other genome sharing its clusters',
    subtitle='Each sector = one genome; each dot = one triad-complete phaC copy. Chords connect copies from\n'
             'different genomes sharing the same 70%-identity paralog cluster; opacity/width = exact pairwise %identity.',
)
