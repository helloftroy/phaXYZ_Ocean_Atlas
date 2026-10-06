"""One identity per phaC copy in a multi-copy deep-dive genome: a short label
and a colour, shared by every panel of the deep-dive figure.

The CARD22-1 Desulfoluna panels disagreed with each other about how many phaC
copies the genome has, because each applied a different filter:

  circos          triad-complete AND qTM >= 0.5   5 copies
  identity heatmap  triad-complete (hybrid)       6 copies
  contig map        QC-passing only               7 copies

The unified rule is triad-complete, giving 6. The two filters that were adding
and removing a copy on top of that are both wrong for this figure:

  187467 (541 aa) fails the triad check by both the structural geometry and the
    alignment columns, so it is not a functional synthase and the contig map
    should not have been showing it alongside the others.
  185910 (577 aa) was dropped only from the circos, on qTM 0.16. It is
    triad-complete, class III, and 91.1% identical over 99% coverage to a
    characterised phaC reference. A qTM that low against sequence evidence that
    strong points at the structural comparison, not at the protein; dropping a
    copy on it while keeping 184878, whose phaC-vs-phaZ assignment is genuinely
    unresolved, inverted the strength of the evidence.

Copies are ordered by scaffold then gene index, so the label order follows the
genome rather than the target id, and coloured from the sequential blue ramp
used across the manuscript figures.

Usage:
    from _deep_dive_copies import copies_for, CopyStyle
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

# Same ramp as plot_phac_cluster_occupancy_spectrum.py, extended to 6 steps.
COPY_COLORS = ['#9ec5f4', '#6da7ec', '#2a78d6', '#1c5cab', '#0d366b', '#C2622D']

PARALOG_TABLES = {
    'CARD22-1_SAMN24292811_MAG_00000010': ROOT / 'figures/card22_desulfoluna_phac_paralogs.tsv',
}


class CopyStyle:
    __slots__ = ('target_id', 'short', 'color', 'length_aa', 'scaffold', 'gene_index', 'cls', 'triad')

    def __init__(self, target_id, short, color, length_aa, scaffold, gene_index, cls, triad):
        self.target_id = target_id
        self.short = short
        self.color = color
        self.length_aa = length_aa
        self.scaffold = scaffold
        self.gene_index = gene_index
        self.cls = cls
        self.triad = triad

    def __repr__(self):
        return f'<CopyStyle {self.short} {self.length_aa}aa {self.scaffold}:{self.gene_index}>'


def has_copies(genome):
    """True when this genome has a registered paralog table to style from."""
    table = PARALOG_TABLES.get(genome)
    return table is not None and table.exists()


def copies_for(genome, triad_complete_only=True):
    """Ordered CopyStyle list for one deep-dive genome's phaC copies, or an
    empty list for a genome with no registered paralog table -- the other
    deep-dive genomes reuse these scripts and simply fall back to plain
    labels."""
    import _triad_filter

    if not has_copies(genome):
        return []
    table = PARALOG_TABLES[genome]
    rows = list(csv.DictReader(open(table, newline=''), delimiter='\t'))
    keep = _triad_filter.load_triad_complete_set([r['target_id'] for r in rows])
    if triad_complete_only:
        rows = [r for r in rows if r['target_id'] in keep]
    rows.sort(key=lambda r: (int(r['own_scaffold'].split('_')[-1]), int(r['own_gene_index'])))
    out = []
    for i, r in enumerate(rows):
        out.append(CopyStyle(
            target_id=r['target_id'],
            short=r['target_id'][-6:],
            color=COPY_COLORS[i % len(COPY_COLORS)],
            length_aa=int(r['protein_length_aa']),
            scaffold=r['own_scaffold'],
            gene_index=int(r['own_gene_index']),
            cls=r['best_class_relative'],
            triad=r['target_id'] in keep,
        ))
    return out


if __name__ == '__main__':
    for genome in PARALOG_TABLES:
        print(genome)
        for c in copies_for(genome):
            print(f'  {c.short}  {c.length_aa:>4}aa  {c.scaffold}:{c.gene_index:<5} {c.cls:<9} {c.color}')
