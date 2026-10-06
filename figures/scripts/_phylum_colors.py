"""One phylum colour per phylum, shared by every figure in the taxonomy panel.

The radial tree, the prevalence-vs-contribution bubble plot and the synthase
class bars are assembled into a single manuscript figure, so a phylum has to
look the same in all three. The tree already generated colours by golden-angle
hue stepping, but it did so inline from its own ranked phylum list, which made
the mapping impossible to reuse and liable to shift if that ranking changed.
That generator now lives here, the ranking it keys off is written out as a TSV
so the assignment is inspectable, and the other two figures import it.

Ranking is by how many single-phylum phaC 70%-identity clusters a phylum holds,
which is what the tree is drawn from, so the tree's own appearance is unchanged
by the move.

Usage:
    from _phylum_colors import color_for, phylum_order
"""
import colorsys
import csv
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
TABLE = ROOT / 'figures/phylum_color_assignment.tsv'

OTHER_COLOR = '#7A7A7A'
GOLDEN_CONJUGATE = 0.6180339887498949


@lru_cache(maxsize=1)
def _ranked_phyla():
    """Phyla ordered by single-phylum phaC cluster count, descending."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _phac_qc

    assignments = {}
    with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
        for row in csv.reader(f, delimiter='\t'):
            assignments[row[1]] = row[0]

    cluster_phyla = {}
    with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if _phac_qc.is_bad(row['best_query']):
                continue
            cluster_id = assignments.get(row['target_id'])
            if cluster_id is None or not row.get('gtdb_phylum'):
                continue
            cluster_phyla.setdefault(cluster_id, set()).add(row['gtdb_phylum'])

    totals = {}
    for phyla in cluster_phyla.values():
        if len(phyla) == 1:
            phylum = next(iter(phyla))
            totals[phylum] = totals.get(phylum, 0) + 1
    # count descending, then name, so ties cannot reorder between runs
    return [p for p, _ in sorted(totals.items(), key=lambda kv: (-kv[1], kv[0]))], totals


@lru_cache(maxsize=1)
def _palette():
    """Golden-angle hue stepping, so each phylum differs from its NEIGHBOURS
    rather than merely being unique. Evenly spaced hues (hue = i/n) read as a
    gradient instead of a categorical palette here, because rank-adjacent phyla
    are usually also spatially adjacent in the tree -- both are sorted by the
    same count -- so neighbouring wedges came out nearly the same hue. Stepping
    by the golden ratio conjugate scatters hues around the wheel regardless of
    rank while staying fully deterministic."""
    ranked, _ = _ranked_phyla()
    colors, hue = {}, 0.15
    for i, phylum in enumerate(ranked):
        hue = (hue + GOLDEN_CONJUGATE) % 1.0
        saturation = 0.58 if i % 2 == 0 else 0.72
        value = 0.72 if i % 3 != 0 else 0.62
        colors[phylum] = colorsys.hsv_to_rgb(hue, saturation, value)
    return colors


def color_for(phylum):
    """RGB tuple for a phylum; grey for anything not in the tree."""
    return _palette().get(phylum, OTHER_COLOR)


def hex_for(phylum):
    rgb = color_for(phylum)
    if isinstance(rgb, str):
        return rgb
    return '#%02X%02X%02X' % tuple(int(round(255 * c)) for c in rgb)


def phylum_order():
    """Phyla ranked by single-phylum cluster count, descending."""
    return list(_ranked_phyla()[0])


def write_table():
    ranked, totals = _ranked_phyla()
    with open(TABLE, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['rank', 'phylum', 'n_single_phylum_clusters', 'hex'])
        for i, phylum in enumerate(ranked, start=1):
            writer.writerow([i, phylum, totals[phylum], hex_for(phylum)])
    return TABLE


if __name__ == '__main__':
    ranked, totals = _ranked_phyla()
    print(f'{len(ranked)} phyla')
    for phylum in ranked[:12]:
        print(f'  {phylum:24s} {totals[phylum]:5,} clusters  {hex_for(phylum)}')
    print('wrote', write_table())
