"""Generate the data files behind the PHA Ocean Atlas explorer site.

Everything the site displays is produced here from the canonical tables, so no
number is ever typed into the HTML by hand and a rebuild cannot drift from the
pipeline. The site is static: these files are plain CSV/JSON, GitHub Pages
gzips them in transit, and all filtering happens in the browser.

Three levels of granularity, because the site explores at three levels:

  genomes.csv    one row per PHA-carrying genome (the atlas backbone) -- carries
                 coordinates, taxonomy, habitat, depth and per-family copy number,
                 so "show me phaJ in Myxococcota below 1000 m" is answerable.
  proteins.csv   one row per QC-passing phaC protein, with the annotations that
                 only phaC has (class, triad, structure, thermostability).
  clusters.csv   one row per 70%-identity phaC cluster.
  occurrences.csv  protein x genome, as (protein, cluster, genome row index).
                 This is what lets the map plot a protein's or a cluster's
                 occurrences: coordinates are read from the genome row rather
                 than repeated per occurrence, which keeps the file at three
                 integers a row.

Repeated strings (phylum, genus, habitat, study) are replaced by integer codes
with lookup tables in lookups.json; this roughly halves the payload and costs
one array index in the browser.

A deliberate correction is applied here: the genome x family matrix is built
from a metadata join capped at 5 genomes per identical protein, which records
7,298 real phaC carriers as phaC-negative. The phaC presence column is therefore
taken from the uncapped membership list (38,744 genomes), not the matrix. Every
other family keeps the matrix value and is a lower bound -- stated in meta.json
and surfaced in the site's UI.

Usage:
    python site_build/build_site_data.py
"""
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'PHA_bioprospecting/omdb_search/results'
OUT = ROOT / 'docs/data'
sys.path.insert(0, str(ROOT / 'figures/scripts'))
import _phac_qc  # noqa: E402

csv.field_size_limit(1 << 30)
OUT.mkdir(parents=True, exist_ok=True)

TARGET_PREFIX = 'OMDBv2.0_AA_G_NR100_'
FAMILIES = ['phaA', 'phaB', 'phaC', 'phaE', 'phaR_synthase', 'phaJ', 'phaG',
            'phaZ', 'phaY', 'phaP', 'phaF', 'phaI', 'phaM', 'phaD',
            'phaR_regulator', 'phaQ']
TEMP_THRESHOLDS = [40, 45, 50, 55, 60, 65, 70, 75, 80]


def short(target_id):
    """Drop the catalogue prefix; the numeric tail is unique and much shorter."""
    return target_id[len(TARGET_PREFIX):].lstrip('0') or '0' \
        if target_id.startswith(TARGET_PREFIX) else target_id


class Codes:
    """Interns repeated strings as integers, preserving first-seen order."""

    def __init__(self):
        self.index, self.labels = {}, []

    def __call__(self, value):
        value = (value or '').strip()
        if value in ('', 'NA', 'nan', 'None', 'Unknown'):
            return ''
        code = self.index.get(value)
        if code is None:
            code = self.index[value] = len(self.labels)
            self.labels.append(value)
        return str(code)


def num(value, places=None):
    """Format a numeric field, or '' when absent -- never a literal 'nan'."""
    value = (value or '').strip()
    if value in ('', 'NA', 'nan', 'None'):
        return ''
    try:
        f = float(value)
    except ValueError:
        return ''
    if f != f:
        return ''
    if places is None:
        return str(int(f)) if f == int(f) else str(f)
    return f'{f:.{places}f}'.rstrip('0').rstrip('.')


phylum, genus, habitat, ecosystem, study, zone, arch, reference = (Codes() for _ in range(8))

print('reading QC exclusions')
bad = _phac_qc.load_bad_targets()
print(f'  {len(bad):,} excluded targets')

# ---------------------------------------------------------------- phaC truth
print('reading uncapped phaC membership')
phac_genomes = set()
protein_genomes = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r)
    for target_id, genome in r:
        if target_id in bad:
            continue
        phac_genomes.add(genome)
        protein_genomes[target_id].add(genome)
print(f'  {len(phac_genomes):,} phaC-positive genomes, {len(protein_genomes):,} proteins')

# ---------------------------------------------------------------- genomes
print('writing genomes.csv')
GENOME_COLS = ['id', 'lat', 'lon', 'ph', 'gn', 'hb', 'eco', 'st', 'dep', 'dz', 'fam', 'arch', 'phac']
names, n_genomes, with_coords = [], 0, 0
family_genomes = Counter()
with open(RESULTS / 'genome_family_matrix.tsv', newline='') as f, \
        open(OUT / 'genomes.csv', 'w', newline='') as out:
    w = csv.writer(out, lineterminator='\n')
    w.writerow(GENOME_COLS)
    for row in csv.DictReader(f, delimiter='\t'):
        genome = row['genome']
        counts = []
        for fam in FAMILIES:
            n = row.get(f'n_{fam}', '0') or '0'
            try:
                n = int(float(n))
            except ValueError:
                n = 0
            if fam == 'phaC':
                n = max(n, 1) if genome in phac_genomes else 0
            if n:
                family_genomes[fam] += 1
            counts.append(str(min(n, 9)))
        lat, lon = num(row['latitude_degN'], 3), num(row['longitude_degE'], 3)
        with_coords += bool(lat and lon)
        w.writerow([n_genomes, lat, lon, phylum(row['gtdb_phylum']), genus(row['gtdb_genus']),
                    habitat(row['ecosystem_compartment']), ecosystem(row['ecosystem_type']),
                    study(row['study_id']), num(row.get('depth_m'), 1),
                    zone(row.get('depth_zone')), ''.join(counts),
                    arch(row.get('architecture')), '1' if genome in phac_genomes else '0'])
        names.append(genome)
        n_genomes += 1
(OUT / 'genome_names.csv').write_text('\n'.join(names) + '\n')
print(f'  {n_genomes:,} genomes, {with_coords:,} with coordinates')

# phaC carriers absent from the matrix entirely: real, but no metadata row
matrix_genomes = set(names)
unlisted = len(phac_genomes - matrix_genomes)
print(f'  phaC-positive genomes with no matrix row: {unlisted:,}')

# ---------------------------------------------------------------- phaC proteins
print('reading phaC protein annotations')
triad, cls, cls_ok, prot_len = {}, {}, {}, {}
with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        triad[row['target_id']] = '1' if row['triad_complete'] == 'True' else '0'
with open(ROOT / 'catalytic_domain/phac_synthase_class.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        cls[row['target_id']] = row['best_model_class']
        cls_ok[row['target_id']] = '1' if row['confident'] == 'True' else '0'
        prot_len[row['target_id']] = row['length_aa']

qtm = {}
with open(ROOT / 'figures/seq_vs_structure_hero.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        qtm[row['target_id']] = num(row['qtmscore'], 3)
print(f'  qTM available for {len(qtm):,} representatives')

print('reading TemStaPro predictions')
thermo = {}
with open(ROOT / 'temstapro/phaC_temstapro_predictions_with_metadata.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        best = 0
        for t in TEMP_THRESHOLDS:
            if (row.get(f't{t}_binary') or '').strip() in ('1', 'True', 'true'):
                best = t
        thermo[row['target_id']] = str(best)
print(f'  {len(thermo):,} proteins scored')

cluster_of = {}
with open(RESULTS / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for cluster_id, target_id in csv.reader(f, delimiter='\t'):
        cluster_of[target_id] = cluster_id

print('aggregating per-protein taxonomy and best hit')
prot_phylum, prot_genus, best = defaultdict(Counter), defaultdict(Counter), {}
with open(RESULTS / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        tid = row['target_id']
        if tid in bad:
            continue
        prot_phylum[tid][row['gtdb_phylum']] += 1
        prot_genus[tid][row['gtdb_genus']] += 1
        if tid not in best:
            best[tid] = (num(row['best_pident'], 1), num(row['best_tcov'], 3), row['best_query'])

PROTEIN_COLS = ['id', 'len', 'cl', 'cls', 'clsok', 'tr', 'pid', 'tcov', 'qtm', 'th', 'ng', 'ph', 'gn', 'ref']
n_proteins = 0
with open(OUT / 'proteins.csv', 'w', newline='') as out:
    w = csv.writer(out, lineterminator='\n')
    w.writerow(PROTEIN_COLS)
    for tid in sorted(protein_genomes):
        pid, tcov, ref = best.get(tid, ('', '', ''))
        dom_ph = prot_phylum[tid].most_common(1)
        dom_gn = prot_genus[tid].most_common(1)
        w.writerow([short(tid), prot_len.get(tid, ''), short(cluster_of.get(tid, '')),
                    cls.get(tid, ''), cls_ok.get(tid, ''), triad.get(tid, ''),
                    pid, tcov, qtm.get(tid, ''), thermo.get(tid, ''),
                    len(protein_genomes[tid]),
                    phylum(dom_ph[0][0]) if dom_ph else '',
                    genus(dom_gn[0][0]) if dom_gn else '',
                    reference(ref)])
        n_proteins += 1
print(f'  {n_proteins:,} phaC proteins')

# ---------------------------------------------------------------- clusters
print('writing occurrences.csv')
genome_index = {name: i for i, name in enumerate(names)}
OCCURRENCE_COLS = ['p', 'c', 'g']
n_occurrences = n_unmapped = 0
with open(OUT / 'occurrences.csv', 'w', newline='') as out:
    w = csv.writer(out, lineterminator='\n')
    w.writerow(OCCURRENCE_COLS)
    for tid in sorted(protein_genomes):
        cid = short(cluster_of.get(tid, ''))
        pid = short(tid)
        for genome in sorted(protein_genomes[tid]):
            row = genome_index.get(genome)
            if row is None:
                n_unmapped += 1   # phaC carrier with no metadata row; cannot be mapped
                continue
            w.writerow([pid, cid, row])
            n_occurrences += 1
print(f'  {n_occurrences:,} mappable occurrences, {n_unmapped:,} without a metadata row')

print('writing clusters.csv')
cluster_genomes, cluster_proteins, cluster_phyla, cluster_classes = (
    defaultdict(set), Counter(), defaultdict(Counter), defaultdict(set))
for tid, genomes in protein_genomes.items():
    cid = cluster_of.get(tid)
    if cid is None:
        continue
    cluster_genomes[cid] |= genomes
    cluster_proteins[cid] += 1
    if cls.get(tid):
        cluster_classes[cid].add(cls[tid])
    for p, n in prot_phylum[tid].items():
        cluster_phyla[cid][p] += n

landscape = {}
with open(ROOT / 'figures/phac_cluster_landscape.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        landscape[row['cluster_id']] = row

CLUSTER_COLS = ['id', 'ng', 'np', 'nloc', 'ncls', 'nph', 'ph', 'gen']
with open(OUT / 'clusters.csv', 'w', newline='') as out:
    w = csv.writer(out, lineterminator='\n')
    w.writerow(CLUSTER_COLS)
    for cid, genomes in sorted(cluster_genomes.items(), key=lambda kv: -len(kv[1])):
        meta = landscape.get(cid, {})
        dom = cluster_phyla[cid].most_common(1)
        w.writerow([short(cid), len(genomes), cluster_proteins[cid],
                    meta.get('n_locations', ''), len(cluster_classes[cid]),
                    len(cluster_phyla[cid]),
                    phylum(dom[0][0]) if dom else '',
                    genus(meta.get('top_genus', ''))])
print(f'  {len(cluster_genomes):,} clusters')

# ---------------------------------------------------------------- lookups
(OUT / 'lookups.json').write_text(json.dumps({
    'phylum': phylum.labels, 'genus': genus.labels, 'habitat': habitat.labels,
    'ecosystem': ecosystem.labels, 'study': study.labels, 'depth_zone': zone.labels,
    'architecture': arch.labels, 'reference': reference.labels,
    'families': FAMILIES, 'thresholds': TEMP_THRESHOLDS,
    'genome_columns': GENOME_COLS, 'protein_columns': PROTEIN_COLS,
    'cluster_columns': CLUSTER_COLS, 'occurrence_columns': OCCURRENCE_COLS,
}, separators=(',', ':')))

(OUT / 'meta.json').write_text(json.dumps({
    'built': date.today().isoformat(),
    'target_prefix': TARGET_PREFIX,
    'counts': {
        'genomes_in_matrix': n_genomes,
        'genomes_with_coordinates': with_coords,
        'phac_positive_genomes': len(phac_genomes),
        'phac_positive_genomes_without_metadata': unlisted,
        'phac_proteins': n_proteins,
        'phac_clusters': len(cluster_genomes),
        'mappable_occurrences': n_occurrences,
        'unmappable_occurrences': n_unmapped,
        'proteins_with_qtm': len(qtm),
        'family_genome_counts': dict(family_genomes),
    },
    'caveats': [
        'phaC presence comes from the uncapped protein-to-genome membership list '
        f'({len(phac_genomes):,} genomes). Every other family comes from a metadata '
        'join capped at 5 genomes per identical protein, so those counts are lower '
        'bounds, not exact.',
        f'{unlisted:,} phaC-positive genomes have no row in the metadata matrix, so '
        'they carry no taxonomy, habitat or coordinates and cannot appear on the map.',
        'qTM structural scores exist only for the 95%-dereplication representatives, '
        f'{len(qtm):,} of {n_proteins:,} proteins. Blank means not measured, not low.',
        'Catalytic-triad completeness is scored by HMM alignment columns, which '
        'produces false negatives on divergent sequences; treat it as a floor.',
        'ESMFold models and TemStaPro thermostability are predictions, not '
        'experimental measurements.',
    ],
}, indent=2))

print('\nwrote:')
for p in sorted(OUT.iterdir()):
    print(f'  {p.stat().st_size / 1e6:8.2f} MB  {p.name}')
