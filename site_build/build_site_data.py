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

Where each field comes from, and why not all from one table:

  family counts, architecture   genome_family_matrix.tsv -- the only source, and
                                capped at 5 genomes per identical protein, so
                                every count except phaC is a lower bound.
  phaC presence                 the uncapped membership list, which is 38,744
                                genomes against the matrix's 31,464.
  taxonomy                      omdb_all_genome_taxonomy.tsv, which covers all
                                274,282 genomes rather than the matrix's subset.
  habitat, study, coordinates   omdb_all_genomes_with_locations.tsv, same full
                                coverage. Using this instead of the matrix adds
                                habitat for about 1,440 genomes and recovers the
                                1,763 phaC carriers the matrix has no row for.
  depth                         phaC_unique_targets_with_metadata_depth.tsv.
                                The matrix has depth_m and depth_zone columns
                                but they are empty in every one of its 241,864
                                rows, so reading depth from it silently produced
                                an atlas with no depth at all. Depth is a
                                property of the sample rather than the genome,
                                so a genome with no depth of its own inherits
                                it from another genome of the same sample.

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

# Placeholders that mean "not recorded". OMDB retires samples by writing
# DEPRECATED into the environment field rather than removing the row.
BLANK_VALUES = {'', 'NA', 'nan', 'None', 'Unknown', 'unknown',
                'Unknown_environment', 'DEPRECATED', 'not applicable', 'missing'}


def short(target_id):
    """Drop the catalogue prefix; the numeric tail is unique and much shorter."""
    return target_id[len(TARGET_PREFIX):].lstrip('0') or '0' \
        if target_id.startswith(TARGET_PREFIX) else target_id


def clean(value):
    value = (value or '').strip()
    return '' if value in BLANK_VALUES else value


class Codes:
    """Interns repeated strings as integers, preserving first-seen order."""

    def __init__(self):
        self.index, self.labels = {}, []

    def __call__(self, value):
        value = clean(value)
        if not value:
            return ''
        code = self.index.get(value)
        if code is None:
            code = self.index[value] = len(self.labels)
            self.labels.append(value)
        return str(code)


def num(value, places=None):
    """Format a numeric field, or '' when absent -- never a literal 'nan'."""
    value = clean(value)
    if not value:
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

# ---------------------------------------------------------------- metadata
print('reading taxonomy for all genomes')
taxonomy = {}
with open(RESULTS / 'omdb_all_genome_taxonomy.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        taxonomy[row['genome']] = (row.get('gtdb_phylum'), row.get('gtdb_genus'))
print(f'  {len(taxonomy):,} genomes')

print('reading sample provenance for all genomes')
provenance = {}
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        provenance[row['genome']] = row
print(f'  {len(provenance):,} genomes')

print('reading depth')
genome_depth, sample_depth = {}, {}
with open(RESULTS / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        metres = clean(row.get('depth_m'))
        if not metres:
            continue
        record = (metres, clean(row.get('depth_zone')))
        genome_depth.setdefault(row['genome'], record)
        sample = clean(row.get('sample_id'))
        if sample:
            sample_depth.setdefault(sample, record)
print(f'  {len(genome_depth):,} genomes with a depth of their own, '
      f'{len(sample_depth):,} samples with a known depth')


def depth_for(genome):
    """A genome's own depth, else the depth of another genome from its sample."""
    record = genome_depth.get(genome)
    if record:
        return record, 'own'
    sample = clean((provenance.get(genome) or {}).get('sample_id'))
    record = sample_depth.get(sample) if sample else None
    return (record, 'sample') if record else (('', ''), '')


# ---------------------------------------------------------------- genomes
print('writing genomes.csv')
GENOME_COLS = ['id', 'lat', 'lon', 'ph', 'gn', 'hb', 'eco', 'st', 'dep', 'dz',
               'fam', 'arch', 'phac', 'prof']
names, with_coords, with_habitat = [], 0, 0
depth_counts = Counter()
family_genomes = Counter()
rows_out = []

matrix_rows = {}
with open(RESULTS / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        matrix_rows[row['genome']] = row

# Genomes the matrix never profiled but that the uncapped list proves carry
# phaC. They are real carriers and belong in the atlas; only their counts for
# the other families are unknown, which `prof` records.
extra = sorted(phac_genomes - set(matrix_rows))
order = list(matrix_rows) + extra
print(f'  {len(matrix_rows):,} profiled genomes + {len(extra):,} phaC carriers '
      f'the matrix has no row for')

with open(OUT / 'genomes.csv', 'w', newline='') as out:
    w = csv.writer(out, lineterminator='\n')
    w.writerow(GENOME_COLS)
    for index, genome in enumerate(order):
        matrix = matrix_rows.get(genome)
        place = provenance.get(genome, {})
        tax_phylum, tax_genus = taxonomy.get(genome, (None, None))
        if tax_phylum is None and matrix:
            tax_phylum, tax_genus = matrix.get('gtdb_phylum'), matrix.get('gtdb_genus')

        counts = []
        for family in FAMILIES:
            n = 0
            if matrix:
                try:
                    n = int(float(matrix.get(f'n_{family}', '0') or '0'))
                except ValueError:
                    n = 0
            if family == 'phaC':
                n = max(n, 1) if genome in phac_genomes else 0
            if n:
                family_genomes[family] += 1
            counts.append(str(min(n, 9)))

        latitude = num(place.get('latitude_degN'), 3)
        longitude = num(place.get('longitude_degE'), 3)
        if not (latitude and longitude) and matrix:
            latitude = num(matrix.get('latitude_degN'), 3)
            longitude = num(matrix.get('longitude_degE'), 3)
        with_coords += bool(latitude and longitude)

        compartment = clean(place.get('ecosystem_compartment')) or (
            clean(matrix.get('ecosystem_compartment')) if matrix else '')
        with_habitat += bool(compartment)

        (metres, depth_zone), source = depth_for(genome)
        if metres:
            depth_counts[source] += 1

        w.writerow([
            index, latitude, longitude, phylum(tax_phylum), genus(tax_genus),
            habitat(compartment),
            ecosystem(clean(place.get('ecosystem_type')) or
                      (clean(matrix.get('ecosystem_type')) if matrix else '')),
            study(clean(place.get('study_id')) or
                  (clean(matrix.get('study_id')) if matrix else '')),
            num(metres, 1), zone(depth_zone), ''.join(counts),
            arch(matrix.get('architecture') if matrix else ''),
            '1' if genome in phac_genomes else '0',
            '1' if matrix else '0',
        ])
        names.append(genome)

(OUT / 'genome_names.csv').write_text('\n'.join(names) + '\n')
n_genomes = len(names)
print(f'  {n_genomes:,} genomes | {with_coords:,} with coordinates '
      f'({100 * with_coords / n_genomes:.1f}%) | {with_habitat:,} with a habitat '
      f'({100 * with_habitat / n_genomes:.1f}%)')
print(f'  depth: {depth_counts["own"]:,} from the genome, '
      f'{depth_counts["sample"]:,} inherited from its sample, '
      f'{sum(depth_counts.values()):,} total')

phac_habitat = sum(1 for g in phac_genomes if clean(
    (provenance.get(g) or {}).get('ecosystem_compartment')))
phac_coords = sum(1 for g in phac_genomes
                  if clean((provenance.get(g) or {}).get('latitude_degN'))
                  and clean((provenance.get(g) or {}).get('longitude_degE')))
phac_depth = sum(1 for g in phac_genomes if depth_for(g)[0][0])
print(f'  of the {len(phac_genomes):,} phaC-positive genomes: '
      f'{phac_habitat:,} have a habitat ({100 * phac_habitat / len(phac_genomes):.1f}%), '
      f'{phac_coords:,} coordinates ({100 * phac_coords / len(phac_genomes):.1f}%), '
      f'{phac_depth:,} depth ({100 * phac_depth / len(phac_genomes):.1f}%)')

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
best = {}
with open(RESULTS / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        tid = row['target_id']
        if tid in bad or tid in best:
            continue
        best[tid] = (num(row['best_pident'], 1), num(row['best_tcov'], 3), row['best_query'])

# Per-protein taxonomy from the genomes that carry it, using the full-coverage
# taxonomy table rather than the capped annotation rows.
prot_phylum, prot_genus = defaultdict(Counter), defaultdict(Counter)
for tid, genomes in protein_genomes.items():
    for genome in genomes:
        tax_phylum, tax_genus = taxonomy.get(genome, (None, None))
        if clean(tax_phylum):
            prot_phylum[tid][clean(tax_phylum)] += 1
        if clean(tax_genus):
            prot_genus[tid][clean(tax_genus)] += 1

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

# ---------------------------------------------------------------- occurrences
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
                n_unmapped += 1
                continue
            w.writerow([pid, cid, row])
            n_occurrences += 1
print(f'  {n_occurrences:,} occurrences, {n_unmapped:,} unplaceable')

# ---------------------------------------------------------------- clusters
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
        'genomes_total': n_genomes,
        'genomes_profiled': len(matrix_rows),
        'genomes_unprofiled': len(extra),
        'genomes_with_coordinates': with_coords,
        'genomes_with_habitat': with_habitat,
        'genomes_with_depth': sum(depth_counts.values()),
        'depth_from_own_record': depth_counts['own'],
        'depth_from_sample': depth_counts['sample'],
        'phac_positive_genomes': len(phac_genomes),
        'phac_with_habitat': phac_habitat,
        'phac_with_coordinates': phac_coords,
        'phac_with_depth': phac_depth,
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
        f'{len(extra):,} genomes are known phaC carriers but were never profiled for '
        'the other PHA families, so their counts for those families read as zero '
        'rather than unknown.',
        f'{len(phac_genomes) - phac_habitat:,} of the {len(phac_genomes):,} '
        'phaC-positive genomes have no habitat recorded, almost all of them samples '
        'OMDB has since marked deprecated; they also carry no coordinates.',
        'Depth is a property of the sample, not the genome. A genome without a depth '
        'of its own inherits the depth of another genome from the same sample, which '
        f'is where {depth_counts["sample"]:,} of the {sum(depth_counts.values()):,} '
        'depth values come from.',
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
