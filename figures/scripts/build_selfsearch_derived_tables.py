"""Rebuilds the tables behind the divergent35 and hadal-trench figures
from the all-vs-all self-search among phaC targets (see
catalytic_domain/selfsearch/hits.tsv, built via mmseqs2 -s 7.0 -e 1e-5
-c 0.0 --max-seqs 20 against catalytic_domain/verified_phac_targets.faa,
originally the 82,846 targets behind the then-68,424 verified phaC
genomes -- see PHA_CLEAN_RESULTS.md section 1.4/2 for how "verified" is
defined, since that definition has since been corrected twice).

The originals (/tmp/all_phac_best_hit.pkl, /tmp/all_depth_phac.pkl,
/tmp/all_phac_genome_target.pkl) were built ad hoc earlier in this
project's history with no reproducing script committed anywhere -- this
replaces them with a documented, rerunnable equivalent.

IMPORTANT: re-filters BOTH the genome->target mapping AND every self-
search hit row against the CURRENT bad-query list (_phac_qc), rather than
trusting /tmp/verified_phac_genome_targets.pkl or the raw hits.tsv to
already reflect it -- this is the same staleness trap fixed in
plot_phac_triad_hmm_pathway_groups.py, now fixed here too. A self-search
hit is only trusted if NEITHER the query nor its matched target is a
now-excluded target_id: the self-search compares candidate phaC sequences
against EACH OTHER, not against references, so "bad" here means "this
target_id's own phaC status is no longer trusted," which invalidates it
on either side of a hit, not just as a query.

Outputs (same filenames/shapes the existing plot scripts expect):
  - figures/phac_divergent35_candidates.tsv
  - figures/deep_phac_hadal_trench_profile.tsv
  - figures/deep_phac_nearest_known_identity.tsv
  - /tmp/all_phac_best_hit.pkl, /tmp/all_depth_phac.pkl, /tmp/all_phac_genome_target.pkl
    (kept under the same names so plot_phac_divergent35_vs_depth.py needs no changes)
"""
import csv
import math
import pickle
import sys
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
OUT = ROOT / 'figures'

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc
bad_targets = _phac_qc.load_bad_targets()

DIVERGENT_THRESHOLD = 35.0
MIN_LEN, MAX_LEN = 150, 700
MIN_COV = 0.6
HADAL_THRESHOLD = 6000.0

# known hadal-trench sample coordinates in this dataset (all from a single
# study, SAWJ20-1 -- 4 distinct lat/lon clusters, matched here by rounding
# to 1 decimal degree; see conversation notes for the coordinate check)
TRENCH_COORDS = {
    (11.4, 142.4): 'Mariana Trench',
    (41.9, 146.3): 'Kuril-Kamchatka Trench',
    (29.2, 142.8): 'Izu-Ogasawara (Bonin) Trench',
    (36.1, 142.7): 'Japan Trench',
}

def trench_name(lat, lon):
    key = (round(lat, 1), round(lon, 1))
    if key in TRENCH_COORDS:
        return TRENCH_COORDS[key]
    # fallback: nearest known trench centroid within ~2 degrees, else generic
    best, best_d = None, 999
    for (clat, clon), name in TRENCH_COORDS.items():
        d = math.hypot(lat - clat, lon - clon)
        if d < best_d:
            best, best_d = name, d
    return best if best_d < 2.0 else 'trench'


print('loading self-search hits...')
best_hit = {}  # target_id -> (pident, hit_target, bits, qcov, tcov)
n_dropped_bad = 0
with open(ROOT / 'catalytic_domain/selfsearch/hits.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        q, t = row['query'], row['target']
        if q == t:
            continue  # self-hit, skip
        if q in bad_targets or t in bad_targets:
            n_dropped_bad += 1
            continue
        bits = float(row['bits'])
        if q not in best_hit or bits > best_hit[q][2]:
            best_hit[q] = (float(row['pident']), t, bits, float(row['qcov']), float(row['tcov']))

print(f'{len(best_hit)} targets with a resolved best non-self hit ({n_dropped_bad} hit rows dropped, bad query or target)')
pickle.dump(best_hit, open('/tmp/all_phac_best_hit.pkl', 'wb'))

genome_targets_raw = pickle.load(open('/tmp/verified_phac_genome_targets.pkl', 'rb'))
genome_targets = {}
for g, ts in genome_targets_raw.items():
    kept = ts - bad_targets
    if kept:
        genome_targets[g] = kept
print(f'{len(genome_targets_raw):,} genomes before re-filtering against the current '
      f'{len(_phac_qc.BAD_QUERIES)}-accession exclusion list -> {len(genome_targets):,} after')

# one representative target per genome (longest sequence, deterministic)
target_lengths = pickle.load(open('/tmp/all_phac_target_lengths.pkl', 'rb'))
all_genomes_targets = {}
for g, ts in genome_targets.items():
    all_genomes_targets[g] = max(ts, key=lambda t: target_lengths.get(t, 0))
pickle.dump(all_genomes_targets, open('/tmp/all_phac_genome_target.pkl', 'wb'))
print(f'{len(all_genomes_targets)} genomes with a representative target_id')

# ---- genome metadata (depth, taxonomy, location) from the depth-enriched file ----
genome_meta = {}
with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        g = row['genome']
        if g not in genome_targets or g in genome_meta:
            continue
        genome_meta[g] = row

depth_rows = {}
for g, row in genome_meta.items():
    if row.get('depth_m'):
        try:
            depth_rows[g] = {'depth_m': float(row['depth_m'])}
        except ValueError:
            pass
pickle.dump(depth_rows, open('/tmp/all_depth_phac.pkl', 'wb'))
print(f'{len(depth_rows)} genomes with resolved depth')

# ---- best-audited-reference identity, for the "nearest known identity" table ----
seq_evidence = {}
with open(ROOT / 'catalytic_domain/phac_sequence_evidence.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        seq_evidence[row['target_id']] = row

# ============================================================
# 1. divergent35 candidates
# ============================================================
LOC_COLS = ['ecosystem_compartment', 'sample_source', 'ecosystem_type', 'study_id',
            'latitude_degN', 'longitude_degE']

candidates = []
qc_fail_n = 0
for genome, t in all_genomes_targets.items():
    if t not in best_hit:
        continue
    pident, hit_t, bits, qcov, tcov = best_hit[t]
    if pident >= DIVERGENT_THRESHOLD:
        continue
    qlen = target_lengths.get(t, 0)
    if not (MIN_LEN <= qlen <= MAX_LEN):
        continue
    if qcov < MIN_COV or tcov < MIN_COV:
        qc_fail_n += 1
        continue
    meta = genome_meta.get(genome, {})
    depth = depth_rows.get(genome, {}).get('depth_m', '')
    candidates.append({
        'genome': genome, 'target_id': t, 'pident_to_nearest': round(pident, 1),
        'nearest_hit_target_id': hit_t, 'qcov': round(qcov, 2), 'tcov': round(tcov, 2),
        'depth_m': depth, 'domain': meta.get('gtdb_domain', ''), 'phylum': meta.get('gtdb_phylum', ''),
        'genus': meta.get('gtdb_genus', ''), 'ecosystem_compartment': meta.get('ecosystem_compartment', ''),
        'sample_source': meta.get('sample_source', ''), 'ecosystem_type': meta.get('ecosystem_type', ''),
        'study_id': meta.get('study_id', ''), 'lat': meta.get('latitude_degN', ''), 'lon': meta.get('longitude_degE', ''),
    })

print(f'{len(candidates)} genome-level divergent35 candidates ({qc_fail_n} failed coverage QC)')
with open(OUT / 'phac_divergent35_candidates.tsv', 'w', newline='') as out:
    w = csv.DictWriter(out, fieldnames=list(candidates[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(candidates)
print('wrote figures/phac_divergent35_candidates.tsv')

# ============================================================
# 2. hadal trench profile + nearest-known-identity
# ============================================================
hadal_rows = []
nearest_rows = []
for genome, t in all_genomes_targets.items():
    d = depth_rows.get(genome, {}).get('depth_m')
    if d is None or d <= HADAL_THRESHOLD:
        continue
    meta = genome_meta.get(genome, {})
    lat = float(meta.get('latitude_degN') or 0)
    lon = float(meta.get('longitude_degE') or 0)
    trench = trench_name(lat, lon)
    hadal_rows.append({
        'depth_m': d, 'genome': genome, 'trench': trench,
        'domain': meta.get('gtdb_domain', ''), 'phylum': meta.get('gtdb_phylum', ''),
        'genus': meta.get('gtdb_genus', ''), 'compartment': meta.get('ecosystem_compartment', ''),
        'architecture': '', 'target_id': t, 'cluster_total_genomes': '',
        'phac_length_aa': target_lengths.get(t, ''), 'biosample': meta.get('biosample', ''),
    })

    hit = best_hit.get(t)
    ev = seq_evidence.get(t, {})
    nearest_rows.append({
        'depth_m': d, 'genome': genome, 'trench': trench, 'genus': meta.get('gtdb_genus', ''),
        'architecture': '',
        'max_id_nonhadal': hit[0] if hit else '', 'nonhadal_hit_target': hit[1] if hit else '',
        'nonhadal_qcov': hit[3] if hit else '', 'nonhadal_tcov': hit[4] if hit else '',
        'max_id_curated': ev.get('pident', ''), 'curated_hit_acc': ev.get('ref_accession', ''),
        'target_id': t,
    })

print(f'{len(hadal_rows)} verified genomes at hadal-trench depth (>{HADAL_THRESHOLD:.0f}m)')
with open(OUT / 'deep_phac_hadal_trench_profile.tsv', 'w', newline='') as out:
    w = csv.DictWriter(out, fieldnames=list(hadal_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(hadal_rows)
print('wrote figures/deep_phac_hadal_trench_profile.tsv')

with open(OUT / 'deep_phac_nearest_known_identity.tsv', 'w', newline='') as out:
    w = csv.DictWriter(out, fieldnames=list(nearest_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(nearest_rows)
print('wrote figures/deep_phac_nearest_known_identity.tsv')
