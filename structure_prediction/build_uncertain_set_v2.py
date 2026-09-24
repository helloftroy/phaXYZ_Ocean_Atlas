"""Rebuilds structure_prediction/uncertain_cluster_representatives.faa with
a hand-curated set, replacing the original 4,992-sequence version from
build_fold_inputs.py. That version was built against the pre-fix tier
definitions (see PHA_CLEAN_RESULTS.md section 6, marked STALE) -- the
"no direct evidence" tiers it sampled from collapsed to 27 genomes total
after the 2026-09-22 phaC reference-query fix, so most of its 4,992
sequences no longer correspond to anything uncertain about this dataset.

This version is picked directly by the user from five sources, each
targeting a specific open question raised during the section 9.1-9.5
deep dives rather than a broad, generic sample:

  1. Modicisalibacter zincidurans -- all 10 distinct phaC target_ids
     across all 6 sister genomes (section 9.5). Includes the paralog that
     is confidently phaC by reference identity but scores weakly on every
     TIGRFAM class HMM (738aa, cluster ...753017) -- exactly the kind of
     case a folded structure can help resolve where sequence HMMs did not.
  2. Desulfoluna (CARD22-1 study) -- all 8 distinct target_ids across both
     CARD22-1 Desulfoluna genomes (section 9.4). Includes the two
     phaC/phaZ-ambiguous copies (708aa/736aa) found there.
  3. "Unknown HK1" -- all 330 distinct target_ids across all 167 phaC-
     positive genomes of this species (considered as a genome-comparison
     candidate in section 9.5 but not chosen; still worth structural
     resolution given its size and sponge-symbiont prevalence).
  4. The 26 distinct target_ids behind section 3's two smallest evidence
     tiers ("No HMM/triad, >=5 other PHA genes" + "1-4 other PHA genes",
     27 genomes total, one shares a cluster) -- genomes with pathway
     context (other verified PHA genes present) but literally zero direct
     sequence evidence for phaC itself. This is the one tier from the
     original build_fold_inputs.py plan that is still exactly what it was:
     a small enough set now to fold every member, not sample from it.
  5. All 273 targets dataset-wide with <50% identity to their nearest
     non-self phaC in the all-vs-all self-search (264 with a resolved hit
     below 50%, plus 9 with literally no self-hit at all even at the loose
     e-value cutoff used -- the most extreme form of "divergent"). Same
     self-search table and re-filtering discipline as
     build_selfsearch_derived_tables.py's divergent35 table, just at 50%
     instead of 35% and applied per-target rather than one representative
     per genome, since the goal here is folding specific proteins, not
     reporting genome-level counts.

Length filter: the project's usual 150-700aa fragment filter (see
plot_phac_divergent35_vs_depth.py) is applied ONLY to sources 4 and 5,
where it serves its original purpose (excluding likely-unfoldable
fragments from an otherwise-generic scan). It is deliberately NOT applied
to sources 1-3 -- those are hand-picked, already-investigated candidates,
and three of them (the CARD22-1 708aa/736aa pair and the Modicisalibater
738aa paralog) are >700aa specifically because their unusual length is
part of why they are worth a fold in the first place; dropping them here
would silently remove the most interesting cases this rebuild exists for.
Below MIN_LEN (150aa) is excluded everywhere, organism-specific sets
included, since a fold from <150 residues is not structurally
interpretable regardless of how the candidate was chosen.

Outputs:
    structure_prediction/uncertain_cluster_representatives.faa
    structure_prediction/fold_manifest.tsv (uncertain rows replaced;
      positive_control/reference rows carried over unchanged)
    structure_prediction/uncertain_set_v2_provenance.tsv (which source(s)
      contributed each target_id, for traceability)

Usage:
    python structure_prediction/build_uncertain_set_v2.py
"""
import csv
import pickle
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
OUT = Path(__file__).resolve().parent

sys.path.insert(0, str(ROOT / 'figures/scripts'))
import _phac_qc

MIN_LEN, MAX_LEN = 150, 700
DIVERGENT_THRESHOLD = 50.0

bad_targets = _phac_qc.load_bad_targets()

# ---------------------------------------------------------------------
# 1-3: organism-specific sets, corrected NR100 target_id<->genome join
# ---------------------------------------------------------------------
genomes = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        genomes[row['genome']] = row

modz_genomes = {g for g, r in genomes.items() if r['gtdb_species'] == 'Modicisalibacter zincidurans'}
card_genomes = {'CARD22-1_SAMN24292811_MAG_00000010', 'CARD22-1_SAMN24292812_MAG_00000014'}
hk1_genomes = {g for g, r in genomes.items() if r['gtdb_species'] == 'Unknown HK1' and int(r['n_phaC']) > 0}

organism_genomes = modz_genomes | card_genomes | hk1_genomes
genome_targets = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['genome'] in organism_genomes and row['target_id'] not in bad_targets:
            genome_targets[row['genome']].add(row['target_id'])

modz_targets, card_targets, hk1_targets = set(), set(), set()
for g in modz_genomes:
    modz_targets |= genome_targets[g]
for g in card_genomes:
    card_targets |= genome_targets[g]
for g in hk1_genomes:
    hk1_targets |= genome_targets[g]

print(f'Modicisalibacter zincidurans: {len(modz_genomes)} genomes, {len(modz_targets)} distinct targets')
print(f'CARD22-1 Desulfoluna: {len(card_genomes)} genomes, {len(card_targets)} distinct targets')
print(f'Unknown HK1: {len(hk1_genomes)} genomes, {len(hk1_targets)} distinct targets')

# ---------------------------------------------------------------------
# 4: the two smallest section-3 evidence tiers (no HMM, no triad)
# ---------------------------------------------------------------------
genome_group = pickle.load(open('/tmp/phac_verified_triad_hmm_group.pkl', 'rb'))
genome_targets_raw = pickle.load(open('/tmp/verified_phac_genome_targets.pkl', 'rb'))
NO_SUPPORT_TIERS = {'No HMM/triad, >=5 other PHA genes', 'No HMM/triad, 1-4 other PHA genes'}
no_support_genomes = [g for g, t in genome_group.items() if t in NO_SUPPORT_TIERS]
no_support_targets = set()
for g in no_support_genomes:
    no_support_targets |= (genome_targets_raw[g] - bad_targets)
print(f'No HMM/triad support tiers: {len(no_support_genomes)} genomes, {len(no_support_targets)} distinct targets')

# ---------------------------------------------------------------------
# 5: dataset-wide, <50% identity to nearest non-self phaC (or no self-hit)
# ---------------------------------------------------------------------
all_candidates = set()
for g, ts in genome_targets_raw.items():
    all_candidates |= (ts - bad_targets)

best_hit = {}
with open(ROOT / 'catalytic_domain/selfsearch/hits.tsv') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        q, t = row['query'], row['target']
        if q == t or q not in all_candidates or t in bad_targets:
            continue
        bits = float(row['bits'])
        if q not in best_hit or bits > best_hit[q][2]:
            best_hit[q] = (float(row['pident']), t, bits)

divergent_targets = {t for t, (pident, hit_t, bits) in best_hit.items() if pident < DIVERGENT_THRESHOLD}
no_selfhit_targets = all_candidates - set(best_hit)
divergent_targets |= no_selfhit_targets
print(f'Divergent (<{DIVERGENT_THRESHOLD:.0f}% to nearest non-self, or no self-hit): {len(divergent_targets)} targets '
      f'({len(divergent_targets) - len(no_selfhit_targets)} resolved + {len(no_selfhit_targets)} no-hit)')

# ---------------------------------------------------------------------
# assemble, length-filter sources 4-5 only, write provenance + FASTA
# ---------------------------------------------------------------------
target_lengths = pickle.load(open('/tmp/all_phac_target_lengths.pkl', 'rb'))

sources = {
    'modicisalibacter': modz_targets, 'card22_desulfoluna': card_targets, 'hk1': hk1_targets,
    'no_hmm_triad_support': {t for t in no_support_targets if target_lengths.get(t) and MIN_LEN <= target_lengths[t] <= MAX_LEN},
    'divergent_lt50pct': {t for t in divergent_targets if target_lengths.get(t) and MIN_LEN <= target_lengths[t] <= MAX_LEN},
}
provenance = defaultdict(list)
for name, tset in sources.items():
    for t in tset:
        provenance[t].append(name)

# still apply the MIN_LEN floor everywhere (unfoldable fragments), just
# not the MAX_LEN ceiling for the organism-specific sets 1-3
final = set()
for t, srcs in provenance.items():
    L = target_lengths.get(t)
    if L is None or L < MIN_LEN:
        continue
    final.add(t)

print(f'\n{len(final)} distinct target_ids in the final union (of {len(provenance)} candidates before the MIN_LEN floor)')

prov_out = OUT / 'uncertain_set_v2_provenance.tsv'
with open(prov_out, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['target_id', 'protein_length_aa', 'sources'])
    for t in sorted(final):
        w.writerow([t, target_lengths.get(t, ''), ','.join(sorted(provenance[t]))])
print(f'wrote {prov_out}')

# ---- pull sequences from the phaC-family sequence catalog ----
sequences = {}
cur_id, cur_seq = None, []
with open(FA / 'phaC_cluster_sequences.faa') as f:
    for line in f:
        if line.startswith('>'):
            if cur_id in final:
                sequences[cur_id] = ''.join(cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_id in final:
        sequences[cur_id] = ''.join(cur_seq)

missing = final - set(sequences)
if missing:
    print(f'WARNING: {len(missing)} target_ids not found in phaC_cluster_sequences.faa: {sorted(missing)[:10]}...')

fasta_out = OUT / 'uncertain_cluster_representatives.faa'
n_written = 0
with open(fasta_out, 'w') as out:
    for t in sorted(final):
        seq = sequences.get(t)
        if not seq:
            continue
        out.write(f'>{t}\n')
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + '\n')
        n_written += 1
print(f'wrote {fasta_out} ({n_written} sequences)')

# ---- update fold_manifest.tsv: replace uncertain rows, keep the rest ----
manifest_path = OUT / 'fold_manifest.tsv'
kept_rows = []
if manifest_path.exists():
    with open(manifest_path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['set'] != 'uncertain':
                kept_rows.append(row)
with open(manifest_path, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['target_id', 'set'])
    for t in sorted(final):
        if t in sequences:
            w.writerow([t, 'uncertain'])
    for row in kept_rows:
        w.writerow([row['target_id'], row['set']])
print(f'wrote {manifest_path} ({n_written} uncertain + {len(kept_rows)} carried-over rows)')

# length distribution sanity check
lens = sorted(target_lengths.get(t, 0) for t in final if t in sequences)
if lens:
    print(f'\nlength range: {lens[0]}-{lens[-1]}aa, median {lens[len(lens)//2]}aa, '
          f'{sum(1 for l in lens if l > MAX_LEN)} sequences >{MAX_LEN}aa (kept deliberately, see docstring)')
