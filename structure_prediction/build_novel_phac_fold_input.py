"""Register the missed-synthase search's genuine finds as a fifth fold set.

PHA_CLEAN_RESULTS.md section 13.3: 255 triad-complete proteins that are
absent from the phaC target set entirely and sit in genomes with no
QC-passing phaC. Sequence evidence already puts them in the phaC family
(53-99% identity over at least 80% of both proteins to a triad-complete
validated synthase, with an intact Cys/Asp/His triad at the expected HMM
columns), so folding is not being asked to decide whether they are
synthases. It is being asked the two questions sequence cannot answer:

  1. Is the ~415 aa median length a real, compact architecture, or are
     these truncated gene calls at contig edges? A complete alpha/beta
     hydrolase core that closes on itself says architecture; a structure
     that stops mid-fold says truncation.
  2. Do the three triad residues actually converge in space? The HMM
     column projection says they are present in the alignment; only
     geometry says they form an active site.

Writes structure_prediction/novel_phac_candidates.faa (plain headers --
run_esmfold.py names each PDB after the first whitespace-delimited token,
and the tier/class annotations in the source FASTA's headers would
otherwise be dropped silently rather than deliberately) and appends the
ids to fold_manifest.tsv under set='novel_phac_candidates', which is what
run_foldseek_search.sh and build_structural_evidence_table.py both
discover query sets from.

Idempotent: re-running neither duplicates manifest rows nor changes the
FASTA.

Usage:
    python structure_prediction/build_novel_phac_fold_input.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'phac_recovery/novel_phac_new_candidates.faa'
FASTA = ROOT / 'novel_phac_candidates.faa'
MANIFEST = ROOT / 'fold_manifest.tsv'
SET_NAME = 'novel_phac_candidates'

if not SOURCE.exists():
    raise SystemExit(f'{SOURCE} is missing -- run phac_recovery/classify_novel_phac_candidates.py first.')

records, name, buf = [], None, []
with open(SOURCE) as f:
    for line in f:
        if line.startswith('>'):
            if name:
                records.append((name, ''.join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        records.append((name, ''.join(buf)))

with open(FASTA, 'w') as f:
    for gene_id, seq in records:
        f.write(f'>{gene_id}\n')
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + '\n')
lengths = sorted(len(s) for _, s in records)
median = lengths[len(lengths) // 2]
print(f'wrote {FASTA} -- {len(records)} sequences, median {median} aa, '
      f'range {lengths[0]}-{lengths[-1]} aa')

existing = set()
with open(MANIFEST, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    for row in reader:
        if row and row[1] == SET_NAME:
            existing.add(row[0])

new = [gene_id for gene_id, _ in records if gene_id not in existing]
# lineterminator='\n' is not optional here: csv.writer defaults to '\r\n', and the
# rest of fold_manifest.tsv uses '\n'. Mixed endings leave a trailing '\r' on the
# set name, which run_foldseek_search.sh's awk discovery turns into a query set
# called "novel_phac_candidates\r" -- pointing it at a directory that does not exist.
with open(MANIFEST, 'a', newline='') as f:
    writer = csv.writer(f, delimiter='\t', lineterminator='\n')
    for gene_id in new:
        writer.writerow([gene_id, SET_NAME])
print(f'{len(new)} ids appended to {MANIFEST.name} under set={SET_NAME} '
      f'({len(existing)} were already registered)')
print('\nNext: sbatch --account=191001-364393 cluster/run_novel_phac_structure_prediction.sbatch')
