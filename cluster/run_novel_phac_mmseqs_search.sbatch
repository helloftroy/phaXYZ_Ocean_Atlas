#!/bin/bash
# Novel-phaC recovery, mmseqs method: searches the 3,478 triad-complete,
# one-per-cluster0.7-representative phaC query sequences (query side,
# built locally: phac_recovery/build_novel_phac_query_set.py) against the
# combined whole-proteome target set for all 7,826 qualifying_genomes.txt
# genomes (target side, cluster-built:
# cluster/run_phac_recovery_download_proteomes.sbatch) -- these genomes
# carry >5 other PHA-pathway genes but no phaC hit from the original
# search, so a real hit here is a candidate that structure/reference
# identity established phaC has diverged too far in sequence for the
# original NR100 seed search to catch, but mmseqs' own local-alignment
# search still can.
#
# CPU, not GPU: this is a small search (3.5k queries x a proteome target
# set several orders of magnitude smaller than the full 249.5M-sequence
# OMDB database run_omdb_search.sbatch handles), so plain `mmseqs
# easy-search` on `morrill` (the confirmed real CPU-compute partition for
# this account -- see cluster/run_foldseek.sbatch's own note; no
# internet needed for this step, unlike the download step before it) is
# the right tool, no GPU target-db build required.
#
# Deliberately permissive SENSITIVITY/EVALUE defaults (higher sensitivity,
# looser e-value than cluster/run_omdb_search.sbatch's SENSITIVITY=7.0/
# EVALUE=1e-5) -- the whole point of this search is catching phaC too
# divergent for the original seed search's default strictness, mirroring
# the same "relaxed --domE 10" philosophy the parallel HMM-based recovery
# attempt used (cluster/run_phac_recovery_hmmsearch.sbatch). Tighten later
# once the raw candidate count is in hand, don't start tight and miss the
# point of the exercise.
#
# One-time prerequisite: cluster/bin/mmseqs (CPU build) -- see
# cluster/install_mmseqs2.sh, already needed by pha-reference cluster95/
# gopc-build-db, so likely already installed.
#
#   sbatch --account=191001-364393 cluster/run_novel_phac_mmseqs_search.sbatch
#SBATCH --job-name=novel_phac_mmseqs_search
#SBATCH --partition=morrill
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/novel_phac_mmseqs_search_%j.out
#SBATCH --error=logs/novel_phac_mmseqs_search_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs
export PYTHONUNBUFFERED=1

MMSEQS_BIN="${MMSEQS_BIN:-$(pwd)/cluster/bin/mmseqs}"
if [ ! -x "${MMSEQS_BIN}" ]; then
  echo "${MMSEQS_BIN} not found/executable -- install it first: ./cluster/install_mmseqs2.sh" >&2
  exit 2
fi
echo "mmseqs (CPU build): $("${MMSEQS_BIN}" version)"
echo "Job started: $(date)"
echo "Host: $(hostname)"

QUERY_FAA="phac_recovery/novel_phac_query_representatives.faa"
TARGET_FAA="phac_recovery/novel_phac_target_proteomes.faa"
OUT_TSV="phac_recovery/novel_phac_mmseqs_hits.tsv"
TMP_DIR="$(mktemp -d -p "${TMPDIR:-/tmp}" novel_phac_mmseqs.XXXXXX)"
trap 'rm -rf "${TMP_DIR}"' EXIT

for required in "${QUERY_FAA}" "${TARGET_FAA}"; do
  if [ ! -f "${required}" ]; then
    echo "${required} not found." >&2
    echo "  Query side (run locally): .venv/bin/python phac_recovery/build_novel_phac_query_set.py" >&2
    echo "  Target side (cluster): sbatch --account=191001-364393 cluster/run_phac_recovery_download_proteomes.sbatch" >&2
    exit 2
  fi
done

SENSITIVITY="${SENSITIVITY:-7.5}"
EVALUE="${EVALUE:-1e-3}"
COVERAGE="${COVERAGE:-0.3}"
COV_MODE="${COV_MODE:-0}"
MAX_SEQS="${MAX_SEQS:-1000}"
FORMAT_OUTPUT="query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen"
echo "SENSITIVITY=${SENSITIVITY} EVALUE=${EVALUE} COVERAGE=${COVERAGE} COV_MODE=${COV_MODE} MAX_SEQS=${MAX_SEQS}"

"${MMSEQS_BIN}" easy-search \
  "${QUERY_FAA}" "${TARGET_FAA}" "${OUT_TSV}" "${TMP_DIR}" \
  -s "${SENSITIVITY}" \
  -e "${EVALUE}" \
  -c "${COVERAGE}" \
  --cov-mode "${COV_MODE}" \
  --max-seqs "${MAX_SEQS}" \
  --format-output "${FORMAT_OUTPUT}" \
  --threads "${SLURM_CPUS_PER_TASK:-16}"

echo ""
echo "=== Candidate summary ==="
N_HITS=$(wc -l < "${OUT_TSV}")
echo "${N_HITS} raw hit rows -> ${OUT_TSV}"
if [ "${N_HITS}" -gt 0 ]; then
  N_GENES=$(cut -f2 "${OUT_TSV}" | sort -u | wc -l)
  N_GENOMES=$(cut -f2 "${OUT_TSV}" | sed -E 's/-scaffold_.*$//' | sort -u | wc -l)
  N_QUERIES_HIT=$(cut -f1 "${OUT_TSV}" | sort -u | wc -l)
  echo "${N_GENES} distinct candidate target genes"
  echo "${N_GENOMES} distinct genomes (of the 7,826 qualifying) with >=1 candidate"
  echo "${N_QUERIES_HIT} of the 3,478 query representatives found >=1 hit"
fi

echo ""
echo "Job finished: $(date)"
echo "Next (not built yet): triad-completeness + structural check of the surfaced candidates."
