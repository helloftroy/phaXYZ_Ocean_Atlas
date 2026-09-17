#!/bin/bash
# Step 2 of "cluster PhaC (or any family) more deeply than the best_query
# proxy": clusters the sequences extracted by
# cluster/run_omdb_extract_cluster_sequences.sbatch at a real identity
# threshold (50-70% is a reasonable starting range for resolving
# subfamilies within a protein family -- defaults to 60%), using
# `mmseqs easy-cluster` -- the same tool/style OMDB's own NR50/NR30
# catalogs were built with (see PHA_bioprospecting/README.md's OMDB
# downloads section for their own documented invocation). CPU-only, no
# GPU needed -- this clusters only OUR discovered target_ids (typically
# thousands to tens of thousands), not the full OMDB catalog, so it is
# far smaller/faster than anything else in this pipeline.
#
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_omdb_cluster_sequences.sbatch
# Override the identity threshold or coverage if 60%/0.8 is not what you want:
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC,MIN_SEQ_ID=0.5 cluster/run_omdb_cluster_sequences.sbatch
#SBATCH --job-name=pha_omdb_cluster_seqs
#SBATCH --partition=service
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/omdb_cluster_seqs_%j.out
#SBATCH --error=logs/omdb_cluster_seqs_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs
export PYTHONUNBUFFERED=1

if [ ! -f cluster/env_activate.sh ]; then
  echo "cluster/env_activate.sh is missing -- run ./cluster/setup_env.sh first." >&2
  exit 2
fi
# shellcheck source=/dev/null
source cluster/env_activate.sh

echo "Job started: $(date)"
echo "Host: $(hostname)"

: "${FAMILY:?Set FAMILY to a single family_id.}"
MMSEQS_BIN="${MMSEQS_BIN:-$(pwd)/cluster/bin/mmseqs}"
RESULTS_DIR="${RESULTS_DIR:-PHA_bioprospecting/omdb_search/results}"
MIN_SEQ_ID="${MIN_SEQ_ID:-0.6}"
COVERAGE="${COVERAGE:-0.8}"
echo "FAMILY=${FAMILY} MIN_SEQ_ID=${MIN_SEQ_ID} COVERAGE=${COVERAGE}"

FASTA_IN="${RESULTS_DIR}/${FAMILY}_cluster_sequences.faa"
if [ ! -f "${FASTA_IN}" ]; then
  echo "${FASTA_IN} not found -- run cluster/run_omdb_extract_cluster_sequences.sbatch first." >&2
  exit 2
fi
if [ ! -x "${MMSEQS_BIN}" ]; then
  echo "${MMSEQS_BIN} not found/executable -- run ./cluster/install_mmseqs2.sh first." >&2
  exit 2
fi

OUT_PREFIX="${RESULTS_DIR}/${FAMILY}_cluster${MIN_SEQ_ID}"
TMP_DIR="PHA_bioprospecting/omdb_search/tmp/${FAMILY}_cluster_seq_tmp"
mkdir -p "${TMP_DIR}"

"${MMSEQS_BIN}" easy-cluster "${FASTA_IN}" "${OUT_PREFIX}" "${TMP_DIR}" \
  --min-seq-id "${MIN_SEQ_ID}" -c "${COVERAGE}" --cov-mode 1 --cluster-mode 2 \
  --threads "${SLURM_CPUS_PER_TASK:-8}"

echo "Job finished: $(date)"
echo "Cluster assignments: ${OUT_PREFIX}_cluster.tsv"
echo "Next: pha-reference cluster-ecology --family ${FAMILY} --cluster-tsv ${OUT_PREFIX}_cluster.tsv"
echo "  and: pha-reference sequence-embedding --family ${FAMILY} --fasta ${FASTA_IN} --cluster-tsv ${OUT_PREFIX}_cluster.tsv"
