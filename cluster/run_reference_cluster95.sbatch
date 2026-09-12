#!/bin/bash
# Exact-dedup + per-family MMseqs2 clustering (95% identity / 90%
# coverage) of PHA_reference/pha_reference.sqlite, then regenerates every
# export (master CSV + all/nr95 FASTA + cluster95 TSV). No internet needed
# -- this only reads the already-built local database. Needs
# cluster/bin/mmseqs -- run ./cluster/install_mmseqs2.sh first if you
# haven't already (cluster/env_activate.sh puts cluster/bin on PATH).
#
# The reference database is small (thousands of sequences per family, not
# millions), so this is a light CPU job -- these resource numbers are
# deliberately modest; raise --cpus-per-task if you profile it as
# CPU-bound on your cluster's node type.
#
# Add --account=<your account> on submission if your cluster requires one
# for this partition/QOS, e.g.:
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_reference_cluster95.sbatch
#SBATCH --job-name=pha_reference_cluster95
#SBATCH --partition=service
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/reference_cluster95_%j.out
#SBATCH --error=logs/reference_cluster95_%j.err

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

if ! command -v mmseqs >/dev/null 2>&1; then
  echo "mmseqs not found on PATH -- run ./cluster/install_mmseqs2.sh first." >&2
  exit 2
fi
echo "mmseqs: $(mmseqs version)"

echo "Job started: $(date)"
echo "Host: $(hostname)"

FAMILY="${FAMILY:-all}"
MIN_SEQ_ID="${MIN_SEQ_ID:-0.95}"
MIN_COV="${MIN_COV:-0.90}"
echo "FAMILY=${FAMILY} MIN_SEQ_ID=${MIN_SEQ_ID} MIN_COV=${MIN_COV}"

pha-reference init-db
pha-reference cluster95 --family "${FAMILY}" --min-seq-id "${MIN_SEQ_ID}" --min-cov "${MIN_COV}"
pha-reference export
pha-reference export-cluster95
pha-reference status

echo "Job finished: $(date)"
echo "Results: PHA_reference/exports/ (pha_reference_master.csv now carries cluster95_*;"
echo "         pha_reference_all.faa, pha_reference_nr95.faa, pha_reference_cluster95.tsv)"
echo "Also reachable via the PHA_bioprospecting/reference symlink."
