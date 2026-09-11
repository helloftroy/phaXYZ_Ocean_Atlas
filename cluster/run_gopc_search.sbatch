#!/bin/bash
# Sensitive MMseqs2 search of the NR95 PHA reference set against GOPC, per
# family: export-queries -> gopc-search -> gopc-combine. Run
# cluster/run_gopc_build_db.sbatch first (once). No internet needed.
#
# Test ONE family first, same philosophy as every other stage in this
# pipeline -- phaC is the cleanest case (a specific EC, not a shared
# superfamily), so it's the best first check that alignment settings and
# --max-seqs are behaving before spending compute on all 16:
#   sbatch --export=ALL,FAMILY=phaC cluster/run_gopc_search.sbatch
# Then everything:
#   sbatch cluster/run_gopc_search.sbatch   # FAMILY defaults to 'all'
#
# --max-seqs is the one setting most likely to need raising after a first
# look at the logs -- watch this job's own .out file (or
# results/<family>_summary.tsv's query_hit_cap_warning/queries_at_cap
# afterward) for the WARNING this script's own `pha-reference gopc-search`
# prints when a query's hit count approaches --max-seqs; that means the
# prefilter cap, not genuine biology, limited what came back for that
# query. Expected most for phaA/phaB/phaJ (shared enzyme superfamilies),
# less likely for phaC. Rerun just the flagged family with e.g.
# MAX_SEQS=50000 or 100000 if you see it.
#SBATCH --job-name=pha_gopc_search
#SBATCH --partition=service
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gopc_search_%j.out
#SBATCH --error=logs/gopc_search_%j.err

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

TARGET_DB="${TARGET_DB:-PHA_bioprospecting/gopc_search/target_db/gopc_db}"
if [ ! -f "${TARGET_DB}" ]; then
  echo "${TARGET_DB} not found -- build it first: cluster/run_gopc_build_db.sbatch" >&2
  exit 2
fi

FAMILY="${FAMILY:-all}"
SENSITIVITY="${SENSITIVITY:-7.0}"
EVALUE="${EVALUE:-1e-5}"
COVERAGE="${COVERAGE:-0.0}"
MAX_SEQS="${MAX_SEQS:-10000}"
CAP_WARNING_FRACTION="${CAP_WARNING_FRACTION:-0.95}"
echo "FAMILY=${FAMILY} SENSITIVITY=${SENSITIVITY} EVALUE=${EVALUE} COVERAGE=${COVERAGE} MAX_SEQS=${MAX_SEQS}"

pha-reference init-db
pha-reference export-queries
pha-reference gopc-search \
  --family "${FAMILY}" \
  --sensitivity "${SENSITIVITY}" \
  --evalue "${EVALUE}" \
  --coverage "${COVERAGE}" \
  --max-seqs "${MAX_SEQS}" \
  --cap-warning-fraction "${CAP_WARNING_FRACTION}" \
  --threads "${SLURM_CPUS_PER_TASK:-16}"
pha-reference gopc-combine

echo "Job finished: $(date)"
echo "Results: PHA_bioprospecting/gopc_search/results/"
echo "  <family>_hits.tsv, <family>_unique_targets.tsv, <family>_summary.tsv per family,"
echo "  plus all_families_unique_targets.tsv / all_families_summary.tsv."
echo "Check this log (or the summary.tsv files' query_hit_cap_warning column) for any"
echo "'query_hit_cap_warning' lines before trusting a family's results as sensitivity-complete."
