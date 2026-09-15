#!/bin/bash
# GPU-accelerated MMseqs2 search of the NR95 PHA reference set against
# OMDBv2.0_AA_G_NR100, per family: export-queries -> omdb-search --gpu ->
# omdb-combine. Mirrors cluster/run_gopc_search.sbatch exactly, just
# pointed at OMDB. No internet needed.
#
# One-time prerequisites (both via cluster/install_mmseqs2.sh, see there --
# the GPU binary is the SAME one gopc search already uses if you've set
# that up, no need to install it twice):
#   1. A GPU-enabled mmseqs binary:
#        MMSEQS_VARIANT=linux-gpu MMSEQS_OUTPUT_DIR=cluster/bin-gpu ./cluster/install_mmseqs2.sh
#   2. A GPU-compatible OMDB target database:
#        sbatch --account=191001-364393 cluster/run_omdb_build_db_gpu.sbatch
#
# Test ONE family first, same philosophy as every other stage in this
# pipeline:
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_omdb_search.sbatch
# Then everything:
#   sbatch --account=191001-364393 cluster/run_omdb_search.sbatch   # FAMILY defaults to 'all'
#
# EXPECTED faster/lighter than GOPC: OMDBv2.0_AA_G_NR100 is ~10x fewer
# sequences (249.5M vs. GOPC's ~2.46B), so hopefully no sharding/chunking
# bottleneck like GOPC needed. Not confirmed live yet, though -- watch the
# first real run's timing and GPU memory/utilization
# (`nvidia-smi`/`srun --jobid=<jobid> --pty nvidia-smi`) before assuming
# it behaves the same as GOPC did. --mem/--cpus-per-task below are left at
# the SAME values as run_gopc_search.sbatch for now (see that script's own
# header for the sinfo-derived ~1/8-share reasoning) -- safer to
# over-provision once more than guess a smaller number down and hit
# another OOM; tune down once real OMDB numbers are in hand.
#
# --time=48:00:00: same account-level GPU cap confirmed live for GOPC
# applies here too.
#SBATCH --job-name=pha_omdb_search
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/omdb_search_%j.out
#SBATCH --error=logs/omdb_search_%j.err

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

MMSEQS_GPU_BIN="${MMSEQS_GPU_BIN:-$(pwd)/cluster/bin-gpu/mmseqs}"
if [ ! -x "${MMSEQS_GPU_BIN}" ]; then
  echo "${MMSEQS_GPU_BIN} not found/executable -- install it first:" >&2
  echo "  MMSEQS_VARIANT=linux-gpu MMSEQS_OUTPUT_DIR=cluster/bin-gpu ./cluster/install_mmseqs2.sh" >&2
  exit 2
fi
echo "mmseqs (GPU build): $("${MMSEQS_GPU_BIN}" version)"
echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi || echo "nvidia-smi not available/failed -- is a GPU actually attached to this job? (--gres=gpu:1 above)" >&2

TARGET_DB="${TARGET_DB:-PHA_bioprospecting/omdb_search/target_db_gpu/omdb_db}"
if [ ! -f "${TARGET_DB}" ]; then
  echo "${TARGET_DB} not found -- build a GPU-compatible target database first:" >&2
  echo "  sbatch --account=191001-364393 --export=ALL,TARGET_DB=${TARGET_DB} cluster/run_omdb_build_db_gpu.sbatch" >&2
  exit 2
fi

FAMILY="${FAMILY:-all}"
SENSITIVITY="${SENSITIVITY:-7.0}"
EVALUE="${EVALUE:-1e-5}"
COVERAGE="${COVERAGE:-0.0}"
MAX_SEQS="${MAX_SEQS:-10000}"
CAP_WARNING_FRACTION="${CAP_WARNING_FRACTION:-0.95}"
echo "FAMILY=${FAMILY} SENSITIVITY=${SENSITIVITY} EVALUE=${EVALUE} COVERAGE=${COVERAGE} MAX_SEQS=${MAX_SEQS} TARGET_DB=${TARGET_DB}"

pha-reference init-db
pha-reference export-queries
pha-reference omdb-search \
  --family "${FAMILY}" \
  --target-db "${TARGET_DB}" \
  --sensitivity "${SENSITIVITY}" \
  --evalue "${EVALUE}" \
  --coverage "${COVERAGE}" \
  --max-seqs "${MAX_SEQS}" \
  --cap-warning-fraction "${CAP_WARNING_FRACTION}" \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  --mmseqs-bin "${MMSEQS_GPU_BIN}" \
  --gpu
pha-reference omdb-combine

echo "Job finished: $(date)"
echo "Results: PHA_bioprospecting/omdb_search/results/"
echo "  <family>_hits.tsv, <family>_unique_targets.tsv, <family>_summary.tsv per family,"
echo "  plus all_families_unique_targets.tsv / all_families_summary.tsv."
echo "Check this log (or the summary.tsv files' query_hit_cap_warning column) for any"
echo "'query_hit_cap_warning' lines before trusting a family's results as sensitivity-complete."
