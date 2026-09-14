#!/bin/bash
# GPU-accelerated MMseqs2 search of the NR95 PHA reference set against
# GOPC, per family: export-queries -> gopc-search --gpu -> gopc-combine.
# No internet needed -- everything here reads local files.
#
# Runs on gpu-a100 (mirroring ../fair_ocean_agent/cluster/run_extraction.sbatch's
# own GPU job), not `service` -- `service` in this project is specifically
# for the stages that need internet (see cluster/run_download_databases.sbatch,
# run_pha_reference_build.sbatch), and this search needs none.
#
# One-time prerequisites (both via cluster/install_mmseqs2.sh, see there):
#   1. A GPU-enabled mmseqs binary -- installed separately from the plain
#      CPU one used for cluster/run_gopc_build_db.sbatch:
#        MMSEQS_VARIANT=linux-gpu MMSEQS_OUTPUT_DIR=cluster/bin-gpu ./cluster/install_mmseqs2.sh
#   2. A GPU-COMPATIBLE target database -- a database built the normal way
#      (without --gpu-compatible) cannot be searched with --gpu; build one once:
#        sbatch --account=191001-364393 cluster/run_gopc_build_db_gpu.sbatch
#
# Test ONE family first, same philosophy as every other stage in this
# pipeline -- phaC is the cleanest case (a specific EC, not a shared
# superfamily), so it's the best first check that alignment settings and
# --max-seqs are behaving before spending compute on all 16:
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_gopc_search.sbatch
# Then everything:
#   sbatch --account=191001-364393 cluster/run_gopc_search.sbatch   # FAMILY defaults to 'all'
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
#
# IMPORTANT, please read: this GPU path (mmseqs' own --gpu 1, described in
# Kallenborn et al., "GPU-accelerated homology search with MMseqs2,"
# bioRxiv 2024.11.13.623350) has NOT been run end-to-end by us before --
# there's no way to test it without an actual CUDA GPU, which this
# development environment doesn't have. Everything here is based on the
# exact flags confirmed present in `mmseqs --help` for this specific
# release (18-8cc5c) plus that paper's description, not a working run we
# watched succeed. If this errors out in a way that looks resource-related,
# --mem is the first thing to reconsider, and please report back exactly
# what happens so this can be corrected quickly.
#
# --mem/--cpus-per-task: `sinfo -p gpu-a100 -o "%P %m %c %l %G"` reports
# each gpu-a100 node as ~1031735MB (~1007GiB) RAM, 128 CPUs, 8 A100 GPUs
# (gpu:a100:8) -- shared across up to 8 concurrent single-GPU jobs. Since
# we only request --gres=gpu:1, below is a considerate ~1/8 share
# (1031735/8 ~ 129000MB, 128/8 = 16 CPUs) rather than a full node's worth
# for one GPU -- raise --mem if the search itself (not just the build
# step, which has its own separate job) turns out to need more.
#
# --time=48:00:00: confirmed live that GPU jobs on this account are capped
# at 48h regardless of gpu-a100's own nominal "infinite" TIMELIMIT (a
# QOS-level restriction tied to the account/allocation, not the raw
# partition config sinfo reports) -- `scontrol update ... TimeLimit=`
# past 48h was rejected. If a family's search doesn't finish in one 48h
# job, resubmit the exact same command -- it reuses the same tmp/<family>/
# directory, so at minimum query-db creation is skipped; whether mmseqs'
# own GPU prefilter resumes partway through a single target-database scan
# (vs. restarting that scan from zero) has NOT been confirmed by us, so
# watch whether a resubmitted run reaches further than the last one did
# before assuming it's making incremental progress. If a family is taking
# MUCH longer than the others (more than ~a day for phaC specifically, the
# cleanest/smallest family, would be worth investigating) it's worth
# checking whether the GPU is actually engaged rather than assuming it's
# just slow: `squeue -u $USER` for the node, then e.g. `srun
# --jobid=<jobid> --pty nvidia-smi` -- near-0% GPU utilization AND
# near-empty GPU memory would mean --gpu 1 isn't actually engaging; high
# memory usage with low utilization (confirmed live for phaC: ~80GB/80GB
# VRAM used at 1% util) more likely means it's data-movement-bound loading
# GOPC in large chunks, not stuck.
#SBATCH --job-name=pha_gopc_search
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
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

# Deliberately a SEPARATE path from the CPU target_db/ -- a GPU-compatible
# database (--createdb-mode 2) is a different on-disk format and can't be
# used interchangeably with the plain one cluster/run_gopc_build_db.sbatch
# builds by default.
TARGET_DB="${TARGET_DB:-PHA_bioprospecting/gopc_search/target_db_gpu/gopc_db}"
if [ ! -f "${TARGET_DB}" ]; then
  echo "${TARGET_DB} not found -- build a GPU-compatible target database first:" >&2
  echo "  sbatch --account=191001-364393 --export=ALL,TARGET_DB=${TARGET_DB} cluster/run_gopc_build_db_gpu.sbatch" >&2
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
pha-reference gopc-search \
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
pha-reference gopc-combine

echo "Job finished: $(date)"
echo "Results: PHA_bioprospecting/gopc_search/results/"
echo "  <family>_hits.tsv, <family>_unique_targets.tsv, <family>_summary.tsv per family,"
echo "  plus all_families_unique_targets.tsv / all_families_summary.tsv."
echo "Check this log (or the summary.tsv files' query_hit_cap_warning column) for any"
echo "'query_hit_cap_warning' lines before trusting a family's results as sensitivity-complete."
