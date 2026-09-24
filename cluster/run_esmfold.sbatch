#!/bin/bash
# ESMFold structure prediction, batched array job. Run once per input
# FASTA (uncertain_cluster_representatives.faa, positive_control_
# representatives.faa, audited_reference_set.faa -- see
# structure_prediction/build_fold_inputs.py), pointing FASTA/OUT_DIR at
# each in turn.
#
# --array MUST cover exactly 0..N_BATCHES-1. Sizing guide -- the
# uncertain/positive_control pair below has now actually been run
# end-to-end (confirmed live: completed without crashing, and ran
# noticeably faster than the Foldseek step that follows it -- see
# cluster/run_foldseek.sbatch), so N_BATCHES for those two is closer to
# a real number than a guess; audited_reference_set.faa and
# all_phac_dedup_representatives.faa are still sized by extrapolation:
#   uncertain_cluster_representatives.faa    (605 seqs, median 475aa, 7 >700aa) -> N_BATCHES=6 (confirmed working)
#   positive_control_representatives.faa     (428 seqs, median 420aa) -> N_BATCHES=4 (confirmed working)
#   audited_reference_set.faa              (1,875 seqs, median 578aa) -> try N_BATCHES=20
#   all_phac_dedup_representatives.faa    (15,880 seqs, median 444aa) -> N_BATCHES=20 (see below)
#   legionella_top2_representatives.faa      (18 seqs, the W1046/FDAARGOS_200 deep-dive
#     pair's own copies, section 9.2) -> N_BATCHES=1, single non-array job is fine
# uncertain_cluster_representatives.faa was rebuilt 2026-09-24 by
# structure_prediction/build_uncertain_set_v2.py -- a hand-curated 605-
# sequence set (down from the original 4,992), see that script's docstring
# and PHA_CLEAN_RESULTS.md section 9.6 for what it contains and why.
# all_phac_dedup_representatives.faa (structure_prediction/
# build_all_phac_dedup_set.py) is a separate, broader run covering nearly
# all verified phaC diversity (short sequences and >=95%-identical near-
# duplicates collapsed out) -- meant as a comprehensive structural atlas
# for future figures, not a hand-picked subset; point OUT_DIR at
# structure_prediction/esmfold_out/all_phac_dedup for it.
#
# N_BATCHES=20 for this one is deliberately NOT the naive scale-up from the
# confirmed 605-sequence run (~101 seqs/task, ~10min/task there -> a naive
# 26x scale-up to 15,880 seqs would keep ~10min/task but need N_BATCHES=160).
# 160 short-lived array tasks is wasteful of SLURM scheduling overhead for
# not much benefit -- confirmed live that ~10min/task is too fine-grained,
# so this instead targets ~1-2hr/task: 20 batches -> ~794 seqs/task -> ~79
# minutes/task at the same ~5.9sec/seq rate, comfortably under the 12-hour
# limit even accounting for this set's slightly longer tail (max 1024aa vs.
# 741aa before). Re-derive if the confirmed per-sequence rate turns out
# different here (longer median could shift this either way) -- watch the
# first array task's own stdout (logs per-sequence time) before assuming
# 20 is still right at a much larger scale.
#
# Usage:
#   sbatch --account=191001-364393 --array=0-5 \
#     --export=ALL,FASTA=structure_prediction/uncertain_cluster_representatives.faa,OUT_DIR=structure_prediction/esmfold_out/uncertain,N_BATCHES=6 \
#     cluster/run_esmfold.sbatch
#   # the comprehensive all_phac_dedup run:
#   sbatch --account=191001-364393 --array=0-19 \
#     --export=ALL,FASTA=structure_prediction/all_phac_dedup_representatives.faa,OUT_DIR=structure_prediction/esmfold_out/all_phac_dedup,N_BATCHES=20 \
#     cluster/run_esmfold.sbatch
#   # the tiny Legionella top-2-genomes run (18 seqs, a couple minutes):
#   sbatch --account=191001-364393 --array=0-0 \
#     --export=ALL,FASTA=structure_prediction/legionella_top2_representatives.faa,OUT_DIR=structure_prediction/esmfold_out/legionella_top2,N_BATCHES=1 \
#     cluster/run_esmfold.sbatch
#
# Resumable: run_esmfold.py skips any target_id whose .pdb already exists
# in OUT_DIR, so re-submitting the same array after a partial failure (OOM
# on a few long sequences, a node preemption, etc.) only redoes the gaps.
#
#SBATCH --job-name=pha_esmfold
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/esmfold_%A_%a.out
#SBATCH --error=logs/esmfold_%A_%a.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is unset -- submit with --array=0-<N_BATCHES-1>, not a plain sbatch.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs
export PYTHONUNBUFFERED=1

: "${FASTA:?Set FASTA=<path to one of the three structure_prediction/*.faa input files>.}"
: "${OUT_DIR:?Set OUT_DIR=<output directory for the PDB files from this FASTA>.}"
: "${N_BATCHES:?Set N_BATCHES=<total batch count> -- must match the --array range used (0..N_BATCHES-1).}"

if [ ! -f cluster/env_activate.sh ]; then
  echo "cluster/env_activate.sh is missing -- run ./cluster/setup_env.sh first, then" >&2
  echo "./cluster/install_esmfold.sh to add ESMFold's packages to that same env." >&2
  exit 2
fi
# shellcheck source=/dev/null
source cluster/env_activate.sh

if [ ! -f cluster/esmfold_pythonpath.sh ]; then
  echo "cluster/esmfold_pythonpath.sh is missing -- run ./cluster/install_esmfold.sh first" >&2
  echo "(it writes this file pointing PYTHONPATH at the scratch package install)." >&2
  exit 2
fi
# shellcheck source=/dev/null
source cluster/esmfold_pythonpath.sh   # sets PYTHONPATH, HF_HOME, PIP_CACHE_DIR -- all on scratch, not $HOME

# GPU compute nodes have no internet access. install_esmfold.sh pre-downloads
# the model weights into $HF_HOME on the login node for exactly this reason --
# HF_HUB_OFFLINE=1 here makes transformers use only that local cache and fail
# immediately with a clear "not found in cache" error if it's missing, rather
# than this job hanging (or timing out slowly) trying to reach huggingface.co
# from a node that can't. If you see that error: run ./cluster/install_esmfold.sh
# again from a login node first.
export HF_HUB_OFFLINE=1

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Array task: ${SLURM_ARRAY_TASK_ID}/${N_BATCHES}"
echo "FASTA=${FASTA} OUT_DIR=${OUT_DIR}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

python structure_prediction/run_esmfold.py \
  --fasta "${FASTA}" \
  --out-dir "${OUT_DIR}" \
  --n-batches "${N_BATCHES}" \
  --batch-index "${SLURM_ARRAY_TASK_ID}"

echo "Job finished: $(date)"
echo "Batch ${SLURM_ARRAY_TASK_ID}/${N_BATCHES} done for ${FASTA}."
