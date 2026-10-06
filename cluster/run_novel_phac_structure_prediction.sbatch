#!/bin/bash
# ESMFold the 255 genuine missed-synthase finds (PHA_CLEAN_RESULTS.md
# section 13.3). A thin wrapper over the generic cluster/run_esmfold.sbatch
# path so this set can be submitted without remembering the FASTA/OUT_DIR/
# N_BATCHES export line -- and, unlike that script, as a plain single job
# rather than an array, because 255 sequences at a median 418 aa is roughly
# 25 minutes of GPU time at the ~5.9 sec/sequence rate confirmed on earlier
# runs. Splitting that across array tasks would cost more in scheduling and
# model-load overhead than it saves.
#
# On the cluster there is nothing to prepare: git pull is enough. Both inputs --
# structure_prediction/novel_phac_candidates.faa and its fold_manifest.tsv rows --
# are committed. Do NOT run build_novel_phac_fold_input.py here; it reads
# phac_recovery/novel_phac_new_candidates.faa, which .gitignore keeps out of the
# repo, so it will only report that file as missing. That script is a local step,
# re-run only when the candidate set itself changes.
#
# Prerequisite, once per checkout:
#   ./cluster/setup_env.sh and ./cluster/install_esmfold.sh, if this checkout has
#   not folded anything before. install_esmfold.sh must run from a LOGIN node --
#   it pre-downloads the model weights, and the GPU nodes have no internet access.
#
# Usage:
#   sbatch --account=191001-364393 cluster/run_novel_phac_structure_prediction.sbatch
#
# Then, once this finishes, the structural comparison against the folded
# audited reference set -- a separate CPU job, see that script's header:
#   sbatch --account=191001-364393 \
#     --export=ALL,QUERY_SETS=novel_phac_candidates \
#     cluster/run_foldseek.sbatch
#
# Resumable: run_esmfold.py skips any id whose .pdb already exists in
# OUT_DIR, so re-submitting after a partial failure only redoes the gaps.
#
#SBATCH --job-name=novel_phac_esmfold
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/novel_phac_esmfold_%j.out
#SBATCH --error=logs/novel_phac_esmfold_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs
export PYTHONUNBUFFERED=1

FASTA="${FASTA:-structure_prediction/novel_phac_candidates.faa}"
OUT_DIR="${OUT_DIR:-structure_prediction/esmfold_out/novel_phac_candidates}"

if [ ! -f "${FASTA}" ]; then
  echo "${FASTA} is missing -- run python structure_prediction/build_novel_phac_fold_input.py first." >&2
  exit 2
fi

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
source cluster/esmfold_pythonpath.sh   # PYTHONPATH, HF_HOME, PIP_CACHE_DIR -- all on scratch, not $HOME

# GPU nodes have no internet; use only the weights install_esmfold.sh cached.
export HF_HUB_OFFLINE=1

# Preflight: confirm the scratch package install matches THIS interpreter before
# burning a GPU allocation on it. torch's own failure mode here is actively
# misleading -- a Python-version mismatch makes its ABI-tagged _C extension
# invisible, Python falls back to the torch/_C/ stub directory, and torch reports
# "It appears that PyTorch has loaded the torch/_C folder of the PyTorch
# repository", which reads as a source-checkout problem rather than the env
# problem it is. Say what is actually wrong, and how to fix it.
if ! python -c 'import torch' 2>/dev/null; then
  echo "ERROR: torch does not import in this environment." >&2
  echo "  interpreter : $(python -c 'import sys; print(sys.executable)')" >&2
  echo "  version     : $(python -V 2>&1)" >&2
  echo "  packages    : ${PYTHONPATH%%:*}" >&2
  echo "  built for   : $(ls "${PYTHONPATH%%:*}"/torch/_C.cpython-*.so 2>/dev/null | xargs -n1 basename 2>/dev/null || echo 'no _C extension found at all')" >&2
  echo "If those Python versions disagree, the packages were installed by a different" >&2
  echo "interpreter than the one cluster/env_activate.sh now provides. Fix from a LOGIN" >&2
  echo "node: source cluster/env_activate.sh && ./cluster/install_esmfold.sh" >&2
  exit 3
fi

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "python: $(python -c 'import sys; print(sys.executable, sys.version.split()[0])')"
echo "torch:  $(python -c 'import torch; print(torch.__version__, "CUDA", torch.cuda.is_available())')"
echo "FASTA=${FASTA} OUT_DIR=${OUT_DIR}"
echo "sequences: $(grep -c '^>' "${FASTA}")"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

python structure_prediction/run_esmfold.py \
  --fasta "${FASTA}" \
  --out-dir "${OUT_DIR}" \
  --n-batches 1 \
  --batch-index 0

echo "Structures written: $(ls -1 "${OUT_DIR}"/*.pdb 2>/dev/null | wc -l)"
echo "Job finished: $(date)"
echo "Next: sbatch --account=191001-364393 --export=ALL,QUERY_SETS=novel_phac_candidates cluster/run_foldseek.sbatch"
