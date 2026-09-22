#!/usr/bin/env bash
# One-time ESMFold environment setup on the cluster -- separate from
# cluster/setup_env.sh's main "pha-reference" conda env deliberately: that
# env is explicitly kept free of GPU/heavy-ML dependencies (BRENDA/UniProt
# retrieval doesn't need them), and ESMFold's dependency footprint (torch,
# transformers, a CUDA build) is large enough to deserve its own env
# rather than risk breaking the main one.
#
# Uses HuggingFace transformers' EsmForProteinFolding rather than Meta's
# original fair-esm[esmfold] package -- the original needs openfold's
# custom CUDA attention kernels compiled against the exact local
# CUDA/PyTorch version, a common source of cluster-install failures never
# tested end-to-end here; the transformers port reimplements the same
# model in plain PyTorch with no custom kernel compilation step. See
# structure_prediction/run_esmfold.py's own docstring for the same note.
#
# Usage:
#   ./cluster/install_esmfold.sh
#   # or, to reuse existing scratch space for the (large, ~5GB+ with model
#   # weights cached) env / HF cache:
#   CONDA_ENV_PREFIX=/scratch/$USER/conda_envs/esmfold ./cluster/install_esmfold.sh
#
# NOT tested end-to-end on an actual cluster GPU node yet -- written from
# documented transformers/ESMFold usage patterns. If a step here fails,
# that's expected to happen at least once; paste the error back rather
# than assuming this whole approach is wrong first.
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root (PHA_Ocean_Atlas/)

CONDA_ENV_NAME="${CONDA_ENV_NAME:-esmfold}"
CONDA_ENV_PREFIX="${CONDA_ENV_PREFIX:-}"
MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
CONDA_SH="${MINIFORGE_HOME}/etc/profile.d/conda.sh"

if [ ! -f "${CONDA_SH}" ]; then
  echo "conda.sh not found at ${CONDA_SH}." >&2
  echo "Set MINIFORGE_HOME if your (mini)conda/anaconda install lives elsewhere." >&2
  exit 2
fi
# shellcheck source=/dev/null
source "${CONDA_SH}"

if [ -n "${CONDA_ENV_PREFIX}" ]; then
  CONDA_CREATE_FLAG=(-p "${CONDA_ENV_PREFIX}")
  CONDA_ACTIVATE_TARGET="${CONDA_ENV_PREFIX}"
else
  CONDA_CREATE_FLAG=(-n "${CONDA_ENV_NAME}")
  CONDA_ACTIVATE_TARGET="${CONDA_ENV_NAME}"
fi

if conda env list | grep -qE "^\S*${CONDA_ACTIVATE_TARGET//\//\\/}\s"; then
  echo "Conda env at/named '${CONDA_ACTIVATE_TARGET}' already exists -- reusing it."
else
  echo "Creating conda env: ${CONDA_ACTIVATE_TARGET} (python 3.11) ..."
  conda create -y "${CONDA_CREATE_FLAG[@]}" python=3.11
fi

conda activate "${CONDA_ACTIVATE_TARGET}"
echo "Active env: $(python -c 'import sys; print(sys.prefix)')"

# CUDA 12.1 build -- matches what's been used elsewhere on this project's
# gpu-a100 partition (see cluster/install_mmseqs2.sh's GPU-variant notes);
# adjust the index-url below if the cluster's actual CUDA toolkit differs.
echo "Installing PyTorch (CUDA 12.1 build) ..."
pip install --index-url https://download.pytorch.org/whl/cu121 torch

echo "Installing transformers + supporting packages ..."
pip install "transformers>=4.35" accelerate einops

echo
echo "Verifying import + CUDA visibility (no model download yet) ..."
python - <<'PYEOF'
import torch
print("torch:", torch.__version__, "CUDA available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device:", torch.cuda.get_device_name(0))
import transformers
print("transformers:", transformers.__version__)
from transformers import EsmForProteinFolding  # import-only check, no download
print("EsmForProteinFolding import OK")
PYEOF

echo
echo "Environment ready. The ESMFold model weights (~2.7GB) download on first"
echo "use via structure_prediction/run_esmfold.py, cached under \$HF_HOME (default"
echo "~/.cache/huggingface) -- set HF_HOME=/scratch/\$USER/hf_cache first if \$HOME"
echo "has a small quota. Consider pre-warming the cache once interactively before"
echo "the first real sbatch array job, so 20+ concurrent array tasks don't all"
echo "try to download the same weights at once:"
echo "  python -c \"from transformers import EsmForProteinFolding; EsmForProteinFolding.from_pretrained('facebook/esmfold_v1')\""
