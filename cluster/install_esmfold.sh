#!/usr/bin/env bash
# ESMFold environment setup on the cluster -- adds to the EXISTING
# pha-reference env (the one `source cluster/env_activate.sh` activates)
# rather than creating a separate one: ESMFold only needs torch +
# transformers + accelerate + einops on top of what's already there
# (SQLAlchemy, httpx, typer, PyYAML, tenacity, rich) -- not enough new
# packages to justify a second env and a second activation step in every
# workflow. If that ever changes (a real dependency conflict surfaces),
# fall back to a dedicated env -- see the bottom of this file -- but
# don't default to that; try the shared env first.
#
# Uses HuggingFace transformers' EsmForProteinFolding rather than Meta's
# original fair-esm[esmfold] package deliberately: the original needs
# openfold's custom CUDA attention kernels compiled against the exact
# local CUDA/PyTorch version, a common source of cluster-install
# failures; the transformers port reimplements the same model in plain
# PyTorch with no custom kernel compilation step. See
# structure_prediction/run_esmfold.py's own docstring for the same note.
#
# Usage (default -- adds to the existing pha-reference env):
#   source cluster/env_activate.sh
#   ./cluster/install_esmfold.sh
#
# IMPORTANT -- scratch space, not $HOME: the ESMFold model weights
# (~2.7GB) download on first use into $HF_HOME (default ~/.cache/
# huggingface), which will blow a small $HOME quota. Point it at scratch
# BEFORE running this script and before every later run_esmfold.py call
# (interactive or via sbatch):
#   export HF_HOME=/scratch/morrill/users/hmp278/hf_cache
# cluster/run_esmfold.sbatch already sets this default itself, but set it
# here too if you pre-warm the cache interactively (see the end of this
# script) so the same cache gets reused instead of a second copy landing
# in $HOME.
#
# NOT tested end-to-end on an actual cluster GPU node yet -- written from
# documented transformers/ESMFold usage patterns. If a step here fails,
# that's expected to happen at least once; paste the error back rather
# than assuming this whole approach is wrong first.
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root (PHA_Ocean_Atlas/)

export HF_HOME="${HF_HOME:-/scratch/morrill/users/hmp278/hf_cache}"
mkdir -p "${HF_HOME}"

if [ -z "${VIRTUAL_ENV:-}" ] && [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
  echo "No active env detected -- run 'source cluster/env_activate.sh' first" >&2
  echo "(from cluster/setup_env.sh's pha-reference env), then re-run this script." >&2
  exit 2
fi
echo "Installing into active env: ${CONDA_DEFAULT_ENV:-${VIRTUAL_ENV}}"

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
echo "Environment ready (added to the existing pha-reference env, no new env created)."
echo "HF_HOME is set to ${HF_HOME} for this shell -- re-export it in any new shell/sbatch"
echo "job before running run_esmfold.py (run_esmfold.sbatch already does this)."
echo
echo "Consider pre-warming the model-weight cache once interactively before the first real"
echo "sbatch array job, so 20+ concurrent array tasks don't all try to download the same"
echo "~2.7GB of weights at once:"
echo "  HF_HOME=${HF_HOME} python -c \"from transformers import EsmForProteinFolding; EsmForProteinFolding.from_pretrained('facebook/esmfold_v1')\""
echo
echo "--- Fallback: dedicated env instead (only if the shared env approach hits a real"
echo "    conflict) -- put it on scratch, NOT \$HOME, same quota reasoning as above:"
echo "  conda create -p /scratch/morrill/users/hmp278/conda_envs/esmfold python=3.11"
echo "  conda activate /scratch/morrill/users/hmp278/conda_envs/esmfold"
echo "  <then re-run the pip installs above in that env>"
