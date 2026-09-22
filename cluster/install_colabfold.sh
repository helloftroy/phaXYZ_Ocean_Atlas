#!/usr/bin/env bash
# ColabFold environment setup -- a DEDICATED conda env on scratch, not an
# addition to the shared pha-reference/faire-agent env (unlike ESMFold's
# install_esmfold.sh, which added its ~4 packages to the existing env via
# PYTHONPATH). Different call here on purpose: ColabFold pulls in JAX +
# a whole pinned fork of AlphaFold (haiku, dm-tree, ml-collections, a
# specific jax/jaxlib combo tied to the local CUDA version) -- a much
# larger and more version-sensitive dependency tree than ESMFold's, and
# one that's genuinely likely to conflict with something the shared env's
# other users depend on. A dedicated env avoids that risk entirely, at
# the cost of one extra `conda activate` step when actually running it --
# worth it here. Uses this project's existing miniforge3 install (same
# one cluster/setup_env.sh uses), not a second miniforge -- just a new
# env prefix under it.
#
# Why ColabFold over raw AlphaFold2: AlphaFold2's own database-search
# path needs ~2.2TB of downloaded genetic databases (BFD, MGnify, PDB70,
# Uniclust30, UniRef90...) -- not remotely proportional to this project's
# actual need (a few hundred to a few thousand "hard" structures ESMFold
# handled poorly). ColabFold's MMseqs2-based MSA generation defaults to
# a free remote API (api.colabfold.com) -- no local database download at
# all for a first pass. Structure inference itself still runs locally on
# GPU either way. Fall back to local MSA databases only if the remote
# server's rate limits become a real bottleneck at this dataset's scale.
#
# NOT tested end-to-end on an actual cluster GPU node yet -- written from
# documented ColabFold/JAX usage patterns, same honesty as
# install_esmfold.sh. JAX's CUDA-version sensitivity is a well-known
# common source of install friction for this tool specifically (not just
# this script's own uncertainty) -- expect at least one round of pasting
# an error back.
#
# Usage:
#   ./cluster/install_colabfold.sh
#   # override the env location if scratch isn't at this path for you:
#   SCRATCH_BASE=/scratch/morrill/users/hmp278 ./cluster/install_colabfold.sh
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root (PHA_Ocean_Atlas/)

SCRATCH_BASE="${SCRATCH_BASE:-/scratch/morrill/users/hmp278}"
CONDA_ENV_PREFIX="${CONDA_ENV_PREFIX:-${SCRATCH_BASE}/conda_envs/colabfold}"
export COLABFOLD_DATA_DIR="${COLABFOLD_DATA_DIR:-${SCRATCH_BASE}/colabfold_weights}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-${SCRATCH_BASE}/pip_cache}"
MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
CONDA_SH="${MINIFORGE_HOME}/etc/profile.d/conda.sh"

if [ ! -f "${CONDA_SH}" ]; then
  echo "conda.sh not found at ${CONDA_SH}." >&2
  echo "Set MINIFORGE_HOME if your (mini)conda/anaconda install lives elsewhere." >&2
  exit 2
fi
# shellcheck source=/dev/null
source "${CONDA_SH}"

mkdir -p "${COLABFOLD_DATA_DIR}" "${PIP_CACHE_DIR}" "$(dirname "${CONDA_ENV_PREFIX}")"

if conda env list | grep -qF "${CONDA_ENV_PREFIX}"; then
  echo "Conda env already exists at ${CONDA_ENV_PREFIX} -- reusing it."
else
  echo "Creating conda env at ${CONDA_ENV_PREFIX} (python 3.10 -- ColabFold's own"
  echo "documented recommendation as of this writing; newer Python has caused"
  echo "resolver issues for some of its pinned deps in the past) ..."
  conda create -y -p "${CONDA_ENV_PREFIX}" python=3.10
fi
conda activate "${CONDA_ENV_PREFIX}"
echo "Active env: $(python -c 'import sys; print(sys.prefix)')"

echo "Installing ColabFold (includes the AlphaFold fork it depends on) ..."
pip install "colabfold[alphafold]"

# Modern extras-based CUDA install (bundles its own matching CUDA runtime
# via pip deps), same lesson just learned installing ESMFold's torch: don't
# force an old hardcoded CUDA wheel-index URL, let pip resolve a build that
# actually matches what's available. Fall back to the older
# `-f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`
# form only if this fails to find a CUDA-enabled build.
echo "Installing JAX (CUDA-enabled) ..."
pip install -U "jax[cuda12]"

echo
echo "Verifying import (no model download yet) ..."
python - <<'PYEOF'
import jax
print("jax:", jax.__version__, "devices:", jax.devices())
import colabfold
print("colabfold import OK")
PYEOF

echo
echo "Environment ready at ${CONDA_ENV_PREFIX}."
echo "COLABFOLD_DATA_DIR=${COLABFOLD_DATA_DIR} -- model weights (~a few GB) download here on"
echo "first real run, not into \$HOME. Re-export this same COLABFOLD_DATA_DIR (and re-activate"
echo "this env) in every later shell/sbatch job that calls colabfold_batch."
echo
echo "If jax.devices() above printed CpuDevice instead of a GPU device: expected if you're on"
echo "a login node (no GPU attached there); only a real problem if still CPU-only inside an"
echo "actual --gres=gpu:1 job."
echo
echo "Quick smoke test once you have a GPU allocation (uses the free remote MSA server, no"
echo "local database needed):"
echo "  conda activate ${CONDA_ENV_PREFIX}"
echo "  export COLABFOLD_DATA_DIR=${COLABFOLD_DATA_DIR}"
echo "  echo '>test' > /tmp/test.fasta && echo 'MSEQVENCEHERE...' >> /tmp/test.fasta"
echo "  colabfold_batch /tmp/test.fasta /tmp/test_out"
