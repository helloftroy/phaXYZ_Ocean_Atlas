#!/bin/bash
# Lightweight figure-stage structure predictions for the 12 PhaC regional
# specialist representatives. This uses ColabFold in single-sequence mode so
# it does not need internet/MMseqs access from the GPU node.
#
# Usage from repo root:
#   sbatch --account=191001-364393 cluster/run_phac_regional_structures_colabfold.sbatch
#
# If your cluster exposes ColabFold as a module, set COLABFOLD_MODULE:
#   sbatch --export=ALL,COLABFOLD_MODULE=colabfold cluster/run_phac_regional_structures_colabfold.sbatch
#
# If colabfold_batch lives somewhere non-standard, set COLABFOLD_BIN:
#   sbatch --export=ALL,COLABFOLD_BIN=/path/to/colabfold_batch cluster/run_phac_regional_structures_colabfold.sbatch
#SBATCH --job-name=phac_regional_cf
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phac_regional_structures_%j.out
#SBATCH --error=logs/phac_regional_structures_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- run with sbatch from the repo root.}"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

if [ -n "${COLABFOLD_MODULE:-}" ]; then
  module load "${COLABFOLD_MODULE}"
fi

COLABFOLD_BIN="${COLABFOLD_BIN:-colabfold_batch}"
if ! command -v "${COLABFOLD_BIN}" >/dev/null 2>&1; then
  echo "Could not find colabfold_batch. Load a module or set COLABFOLD_BIN=/path/to/colabfold_batch." >&2
  exit 2
fi

python3 figures/scripts/prepare_regional_specialist_structure_inputs.py

INPUT="figures/structures/regional_specialists/phaC_regional_specialists_for_structure.faa"
OUT="figures/structures/regional_specialists/colabfold_single_sequence"
mkdir -p "${OUT}"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Using: $("${COLABFOLD_BIN}" --help 2>/dev/null | head -1 || true)"

"${COLABFOLD_BIN}" \
  --msa-mode single_sequence \
  --model-type alphafold2_ptm \
  --num-recycle 3 \
  --num-models 1 \
  --rank plddt \
  "${INPUT}" \
  "${OUT}"

echo "Job finished: $(date)"
find "${OUT}" -maxdepth 1 -type f \( -name '*.pdb' -o -name '*.cif' \) -print | sort
