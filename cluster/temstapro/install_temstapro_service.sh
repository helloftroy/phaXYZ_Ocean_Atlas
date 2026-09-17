#!/usr/bin/env bash
# Install TemStaPro and pre-download the ProtTrans model on an internet-capable
# service/login node. Run this before submitting GPU array jobs.
set -euo pipefail

PREFIX="${1:-$PWD/.external/TemStaPro}"
ENV_PREFIX="${2:-${TEMSTAPRO_ENV_PREFIX:-${PREFIX}/conda_env}}"
CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$(dirname "${PREFIX}")/conda_pkgs}"
export CONDA_PKGS_DIRS

mkdir -p "$(dirname "${PREFIX}")"
mkdir -p "${CONDA_PKGS_DIRS}"
if [ ! -d "${PREFIX}/.git" ]; then
  git clone https://github.com/ievapudz/TemStaPro.git "${PREFIX}"
else
  git -C "${PREFIX}" pull --ff-only
fi

# TemStaPro currently calls transformers.utils.logging.disable_progress_bar(),
# which is absent in some Python-3.7-compatible transformers builds. Make that
# call conditional instead of pinning an awkward transformer/Python matrix.
python - "${PREFIX}/prottrans_models.py" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
if path.exists():
    text = path.read_text()
    old = "hf_logging.disable_progress_bar()"
    new = "getattr(hf_logging, 'disable_progress_bar', lambda: None)()"
    if old in text and new not in text:
        path.write_text(text.replace(old, new))
        print(f"patched {path}")
PY

if [ ! -d "${ENV_PREFIX}/conda-meta" ]; then
  conda create -y -p "${ENV_PREFIX}" python=3.7
fi

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_PREFIX}"
set -u

conda install -y -c conda-forge transformers sentencepiece matplotlib
if [ "${TEMSTAPRO_GPU:-1}" = "1" ]; then
  conda install -y pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
else
  conda install -y pytorch torchvision torchaudio -c pytorch
fi
# Avoid PyTorch/MKL runtime mismatch seen as:
# libtorch_cpu.so: undefined symbol: iJIT_NotifyEvent
conda install -y "mkl<2024.1"

python - <<'PY'
import torch
print("torch import ok", torch.__version__, "cuda_available", torch.cuda.is_available())
PY

cd "${PREFIX}"
chmod +x ./temstapro

# Download the HuggingFace ProtTrans model into the TemStaPro checkout so GPU
# nodes do not need internet access.
python - <<'PY'
from pathlib import Path
from transformers import T5Tokenizer, T5EncoderModel
name = "Rostlab/prot_t5_xl_half_uniref50-enc"
out = Path("ProtTrans")
print("downloading", name, "->", out)
tokenizer = T5Tokenizer.from_pretrained(name, do_lower_case=False)
model = T5EncoderModel.from_pretrained(name)
tokenizer.save_pretrained(out)
model.save_pretrained(out)
print("done")
PY

echo "TemStaPro installed at: ${PREFIX}"
echo "Conda env: ${ENV_PREFIX}"
echo "Conda package cache: ${CONDA_PKGS_DIRS}"
echo "Submit jobs with:"
echo "  sbatch --export=ALL,TEMSTAPRO_DIR=${PREFIX},TEMSTAPRO_ENV=${ENV_PREFIX} cluster/temstapro/run_phac_temstapro_array.sbatch"
