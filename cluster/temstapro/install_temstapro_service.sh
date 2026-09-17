#!/usr/bin/env bash
# Install TemStaPro and pre-download the ProtTrans model on an internet-capable
# service/login node. Run this before submitting GPU array jobs.
set -euo pipefail

PREFIX="${1:-$PWD/.external/TemStaPro}"
ENV_NAME="${TEMSTAPRO_ENV_NAME:-temstapro_env}"

mkdir -p "$(dirname "${PREFIX}")"
if [ ! -d "${PREFIX}/.git" ]; then
  git clone https://github.com/ievapudz/TemStaPro.git "${PREFIX}"
else
  git -C "${PREFIX}" pull --ff-only
fi

if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  conda create -y -n "${ENV_NAME}" python=3.7
fi

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"
set -u

conda install -y -c conda-forge transformers sentencepiece matplotlib
if [ "${TEMSTAPRO_GPU:-1}" = "1" ]; then
  conda install -y pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
else
  conda install -y pytorch torchvision torchaudio -c pytorch
fi

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
echo "Conda env: ${ENV_NAME}"
echo "Submit jobs with:"
echo "  sbatch --export=ALL,TEMSTAPRO_DIR=${PREFIX},TEMSTAPRO_ENV=${ENV_NAME} cluster/temstapro/run_phac_temstapro_array.sbatch"
