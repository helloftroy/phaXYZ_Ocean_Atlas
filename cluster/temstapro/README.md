# phaC TemStaPro Cluster Workflow

Run these commands from the `PHA_Ocean_Atlas` repository root on the cluster.

## 1. Put required inputs inside this repo

The large phaC data files are not committed to git. Copy or symlink them into:

```bash
mkdir -p data/temstapro_inputs
```

Required for FASTA chunk export:

```bash
data/temstapro_inputs/phaC_cluster_sequences.faa
```

Required later for metadata merge:

```bash
data/temstapro_inputs/phaC_unique_targets_with_metadata_depth.tsv
data/temstapro_inputs/phaC_genomes_woa23_annual_temperature.tsv
```

If the files already live elsewhere on scratch, symlinks are fine:

```bash
ln -s /path/to/phaC_cluster_sequences.faa data/temstapro_inputs/phaC_cluster_sequences.faa
ln -s /path/to/phaC_unique_targets_with_metadata_depth.tsv data/temstapro_inputs/phaC_unique_targets_with_metadata_depth.tsv
ln -s /path/to/phaC_genomes_woa23_annual_temperature.tsv data/temstapro_inputs/phaC_genomes_woa23_annual_temperature.tsv
```

## 2. Prepare TemStaPro inputs

```bash
python3 figures/scripts/export_phac_temstapro_inputs.py
ls temstapro/chunks/*.faa | wc -l
```

This writes chunked FASTA files to `temstapro/chunks/`.

## 3. Install TemStaPro on an internet-capable service node

```bash
bash cluster/temstapro/install_temstapro_service.sh /scratch/morrill/users/hmp278/TemStaPro
```

This clones TemStaPro, creates the conda environment under scratch, and
downloads ProtTrans under scratch:

```bash
/scratch/morrill/users/hmp278/TemStaPro/conda_env
/scratch/morrill/users/hmp278/TemStaPro/ProtTrans
```

The installer also keeps the conda package cache in scratch by default:

```bash
/scratch/morrill/users/hmp278/conda_pkgs
```

If PyTorch fails with `libtorch_cpu.so: undefined symbol: iJIT_NotifyEvent`,
repair the existing scratch environment on a service/login node:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/morrill/users/hmp278/TemStaPro/conda_env
conda install -y "mkl<2024.1"
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
```

If TemStaPro fails with
`AttributeError: module 'transformers.utils.logging' has no attribute 'disable_progress_bar'`,
patch the existing checkout:

```bash
python - <<'PY'
from pathlib import Path
p = Path("/scratch/morrill/users/hmp278/TemStaPro/prottrans_models.py")
text = p.read_text()
p.write_text(text.replace(
    "hf_logging.disable_progress_bar()",
    "getattr(hf_logging, 'disable_progress_bar', lambda: None)()",
))
print(f"patched {p}")
PY
```

If GPU jobs fail with `Connection error, and we cannot find the requested
files in the cached path`, pre-populate the HF cache on a service/login node:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/morrill/users/hmp278/TemStaPro/conda_env
cd /scratch/morrill/users/hmp278/TemStaPro
python - <<'PY'
from pathlib import Path
from transformers import T5Tokenizer, T5EncoderModel
name = "Rostlab/prot_t5_xl_half_uniref50-enc"
cache = Path("ProtTrans_cache")
tok = T5Tokenizer.from_pretrained(name, do_lower_case=False, cache_dir=str(cache))
model = T5EncoderModel.from_pretrained(name, cache_dir=str(cache))
T5Tokenizer.from_pretrained(name, do_lower_case=False, cache_dir=str(cache), local_files_only=True)
T5EncoderModel.from_pretrained(name, cache_dir=str(cache), local_files_only=True)
print(f"cached and offline-verified {name} in {cache.resolve()}")
PY
```

Then submit with `TEMSTAPRO_PT_CACHE`:

```bash
sbatch --array=0-128 \
  --export=ALL,TEMSTAPRO_DIR=/scratch/morrill/users/hmp278/TemStaPro,TEMSTAPRO_ENV=/scratch/morrill/users/hmp278/TemStaPro/conda_env,TEMSTAPRO_PT_CACHE=/scratch/morrill/users/hmp278/TemStaPro/ProtTrans_cache \
  cluster/temstapro/run_phac_temstapro_array.sbatch
```

To choose a different env path, pass it as the second argument:

```bash
bash cluster/temstapro/install_temstapro_service.sh \
  /scratch/morrill/users/hmp278/TemStaPro \
  /scratch/morrill/users/hmp278/conda_envs/temstapro
```

## 4. Submit the GPU array

Adjust the array range to match the chunk count from step 2.

```bash
sbatch --array=0-128 \
  --export=ALL,TEMSTAPRO_DIR=/scratch/morrill/users/hmp278/TemStaPro,TEMSTAPRO_ENV=/scratch/morrill/users/hmp278/TemStaPro/conda_env \
  cluster/temstapro/run_phac_temstapro_array.sbatch
```

Outputs go to `temstapro/temstapro_outputs/`.

## 5. Merge predictions with metadata

```bash
python3 figures/scripts/merge_phac_temstapro_outputs.py
```

Final merged output:

```bash
temstapro/phaC_temstapro_predictions_with_metadata.tsv
```
