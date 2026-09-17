# PhaC Regional Specialist Structure Cartoons

Inputs come from `figures/phaC_regional_specialists_sequences.csv`.

Prepared files:

- `phaC_regional_specialists_for_structure.faa` - combined FASTA for all 12 representative proteins.
- `inputs/*.faa` - one FASTA per cluster.
- `phaC_regional_specialists_structure_metadata.tsv` - cluster/ecology metadata plus raw and cleaned sequence lengths.

Sequence cleanup:

- Each exported sequence ended with one terminal `X`.
- The structure input script trims terminal `X` only and records `trimmed_terminal_X:1`.
- No other noncanonical residues remain after cleanup.

Cluster prediction:

```bash
cd /path/to/PHA_Ocean_Atlas
git pull
sbatch --account=191001-364393 cluster/run_phac_regional_structures_colabfold.sbatch
```

The SLURM script runs ColabFold with:

```bash
colabfold_batch --msa-mode single_sequence --model-type alphafold2_ptm --num-recycle 3 --num-models 1
```

This is intentionally lightweight and does not require internet from the GPU node. It is meant for small figure cartoons, not final structural analysis.

If the cluster uses modules:

```bash
sbatch --account=191001-364393 \
  --export=ALL,COLABFOLD_MODULE=colabfold \
  cluster/run_phac_regional_structures_colabfold.sbatch
```

If `colabfold_batch` is installed at a specific path:

```bash
sbatch --account=191001-364393 \
  --export=ALL,COLABFOLD_BIN=/path/to/colabfold_batch \
  cluster/run_phac_regional_structures_colabfold.sbatch
```

Rendering cartoons:

After predictions finish, render transparent PNG cartoons with PyMOL:

```bash
python3 figures/scripts/render_regional_specialist_structure_cartoons.py
```

If PyMOL is a module, load it first. The renderer writes individual PNGs and a contact sheet under `cartoons/`.
