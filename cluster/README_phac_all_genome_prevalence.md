# phaC Prevalence With All-Genome Denominator

This workflow estimates:

```text
% phaC = distinct phaC-positive OMDB genomes / all OMDB genomes
```

by geographic grid cell. It replaces the older denominator in
`figures/scripts/plot_phac_pct_by_region.py`, which only counted genomes
that had at least one searched PHA-pathway gene.

Run from the `PHA_Ocean_Atlas` repository root on a service/login node with
internet access for the OMDB metadata API step.

## 1. Make sure the small OMDB catalog is present

```bash
cd PHA_bioprospecting/scripts
./download_omdb.sh data
cd ../..
```

Expected file:

```bash
PHA_bioprospecting/databases/OMDBv2/OMDBv2.0_data.tsv.gz
```

## 2. Export all OMDB genomes with locations

```bash
python3 figures/scripts/export_omdb_all_genome_locations.py
```

Output:

```bash
data/all_genomes/omdb_all_genomes_with_locations.tsv
```

This queries OMDB sample metadata to attach latitude/longitude to every genome
in the OMDB catalog.

## 3. Export all phaC-positive genomes from NR100 clusters

This builds an uncapped numerator by streaming the full OMDB NR100 cluster
membership table once.

Make sure the cluster-membership file exists:

```bash
cd PHA_bioprospecting/scripts
./download_omdb.sh nr100-clusters
cd ../..
```

Then run:

```bash
python3 figures/scripts/export_phac_genomes_from_omdb_clusters.py
```

Output:

```bash
data/all_genomes/phaC_all_genomes_from_nr100_clusters.tsv
```

## 4. Plot phaC prevalence using all genomes as denominator

The plot script uses the uncapped numerator table by default:

```bash
data/all_genomes/phaC_all_genomes_from_nr100_clusters.tsv
```

Run:

```bash
python3 figures/scripts/plot_phac_pct_all_genomes_by_region.py
```

Outputs:

```bash
data/all_genomes/phaC_prevalence_all_genomes_by_region.tsv
figures/phaC_pct_all_genomes_by_region.png
figures/phaC_pct_all_genomes_by_region.pdf
```

If you cannot stream the NR100 cluster file yet, the plot can fall back to a
metadata table, but that may undercount genomes if metadata enrichment was run
with a per-target cap:

```bash
python3 figures/scripts/plot_phac_pct_all_genomes_by_region.py \
  --phac-metadata /path/to/phaC_unique_targets_with_metadata_depth.tsv
```
