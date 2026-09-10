# PHA_bioprospecting

Phase 2: cluster the PHA reference database at 95% identity, and stage the
two discovery databases (GOPC, OMDBv2) it will later be searched against.

## Layout

```
PHA_bioprospecting/
├── reference -> ../PHA_reference        # symlink, NOT a copy (see below)
├── databases/
│   ├── GOPC/
│   │   ├── raw/                          # gitignored -- re-downloadable
│   │   └── download_manifest.tsv         # committed
│   └── OMDBv2/
│       ├── OMDBv2.0_*                    # gitignored -- re-downloadable
│       └── download_manifest.tsv         # committed
└── scripts/
    ├── lib_download.sh                   # shared: resumable fetch, checksum, manifest
    ├── download_gopc.sh
    └── download_omdb.sh
```

`reference/` is a symlink to `../PHA_reference` (the same canonical
database and exports Phase 1 built), not a copy -- there is exactly one
`pha_reference.sqlite` in this repo. Clustering writes back into that same
file; nothing here forks it.

## Reference clustering (95% identity / 90% coverage, per family)

```bash
./cluster/install_mmseqs2.sh          # once, if not already done
pha-reference cluster95 --family phaC # test one family first, per the same
                                       # philosophy as Phase 1's ingest step
pha-reference export                  # master CSV now carries cluster95_*
pha-reference export-cluster95        # pha_reference_all.faa / _nr95.faa / _cluster95.tsv

pha-reference cluster95 --family all  # then everything
pha-reference export && pha-reference export-cluster95
```

Or as a cluster job: `sbatch cluster/run_reference_cluster95.sbatch`
(`FAMILY=phaC` to test one family first, defaults to `all`).

**Exact-dedup, then cluster.** Sequences are grouped by `sequence_sha256`
first; MMseqs2 only ever sees one representative per unique sequence
within a family (`easy-cluster --min-seq-id 0.95 -c 0.90 --cov-mode 0`).
Every actual protein/family row sharing that hash is then stamped with the
same `cluster95_id` afterward, so `cluster95_size` reflects the real
number of database entries in the cluster, not just how many distinct
sequences MMseqs2 saw. **Clustering is scoped within each `pha_family`
independently** -- two different families' sequences never merge into the
same cluster even if similar. **Nothing is deleted**: `cluster95_id`/
`cluster95_representative`/`cluster95_size` are three new columns on
`family_assignment` (added via an idempotent `ALTER TABLE`, run
automatically by `init-db`), populated in place.

## GOPC and OMDBv2 downloads

Both scripts are resumable (`curl -C -`), verify against a published
checksum where one exists, and always append a row to their
`download_manifest.tsv` (filename, source URL, size, checksum, checksum
status, download timestamp) -- manifests are never overwritten, only
appended to, so a re-download after an upstream file changes still leaves
the old record for comparison.

```bash
./scripts/download_gopc.sh md5        # tiny -- safe to run anywhere, fetches the published checksums
./scripts/download_omdb.sh data       # ~7.5MB genome->sample->study link table -- safe anywhere

# The following are genuinely large -- run via cluster/run_download_databases.sbatch instead
# of interactively, and confirm free scratch space first:
./scripts/download_gopc.sh geneset          # GOPC.geneset.pep.fa.gz, ~184GB
./scripts/download_omdb.sh nr100            # OMDBv2.0_AA_G_NR100.faa.gz, ~49GB
./scripts/download_omdb.sh nr100-clusters   # OMDBv2.0_AA_G_NR100.cluster.tsv.gz, ~4.2GB
```

As a cluster job:

```bash
sbatch --export=ALL,DATASET=gopc,TARGET=geneset cluster/run_download_databases.sbatch
sbatch --export=ALL,DATASET=omdb,TARGET=nr100 cluster/run_download_databases.sbatch
sbatch --export=ALL,DATASET=omdb,TARGET=nr100-clusters cluster/run_download_databases.sbatch
```

Each is independent and resumable, so they're safe to submit in parallel
or resubmit after a timeout (walltime is generous but these are real
multi-hour transfers over the network, not compute).

**Sources** (both confirmed live with real `HTTP 206 Partial Content`
responses to range requests, i.e. genuinely resumable, not just
theoretically):

- **GOPC** -- the gene/protein-catalog FASTA output of the Global Ocean
  Microbiome Catalogue (GOMC). Chen, Jia, Sun et al., "Global marine
  microbial diversity and its potential in bioprospecting," *Nature* 633,
  371-379 (2024), doi:[10.1038/s41586-024-07891-2](https://doi.org/10.1038/s41586-024-07891-2).
  Hosted as CNGBdb dataset `MDB0000002`; files served from CNGB's own
  public FTP-over-HTTPS mirror (`ftp.cngb.org`), no login required.
  `GOPC.geneset.pep.fa.gz` is ~184GB (183,959,496,042 bytes) -- larger
  than the earlier working assumption of "small," confirmed via `HEAD`.
  MD5 published in the dataset's own `md5.txt`.
- **OMDBv2** -- Ocean Microbiomics Database v2, Sunagawa Lab / ETH Zurich
  (same long-running project as Paoli, Ruscheweyh, Forneris et al.,
  "Biosynthetic potential of the global ocean microbiome," *Nature* 607,
  111-118 (2022), doi:[10.1038/s41586-022-04862-3](https://doi.org/10.1038/s41586-022-04862-3);
  v2.0 is a later data release of the same portal -- no separate
  publication specific to v2 could be confirmed). Files served from ETH's
  own file server (`sunagawalab.ethz.ch`), no login required. Only
  `OMDBv2.0_data.tsv.gz` has a published checksum (MD5, from the portal's
  suppl_info page); the two large `NR100` files have none published --
  `download_omdb.sh` computes and records our own SHA256 for those
  instead, so at least a later re-download of *our* copy can be checked
  for corruption even without an upstream hash to compare against.

Deliberately not decompressed on download -- `.gz` stays `.gz` until a
downstream search tool (MMseqs2/DIAMOND/etc. against these catalogs, not
part of this phase yet) actually needs it unpacked, per the instruction
that motivated this whole layout: don't duplicate enormous files onto disk
before something actually needs the duplicate.

## What's committed vs. what's not

Committed: this README, both download scripts + their shared library, the
`reference` symlink, and both `download_manifest.tsv` files (small,
they're the provenance record). Never committed: anything under
`databases/GOPC/raw/` or the downloaded `OMDBv2.0_*` files themselves --
all multi-GB to multi-hundred-GB, and fully reproducible from the
manifests' own URLs.
