# PHA_bioprospecting

Phase 2: cluster the PHA reference database at 95% identity, stage the two
discovery databases (GOPC, OMDBv2) it's searched against, and search the
NR95 reference set against GOPC. Does NOT classify hits as true PHA
proteins yet -- that (and OMDB's genome-resolved search) is later work.

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
├── gopc_search/                          # gitignored -- fully regenerable, see below
│   ├── queries/<family>.faa
│   ├── target_db/gopc_db*
│   └── results/<family>_{hits,unique_targets,summary}.tsv + all_families_*.tsv
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

Deliberately not decompressed on download -- `.gz` stays `.gz`; `mmseqs
createdb` (below) reads gzipped FASTA directly, so it's never unpacked to
disk at all.

## Searching GOPC (sensitive MMseqs2, per family)

Does NOT classify hits as true PHA proteins -- alignment statistics only,
retained in full so filtering/classification can happen later without
rerunning the (expensive) search itself.

```bash
pha-reference export-queries              # one queries/<family>.faa per family, headers = protein_id only
pha-reference gopc-build-db PHA_bioprospecting/databases/GOPC/raw/GOPC.geneset.pep.fa.gz   # once

pha-reference gopc-search --family phaC   # test one family first
pha-reference gopc-search --family all    # then everything
pha-reference gopc-combine                # all_families_unique_targets.tsv / all_families_summary.tsv
```

**`gopc-build-db` memory**: confirmed live, `mmseqs createdb`'s own
`--shuffle` default (on) got this step OOM-killed at 64GB given GOPC's
real scale (hundreds of millions of sequences) -- `--no-shuffle` is now
the `phaatlas` default (pass `--shuffle` to re-enable it, trading memory
for better target-split load-balancing in later searches), and
`cluster/run_gopc_build_db.sbatch` requests 256GB. Raise further if that
still isn't enough on your cluster.

Or as cluster jobs (build once, then search):

```bash
sbatch cluster/run_gopc_build_db.sbatch
sbatch --export=ALL,FAMILY=phaC cluster/run_gopc_search.sbatch   # test one family first
sbatch cluster/run_gopc_search.sbatch                             # FAMILY defaults to 'all'
```

**Search settings** (`pipeline/gopc_search.py`, all configurable via CLI
flags / the sbatch job's own env vars): `-s 7.0` (MMseqs2's very-sensitive
end), `-e 1e-5`, `-c 0.0` -- deliberately no coverage floor at search time,
so a fragmented metagenomic protein isn't lost before `qcov`/`tcov` are
even recorded; filter on those afterward (e.g. "≥50% query coverage" or
"≥70% bidirectional") without rerunning the search. `--alignment-mode 3`
so reported `pident` comes from the real alignment, not MMseqs2's faster
estimate.

**`--max-seqs` (default 10000, configurable) is a prefilter cap, not a
biological result** -- if a query's hit count lands at/near that cap, its
true hit set may have been truncated rather than genuinely exhausted. Each
family's own `run_family_search` checks this automatically: if any query's
hit count reaches ≥95% of `--max-seqs` (also configurable,
`--cap-warning-fraction`), the job logs a `WARNING: query_hit_cap_warning`
line naming that family and how many queries hit the cap, and
`<family>_summary.tsv` records `query_hit_cap_warning=True` plus the exact
`queries_at_cap` list. Expected most for phaA/phaB/phaJ (embedded in
larger enzyme superfamilies); phaC should behave more cleanly. Rerun just
the flagged family at a higher `--max-seqs` (50000/100000) rather than
redoing everything.

**Per-family outputs**, all under `gopc_search/results/`:

- `<family>_hits.tsv` -- every query-target alignment, all 14 columns
  (`query,target,evalue,bits,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov`).
- `<family>_unique_targets.tsv` -- one row per unique GOPC protein hit,
  its single best alignment (lowest e-value, ties broken by bitscore)
  retained, plus `n_reference_queries_matching` (how many distinct PHA
  reference queries hit this same target).
- `<family>_summary.tsv` -- `n_reference_queries`, `n_alignments`,
  `n_unique_GOPC_targets`, counts at e-value ≤1e-5/1e-10/1e-20, counts at
  query-coverage ≥0.5 and bidirectional-coverage ≥0.5, identity quantiles
  (min/p10/p25/median/p75/p90/max), and the cap-warning fields above.

Then combined across every family searched so far: `results/
all_families_unique_targets.tsv` and `results/all_families_summary.tsv`.
**A GOPC target hit by more than one family is deliberately NOT
deduplicated** -- both families' rows are kept in
`all_families_unique_targets.tsv` rather than picking one; resolving that
conflict is later work, not this stage's job.

## What's committed vs. what's not

Committed: this README, both download scripts + their shared library, the
`reference` symlink, and both `download_manifest.tsv` files (small,
they're the provenance record). Never committed: anything under
`databases/GOPC/raw/`, the downloaded `OMDBv2.0_*` files, or anything
under `gopc_search/` (queries/target_db/results -- all regenerated from
`pha_reference.sqlite` + the downloaded GOPC target by the commands
above). The search *code* lives in `phaatlas/pipeline/gopc_search.py`
(committed, same repo as everything else) --
all multi-GB to multi-hundred-GB, and fully reproducible from the
manifests' own URLs.
