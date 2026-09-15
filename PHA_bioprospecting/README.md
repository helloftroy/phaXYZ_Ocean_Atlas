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

Or as a cluster job: `sbatch --account=191001-364393 cluster/run_reference_cluster95.sbatch`
(`--export=ALL,FAMILY=phaC` to test one family first, defaults to `all`).
Add `--account=<your account>` to every `sbatch` command on this page if
your cluster requires one for this partition/QOS -- omitting it can
silently route the job through a different default QOS with tighter
resource enforcement than what you actually requested.

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
sbatch --account=191001-364393 --export=ALL,DATASET=gopc,TARGET=geneset cluster/run_download_databases.sbatch
sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100 cluster/run_download_databases.sbatch
sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100-clusters cluster/run_download_databases.sbatch
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
pha-reference gopc-build-db PHA_bioprospecting/databases/GOPC/raw/GOPC.geneset.pep.fa.gz --gpu-compatible   # once

pha-reference gopc-search --family phaC --gpu --mmseqs-bin cluster/bin-gpu/mmseqs   # test one family first
pha-reference gopc-search --family all --gpu --mmseqs-bin cluster/bin-gpu/mmseqs    # then everything
pha-reference gopc-combine                # all_families_unique_targets.tsv / all_families_summary.tsv
```

**GPU vs. CPU search.** The actual search (`gopc-search`) runs on
`gpu-a100` via mmseqs' own `--gpu 1` (GPU-accelerated prefilter, see
Kallenborn et al., "GPU-accelerated homology search with MMseqs2,"
bioRxiv 2024.11.13.623350) -- `service` in this project is specifically
for stages that need internet, and this search needs none. This requires
its own GPU-enabled mmseqs binary (separate from the plain CPU one) and a
target database built with `--gpu-compatible` (`--createdb-mode 2`,
mmseqs' GPU-compatible storage format) -- a plain database can't be
searched with `--gpu`, and vice versa. **We have not been able to run
this GPU path ourselves** (no CUDA GPU in this development environment) --
it's built from the exact flags confirmed present in `mmseqs --help` for
this release plus that paper's description, not a working run we watched
succeed, so treat the resource numbers in `cluster/run_gopc_search.sbatch`
as a first guess and report back what actually happens. The plain CPU
path (`pha-reference gopc-search` without `--gpu`, `--split-memory-limit`
instead) still works as a fallback if GPU search doesn't pan out -- see
the git history for that version if needed.

**`gopc-build-db` memory**: `mmseqs createdb`'s own `--shuffle` default
(on) reorders the whole input up front, which costs extra memory at
GOPC's real scale (hundreds of millions of sequences) for no benefit to
a one-time build step -- `--no-shuffle` is the `phaatlas` default for a
plain build. Confirmed live, though: mmseqs **refuses** `--no-shuffle`
together with `--gpu-compatible` ("Shuffle database cannot be turned off
for --createdb-mode 2") and silently re-enables shuffle regardless -- so
the GPU-compatible build needs more memory than the plain build does at
the same `--mem`, not the same. `sinfo -p gpu-a100 -o "%P %m %c %l %G"`
confirmed each `gpu-a100` node has ~1031735MB (~1007GiB) RAM, 128 CPUs,
and 8 A100 GPUs (`gpu:a100:8`) -- shared across up to 8 concurrent
single-GPU jobs, so `cluster/run_gopc_build_db_gpu.sbatch` and
`run_gopc_search.sbatch` both request a considerate ~1/8 share
(~120-129000MB, 16 CPUs) for their single `--gres=gpu:1` rather than a
full node's worth; there's room to raise `--mem` further (up to the full
~1TB) if a specific run needs it and the cluster isn't busy. Separately, a
job OOM-killed or failing to schedule at a given `--mem` is also worth
double-checking against missing `--account` first (confirmed live:
omitting it can silently route the job through a different default QOS
with its own, possibly tighter, resource enforcement) before assuming
`--mem` itself needs to go up.

**Why the GPU-compatible build also runs on `gpu-a100`, not `service`,
even though `--createdb-mode 2` itself doesn't strictly need CUDA to
write** (confirmed live: it completes fine with no GPU present) -- there
was no way to confirm the resulting database is actually fully correct
for a later `--gpu` search without a real GPU to test against, and
getting that wrong after a ~184GB rebuild would cost far more than
running this one step on a GPU node it may not have strictly needed. When
unsure, match the resource type the consuming step needs.

Or as cluster jobs (build once, then search) -- add `--account=<your
account>` if your cluster requires one for this partition/QOS:

```bash
# One-time: GPU-enabled mmseqs binary + GPU-compatible target database
MMSEQS_VARIANT=linux-gpu MMSEQS_OUTPUT_DIR=cluster/bin-gpu ./cluster/install_mmseqs2.sh
sbatch --account=191001-364393 cluster/run_gopc_build_db_gpu.sbatch

sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_gopc_search.sbatch   # test one family first
sbatch --account=191001-364393 cluster/run_gopc_search.sbatch                             # FAMILY defaults to 'all'
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

**If you fall back to the CPU search path** (`gopc-search` without
`--gpu`, against a plain non-`--gpu-compatible` target database),
`--split-memory-limit` matters a lot at GOPC's scale, under SLURM.
Confirmed live on `service` before this moved to GPU: the search died
right after loading the (tiny) query database -- `Error: Prefilter died` /
`Error: Search died` -- because without this flag, mmseqs sizes its
target-database split against the *node's* total physical memory, not the
job's actual cgroup allocation. Also confirmed live: even WITH it set (80%
of a 64G job), it still died -- past mmseqs' own "can this fit" estimate
check, i.e. a genuine crash during execution, meaning that estimate
doesn't cover everything mmseqs actually uses at runtime; a 100G job with
a 60% split limit (~60G) is the last CPU configuration tried, itself
unconfirmed since the GPU path superseded it before a retry. Pass
`--split-memory-limit` comfortably below whatever `--mem` you actually
request, and check your own cluster's real ceiling with `sinfo -p
<partition> -o "%P %m %c %l"` rather than assuming these numbers.

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

## Resuming after a late-stage failure (e.g. convertalis OOM)

**Confirmed live: resubmitting `gopc-search`/`omdb-search` now reuses a
prior search that already finished, and only redoes the final
formatting step.** `run_mmseqs_only` no longer calls mmseqs' `easy-search`
wrapper (which creates a fresh, randomly-named tmp subdirectory on every
invocation, making resume impossible) -- it calls the lower-level
`createdb`/`search`/`convertalis` modules directly against a FIXED path
under `tmp/<family>/`. This is what fixed phaA: it was OOM-killed in
convertalis (the final TSV-formatting step) well after the expensive
prefilter+alignment work had already finished successfully, and the old
code lost the entire ~2-day run. Now, if `tmp/<family>/result.dbtype`
already exists from a completed search, a resubmit skips straight to
convertalis (seconds to minutes, not hours/days) instead of rescanning
the target from scratch. Verified with a dedicated test
(`test_run_mmseqs_only_resumes_convertalis_without_rerunning_search`)
that this reproduces byte-identical output to an uninterrupted run.

**This does NOT help if the job is killed while `search` itself is still
scanning the target** (as opposed to having already finished) -- that
work is not preserved, and a resubmit restarts the scan from zero, same
as before. If a family's full query set can't realistically finish a
single search inside your cluster's job time limit at all, that's what
batched search (below) is for.

## When one family's search doesn't fit in one job (batched search)

If a family's full query set can't realistically finish inside your
cluster's job time limit in one shot even with the resume behavior above,
split it into independent batches instead -- each batch still scans all
of GOPC once (so total GPU-time goes up, not down), but each is short
enough to actually complete rather than looping forever without progress
(and each batch's own search is separately resumable the same
convertalis-only way once it completes):

```bash
pha-reference export-queries
for i in 0 1 2 3; do   # N_BATCHES=4 example
  pha-reference gopc-search-batch --family phaC --n-batches 4 --batch-index "$i" --gpu --mmseqs-bin cluster/bin-gpu/mmseqs
done
pha-reference gopc-finalize-batches --family phaC --n-batches 4
```

Or as cluster jobs -- a SLURM job array for the batches, then a finalize
job chained to run once the whole array completes:

```bash
sbatch --account=191001-364393 --array=0-3 --export=ALL,FAMILY=phaC,N_BATCHES=4 \
  cluster/run_gopc_search_batched.sbatch
# note the printed array job id, e.g. 12345:
sbatch --account=191001-364393 --dependency=afterok:12345 \
  --export=ALL,FAMILY=phaC,N_BATCHES=4 cluster/run_gopc_finalize_batches.sbatch
```

`N_BATCHES` must match on both submissions and must equal the `--array`
range size (`--array=0-3` means 4 batches, indices 0..3). Output format
(`<family>_hits.tsv`/`_unique_targets.tsv`/`_summary.tsv`) is identical to
an unbatched run -- batching is purely how the compute gets scheduled, not
a change in what's found (verified: a batched run against a small test
target produced byte-identical `unique_targets.tsv` rows to an unbatched
run over the same data).

**How many batches?** No formula yet -- there isn't enough real timing
data to derive one. Pick `N_BATCHES` so that (this family's query count /
`N_BATCHES`) is in the same ballpark as a family that already finished
successfully in one 48h shot, and adjust from what you actually observe;
each batch's own log shows `Query database size: N` near the top, same as
an unbatched run, so you can compare directly.

The same batching applies to OMDB -- swap `gopc-search-batch`/
`gopc-finalize-batches` for `omdb-search-batch`/`omdb-finalize-batches`
(and the corresponding `run_omdb_search_batched.sbatch`/
`run_omdb_finalize_batches.sbatch`), though it's less likely to be needed
given OMDB's ~10x smaller size.

## Searching OMDB (same pipeline, against OMDBv2.0_AA_G_NR100 instead of GOPC)

**Why OMDB, not just GOPC**: GOPC's gene catalog is built by predicting
genes from each sample's assembled contigs and concatenating them --
there's no clean per-protein link back to which genome/sample/location it
came from. OMDB keeps that chain intact (`OMDBv2.0_data.tsv.gz` maps
GENOME → SAMPLE → STUDY, and `OMDBv2.0_AA_G_NR100.cluster.tsv.gz` maps an
NR100 representative back to every genome carrying that exact sequence),
so it's the right catalog to search when the actual question is "which
organisms/locations/depths carry this" rather than just "does a similar
protein exist somewhere."

Identical commands to the `gopc-*` ones throughout this doc, with `gopc`
swapped for `omdb` (`omdb-build-db`, `omdb-search`, `omdb-search-batch`,
`omdb-finalize-batches`, `omdb-combine`), reusing the exact same
underlying pipeline code (`pipeline/gopc_search.py`'s functions all take
`target_db_path` as a plain parameter -- nothing GOPC-specific in the
logic itself) and the exact same `queries/` directory (the NR95 reference
FASTAs don't depend on which catalog they're searched against). Results
land in a separate `PHA_bioprospecting/omdb_search/` so they never
collide with GOPC's:

```bash
# One-time: download + build (target file is OMDBv2.0_AA_G_NR100.faa.gz, not the geneset)
sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100 cluster/run_download_databases.sbatch
sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100-clusters cluster/run_download_databases.sbatch
sbatch --account=191001-364393 cluster/run_omdb_build_db_gpu.sbatch

sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_omdb_search.sbatch   # test one family first
sbatch --account=191001-364393 cluster/run_omdb_search.sbatch                             # FAMILY defaults to 'all'
```

**Expected to be faster/lighter than GOPC** (~10x fewer sequences: 249.5M
vs. GOPC's ~2.46B, so hopefully no sharding/chunking bottleneck) but this
has not been confirmed live yet -- `run_omdb_build_db_gpu.sbatch`/
`run_omdb_search.sbatch` currently request the exact same `--mem`/
`--cpus-per-task` as their GOPC counterparts (safer to over-provision once
more than guess a smaller number down and repeat the OOM churn already
hit twice on GOPC); tune down once real OMDB numbers are in hand. Watch
the first real run's timing and `nvidia-smi` GPU memory/utilization
before assuming it behaves the same as GOPC did.

Once you have `<family>_unique_targets.tsv` from an OMDB search, joining
back to `OMDBv2.0_data.tsv.gz` (by whatever genome/sample identifier the
target IDs carry -- not yet confirmed exactly what format OMDB's target
IDs take, unlike GOPC's, since no OMDB search has been run yet) is what
actually answers "where/what organism was this found in" -- that join
itself isn't built yet, worth doing once real OMDB results exist to
design it against actual ID formats rather than guessing.

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
