# Running the PHA reference pipeline on a SLURM cluster

This pipeline needs internet (BRENDA's download endpoint + the UniProt REST
API) but **no GPU** -- Phase 1 is pure retrieval and normalization, no LLM
calls. `cluster/run_pha_reference_build.sbatch` runs on the same kind of
CPU/HTTPS-capable partition as
[`../fair_ocean_agent`](https://github.com/helloftroy/FAIRe_Ocean_Agent)'s
own `run_discovery.sbatch` (`service` in these scripts -- confirm the real
partition name for your account with `sinfo`).

## One-time setup

```bash
git clone https://github.com/helloftroy/phaXYZ_Ocean_Atlas.git
cd phaXYZ_Ocean_Atlas
./cluster/setup_env.sh
```

If you already have `fair_ocean_agent`'s own `faire-agent` conda env set up
on this cluster, it already carries most of this package's dependencies
(SQLAlchemy, httpx, typer, PyYAML, tenacity, rich) -- reuse it instead of
creating a new env:

```bash
CONDA_ENV_NAME=faire-agent ./cluster/setup_env.sh
```

`setup_env.sh` runs `pip install -e ".[dev]"` either way, which adds
`phaatlas`'s own remaining dependencies (ijson, python-dotenv) into
whichever env you point it at.

**Contact email** (sent as a polite User-Agent string on every BRENDA/
UniProt request -- unset falls back to a placeholder, fine for testing but
not for a real run):

```bash
cp .env.example .env
# edit .env: PHA_REFERENCE_CONTACT_EMAIL=you@example.org
```

## Test with phaC first, then run everything

Per the spec: pressure-test the whole workflow against one family before
running all sixteen.

Add `--account=<your account>` to every `sbatch` command below if your
cluster requires one for this partition/QOS (confirmed live: omitting it
can silently route the job through a different default QOS with its own,
possibly tighter, resource enforcement -- a job that OOM-kills or fails to
schedule without `--account` may work fine at the exact same resource
request once submitted correctly).

```bash
mkdir -p logs
sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_pha_reference_build.sbatch
```

Check progress:

```bash
tail -f logs/pha_reference_build_<job_id>.out
```

Once that looks right (check `PHA_reference/exports/pha_reference_master.csv`
and the per-tier counts `pha-reference status` prints at the end of the
job), run every family:

```bash
sbatch --account=191001-364393 cluster/run_pha_reference_build.sbatch   # FAMILY defaults to 'all'
```

Re-running is safe: every `pha-reference ingest-*` command upserts against
the same database (protein rows are keyed by UniProt accession or BRENDA
local id; family_assignment is keyed by protein+family and only ever
widens, never narrows, its evidence tier -- see `pipeline/upsert.py`), so a
failed or interrupted job can just be resubmitted, and re-running with a
newer BRENDA release or an edited `family_definitions.yaml` re-syncs
in place without duplicating rows.

## Pulling results back

```bash
scp -r <cluster>:<path-to-repo>/PHA_reference ./
```

`PHA_reference/pha_reference.sqlite` is the canonical database;
`PHA_reference/exports/pha_reference_master.csv` and `pha_reference.faa`
are both regenerated from it by `pha-reference export` and never hand-edited.

## Notes

- **BRENDA download**: `pha-reference download-brenda` (or `ingest-brenda`,
  which calls it automatically) POSTs the same license-acceptance form the
  website itself uses (`accept-license=1` + the JSON-file button) to fetch
  `brenda_<release>.json.tar.gz` -- no account/login needed. The extracted
  JSON (~700MB) is cached under `data/brenda/`; re-running reuses it unless
  you pass `--force`.
- **Families without `ec_brenda`** (phasins, regulators, synthase
  partners) are silently skipped by `ingest-brenda` -- that's expected,
  per spec; `ingest-uniprot` is their only route.
- **phaR disambiguation**: requesting `phaR_synthase` or `phaR_regulator`
  alone still automatically pulls in its partner family for the joint
  `gene:phaR` disambiguation pass (see `pipeline/ingest_uniprot.py`) -- you
  never need to remember to request both yourself.
- **"database is locked"** after an interrupted job: same SQLite recovery
  steps as `../fair_ocean_agent/cluster/README.md`'s own Troubleshooting
  section apply here (`PRAGMA quick_check` after backing up the file) --
  more likely if `PHA_reference/` sits on NFS-mounted cluster storage.
