# PHA_Ocean_Atlas -- PHA reference database (Phase 1)

A reproducible reference database of known/annotated PHA-related proteins
(BRENDA + UniProt), family-config driven. Built to later identify PHA
proteins in GOPC and OMDB.

## Canonical storage

`PHA_reference/pha_reference.sqlite` is the source of truth. Everything
else is generated from it and never hand-edited:

- `PHA_reference/exports/pha_reference_master.csv` -- one flattened row per
  (protein, PHA family), from the `protein_master_export` SQL view.
- `PHA_reference/exports/pha_reference.faa` -- FASTA, generated straight
  from `protein.sequence`.

Neither the database nor the exports are committed to git (see
`.gitignore`) -- both are fully reproducible from BRENDA + UniProt by
re-running the pipeline below. Same for `data/brenda/` (the ~700MB
extracted BRENDA JSON dump, cached, re-downloadable).

## Quickstart (local pressure test)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

cp .env.example .env   # optional locally; set PHA_REFERENCE_CONTACT_EMAIL for a real cluster run

pha-reference init-db

# Test the whole workflow against phaC first, per spec:
pha-reference ingest-uniprot --family phaC
pha-reference download-brenda          # one-time, ~80MB compressed / ~700MB extracted, cached under data/brenda/
pha-reference ingest-brenda --family phaC

pha-reference export
pha-reference status
```

Then run every family:

```bash
pha-reference ingest-brenda --family all
pha-reference ingest-uniprot --family all
pha-reference export
```

Offline unit tests (no network calls):

```bash
python -m pytest tests/unit -q
```

## Family configuration

Nothing family-specific (EC numbers, gene aliases, UniProt search terms) is
hardcoded in `phaatlas/*.py` -- it all comes from
[`config/family_definitions.yaml`](config/family_definitions.yaml). Add or
edit a family there; no code changes needed. All 16 families from the
spec are defined: `phaA, phaB, phaC, phaE, phaR_synthase, phaJ, phaG, phaZ,
phaY, phaP, phaF, phaI, phaM, phaD, phaR_regulator, phaQ`.

`phaR_synthase` and `phaR_regulator` are deliberately separate families
that happen to share the exact gene name `phaR` -- see
`pipeline/ingest_uniprot.py`'s `ingest_phaR_disambiguation`, which resolves
every `gene:phaR` hit against both families' `disambiguation_terms` rather
than assuming they're the same protein.

Monomer/polymer controlled vocabularies live in
[`config/vocab/monomers.yaml`](config/vocab/monomers.yaml) and
[`polymers.yaml`](config/vocab/polymers.yaml); grow them the same way.

## Evidence tiers and the "demonstrated" evidence rule

Every family_assignment gets one of:

- `GOLD_EXPERIMENTAL` -- literature/BRENDA-backed or a highly specific
  experimental UniProt annotation.
- `CURATED_REVIEWED` -- reviewed Swiss-Prot PHA annotation, no
  protein-specific experiment found.
- `ANNOTATED_UNREVIEWED` -- TrEMBL/unreviewed PHA annotation.
- `CANDIDATE_AMBIGUOUS` -- matched only through a broad/shared-chemistry
  signal (an EC class or gene name that also covers non-PHA proteins) with
  no confirming PHA-context text found. Never forced into a higher tier.

Separately, `phenotype.demonstrated_monomers`/`demonstrated_polymers` are
populated **only** when `phenotype_evidence_level` is
`experimental_literature` or `curated_specific_function` -- a generic
Rhea/catalytic-activity annotation or a family-membership-only inference
never populates those two fields, even if the raw text happens to name a
specific monomer. See `phaatlas/phenotype/normalize.py`.

## Database schema

- `family_definition` -- mirrors `config/family_definitions.yaml`.
- `protein` -- one row per UniProt accession, or per BRENDA protein record
  with no accession at all (kept, not discarded -- `BRENDA:<ec>:<local_id>`
  as its protein_id).
- `family_assignment` -- one row per (protein, family); `retrieval_reasons`
  widens (pipe-separated) and `evidence_tier` only ever improves as
  independent searches re-confirm the same protein.
- `source_evidence` -- one row per BRENDA/UniProt source record, raw
  payload retained in full (`raw_payload_json`).
- `phenotype` -- one row per literature/annotation statement (polymer/
  monomer composition, substrate specificity, temperature/pH/kinetics),
  original sentence always preserved in `phenotype_raw`.
- `retrieval_run` -- append-only provenance: source, release, query,
  count, timestamp, software version. Never overwritten.
- `protein_master_export` (view) -- the flattened join behind the CSV
  export.

## A note on `PHA_reference/` vs. the `phaatlas` package, on macOS

The importable Python package here is `phaatlas`, not `pha_reference` --
deliberately. macOS's default APFS volume is case-insensitive, and a
package literally named `pha_reference` would resolve to the *same
directory entry* as the canonical `PHA_reference/` data directory required
by spec, differing only in case. That's not just a naming clash: it broke
`pip install -e .`'s package discovery outright (confirmed live -- even
`setuptools.packages.find`'s glob patterns matched `PHA_reference/`
case-insensitively at the OS level). Real HPC storage (ext4/Lustre/GPFS/
NFS, i.e. the cluster this is meant to run on) is case-sensitive and
wouldn't hit this at all, but the package still keeps the distinct name so
local development on macOS isn't fighting the filesystem.

## Running on the cluster

See [`cluster/README.md`](cluster/README.md). No GPU needed -- this stage
is pure BRENDA/UniProt retrieval, so it runs on the same kind of
CPU/HTTPS-capable partition as
[`fair_ocean_agent`](https://github.com/helloftroy/FAIRe_Ocean_Agent)'s own
`run_discovery.sbatch`, and can reuse that project's conda env if you
already have one set up.
