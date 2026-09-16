"""Genome x PHA-family presence/count matrix, collapsed into short
"pathway architecture" labels (e.g. "ABC", "ABCJ") -- which combinations
of PHA genes co-occur in the same genome, how common each combination is,
and which taxa carry it.

Input: one or more <family>_unique_targets_with_metadata.tsv files (the
output of pipeline/omdb_metadata.py's enrich_unique_targets, one row per
target_id x genome, joined to that genome's GTDB taxonomy and sample
location/ecosystem). This module only reads those files -- no new network
calls or database access.

Scope note: a genome's architecture label only reflects the families
whose <family>_unique_targets_with_metadata.tsv has actually been built
(via `pha-reference omdb-enrich-metadata`) and fed into this module. A
family missing from the input set is indistinguishable here from a
genome genuinely lacking that gene -- both read as "absent". Common/rare
frequencies and "which taxa carry them" are only meaningful once most/all
families of interest have been enriched and included; with only one
family enriched so far, every genome's architecture is trivially just
that one family's own code, by construction, not a real finding.

This module does NOT infer gene order/operon structure (e.g. the
"C1-Z-C2" kind of label the pipeline may eventually want) -- only
presence/absence + counts per family per genome. Gene order would need
each hit's scaffold and within-scaffold gene index, which IS present in
the underlying NR100 cluster member IDs (<GENOME>-scaffold_<N>_<gene#>,
see pipeline/omdb_metadata.py) but is not threaded through to
unique_targets_with_metadata.tsv today -- a deliberately separate,
later piece of work.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from phaatlas.config_loader import DEFAULT_FAMILY_CONFIG_PATH, load_family_definitions

TAXONOMY_COLUMNS = [
    "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order", "gtdb_family", "gtdb_genus", "gtdb_species",
]
LOCATION_COLUMNS = [
    "sample_id", "study_id", "latitude_degN", "longitude_degE",
    "ecosystem_type", "ecosystem_name", "ecosystem_compartment", "sample_source",
    # Only present in a *_with_metadata_depth.tsv (see pipeline/ncbi_depth.py)
    # -- a plain *_with_metadata.tsv simply has none of these three columns
    # in its header, so row.get(...) below reads them as "" gracefully.
    "depth_raw", "depth_m", "depth_zone",
]


def short_code(family_id: str) -> str:
    """phaA -> A, phaR_regulator -> RReg, phaR_synthase -> RSyn. Derived
    from family_id rather than a hardcoded table so it can never drift
    from config/family_definitions.yaml."""
    name = family_id.removeprefix("pha")
    if "_" in name:
        base, suffix = name.split("_", 1)
        return base + suffix[:3].capitalize()
    return name


def default_family_order(config_path: Path = DEFAULT_FAMILY_CONFIG_PATH) -> list[str]:
    """Canonical family ordering for architecture labels -- config file's
    own order, so labels are deterministic and match this project's
    existing family_id ordering elsewhere (e.g. `status` command output)."""
    return [f.family_id for f in load_family_definitions(config_path)]


def discover_metadata_files(results_dir: Path) -> list[Path]:
    return sorted(results_dir.glob("*_unique_targets_with_metadata.tsv"))


@dataclass
class GenomeRecord:
    genome: str
    taxonomy: dict[str, str] = field(default_factory=dict)
    location: dict[str, str] = field(default_factory=dict)
    family_target_ids: dict[str, set[str]] = field(default_factory=dict)  # family_id -> distinct target_ids hitting this genome

    def family_counts(self) -> dict[str, int]:
        return {fam: len(ids) for fam, ids in self.family_target_ids.items()}


def load_genome_records(metadata_paths: list[Path]) -> dict[str, GenomeRecord]:
    """Reads every given <family>_unique_targets_with_metadata.tsv and
    aggregates per-genome family hit counts + taxonomy/location (taken
    from whichever row is seen first for that genome -- these are
    genome/sample facts, not family-specific, so any row agrees)."""
    genomes: dict[str, GenomeRecord] = {}
    for path in metadata_paths:
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                genome = row.get("genome", "")
                if not genome:
                    continue  # no genome resolved for this target_id -- nothing to attribute
                family_id = row["pha_family"]
                target_id = row["target_id"]

                rec = genomes.get(genome)
                if rec is None:
                    rec = GenomeRecord(
                        genome=genome,
                        taxonomy={c: row.get(c, "") for c in TAXONOMY_COLUMNS},
                        location={c: row.get(c, "") for c in LOCATION_COLUMNS},
                    )
                    genomes[genome] = rec
                rec.family_target_ids.setdefault(family_id, set()).add(target_id)
    return genomes


def architecture_label(counts: dict[str, int], family_order: list[str]) -> str:
    return "".join(short_code(fam) for fam in family_order if counts.get(fam, 0) > 0)


def write_genome_family_matrix(
    genomes: dict[str, GenomeRecord],
    family_order: list[str],
    out_path: Path,
) -> int:
    """One row per genome: taxonomy/location + one count column per
    family in family_order + the collapsed architecture label. Returns
    the number of rows written."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = (
        ["genome"] + TAXONOMY_COLUMNS + LOCATION_COLUMNS
        + [f"n_{fam}" for fam in family_order] + ["n_families_present", "architecture"]
    )
    n = 0
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for genome in sorted(genomes):
            rec = genomes[genome]
            counts = rec.family_counts()
            row = {"genome": genome, **rec.taxonomy, **rec.location}
            for fam in family_order:
                row[f"n_{fam}"] = counts.get(fam, 0)
            row["n_families_present"] = sum(1 for fam in family_order if counts.get(fam, 0) > 0)
            row["architecture"] = architecture_label(counts, family_order)
            writer.writerow(row)
            n += 1
    return n


@dataclass
class ArchitectureStats:
    architecture: str
    families_present: list[str]
    n_genomes: int
    pct_of_genomes: float
    top_genera: list[tuple[str, int]]
    top_species: list[tuple[str, int]]
    n_distinct_studies: int
    top_studies: list[tuple[str, int]]


def summarize_architectures(
    genomes: dict[str, GenomeRecord],
    family_order: list[str],
    top_n: int = 5,
) -> list[ArchitectureStats]:
    """One ArchitectureStats per distinct architecture found, sorted most
    to least common (ties broken alphabetically for determinism)."""
    total = len(genomes)
    by_arch: dict[str, list[GenomeRecord]] = {}
    for rec in genomes.values():
        arch = architecture_label(rec.family_counts(), family_order)
        by_arch.setdefault(arch, []).append(rec)

    results: list[ArchitectureStats] = []
    for arch, recs in by_arch.items():
        genera = Counter(r.taxonomy.get("gtdb_genus", "") for r in recs if r.taxonomy.get("gtdb_genus"))
        species = Counter(r.taxonomy.get("gtdb_species", "") for r in recs if r.taxonomy.get("gtdb_species"))
        studies = Counter(r.location.get("study_id", "") for r in recs if r.location.get("study_id"))
        families_present = [fam for fam in family_order if short_code(fam) in _split_codes(arch, family_order)]
        results.append(ArchitectureStats(
            architecture=arch,
            families_present=families_present,
            n_genomes=len(recs),
            pct_of_genomes=100.0 * len(recs) / total if total else 0.0,
            top_genera=genera.most_common(top_n),
            top_species=species.most_common(top_n),
            n_distinct_studies=len(studies),
            top_studies=studies.most_common(top_n),
        ))
    results.sort(key=lambda s: (-s.n_genomes, s.architecture))
    return results


def _split_codes(arch: str, family_order: list[str]) -> set[str]:
    """Architecture labels are a straight concatenation of variable-length
    short codes (e.g. RReg is 4 chars) -- reconstructing which families
    are present from the label alone means matching against the known
    code set greedily longest-first, rather than assuming a fixed width."""
    codes_by_length = sorted({short_code(f) for f in family_order}, key=len, reverse=True)
    present = set()
    i = 0
    while i < len(arch):
        for code in codes_by_length:
            if arch.startswith(code, i):
                present.add(code)
                i += len(code)
                break
        else:
            raise ValueError(f"Could not decompose architecture label {arch!r} using known family codes")
    return present


def write_architecture_summary(stats: list[ArchitectureStats], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "architecture", "families_present", "n_genomes", "pct_of_genomes",
            "n_distinct_studies", "top_genera", "top_species", "top_studies",
        ])
        for s in stats:
            writer.writerow([
                s.architecture,
                ",".join(s.families_present),
                s.n_genomes,
                f"{s.pct_of_genomes:.2f}",
                s.n_distinct_studies,
                "; ".join(f"{g} ({n})" for g, n in s.top_genera),
                "; ".join(f"{sp} ({n})" for sp, n in s.top_species),
                "; ".join(f"{st} ({n})" for st, n in s.top_studies),
            ])
