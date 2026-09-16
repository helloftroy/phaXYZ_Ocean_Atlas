"""Builds a "protein clade x depth zone" abundance matrix from one
family's depth-enriched search results (pipeline/ncbi_depth.py's output),
for exactly the kind of question this was built for: is one PhaC clade
almost exclusively deep-ocean while another dominates the photic zone?

"Clade" here means best_query -- the single closest NR95 reference
protein each OMDB hit matched (already present on every row from the
search stage, see pipeline/gopc_search.py's aggregate_hits_file). It is
a coarse but immediately available proxy for a hit's likely
subfamily/class -- not a phylogenetic tree. Genuine clade assignment
(e.g. actual PhaC class I-IV) would need a real tree/alignment of the
hits themselves, deliberately out of scope here.

Depth zone comes from pipeline/ncbi_depth.py's depth_zone column; rows
with no resolvable depth (no BioSample record, or a record with no
reported depth) are kept in an explicit "(unknown depth)" column rather
than silently dropped, since that coverage gap is itself worth seeing
next to the real bins.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from phaatlas.pipeline.ncbi_depth import DEPTH_BINS

UNKNOWN_DEPTH_LABEL = "(unknown depth)"
DEPTH_ZONE_ORDER = [label for label, _, _ in DEPTH_BINS] + [UNKNOWN_DEPTH_LABEL]


def build_clade_depth_matrix(
    metadata_depth_path: Path,
    clade_column: str = "best_query",
) -> tuple[list[str], dict[str, dict[str, int]], dict[str, dict[str, int]]]:
    """Reads a <family>_unique_targets_with_metadata_depth.tsv. Returns
    (clades_ranked_by_total_abundance_desc, genome_counts, target_id_counts)
    where both count dicts are {clade: {depth_zone_label: count}} --
    genome_counts counts DISTINCT genomes per cell (how many organisms
    carry this clade in this depth zone), target_id_counts counts DISTINCT
    target_ids (how many distinct NR100 hits, i.e. exact protein
    sequences). Missing depth reads as UNKNOWN_DEPTH_LABEL, never dropped.
    """
    genomes_seen: dict[str, dict[str, set[str]]] = {}
    target_ids_seen: dict[str, dict[str, set[str]]] = {}

    with open(metadata_depth_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            clade = row.get(clade_column, "") or "(no clade)"
            zone = row.get("depth_zone", "") or UNKNOWN_DEPTH_LABEL
            genome = row.get("genome", "")
            target_id = row.get("target_id", "")

            genomes_seen.setdefault(clade, {}).setdefault(zone, set())
            target_ids_seen.setdefault(clade, {}).setdefault(zone, set())
            if genome:
                genomes_seen[clade][zone].add(genome)
            if target_id:
                target_ids_seen[clade][zone].add(target_id)

    genome_counts = {c: {z: len(s) for z, s in zones.items()} for c, zones in genomes_seen.items()}
    target_id_counts = {c: {z: len(s) for z, s in zones.items()} for c, zones in target_ids_seen.items()}

    totals = Counter({c: sum(zones.values()) for c, zones in genome_counts.items()})
    clades_ranked = [c for c, _ in totals.most_common()]
    return clades_ranked, genome_counts, target_id_counts


def write_long_format(
    clades_ranked: list[str],
    genome_counts: dict[str, dict[str, int]],
    target_id_counts: dict[str, dict[str, int]],
    out_path: Path,
) -> int:
    """One row per (clade, depth_zone) that has any data at all --
    convenient for feeding straight into pandas/seaborn without a pivot
    step. Zones always appear in DEPTH_ZONE_ORDER, not alphabetically."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["clade", "depth_zone", "n_distinct_genomes", "n_distinct_target_ids"])
        for clade in clades_ranked:
            for zone in DEPTH_ZONE_ORDER:
                n_g = genome_counts.get(clade, {}).get(zone, 0)
                n_t = target_id_counts.get(clade, {}).get(zone, 0)
                if n_g == 0 and n_t == 0:
                    continue
                writer.writerow([clade, zone, n_g, n_t])
                n += 1
    return n


def write_wide_matrix(
    clades_ranked: list[str],
    counts: dict[str, dict[str, int]],
    out_path: Path,
) -> None:
    """clade rows x depth-zone columns, ready to paste straight into a
    heatmap tool (values = whichever count dict is passed in -- genome
    counts or target_id counts)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["clade"] + DEPTH_ZONE_ORDER)
        for clade in clades_ranked:
            writer.writerow([clade] + [counts.get(clade, {}).get(zone, 0) for zone in DEPTH_ZONE_ORDER])
