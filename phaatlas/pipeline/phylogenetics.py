"""Builds a real phylogenetic tree (MAFFT alignment + FastTree, not a
k-mer-distance approximation -- pipeline/sequence_embedding.py already
showed that composition-based clustering doesn't resolve much for a
protein this conserved) of a family's major sequence-cluster
representatives, annotated with depth, ocean basin, taxonomy, and
pathway architecture per tip -- an iTOL-style figure asking whether
clades are ecologically structured, whether lineages track pathway type,
and whether host taxonomy and protein phylogeny agree (a mismatch there
is exactly the kind of hint that points at horizontal gene transfer).

Neither MAFFT nor FastTree are Python packages -- both are shelled out to
(mirroring how pipeline/gopc_search.py and pipeline/sequence_clustering.py
already shell out to mmseqs), so a working `mafft_bin`/`fasttree_bin` is
required. Confirmed live: neither is available via Homebrew on this
machine (broken Xcode Command Line Tools, and brewsci/bio's FastTree tap
needs manual trust) -- both installed cleanly and fast via
`conda create -n pha_phylo -c bioconda -c conda-forge mafft fasttree`
instead, which is precompiled and needs no local compiler at all.
"""

from __future__ import annotations

import csv
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from phaatlas.pipeline.community_ordination import classify_ocean_basin
from phaatlas.pipeline.sequence_embedding import read_fasta


class ExternalToolNotFoundError(RuntimeError):
    pass


def _require_tool(bin_path: str, install_hint: str) -> None:
    import shutil

    if shutil.which(bin_path) is None:
        raise ExternalToolNotFoundError(f"'{bin_path}' not found on PATH. {install_hint}")


def select_representatives(
    ecology_tsv_path: Path,
    top_n: int = 80,
    extra_cluster_ids: list[str] | None = None,
) -> list[str]:
    """Top top_n clusters by genome count (the "major" lineages), plus
    any extra_cluster_ids explicitly requested (e.g. the same "regional"
    clusters already featured in an earlier map figure, for continuity
    across figures) -- deduplicated, original ranking order preserved for
    the top_n, extras appended after."""
    rows = list(csv.DictReader(open(ecology_tsv_path, newline=""), delimiter="\t"))
    rows.sort(key=lambda r: -int(r["n_genomes"]))
    top_ids = [r["cluster_id"] for r in rows[:top_n]]
    seen = set(top_ids)
    result = list(top_ids)
    for cid in extra_cluster_ids or []:
        if cid not in seen:
            result.append(cid)
            seen.add(cid)
    return result


def extract_sequences_subset(fasta_path: Path, wanted_ids: list[str], out_path: Path) -> int:
    """Pulls just the wanted headers' sequences out of an already-local
    FASTA (e.g. the one pipeline/sequence_clustering.py's cluster
    extraction step already produced) -- no mmseqs/network involved,
    just a plain-text filter."""
    records = read_fasta(fasta_path)
    by_id = {h: s for h, s in records}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out_path, "w") as f:
        for wid in wanted_ids:
            seq = by_id.get(wid)
            if seq is None:
                continue
            f.write(f">{wid}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")
            n += 1
    return n


def run_mafft(input_fasta: Path, output_fasta: Path, mafft_bin: str = "mafft", threads: int | None = None) -> None:
    _require_tool(mafft_bin, "Install via: conda create -n pha_phylo -c bioconda -c conda-forge mafft fasttree")
    cmd = [mafft_bin, "--auto"]
    if threads:
        cmd += ["--thread", str(threads)]
    cmd.append(str(input_fasta))
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)


def run_fasttree(aligned_fasta: Path, output_newick: Path, fasttree_bin: str = "FastTree") -> None:
    _require_tool(fasttree_bin, "Install via: conda create -n pha_phylo -c bioconda -c conda-forge mafft fasttree")
    output_newick.parent.mkdir(parents=True, exist_ok=True)
    with open(aligned_fasta) as in_f, open(output_newick, "w") as out_f:
        subprocess.run([fasttree_bin], stdin=in_f, stdout=out_f, check=True)


@dataclass
class ClusterAnnotation:
    cluster_id: str
    n_genomes: int
    dominant_genus: str
    dominant_phylum: str
    median_depth_m: float | None
    dominant_depth_zone: str
    dominant_ocean_basin: str
    dominant_architecture: str
    top_study: str


def annotate_clusters(
    metadata_depth_path: Path,
    cluster_assignments: dict[str, str],
    genome_architecture: dict[str, str],
    wanted_cluster_ids: list[str],
) -> dict[str, ClusterAnnotation]:
    """One pass over metadata_depth_path, gathering per-cluster (only for
    wanted_cluster_ids, to keep this cheap) genus/phylum/depth/basin/
    architecture -- each taken as the modal (most common) value across
    that cluster's distinct host genomes, since a tree tip needs one
    label per gene cluster, not a distribution."""
    wanted = set(wanted_cluster_ids)
    genus_counts: dict[str, Counter] = {c: Counter() for c in wanted}
    phylum_counts: dict[str, Counter] = {c: Counter() for c in wanted}
    depth_zone_counts: dict[str, Counter] = {c: Counter() for c in wanted}
    basin_counts: dict[str, Counter] = {c: Counter() for c in wanted}
    arch_counts: dict[str, Counter] = {c: Counter() for c in wanted}
    study_counts: dict[str, Counter] = {c: Counter() for c in wanted}
    depths: dict[str, list[float]] = {c: [] for c in wanted}
    seen_genomes: dict[str, set[str]] = {c: set() for c in wanted}

    with open(metadata_depth_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cluster_id = cluster_assignments.get(row.get("target_id", ""))
            if cluster_id not in wanted:
                continue
            genome = row.get("genome", "")
            if not genome or genome in seen_genomes[cluster_id]:
                continue
            seen_genomes[cluster_id].add(genome)

            genus = row.get("gtdb_genus", "")
            if genus:
                genus_counts[cluster_id][genus] += 1
            phylum = row.get("gtdb_phylum", "")
            if phylum:
                phylum_counts[cluster_id][phylum] += 1
            depth_zone = row.get("depth_zone", "")
            if depth_zone:
                depth_zone_counts[cluster_id][depth_zone] += 1
            depth_m = row.get("depth_m", "")
            if depth_m:
                depths[cluster_id].append(float(depth_m))
            study = row.get("study_id", "")
            if study:
                study_counts[cluster_id][study] += 1

            lat, lon = row.get("latitude_degN", ""), row.get("longitude_degE", "")
            if lat and lon:
                basin_counts[cluster_id][classify_ocean_basin(float(lat), float(lon))] += 1

            arch = genome_architecture.get(genome, "")
            if arch:
                arch_counts[cluster_id][arch] += 1

    from statistics import median

    def top1(counter: Counter) -> str:
        return counter.most_common(1)[0][0] if counter else ""

    result: dict[str, ClusterAnnotation] = {}
    for cid in wanted_cluster_ids:
        result[cid] = ClusterAnnotation(
            cluster_id=cid,
            n_genomes=len(seen_genomes[cid]),
            dominant_genus=top1(genus_counts[cid]),
            dominant_phylum=top1(phylum_counts[cid]),
            median_depth_m=median(depths[cid]) if depths[cid] else None,
            dominant_depth_zone=top1(depth_zone_counts[cid]),
            dominant_ocean_basin=top1(basin_counts[cid]),
            dominant_architecture=top1(arch_counts[cid]),
            top_study=top1(study_counts[cid]),
        )
    return result


def write_cluster_annotations(annotations: dict[str, ClusterAnnotation], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "cluster_id", "n_genomes", "dominant_genus", "dominant_phylum", "median_depth_m",
            "dominant_depth_zone", "dominant_ocean_basin", "dominant_architecture", "top_study",
        ])
        for a in annotations.values():
            writer.writerow([
                a.cluster_id, a.n_genomes, a.dominant_genus, a.dominant_phylum,
                "" if a.median_depth_m is None else f"{a.median_depth_m:.2f}",
                a.dominant_depth_zone, a.dominant_ocean_basin, a.dominant_architecture, a.top_study,
            ])


@dataclass
class TreeLayoutNode:
    id: int
    name: str  # empty for internal nodes
    x: float  # cumulative branch length from the root
    y: float  # tip display order for tips; mean of immediate children's y for internal nodes
    parent_id: int | None
    support: float | None
    is_tip: bool


def compute_tree_layout(newick_path: Path) -> list[TreeLayoutNode]:
    """Parses a Newick tree (via Biopython -- not a hand-rolled parser,
    since a tree-topology parsing bug would be a serious, hard-to-notice
    error) and computes the standard rectangular phylogram layout every
    major tree viewer (iTOL, FigTree, ETE) uses: x = cumulative branch
    length from the root; y = tip display order (0..n_tips-1, in
    Biopython's own stable left-to-right traversal order) for tips, and
    the mean of immediate children's y for internal nodes -- giving the
    familiar "ladder" look. Returns a flat list with each node's
    parent_id, so a renderer can draw one line segment per node without
    re-walking the tree structure."""
    from Bio import Phylo

    tree = Phylo.read(str(newick_path), "newick")
    root = tree.root

    tip_order = tree.get_terminals()
    y_for_tip = {id(t): i for i, t in enumerate(tip_order)}

    x_of: dict[int, float] = {id(root): 0.0}

    def assign_x(node, parent_x: float) -> None:
        x_of[id(node)] = parent_x + (node.branch_length or 0.0)
        for child in node.clades:
            assign_x(child, x_of[id(node)])

    for child in root.clades:
        assign_x(child, 0.0)

    y_of: dict[int, float] = {}

    def assign_y(node) -> float:
        if not node.clades:
            y_of[id(node)] = float(y_for_tip[id(node)])
        else:
            child_ys = [assign_y(c) for c in node.clades]
            y_of[id(node)] = sum(child_ys) / len(child_ys)
        return y_of[id(node)]

    assign_y(root)

    nodes: list[TreeLayoutNode] = []

    def flatten(node, parent_seq_id: int | None) -> None:
        seq_id = len(nodes)
        nodes.append(TreeLayoutNode(
            id=seq_id, name=node.name or "", x=x_of[id(node)], y=y_of[id(node)],
            parent_id=parent_seq_id, support=node.confidence, is_tip=not node.clades,
        ))
        this_id = seq_id
        for child in node.clades:
            flatten(child, this_id)

    flatten(root, None)
    return nodes


def write_tree_layout(nodes: list[TreeLayoutNode], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["id", "name", "x", "y", "parent_id", "support", "is_tip"])
        for n in nodes:
            writer.writerow([
                n.id, n.name, n.x, n.y,
                "" if n.parent_id is None else n.parent_id,
                "" if n.support is None else n.support,
                n.is_tip,
            ])
