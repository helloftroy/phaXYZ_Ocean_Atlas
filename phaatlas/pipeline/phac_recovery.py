"""Recovering a divergent/missed phaC in genomes that carry >5 OTHER PHA
pathway genes but no phaC hit -- see figures/scripts/plot_phac_negative_taxonomy.py
and figures/scripts/plot_phac_negative_treemap.py for how this population
(13,522 genomes as of this writing) was identified and characterized
taxonomically.

Why a genomic-neighborhood search, not a whole-genome one: PHA pathway
genes are typically operonic/co-located, so restricting the search to the
flanking ORFs around a genome's OWN other-pha-gene hits (rather than every
ORF in the genome) both keeps the search tractable across 13k+ genomes and
raises confidence that any hit found is genuinely part of THIS pathway,
not an unrelated hydrolase picked up genome-wide.

Why exact sequence matching to locate anchors, not the NR100 cluster.tsv:
target_id in our metadata is an NR100 CLUSTER id (a 100%-identity
dereplication), not a per-genome gene id -- translating "genome G carries
cluster T" into "which of G's own genes that is" would normally need
OMDBv2.0_AA_G_NR100.cluster.tsv.gz's full member list (4.2GB, not
downloaded locally, mapping representative_id -> member_id pairs for
every protein in every OMDB genome, not just PHA-relevant ones). Since
NR100 dereplication is exact-identity, this is unnecessary: a genome that
truly carries cluster T has a byte-identical copy of T's own sequence
somewhere in its own gene-calling output, so a straight dict lookup by
sequence string resolves it without touching the cluster file at all.
T's own sequence is pulled from the mmseqs target_db already built on the
cluster for the original search (see cluster/run_phac_recovery_extract_targets.sbatch),
the same way cluster/run_omdb_extract_cluster_sequences.sbatch already
pulls sequences for a whole family's target_ids.

Genome ORF order: OMDB's per-genome `<genome>.genes.faa.gz` (prodigal
calls) headers are `<genome>-scaffold_<N>_<gene#>`, already listed in
scaffold-then-gene-index order in the file itself -- no GFF parsing
needed to determine gene order/neighbors, just the FASTA.
"""
from __future__ import annotations

import csv
import gzip
import subprocess
import time
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq)))
    return records


def identify_qualifying_genomes(genome_family_matrix_path: Path, min_other_genes: int = 5) -> set[str]:
    """Genomes with zero phaC hits but more than min_other_genes other PHA
    families present -- the exact criterion used throughout this analysis."""
    qualifying = set()
    with open(genome_family_matrix_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if int(row["n_phaC"]) == 0 and int(row["n_families_present"]) > min_other_genes:
                qualifying.add(row["genome"])
    return qualifying


@dataclass
class GenomeHit:
    genome: str
    family: str
    target_id: str


def collect_other_pha_hits(qualifying_genomes: set[str], family_metadata_paths: dict[str, Path]) -> list[GenomeHit]:
    """Scans each non-phaC family's <family>_unique_targets_with_metadata.tsv
    for rows belonging to a qualifying genome, recording (genome, family,
    target_id) -- these target_ids are the "other pha gene" anchors whose
    genomic neighborhood we want to search for a missed phaC. family_id
    'phaC' should not be included in family_metadata_paths (nothing to
    anchor around a synthase these genomes don't have)."""
    assert "phaC" not in family_metadata_paths, "phaC is the gene we're looking for, not an anchor"
    hits: list[GenomeHit] = []
    for family, path in family_metadata_paths.items():
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                genome = row.get("genome", "")
                if genome in qualifying_genomes:
                    hits.append(GenomeHit(genome=genome, family=family, target_id=row["target_id"]))
    return hits


def write_genome_download_manifest(qualifying_genomes: set[str], omdb_catalog_path: Path, out_path: Path) -> int:
    """genome -> its GENES_AA_FILE URL, from OMDBv2.0_data.tsv (the small
    catalog listing every genome's per-genome file URLs -- not the genome
    sequences themselves). Returns the number of genomes matched."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(omdb_catalog_path, newline="") as f_in, open(out_path, "w", newline="") as f_out:
        reader = csv.DictReader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["genome", "genes_aa_url"])
        for row in reader:
            if row["GENOME"] in qualifying_genomes:
                writer.writerow([row["GENOME"], row["GENES_AA_FILE"]])
                n += 1
    return n


def parse_prodigal_faa(records: list[tuple[str, str]]) -> dict[str, list[tuple[str, str]]]:
    """Groups a genome's own gene-call FASTA records by scaffold, preserving
    file order within each scaffold (== gene index order, since prodigal
    numbers genes sequentially along each scaffold). Returns
    {scaffold_id: [(gene_id, seq), ...]}."""
    by_scaffold: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for gene_id, seq in records:
        # gene_id = "<genome>-scaffold_<N>_<gene#>" -- scaffold_id is everything before the final "_<gene#>"
        scaffold_id = gene_id.rsplit("_", 1)[0]
        by_scaffold[scaffold_id].append((gene_id, seq))
    return dict(by_scaffold)


def find_anchor_gene_ids(genome_genes: dict[str, list[tuple[str, str]]], anchor_sequences: set[str]) -> list[str]:
    """Which of this genome's own genes are an exact-sequence match to one
    of the known anchor (other-pha-family target) sequences."""
    found = []
    for gene_list in genome_genes.values():
        for gene_id, seq in gene_list:
            if seq in anchor_sequences:
                found.append(gene_id)
    return found


def extract_neighborhood_windows(
    genome_genes: dict[str, list[tuple[str, str]]],
    anchor_gene_ids: list[str],
    window: int = 10,
) -> dict[str, str]:
    """For each anchor gene, takes the +/-window flanking genes on its own
    scaffold (naturally the whole scaffold if it has <= 2*window+1 genes
    total). Returns the union across all anchors as {gene_id: seq},
    deduplicated (a scaffold with several anchors close together
    contributes one merged neighborhood, not overlapping copies)."""
    gene_to_scaffold_and_index: dict[str, tuple[str, int]] = {}
    for scaffold_id, gene_list in genome_genes.items():
        for i, (gene_id, _seq) in enumerate(gene_list):
            gene_to_scaffold_and_index[gene_id] = (scaffold_id, i)

    result: dict[str, str] = {}
    for anchor_id in anchor_gene_ids:
        if anchor_id not in gene_to_scaffold_and_index:
            continue
        scaffold_id, idx = gene_to_scaffold_and_index[anchor_id]
        gene_list = genome_genes[scaffold_id]
        lo, hi = max(0, idx - window), min(len(gene_list), idx + window + 1)
        for gene_id, seq in gene_list[lo:hi]:
            result[gene_id] = seq
    return result


def write_neighborhood_fasta(neighborhoods: dict[str, str], out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for gene_id, seq in neighborhoods.items():
            f.write(f">{gene_id}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")
    return len(neighborhoods)


def default_fetch_gzipped_fasta(url: str, timeout: int = 60) -> bytes:
    """Real network fetch of a .gz FASTA (OMDB's GENES_AA_FILE URLs),
    returning the decompressed bytes. Kept as a standalone function
    (rather than inlined) so run_neighborhood_extraction_batch can take a
    fake in tests instead of hitting the network."""
    req = urllib.request.Request(url, headers={"User-Agent": "phac_recovery/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        compressed = resp.read()
    return gzip.decompress(compressed)


@dataclass
class BatchResult:
    genome: str
    status: str  # "ok", "download_failed", "no_anchors_found"
    n_neighborhood_genes: int = 0


def run_neighborhood_extraction_batch(
    genome_urls: dict[str, str],
    genome_anchor_sequences: dict[str, set[str]],
    out_path: Path,
    window: int = 10,
    fetch: Callable[[str], bytes] = default_fetch_gzipped_fasta,
    retries: int = 3,
    retry_delay_s: float = 2.0,
    log: Callable[[str], None] = print,
) -> list[BatchResult]:
    """The heavy cluster-side step: for every genome in genome_urls,
    downloads its own gene calls, locates the anchors it's known to carry
    (genome_anchor_sequences[genome], the actual sequences of that
    genome's other-pha-family hits -- see collect_other_pha_hits +
    a target_id->sequence extraction from the cluster's target_db),
    extracts the neighborhood windows, and appends them to one combined
    output FASTA -- so 13k+ raw per-genome downloads don't all need to be
    kept on disk at once, only the much smaller neighborhood subset.
    Per-genome failures (network hiccups, unresolved anchors) are
    recorded and skipped rather than aborting the whole batch, since with
    13k+ individual downloads some failures are expected."""
    results: list[BatchResult] = []
    seen_gene_ids: set[str] = set()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as out_f:
        for i, (genome, url) in enumerate(genome_urls.items()):
            anchors = genome_anchor_sequences.get(genome, set())
            if not anchors:
                results.append(BatchResult(genome=genome, status="no_anchors_found"))
                continue

            raw = None
            for attempt in range(retries):
                try:
                    raw = fetch(url)
                    break
                except Exception as e:  # noqa: BLE001 -- deliberately broad: any network hiccup should retry, not abort the batch
                    log(f"[{i+1}/{len(genome_urls)}] {genome}: fetch attempt {attempt+1} failed ({e})")
                    time.sleep(retry_delay_s)
            if raw is None:
                results.append(BatchResult(genome=genome, status="download_failed"))
                continue

            records = read_fasta_text(raw.decode("utf-8", errors="replace"))
            genome_genes = parse_prodigal_faa(records)
            anchor_gene_ids = find_anchor_gene_ids(genome_genes, anchors)
            if not anchor_gene_ids:
                results.append(BatchResult(genome=genome, status="no_anchors_found"))
                continue

            windows = extract_neighborhood_windows(genome_genes, anchor_gene_ids, window=window)
            n_written = 0
            for gene_id, seq in windows.items():
                if gene_id in seen_gene_ids:
                    continue
                seen_gene_ids.add(gene_id)
                out_f.write(f">{gene_id}\n")
                for j in range(0, len(seq), 60):
                    out_f.write(seq[j : j + 60] + "\n")
                n_written += 1
            results.append(BatchResult(genome=genome, status="ok", n_neighborhood_genes=n_written))

            if (i + 1) % 500 == 0:
                log(f"[{i+1}/{len(genome_urls)}] genomes processed")

    return results


def read_fasta_text(text: str) -> list[tuple[str, str]]:
    """Same as read_fasta but from an in-memory string, for FASTA bytes
    fetched over the network rather than read from a local file."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def write_batch_report(results: list[BatchResult], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["genome", "status", "n_neighborhood_genes"])
        for r in results:
            writer.writerow([r.genome, r.status, r.n_neighborhood_genes])


def load_anchor_sequences_by_genome(other_pha_hits_path: Path, target_id_sequences: dict[str, str]) -> dict[str, set[str]]:
    """Joins other_pha_hits.tsv (genome, family, target_id) with a
    target_id->sequence lookup (from the cluster target_db extraction) to
    build {genome: {anchor_seq, ...}} -- exactly what
    run_neighborhood_extraction_batch needs as genome_anchor_sequences."""
    by_genome: dict[str, set[str]] = defaultdict(set)
    with open(other_pha_hits_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq = target_id_sequences.get(row["target_id"])
            if seq:
                by_genome[row["genome"]].add(seq)
    return dict(by_genome)


@dataclass
class HmmHit:
    gene_id: str
    profile: str
    score: float
    evalue: float


def run_hmmsearch(hmm_path: Path, fasta_path: Path, out_domtbl_path: Path, dom_e: float = 10.0, hmmsearch_bin: str = "hmmsearch") -> None:
    """Relaxed-threshold search (default domE=10, far more permissive than
    Pfam's own curated cutoffs) -- deliberately loose since the point is
    catching divergent phaC our stricter mmseqs2 seed search already
    missed; --domtblout gives one clean row per domain hit for parsing,
    rather than needing to scrape the human-readable report."""
    out_domtbl_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [hmmsearch_bin, "--domtblout", str(out_domtbl_path), "--domE", str(dom_e), str(hmm_path), str(fasta_path)],
        check=True, stdout=subprocess.DEVNULL,
    )


def parse_hmmsearch_domtbl(domtbl_path: Path, profile_name: str) -> list[HmmHit]:
    hits = []
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            gene_id = fields[0]
            evalue = float(fields[12])
            score = float(fields[13])
            hits.append(HmmHit(gene_id=gene_id, profile=profile_name, score=score, evalue=evalue))
    return hits


def write_hmm_hits_report(hits: Iterable[HmmHit], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene_id", "genome", "profile", "score", "evalue"])
        for h in hits:
            genome = h.gene_id.split("-scaffold_")[0]
            writer.writerow([h.gene_id, genome, h.profile, h.score, h.evalue])
