"""Depth vs. phaC cluster identity for Sulfitobacter -- replaces
REDSEA-S09-B13 and DTSX01 as this project's single-organism depth/cluster
case study (both collapsed to near-nothing after the 2026-09-22 phaC
reference-query fix; see PHA_CLEAN_RESULTS.md section 5.4). Sulfitobacter
is also this project's headline "genuinely global cluster" pick (section
5.2) -- picked again here on its own separate merits: 236 phaC-positive
genomes with resolved depth, 17 distinct 70%-identity clusters, spanning
0-4200m, the best combination of genome count, cluster diversity, and
depth range of any genus checked. See _depth_cluster_chart.py for the
shared chart-building logic (also used by plot_vibrio_depth_clusters.py).
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _depth_cluster_chart import build_depth_cluster_chart, FA
from phaatlas.pipeline import sequence_clustering as sc

assignments = sc.load_cluster_assignments(FA / 'phaC_cluster0.7_cluster.tsv')
rows = list(csv.DictReader(open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline=''), delimiter='\t'))

build_depth_cluster_chart(
    genus='Sulfitobacter',
    title='Sulfitobacter: Does Its phaC Cluster Change With Depth?',
    out_name='sulfitobacter_depth_clusters',
    assignments=assignments,
    rows=rows,
)
