"""Depth vs. phaC cluster identity for REDSEA-S09-B13 -- same chart as
plot_vibrio_depth_clusters.py, for the second genus pulled into the
depth/cluster comparison (see plot_three_genera_depth_clusters.py). See
_depth_cluster_chart.py for the shared chart-building logic.
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
    genus='REDSEA-S09-B13',
    title='REDSEA-S09-B13: Does Its phaC Cluster Change With Depth?',
    out_name='redsea_depth_clusters',
    assignments=assignments,
    rows=rows,
)
