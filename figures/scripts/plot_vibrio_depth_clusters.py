"""Depth vs. phaC cluster identity for Vibrio, the genus that motivated
this whole line of inquiry -- "Vibrio C8YNX0 is at all depths." See
_depth_cluster_chart.py for the shared chart-building logic (also used by
plot_redsea_depth_clusters.py).
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
    genus='Vibrio',
    title='Vibrio: Does Its phaC Cluster Change With Depth?',
    out_name='vibrio_depth_clusters',
    assignments=assignments,
    rows=rows,
)
