"""Stage the PHA Ocean Atlas data deposit for Zenodo.

Collects the tables, sequences and structures a reader needs to reproduce or
reuse the paper's results, compresses them, and writes a manifest with byte
sizes, row counts and SHA-256 checksums. The staged tree is gitignored -- it is
an upload artefact, rebuildable from this script; only the script, README and
manifest are version-controlled.

Everything here is already derived data. Raw OMDB genomes and the WOA23 grids
are third-party and are cited in the README rather than redistributed.

Usage:
    python zenodo/build_deposit.py            # stage everything missing
    python zenodo/build_deposit.py --manifest # rewrite the manifest only
"""
import argparse
import csv
import gzip
import hashlib
import os
import shutil
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'zenodo'
RESULTS = ROOT / 'PHA_bioprospecting/omdb_search/results'

PHA_FAMILIES = ['phaA', 'phaB', 'phaC', 'phaD', 'phaE', 'phaF', 'phaG',
                'phaI', 'phaJ', 'phaP', 'phaR_regulator', 'phaY', 'phaZ']

# (destination, source, how, description)
#   how: 'gz' compress a single file, 'copy' verbatim, 'tar' a directory,
#        'tar-glob' a glob of files into one archive
ITEMS = [
    # ---- 01 gene annotations -------------------------------------------
    *[(f'01_gene_annotations/{fam}_unique_targets_with_metadata.tsv.gz',
       RESULTS / f'{fam}_unique_targets_with_metadata.tsv', 'gz',
       f'{fam} hits against the OMDB protein catalogue, one row per protein-genome '
       'pair, with GTDB taxonomy, study/sample metadata and coordinates. '
       'Capped at 5 genomes per identical protein (see README, Caveats).')
      for fam in PHA_FAMILIES],
    ('01_gene_annotations/genome_family_matrix.tsv.gz',
     RESULTS / 'genome_family_matrix.tsv', 'gz',
     'Genome x PHA-family presence/copy-number matrix, built from the capped '
     'annotation tables above.'),
    ('01_gene_annotations/omdb_all_genome_taxonomy.tsv.gz',
     RESULTS / 'omdb_all_genome_taxonomy.tsv', 'gz',
     'GTDB taxonomy for every OMDB genome, including the phaC-negative ones, '
     'used as the denominator for prevalence rates.'),

    # ---- 02 clustering -------------------------------------------------
    *[(f'02_clustering/phaC_cluster{t}_cluster.tsv.gz',
       RESULTS / f'phaC_cluster{t}_cluster.tsv', 'gz',
       f'MMseqs2 cluster assignment of all phaC proteins at {float(t)*100:.0f}% identity '
       '(two columns: representative, member).')
      for t in ('0.5', '0.6', '0.7', '0.9')],
    ('02_clustering/phaC_all_genomes_from_nr100_clusters.tsv.gz',
     ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', 'gz',
     'Uncapped protein-to-genome membership for phaC. This is the authoritative '
     'source for genome counts (38,744 phaC-positive genomes); prefer it over the '
     'capped annotation tables whenever counting genomes.'),
    ('02_clustering/phaC_all_representatives.faa.gz',
     RESULTS / 'phaC_all_representatives.faa', 'gz',
     'Amino-acid sequences of the phaC cluster representatives.'),

    # ---- 03 catalytic site ---------------------------------------------
    ('03_catalytic_site/phac_catalytic_triad.tsv.gz',
     ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', 'gz',
     'Catalytic-triad scoring for every phaC protein, by alignment against a '
     '625-column phaC HMM: the residues at the Cys/Asp/His columns and a '
     'triad_complete flag.'),
    ('03_catalytic_site/phac_synthase_class.tsv.gz',
     ROOT / 'catalytic_domain/phac_synthase_class.tsv', 'gz',
     'Class I/II/III HMM scores per phaC protein, with the best-scoring class and '
     'whether it clears that model\'s trusted cutoff.'),
    ('03_catalytic_site/phac_sequence_evidence.tsv.gz',
     ROOT / 'catalytic_domain/phac_sequence_evidence.tsv', 'gz',
     'Per-protein sequence evidence tier: length, best reference, identity, '
     'coverage and query-to-reference length ratio.'),

    # ---- 04 structure --------------------------------------------------
    ('04_structure/structural_evidence_best_hit.tsv.gz',
     ROOT / 'structural_evidence_best_hit.tsv', 'gz',
     'Best Foldseek hit of each ESMFold model against the folded reference set, '
     'with TM-scores (query/target/alignment), coverages and structural identity.'),
    ('04_structure/novel_phac_structural_verdicts.tsv',
     ROOT / 'structure_prediction/novel_phac_structural_verdicts.tsv', 'copy',
     'Fold verdict for each of the 255 novel phaC candidates. The active-site '
     'columns are intentionally empty: the per-candidate PDBs were not retained, '
     'so triad geometry was never measured for this set.'),
    ('04_structure/phac_substrate_pocket.tsv',
     ROOT / 'structure_prediction/phac_substrate_pocket.tsv', 'copy',
     'Substrate-channel geometry (LIGSITE-style pocket volume and hydropathy by '
     'distance shell) for the Desulfoluna phaC copies.'),
    ('04_structure/phac_pocket_vs_experiment.tsv',
     ROOT / 'structure_prediction/phac_pocket_vs_experiment.tsv', 'copy',
     'The same pocket measurement applied to PhaC crystal structures of known '
     'chain-length preference, used to calibrate the predictions.'),
    ('04_structure/esmfold_models.tar.gz', ROOT / 'PDB_files', 'tar',
     'ESMFold predicted structures (2,108 PDB files) for phaC proteins and '
     'references. Models, not experimental structures.'),

    # ---- 05 thermostability --------------------------------------------
    ('05_thermostability/phaC_temstapro_predictions_with_metadata.tsv.gz',
     ROOT / 'temstapro/phaC_temstapro_predictions_with_metadata.tsv', 'gz',
     'TemStaPro thermostability predictions per phaC protein at eight temperature '
     'thresholds (40-80 C), joined to genome metadata.'),

    # ---- 06 novel candidates -------------------------------------------
    ('06_novel_candidates/novel_phac_new_candidates.tsv',
     ROOT / 'phac_recovery/novel_phac_new_candidates.tsv', 'copy',
     'The 255 phaC genes recovered from genomes that the primary annotation '
     'called phaC-negative despite carrying more than five other pha genes.'),
    ('06_novel_candidates/novel_phac_new_candidates.faa',
     ROOT / 'phac_recovery/novel_phac_new_candidates.faa', 'copy',
     'Amino-acid sequences of those 255 candidates.'),
    ('06_novel_candidates/novel_phac_candidate_verdicts.tsv.gz',
     ROOT / 'phac_recovery/novel_phac_candidate_verdicts.tsv', 'gz',
     'All candidates from that search, each labelled as a re-find of a '
     'QC-removed protein, a re-find of a QC-passing protein, or genuinely new.'),
    ('06_novel_candidates/novel_phac_candidate_triads.tsv.gz',
     ROOT / 'phac_recovery/novel_phac_candidate_triads.tsv', 'gz',
     'Triad scoring for every candidate gene in that search, not only the 255 kept.'),
    ('06_novel_candidates/qualifying_genomes.txt',
     ROOT / 'phac_recovery/qualifying_genomes.txt', 'copy',
     'The 4,170 phaC-negative genomes with more than five other pha genes that '
     'the recovery search targeted.'),

    # ---- 07 environment ------------------------------------------------
    ('07_environment/omdb_all_genomes_with_locations.tsv.gz',
     ROOT / 'omdb_all_genomes_with_locations.tsv', 'gz',
     'Sample provenance and coordinates for every OMDB genome.'),
    ('07_environment/phaC_genomes_woa23_annual_temperature.tsv.gz',
     RESULTS / 'phaC_genomes_woa23_annual_temperature.tsv', 'gz',
     'phaC-positive genomes joined to World Ocean Atlas 2023 annual mean '
     'temperature at their sampling coordinate and depth.'),
    ('07_environment/phaC_genomes_woa23_salinity_oxygen.tsv.gz',
     RESULTS / 'phaC_genomes_woa23_salinity_oxygen.tsv', 'gz',
     'The same join for WOA23 salinity and dissolved oxygen.'),

    # ---- 08 figure data ------------------------------------------------
    ('08_figure_data/figure_source_tables.tar.gz', ROOT / 'figures/*.tsv', 'tar-glob',
     'The companion table behind each published figure, one per figure, so every '
     'plotted number can be read without rerunning the pipeline.'),
]


def sha256(path, blocks=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(blocks), b''):
            h.update(chunk)
    return h.hexdigest()


def rows_in(path):
    """Data rows (excluding one header) in a text table, else blank."""
    name = path.name
    opener = gzip.open if name.endswith('.gz') else open
    inner = name[:-3] if name.endswith('.gz') else name
    if not inner.endswith(('.tsv', '.txt', '.csv')):
        return ''
    try:
        with opener(path, 'rt', errors='replace') as f:
            n = sum(1 for _ in f)
    except OSError:
        return ''
    return str(max(n - 1, 0)) if inner.endswith(('.tsv', '.csv')) else str(n)


def stage(dest, source, how):
    """Produce OUT/dest from source. Returns False if the source is missing."""
    target = OUT / dest
    if target.exists():
        return True
    sources = sorted(ROOT.glob(str(source.relative_to(ROOT)))) if how == 'tar-glob' else [source]
    if not sources or (how != 'tar-glob' and not source.exists()):
        return False
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(target.suffix + '.partial')
    if how == 'copy':
        shutil.copy2(source, tmp)
    elif how == 'gz':
        # mtime=0 so rebuilding the deposit gives byte-identical archives
        with open(source, 'rb') as src, gzip.GzipFile(tmp, 'wb', compresslevel=9, mtime=0) as dst:
            shutil.copyfileobj(src, dst, length=1 << 22)
    elif how in ('tar', 'tar-glob'):
        with tarfile.open(tmp, 'w:gz', compresslevel=9) as tar:
            if how == 'tar':
                for child in sorted(source.iterdir()):
                    tar.add(child, arcname=f'{source.name}/{child.name}')
            else:
                for child in sources:
                    tar.add(child, arcname=f'figure_source_tables/{child.name}')
    else:
        raise ValueError(how)
    os.replace(tmp, target)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', action='store_true',
                        help='rewrite MANIFEST.tsv without staging anything new')
    args = parser.parse_args()

    staged, missing = [], []
    for dest, source, how, description in ITEMS:
        if args.manifest:
            ok = (OUT / dest).exists()
        else:
            ok = stage(dest, source, how)
            print(('  ok   ' if ok else '  MISS ') + dest, flush=True)
        (staged if ok else missing).append((dest, description))

    rows = []
    for dest, description in staged:
        path = OUT / dest
        rows.append((dest, path.stat().st_size, rows_in(path), sha256(path), description))

    with open(OUT / 'MANIFEST.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['file', 'bytes', 'data_rows', 'sha256', 'description'])
        w.writerows(rows)

    total = sum(r[1] for r in rows)
    print(f'\n{len(rows)} files staged, {total / 2**30:.2f} GiB total')
    print(f'manifest: {OUT / "MANIFEST.tsv"}')
    for dest, _ in missing:
        print(f'  missing source, not staged: {dest}')


if __name__ == '__main__':
    main()
