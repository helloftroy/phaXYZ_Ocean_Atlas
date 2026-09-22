"""ESMFold structure prediction driver -- folds one FASTA (or a slice of
one, for the sbatch array job) and writes one PDB file + a summary row
(mean pLDDT, pTM, sequence length) per input sequence.

Uses HuggingFace transformers' EsmForProteinFolding rather than Meta's
original fair-esm[esmfold] package deliberately: the original package
depends on openfold's custom CUDA attention kernels, which need to be
compiled against the exact local CUDA/PyTorch version and are a common
source of cluster-install failures; the transformers port reimplements
the same model in plain PyTorch with no custom kernel compilation step,
at some inference-speed cost that doesn't matter much at this dataset's
scale (~7,300 sequences total, see build_fold_inputs.py).

Usage (single shot, all sequences in one FASTA):
    python run_esmfold.py --fasta uncertain_cluster_representatives.faa --out-dir esmfold_out/uncertain

Usage (array job, one 1/N_BATCHES slice per task -- see ../cluster/run_esmfold.sbatch):
    python run_esmfold.py --fasta uncertain_cluster_representatives.faa --out-dir esmfold_out/uncertain \\
        --n-batches 20 --batch-index $SLURM_ARRAY_TASK_ID

Resumable: skips any target_id whose PDB file already exists in --out-dir,
so a re-submitted/retried array task picks up where a prior run left off
rather than re-folding everything.
"""
import argparse
import csv
import time
from pathlib import Path


def read_fasta(path):
    cur = None
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if cur:
                    yield cur, ''.join(seq)
                cur, seq = line[1:].strip().split()[0], []
            else:
                seq.append(line.strip())
        if cur:
            yield cur, ''.join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', type=Path, required=True)
    ap.add_argument('--out-dir', type=Path, required=True)
    ap.add_argument('--n-batches', type=int, default=1)
    ap.add_argument('--batch-index', type=int, default=0)
    ap.add_argument('--chunk-size', type=int, default=128,
                     help='ESMFold axial-attention chunk size -- lower uses less GPU memory at '
                          'some speed cost; 128 is a safe default for sequences up to ~1200aa on '
                          'a 40GB A100. Raise if you confirm more headroom, lower if you hit OOM.')
    ap.add_argument('--max-length', type=int, default=1500,
                     help='sequences longer than this are skipped and logged, not folded -- a '
                          'safety cap, not expected to trigger at this project\'s current inputs '
                          '(max observed is 1160aa, see build_fold_inputs.py\'s length check).')
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.out_dir / f'summary_batch{args.batch_index}.tsv'

    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    print('loading ESMFold (facebook/esmfold_v1) ...', flush=True)
    tokenizer = AutoTokenizer.from_pretrained('facebook/esmfold_v1')
    model = EsmForProteinFolding.from_pretrained('facebook/esmfold_v1', low_cpu_mem_usage=True)
    model = model.cuda()
    model.esm = model.esm.half()  # halve the ESM-2 language-model backbone's memory; structure module stays fp32
    model.trunk.set_chunk_size(args.chunk_size)
    model.eval()
    print(f'model loaded, chunk_size={args.chunk_size}', flush=True)

    all_seqs = list(read_fasta(args.fasta))
    my_seqs = all_seqs[args.batch_index::args.n_batches]
    print(f'{len(all_seqs)} total sequences in {args.fasta}, {len(my_seqs)} assigned to batch '
          f'{args.batch_index}/{args.n_batches}', flush=True)

    rows = []
    for i, (name, seq) in enumerate(my_seqs):
        pdb_path = args.out_dir / f'{name}.pdb'
        if pdb_path.exists():
            continue
        if len(seq) > args.max_length:
            print(f'[{i}/{len(my_seqs)}] {name}: SKIPPED, length {len(seq)} > --max-length {args.max_length}', flush=True)
            rows.append({'target_id': name, 'length': len(seq), 'mean_plddt': '', 'ptm': '', 'status': 'skipped_too_long'})
            continue

        t0 = time.time()
        try:
            with torch.no_grad():
                tokenized = tokenizer([seq], return_tensors='pt', add_special_tokens=False)
                tokenized = {k: v.cuda() for k, v in tokenized.items()}
                output = model(**tokenized)
                pdb_str = model.output_to_pdb(output)[0]
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            print(f'[{i}/{len(my_seqs)}] {name}: OOM at length {len(seq)}, skipping', flush=True)
            rows.append({'target_id': name, 'length': len(seq), 'mean_plddt': '', 'ptm': '', 'status': 'oom'})
            continue

        pdb_path.write_text(pdb_str)

        # pLDDT read back from the PDB's own B-factor column (ESMFold's documented
        # convention) rather than indexed directly off the output tensor -- more
        # robust to not being 100% sure of that tensor's exact axis layout without
        # a live model to check against.
        plddts = []
        for line in pdb_str.splitlines():
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                plddts.append(float(line[60:66]))
        mean_plddt = sum(plddts) / len(plddts) if plddts else None
        try:
            ptm = float(output.ptm.item())
        except (AttributeError, TypeError):
            ptm = None

        dt = time.time() - t0
        ptm_str = f'{ptm:.3f}' if ptm is not None else 'n/a'
        plddt_str = f'{mean_plddt:.1f}' if mean_plddt is not None else 'n/a'
        print(f'[{i}/{len(my_seqs)}] {name}: length={len(seq)} mean_plddt={plddt_str} '
              f'ptm={ptm_str} time={dt:.1f}s', flush=True)
        rows.append({'target_id': name, 'length': len(seq),
                     'mean_plddt': f'{mean_plddt:.2f}' if mean_plddt is not None else '',
                     'ptm': f'{ptm:.4f}' if ptm is not None else '', 'status': 'ok'})

    if rows:
        write_header = not summary_path.exists()
        with open(summary_path, 'a', newline='') as out:
            w = csv.DictWriter(out, fieldnames=['target_id', 'length', 'mean_plddt', 'ptm', 'status'], delimiter='\t')
            if write_header:
                w.writeheader()
            w.writerows(rows)
    print(f'batch {args.batch_index}/{args.n_batches} done, {len(rows)} sequences processed this run '
          f'(some may have been skipped as already-done)', flush=True)


if __name__ == '__main__':
    main()
