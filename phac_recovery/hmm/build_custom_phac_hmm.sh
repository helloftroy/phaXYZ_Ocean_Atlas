#!/bin/bash
# Rebuilds phaC_custom.hmm from PHA_reference/exports/pha_reference_nr95.faa
# (this repo's own curated PhaC reference diversity, 947 sequences from
# BRENDA/UniProt at 95% clustering -- the same seed set the original
# mmseqs2 OMDB/GOPC search used). Complements Pfam PF07167 (fetched
# separately, see below) which only covers PhaC's 173-residue N-terminal
# domain from just 20 Pfam seed sequences -- this custom profile is
# full-length (625-column model) and reflects our own broader diversity,
# so it should be more sensitive to divergent-but-full-length homologs
# that PF07167's narrow N-terminal signature might miss.
#
# Requires mafft + hmmer (both in the `pha_phylo` conda env: `conda create
# -n pha_phylo -c bioconda -c conda-forge mafft fasttree hmmer`).
#
# Length filter (200-900aa) excludes a handful of pathological outliers
# in the reference set (min was 36aa -- a truncated UniProt fragment; max
# was 5178aa -- almost certainly a fusion/misannotated entry) that would
# otherwise dominate/distort the alignment -- confirmed live: 822 of 947
# sequences retained.
#
# Usage: ./build_custom_phac_hmm.sh
set -euo pipefail
cd "$(dirname "$0")"

REPO_ROOT="$(cd ../.. && pwd)"
REF_FAA="${REPO_ROOT}/PHA_reference/exports/pha_reference_nr95.faa"

python3 - "$REF_FAA" <<'PYEOF'
import sys
def read_fasta(path):
    cur=None; seq=[]
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if cur: yield cur, ''.join(seq)
                cur=line[1:].strip(); seq=[]
            else: seq.append(line.strip())
        if cur: yield cur, ''.join(seq)

kept, excluded = 0, 0
with open('phac_ref_filtered.faa', 'w') as out:
    for h, s in read_fasta(sys.argv[1]):
        if 200 <= len(s) <= 900:
            out.write(f'>{h}\n')
            for i in range(0, len(s), 60):
                out.write(s[i:i+60] + '\n')
            kept += 1
        else:
            excluded += 1
print(f'kept {kept}, excluded {excluded} (outside 200-900aa)')
PYEOF

mafft --auto --thread 8 phac_ref_filtered.faa > phac_ref_aligned.faa 2> mafft.log
hmmbuild --amino phaC_custom.hmm phac_ref_aligned.faa
hmmstat phaC_custom.hmm

echo
echo "Also fetch Pfam PF07167 (PhaC N-terminal domain) directly:"
echo '  curl -sL "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF07167/?annotation=hmm" | gunzip > PF07167.hmm'
