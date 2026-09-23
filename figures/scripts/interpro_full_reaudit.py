"""Full InterPro-domain re-audit of every PHA reference query (not just
the top-10-per-family sample from the earlier pass). Two-directional
check per reference: does it carry a domain signature that's a confirmed
red flag (catches a deceptive name passing the keyword filter, as found
repeatedly: Q8KRE8/Q1JR79 AraC regulators, A0A0C9NDJ2/A0AAI8EP32 antiporter
subunit D, A0A9W5B2L2/A0ABP2BCP7 antiporter subunit F, A0A0F3KLG9 PhaR
domain mislabeled as phaZ, A0A1R4IZE9 ABC-transporter mislabeled as phaA)?
Or does it carry that family's own confirmed-specific marker domain
(resolves a keyword-UNCLEAR reference the text-only pass couldn't place)?

Usage: python interpro_full_reaudit.py <family> <audit_tsv> [phaC]
Prints a revised classification; does not overwrite the audit TSVs
in place (caller decides what to do with the output).
"""
import csv, json, sys
from pathlib import Path
from collections import Counter

CACHE = Path('/tmp/interpro_cache')

def load_domains(acc):
    p = CACHE / f'{acc}.json'
    if not p.exists() or p.stat().st_size == 0:
        return set()
    try:
        d = json.loads(p.read_text())
    except Exception:
        return set()
    return {r['metadata']['accession'] for r in d.get('results', [])}

# universal red-flag domains -- unambiguous regardless of family (context-dependent
# ones like ABC-transporter periplasmic-binding, which is a legitimate phaC FUSION
# case, and IPR007208 Na+/H+ subunit F-like, which overlaps phaF's own Pfam, are
# deliberately excluded here and handled per-family below instead)
BAD_DOMAINS_UNIVERSAL = {
    'IPR018060', 'IPR032687',                       # AraC-type HTH regulator
    'IPR000163', 'IPR001107', 'IPR036013',           # Prohibitin / Band_7
    'IPR050091', 'IPR014030', 'IPR014031',           # PKS/NRPS biosynthesis enzyme
    'IPR009014', 'IPR033248',                        # Transketolase C-term / phenylpyruvate DH-type
    'IPR052703', 'IPR011882', 'IPR007814', 'IPR012347',  # phenylacetate-pathway (paa operon)
    'IPR002758',                                     # Na+/H+ antiporter subunit E
    'IPR005133',                                     # Na+/H+ antiporter subunit G
    'IPR001750', 'IPR050586',                        # Mrp/CPA3 antiporter TM domain / subunit D
    'IPR018170', 'IPR020471', 'IPR036812',           # Aldo-keto reductase (non-PHA)
    'IPR003560',                                     # 2,3-dihydro-dihydroxybenzoate DH (entA/phbA collision)
    'IPR004561', 'IPR005801',                        # isochorismate synthase
    'IPR001017', 'IPR050771',                        # 2-oxoacid dehydrogenase E1 (BCKDH-type)
    'IPR001734', 'IPR038377', 'IPR050277',           # Na+/solute symporter (ActP-like)
    'IPR001647', 'IPR015893', 'IPR025722', 'IPR050109',  # TetR-family regulator
}
# phaF-specific: real phaF/phaI legitimately share IPR007208 (Pfam PF04066 = "MrpF /
# PhaF"), so only flag it there when paired with an explicit antiporter-only NCBIFam
# hit (NF009245 "cation:proton antiporter") -- otherwise leave as an open case.
PHAF_ANTIPORTER_TELL = {'NF009245', 'NF004812', 'PIRSF028784'}

# confirmed family-specific marker domains -- presence is strong positive evidence
#
# phaC was MISSING from this dict for the entire earlier "full" InterPro audit --
# every other family below got a genuine positive-marker check (promote-from-
# UNCLEAR AND, via the dedicated re-verification pass described in
# PHA_CLEAN_RESULTS.md, demote-from-ON_TARGET-if-marker-absent); phaC only ever
# got the negative/BAD_DOMAINS_UNIVERSAL blacklist check below, which a protein
# that is simply the WRONG family (e.g. a misannotated acyl-CoA synthetase, or a
# polysaccharide-biosynthesis protein) sails past easily, since neither of those
# domain families was ever on the blacklist. Found via a full 1,867-query direct
# UniProt/InterPro re-verification (2026-09-22) after 152 previously-ON_TARGET
# phaC references turned out to lack real PHA-synthase domain evidence --
# collectively responsible for a 54% drop in the phaC-positive genome count once
# corrected. Do not let this happen to another family: every family MUST have an
# entry here, and reaudit_family() must demote away from ON_TARGET when a
# family's own marker is absent, not just promote toward it.
GOOD_DOMAINS = {
    'phaA': {'IPR002155', 'IPR016039', 'IPR020610', 'IPR020613', 'IPR020615', 'IPR020616', 'IPR020617'},
    'phaB': {'IPR011283'},
    'phaC': {'IPR051321'},
    'phaE': {'IPR010123', 'PF09712'},
    'phaF': {'TIGR01837'},
    'phaI': {'IPR008769', 'TIGR01837', 'PF05597'},
    'phaP': {'IPR018968', 'PF09361', 'TIGR01841', 'IPR010127', 'IPR014176', 'PIRSF028226'},
    'phaQ': {'IPR014091'},
    'phaR_regulator': {'IPR010134', 'TIGR01848', 'PF05233', 'PF07879', 'IPR007897', 'IPR012909'},
    'phaR_synthase': {'IPR011729', 'TIGR02132'},
    'phaJ': {'IPR002539'},
}

def reaudit_family(family, audit_rows):
    """audit_rows: list of dicts with reference_query, n_targets_recruited, audit_class."""
    revised = []
    n_flipped_to_wrong = 0
    n_flipped_to_ontarget = 0
    n_demoted_no_marker = 0
    family_markers = GOOD_DOMAINS.get(family)
    for row in audit_rows:
        acc = row['reference_query'].replace('UNIPROT:', '')
        doms = load_domains(acc)
        cls = row['audit_class']
        new_cls = cls

        bad_hit = doms & BAD_DOMAINS_UNIVERSAL
        if family == 'phaF' and doms & {'IPR007208'} and doms & PHAF_ANTIPORTER_TELL:
            bad_hit = bad_hit | {'IPR007208(+antiporter NF tell)'}

        good_hit = doms & (family_markers or set())

        if bad_hit and not cls.startswith('WRONG'):
            new_cls = 'WRONG_GENE-interpro_fullaudit'
            n_flipped_to_wrong += 1
        elif good_hit and (cls.startswith('UNCLEAR') or cls.startswith('WRONG_SPECIFIC')):
            # only promote out of UNCLEAR, and only demote WRONG_SPECIFIC (wrong-family
            # guess) back to on-target if a real marker for THIS family is present and
            # no bad flag fired
            if not bad_hit:
                new_cls = 'ON_TARGET-interpro_fullaudit'
                n_flipped_to_ontarget += 1
        elif (cls.startswith('ON_TARGET') and family_markers and doms and not good_hit and not bad_hit):
            # the gap that let 152 phaC references through undetected: an
            # ON_TARGET classification is a claim, not a fact -- demote it if
            # domain data exists (so this isn't just an annotation-coverage gap)
            # but none of the family's OWN confirmed marker domains are present.
            # Requires family_markers to be defined -- a family missing from
            # GOOD_DOMAINS entirely is a bug in this dict, not evidence the
            # query is wrong, so it is left alone rather than silently demoted.
            new_cls = 'WRONG_GENE-interpro_fullaudit_no_marker'
            n_demoted_no_marker += 1

        row2 = dict(row)
        row2['audit_class_revised'] = new_cls
        row2['interpro_domains_checked'] = len(doms) > 0
        revised.append(row2)

    return revised, n_flipped_to_wrong, n_flipped_to_ontarget, n_demoted_no_marker


if __name__ == '__main__':
    family = sys.argv[1]
    audit_tsv = sys.argv[2]
    with open(audit_tsv) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        rows = [row for row in r if row.get('family', family) == family or 'family' not in row]

    revised, n_wrong, n_good, n_demoted = reaudit_family(family, rows)
    n_no_domain_data = sum(1 for r in revised if not r['interpro_domains_checked'])
    print(f'{family}: {len(revised)} refs, {n_no_domain_data} with no InterPro data at all, '
          f'{n_wrong} flipped ON_TARGET->WRONG (blacklist hit), '
          f'{n_demoted} flipped ON_TARGET->WRONG (family marker absent), '
          f'{n_good} flipped UNCLEAR/WRONG_SPECIFIC->ON_TARGET')

    out_path = f'/tmp/{family}_interpro_reaudit.tsv'
    with open(out_path, 'w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=list(revised[0].keys()), delimiter='\t')
        w.writeheader()
        w.writerows(revised)
    print(f'wrote {out_path}')
