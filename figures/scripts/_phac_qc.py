"""Shared QC filter for phaC hit data.

Every reference sequence that ever recruited a phaC hit in the original
OMDB search was individually resolved against UniProt and InterPro and
checked against the expected phaC identity -- not just the earlier,
explicitly partial 19-accession pass this module used to carry. See
figures/PHAC_VALIDATION_FINDINGS.md section 10 and figures/PHA_ALL_
FAMILIES_REFERENCE_AUDIT.md for the full methodology and evidence, and
figures/PHA_CLEAN_RESULTS.md for the resulting verified dataset. The
dominant contamination source is a real, unrelated bacterial gene-naming
collision: a multi-subunit K+/Na+-H+ antiporter system historically also
called "Pha" (for pH adaptation), whose subunit letters collide with this
project's own PHA pathway gene names -- confirmed directly in several
baits' own UniProt names (e.g. "K+/H+ antiporter subunit C, putative
PhaC"). The rest is a mix of unrelated enzymes recruited by the same
gene-name text search (isochorismate synthase, a toxin-antitoxin system
toxin, fatty-acid/polyketide synthase) and a few real PHA-pathway
proteins filed under the wrong family (phaB, phaE, phaZ).

BAD_QUERIES below is the exact same 219-accession set as
`phaatlas.pipeline.pathway_architecture.PHAC_BAD_QUERIES` -- duplicated
here rather than imported, per this project's convention of keeping
figures/scripts self-contained rather than reaching into the phaatlas
package. Keep the two in sync if either changes.

Updated 2026-09-22: 152 accessions added on top of the original 67, after
a full direct InterPro re-verification of every previously ON_TARGET-
classified phaC reference found the original "8-query FUSION" exclusion
was not exhaustive -- see PHA_CLEAN_RESULTS.md and the corresponding
comment in pathway_architecture.py for the full methodology.

Every figure-generating script that reads phaC_unique_targets_with_metadata.tsv
should skip rows where best_query is in BAD_QUERIES (or use
load_bad_targets() when it needs to also filter a joined file like
phaC_cluster0.7_cluster.tsv, which is keyed by target_id rather than
carrying best_query itself).
"""
import csv
from pathlib import Path

FA = Path('/Users/hellpark/multimodal_seusmbol/PHA_Ocean_Atlas/PHA_bioprospecting/omdb_search/results')

BAD_QUERIES = frozenset(f'UNIPROT:{acc}' for acc in (
    'A0A090R7C3', 'A0A097ERZ2', 'A0A0C3DJ16', 'A0A0D6JHM2', 'A0A0D8PXE4', 'A0A0G3BRR5', 'A0A0G3WHX8', 'A0A0H2Z7T5',
    'A0A0H3LUG2', 'A0A0H3LVF9', 'A0A0J1GRI2', 'A0A0J1K4W4', 'A0A0R0DVU2', 'A0A0R4J7U1', 'A0A0T7DNW4', 'A0A150WMA4',
    'A0A150WNK1', 'A0A151AEV0', 'A0A151AIL1', 'A0A151L0P1', 'A0A157R4S9', 'A0A157SVZ8', 'A0A178JFB0', 'A0A1B0ZMT2',
    'A0A1B9QV15', 'A0A1E8EW93', 'A0A1G5SDS9', 'A0A1G8Y6C9', 'A0A1L3GI38', 'A0A1Q9G762', 'A0A1S6R531', 'A0A1T5FPD3',
    'A0A1W6NY09', 'A0A1X7AFX8', 'A0A1Y4DIX0', 'A0A1Y6KRP4', 'A0A1Y6M9P7', 'A0A240EG07', 'A0A2M8H6D7', 'A0A2N4UXW8',
    'A0A2P8GDY2', 'A0A2S7V7Z2', 'A0A2T3NAI5', 'A0A2T3PRJ9', 'A0A330LTS8', 'A0A345S6Y8', 'A0A348FW77', 'A0A3G8JFS6',
    'A0A3Q9K0Y4', 'A0A3S9XZH6', 'A0A444JWC5', 'A0A480AP87', 'A0A498Q0R0', 'A0A4D6LUZ9', 'A0A4R1K402', 'A0A4R2NZ35',
    'A0A4U3FLS4', 'A0A4Y3HY53', 'A0A4Y3IMB5', 'A0A4Y3IRN2', 'A0A4Y8WJ87', 'A0A556QKC2', 'A0A557ST61', 'A0A5J6WQD8',
    'A0A5J6WST3', 'A0A5J6WUA4', 'A0A5Q0TAN5', 'A0A5S9QVF8', 'A0A6A7GGQ7', 'A0A6P1M6T9', 'A0A6S6XVD8', 'A0A7J5AGP2',
    'A0A7U7GEB2', 'A0A7X1E4U9', 'A0A7Z2T3G5', 'A0A7Z7ILZ6', 'A0A812VQC3', 'A0A813BGK3', 'A0A829YA03', 'A0A853R7Z3',
    'A0A8J3E0L9', 'A0A8S2BDR2', 'A0A8S2BLM7', 'A0A916RII9', 'A0A916RIM5', 'A0A916S023', 'A0A916ZXM0', 'A0A917QD60',
    'A0A918ITM8', 'A0A918JWF9', 'A0A918JXY7', 'A0A918W8J3', 'A0A9W6GI39', 'A0A9W6GL72', 'A0A9W6GNE3', 'A0A9W6LLT3',
    'A0A9X1WD73', 'A0A9X3HW65', 'A0AAD1BY69', 'A0AAI7ZCY1', 'A0AAN1CWM3', 'A0AAN4VX63', 'A0AAU8BIW8', 'A0AAU9D2C6',
    'A0AAV5NSL7', 'A0AAX0Z1X0', 'A0AAX2LJW2', 'A0ABM7GZ97', 'A0ABM8ZQR0', 'A0ABN5C053', 'A0ABN5CAC2', 'A0ABN5JGB5',
    'A0ABN6T5N2', 'A0ABN8E4V9', 'A0ABN8EHT7', 'A0ABN8JP91', 'A0ABP2DAC8', 'A0ABP7VE29', 'A0ABP8QGS0', 'A0ABP9S7A1',
    'A0ABQ1I0H1', 'A0ABQ1J139', 'A0ABQ1KAJ0', 'A0ABQ2CQ11', 'A0ABQ2WTK0', 'A0ABQ5TXW2', 'A0ABQ5VKV7', 'A0ABQ5Y5A3',
    'A0ABQ6DXS1', 'A0ABR7J2E5', 'A0ABS2G0X7', 'A0ABS2HP62', 'A0ABT1N2I8', 'A0ABT7QAT8', 'A0ABU4W9W4', 'A0ABU7G4P0',
    'A0ABU9FSZ5', 'A0ABU9HU42', 'A0ABU9JAJ3', 'A0ABW9GBA2', 'A0ABX0DA28', 'A0ABX1U6D0', 'A0ABX3AZ03', 'A0ABX4XAX1',
    'A0ABX5GGB4', 'A0ABX5H8K4', 'A0ABX5H8M3', 'A0ABX9KEK0', 'A0ABX9KJH0', 'A0ABX9KJK0', 'A0ABY1HDR6', 'A0ABY3S282',
    'A0ABY5G6F1', 'A0ABZ0F622', 'A0ACE0BGR7', 'A0ACH1R8B4', 'A0ACH1RK58', 'A0ACH1RMQ7', 'A0ACH3UPZ5', 'A0ACM8RTF1',
    'A0ACM8RWN5', 'A0ACN3LVC5', 'A0ACN3P6Y2', 'A0ACN3VS62', 'A0ACN4YJT2', 'A9CK78', 'A9FZG2', 'A9IP81',
    'B2FKM3', 'B6JCB5', 'C7BNH2', 'C8YNX0', 'D0YZB4', 'E1SWP3', 'E3HAX3', 'E8M187',
    'I3R9Z3', 'I3TJ70', 'K1KB67', 'K6DJS5', 'M0DCW3', 'M0I0U4', 'M0ICI7', 'N6UCP8',
    'N6VCC7', 'N9VHT5', 'P45367', 'P45372', 'P73389', 'P76108', 'Q162R1', 'Q17V39',
    'Q17V43', 'Q2KZX2', 'Q3A8F5', 'Q4KDG6', 'Q52980', 'Q5E3C0', 'Q5P8E8', 'Q5UYM1',
    'Q6FF34', 'Q6RI97', 'Q87KX4', 'Q8GFF3', 'Q8GI83', 'Q8KR79', 'Q8PD95', 'Q92RA1',
    'Q9F5P9', 'Q9KNR9', 'R1F852', 'U3B5V2', 'U4KH66', 'U5NXH1', 'U7D663', 'U7VCC2',
    'V5FG99', 'V5FMS0', 'W5U290',
))

# kept for backwards compat with earlier-written scripts/notes referencing a
# single accession -- prefer BAD_QUERIES / is_bad() in new code
BAD_QUERY = 'UNIPROT:C7BNH2'

# Target-level exclusions, as opposed to BAD_QUERIES above (which excludes by
# best_query, i.e. a bad *reference*). These are individual candidates ruled
# out directly by their own structure, not by which reference recruited them
# -- excluding by best_query doesn't apply here since these 19 don't share
# one bad reference (12/19 best-match Q5P962, which is itself a legitimate,
# if weakly-annotated, reference that also recruits many genuine targets --
# see PHA_CLEAN_RESULTS.md section 9.7).
#
# All 19 are the full no_hmm_triad_support population from section 9.6's
# uncertain-set build (figures/structural_no_hmm_triad_active_site_audit.tsv
# + _no_cys_alternative_triad_audit.tsv, section 9.12): every one of them
# was checked directly against its own ESMFold structure for a real
# catalytic triad, canonical Cys-Asp-His or alternative Ser/Thr-Asp-His, and
# none was found -- 15/19 have no cysteine anywhere in the sequence at all
# (ruling out the canonical mechanism outright), the other 4 have a
# cysteine but 8-17A from the nearest His (an order of magnitude looser
# than the 2.4-3.4A seen in every confirmed-real triad in this project),
# and all 15 no-cysteine ones were also checked for a Ser/Thr alternative
# and found none plausible either. No active site by any mechanism checked
# -- not phaC.
BAD_TARGET_IDS = frozenset((
    'OMDBv2.0_AA_G_NR100_000000013407', 'OMDBv2.0_AA_G_NR100_000007484248', 'OMDBv2.0_AA_G_NR100_000012946226',
    'OMDBv2.0_AA_G_NR100_000018865739', 'OMDBv2.0_AA_G_NR100_000057104467', 'OMDBv2.0_AA_G_NR100_000057745835',
    'OMDBv2.0_AA_G_NR100_000060086046', 'OMDBv2.0_AA_G_NR100_000098889079', 'OMDBv2.0_AA_G_NR100_000101178016',
    'OMDBv2.0_AA_G_NR100_000111006886', 'OMDBv2.0_AA_G_NR100_000111912241', 'OMDBv2.0_AA_G_NR100_000160099602',
    'OMDBv2.0_AA_G_NR100_000187164299', 'OMDBv2.0_AA_G_NR100_000204524282', 'OMDBv2.0_AA_G_NR100_000211778674',
    'OMDBv2.0_AA_G_NR100_000218113895', 'OMDBv2.0_AA_G_NR100_000218537902', 'OMDBv2.0_AA_G_NR100_000225793822',
    'OMDBv2.0_AA_G_NR100_000236086924',
))


def is_bad(best_query):
    return best_query in BAD_QUERIES


def load_bad_targets(path=None):
    """Returns the set of target_ids that should be excluded from any
    figure or analysis: either their best_query is one of the confirmed-bad
    reference proteins (BAD_QUERIES), or the target itself was directly
    ruled out by its own structure (BAD_TARGET_IDS) -- for filtering files
    (like the cluster-assignment TSVs) that don't carry best_query
    themselves."""
    path = path or FA / 'phaC_unique_targets_with_metadata.tsv'
    bad = set(BAD_TARGET_IDS)
    with open(path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['best_query'] in BAD_QUERIES:
                bad.add(row['target_id'])
    return bad
