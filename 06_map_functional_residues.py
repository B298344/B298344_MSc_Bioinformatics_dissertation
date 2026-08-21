#!/usr/bin/env python3
"""Map the 12 experimentally characterised SMG5/SMG6 positions to an alignment."""
import argparse
import pandas as pd
from Bio import SeqIO

GAPS = set("-?.*Xx")
SITES = [
    ("A", "SMG6", "Q86US8", 1251, "D", "Catalytic tetrad", 1246),
    ("B", "SMG6", "Q86US8", 1282, "E", "Catalytic tetrad", 1246),
    ("C", "SMG6", "Q86US8", 1353, "D", "Catalytic tetrad", 1246),
    ("D", "SMG6", "Q86US8", 1392, "D", "Catalytic tetrad", 1246),
    ("E", "SMG5", "Q9UPR3", 893,  "D", "Active-site complement", 855),
    ("F", "SMG5", "Q9UPR3", 896,  "K", "RNA positioning", 855),
    ("G", "SMG5", "Q9UPR3", 897,  "K", "RNA positioning", 855),
    ("H", "SMG5", "Q9UPR3", 889,  "I", "PIN-PIN interface", 855),
    ("I", "SMG5", "Q9UPR3", 906,  "I", "PIN-PIN interface", 855),
    ("J", "SMG6", "Q86US8", 1397, "V", "PIN-PIN interface", 1246),
    ("K", "SMG6", "Q86US8", 1400, "L", "PIN-PIN interface", 1246),
    ("L", "SMG6", "Q86US8", 1401, "T", "PIN-PIN interface", 1246),
]


def alignment_column(sequence, full_position, first_full_position):
    target = full_position - first_full_position + 1
    n = 0
    for col, aa in enumerate(sequence, 1):
        if aa not in GAPS:
            n += 1
            if n == target:
                return col
    raise ValueError(full_position)


ap = argparse.ArgumentParser()
ap.add_argument("alignment")
ap.add_argument("--full-length", action="store_true")
ap.add_argument("--output", default="functional_residue_states.tsv")
args = ap.parse_args()

records = list(SeqIO.parse(args.alignment, "fasta"))
seqs = {r.id: str(r.seq).upper() for r in records}
rows = []

for site, protein, ref_id, pos, expected, function, pin_start in SITES:
    reference = next((sid for sid in seqs if ref_id in sid), None)
    if reference is None:
        raise RuntimeError(f"Reference {ref_id} not found")
    first = 1 if args.full_length else pin_start
    col = alignment_column(seqs[reference], pos, first)
    for record in records:
        aa = str(record.seq).upper()[col - 1]
        rows.append({
            "sequence_id": record.id,
            "site": site,
            "reference_protein": protein,
            "reference_position": pos,
            "functional_class": function,
            "expected_residue": expected,
            "alignment_column": col,
            "observed_residue": aa,
            "exact_match": aa == expected,
        })

pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
