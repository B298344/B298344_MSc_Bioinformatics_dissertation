#!/usr/bin/env python3
"""Construct the final 32-sequence PIN-domain dataset."""
import argparse, subprocess, tempfile
from pathlib import Path
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

GAPS = set("-?.*Xx")


def load_fasta(path):
    return {r.id: r for r in SeqIO.parse(path, "fasta")}


def find_id(records, accession):
    hits = [x for x in records if accession in x]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one sequence containing {accession}: {hits}")
    return hits[0]


def column_for_residue(aligned_sequence, residue_number):
    n = 0
    for col, aa in enumerate(str(aligned_sequence), 1):
        if aa not in GAPS:
            n += 1
            if n == residue_number:
                return col
    raise ValueError(residue_number)


def residue_numbers_in_columns(aligned_sequence, columns):
    n, lookup = 0, {}
    for col, aa in enumerate(str(aligned_sequence), 1):
        if aa not in GAPS:
            n += 1
            lookup[col] = n
    return [lookup[c] for c in columns if c in lookup]


ap = argparse.ArgumentParser()
ap.add_argument("--full90-fasta", required=True)
ap.add_argument("--full90-alignment", required=True, help="Untrimmed E-INS-i alignment")
ap.add_argument("--parnas-tree", required=True, help="PARNAS 25-tip subtree")
ap.add_argument("--nmd4-fasta", required=True, help="Five full-length NMD4 proteins")
ap.add_argument("--swt1-fasta", required=True, help="Two full-length SWT1 proteins")
ap.add_argument("--output", default="Final_PIN_32_unaligned.fasta")
args = ap.parse_args()

raw90 = load_fasta(args.full90_fasta)
aln90 = load_fasta(args.full90_alignment)
selected = [x.name for x in Phylo.read(args.parnas_tree, "newick").get_terminals()]

# Project the human SMG5 (855–1016) and SMG6 (1246–1419) PIN reference
# boundaries through the untrimmed full-length alignment.
smg5, smg6 = find_id(aln90, "Q9UPR3"), find_id(aln90, "Q86US8")
boundary_cols = [
    column_for_residue(aln90[smg5].seq, 855), column_for_residue(aln90[smg5].seq, 1016),
    column_for_residue(aln90[smg6].seq, 1246), column_for_residue(aln90[smg6].seq, 1419)
]
columns = range(min(boundary_cols), max(boundary_cols) + 1)

pin_records = []
for sid in selected:
    nums = residue_numbers_in_columns(aln90[sid].seq, columns)
    start, end = min(nums), max(nums)
    pin_records.append(SeqRecord(Seq(str(raw90[sid].seq)[start-1:end]), id=sid, description=""))

# NMD4: align five proteins and project the alignment columns occupied by
# S. cerevisiae Q12129 residues 1–167.
nmd4 = load_fasta(args.nmd4_fasta)
with tempfile.TemporaryDirectory() as td:
    aln_path = Path(td) / "nmd4_linsi.fasta"
    with open(aln_path, "w") as out:
        subprocess.run([
            "mafft", "--localpair", "--maxiterate", "1000", args.nmd4_fasta
        ], stdout=out, check=True)
    aln = load_fasta(aln_path)
    anchor = find_id(aln, "Q12129")
    anchor_cols = [column_for_residue(aln[anchor].seq, p) for p in range(1, 168)]
    for sid, record in nmd4.items():
        nums = residue_numbers_in_columns(aln[sid].seq, anchor_cols)
        start, end = min(nums), max(nums)
        pin_records.append(SeqRecord(Seq(str(record.seq)[start-1:end]), id=sid, description=""))

# SWT1 PIN boundaries from concordant IPR002716/IPR029060 annotations.
swt1 = load_fasta(args.swt1_fasta)
for accession, start, end in [("Q12104", 130, 275), ("Q9P7J1", 70, 202)]:
    sid = find_id(swt1, accession)
    pin_records.append(SeqRecord(Seq(str(swt1[sid].seq)[start-1:end]), id=sid, description=""))

if len(pin_records) != 32:
    raise RuntimeError(f"Expected 32 PIN sequences, observed {len(pin_records)}")
SeqIO.write(pin_records, args.output, "fasta")
