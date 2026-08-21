#!/usr/bin/env python3
"""PARNAS selection used for the PIN-domain representative subset."""
import argparse, re, subprocess
from pathlib import Path
from Bio import Phylo

FIXED_ACCESSIONS = ["Q9UPR3", "Q86US8", "P40456", "P36168", "EP01134_P006785"]

ap = argparse.ArgumentParser()
ap.add_argument("tree", help="Final 90-protein full-length Newick tree")
ap.add_argument("--outdir", default="PARNAS_final25")
args = ap.parse_args()

out = Path(args.outdir)
out.mkdir(exist_ok=True)
tips = [x.name for x in Phylo.read(args.tree, "newick").get_terminals()]

fixed = []
for accession in FIXED_ACCESSIONS:
    hits = [tip for tip in tips if accession in tip]
    if len(hits) != 1:
        raise RuntimeError(f"Could not uniquely identify {accession}: {hits}")
    fixed.append(hits[0])

# Diversity represented by subset sizes k = 2–35.
subprocess.run([
    "parnas", "-t", args.tree, "-n", "35", "--diversity",
    str(out / "diversity_k2-k35.csv")
], check=True)

# Five fixed representatives + 20 selected by PARNAS = 25 proteins.
prior = "^(" + "|".join(re.escape(x) for x in fixed) + ")$"
subprocess.run([
    "parnas", "-t", args.tree, "-n", "20", "--prior", prior, "--include-prior",
    "--subtree", str(out / "PARNAS_final25.tree"),
    "--clusters", str(out / "PARNAS_final25_clusters.tsv")
], check=True)
