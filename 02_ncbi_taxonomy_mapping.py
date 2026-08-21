#!/usr/bin/env python3
"""Map candidate species/TaxIDs to NCBI Taxonomy.

Cleaned version of the taxonomy-mapping code used during dataset curation.
Taxon identifiers supplied by source databases are used when available; if a
TaxID is absent, the species scientific name is searched in NCBI Taxonomy.
Taxonomic lineages are then retrieved with the NCBI Entrez E-utilities EFetch
service, matching the procedure described in the dissertation Methods.

Input: CSV/TSV containing a species column and/or a TaxID column.
Output: the same records annotated with NCBI TaxID, standard taxonomic ranks
and the full NCBI lineage.
"""

from __future__ import annotations
import argparse
import os
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez

RANKS = [
    "superkingdom", "kingdom", "phylum", "class",
    "order", "family", "genus", "species",
]


def clean_taxid(value) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    # Handles values read by pandas as e.g. 9606.0
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return text if text.isdigit() else ""


def find_taxid(scientific_name: str) -> str:
    """Resolve a scientific name to an NCBI TaxID."""
    scientific_name = str(scientific_name).strip()
    if not scientific_name:
        return ""

    query = f'"{scientific_name}"[Scientific Name]'
    with Entrez.esearch(db="taxonomy", term=query, retmax=5) as handle:
        result = Entrez.read(handle)
    ids = result.get("IdList", [])
    time.sleep(0.34)
    return str(ids[0]) if ids else ""


def parse_taxonomy_record(record: dict) -> dict:
    """Extract standard ranks and lineage from one NCBI taxonomy record."""
    rank_values = {rank: "" for rank in RANKS}

    for ancestor in record.get("LineageEx", []):
        rank = str(ancestor.get("Rank", "")).lower()
        if rank in rank_values:
            rank_values[rank] = str(ancestor.get("ScientificName", ""))

    current_rank = str(record.get("Rank", "")).lower()
    if current_rank in rank_values:
        rank_values[current_rank] = str(record.get("ScientificName", ""))

    return {
        "NCBI_TaxID": str(record.get("TaxId", "")),
        "NCBI_scientific_name": str(record.get("ScientificName", "")),
        **{rank.capitalize(): rank_values[rank] for rank in RANKS},
        "NCBI_lineage": str(record.get("Lineage", "")),
    }


def fetch_taxonomy(taxids: list[str], batch_size: int = 100) -> dict[str, dict]:
    """Retrieve NCBI taxonomy records in batches using Entrez EFetch."""
    resolved: dict[str, dict] = {}
    unique_taxids = sorted(set(t for t in taxids if t), key=int)

    for start in range(0, len(unique_taxids), batch_size):
        batch = unique_taxids[start:start + batch_size]
        with Entrez.efetch(
            db="taxonomy",
            id=",".join(batch),
            retmode="xml",
        ) as handle:
            records = Entrez.read(handle)

        for record in records:
            parsed = parse_taxonomy_record(record)
            resolved[parsed["NCBI_TaxID"]] = parsed

        time.sleep(0.4)

    return resolved


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input .tsv or .csv file")
    parser.add_argument("-o", "--output", default="taxonomy_mapped_to_NCBI.tsv")
    parser.add_argument("--species-column", default="Species")
    parser.add_argument("--taxid-column", default="TaxID")
    parser.add_argument(
        "--email",
        default=os.getenv("NCBI_EMAIL"),
        help="Email required by NCBI Entrez (or set NCBI_EMAIL).",
    )
    args = parser.parse_args()

    if not args.email:
        raise SystemExit("Provide --email or set NCBI_EMAIL for NCBI Entrez.")
    Entrez.email = args.email

    path = Path(args.input)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep, dtype=str).fillna("")

    if args.species_column not in df.columns and args.taxid_column not in df.columns:
        raise SystemExit(
            f"Input must contain '{args.species_column}' and/or '{args.taxid_column}'."
        )

    # Keep source-database TaxIDs where available.
    if args.taxid_column in df.columns:
        taxids = df[args.taxid_column].map(clean_taxid)
    else:
        taxids = pd.Series([""] * len(df), index=df.index)

    # Resolve only records lacking a source TaxID.
    if args.species_column in df.columns:
        cache: dict[str, str] = {}
        for idx in df.index[taxids.eq("")]:
            species = str(df.at[idx, args.species_column]).strip()
            if species not in cache:
                cache[species] = find_taxid(species)
            taxids.at[idx] = cache[species]

    df["NCBI_TaxID"] = taxids
    taxonomy = fetch_taxonomy(df["NCBI_TaxID"].tolist())

    annotation_columns = [
        "NCBI_scientific_name", "Superkingdom", "Kingdom", "Phylum",
        "Class", "Order", "Family", "Genus", "Species", "NCBI_lineage",
    ]
    for column in annotation_columns:
        df[column] = df["NCBI_TaxID"].map(
            lambda tid: taxonomy.get(tid, {}).get(column, "")
        )

    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(df)} records to {args.output}")


if __name__ == "__main__":
    main()
