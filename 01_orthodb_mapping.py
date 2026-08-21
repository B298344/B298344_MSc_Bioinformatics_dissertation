#!/usr/bin/env python3
"""Map the three reviewed human references to Eukaryota-level OrthoDB v12 OGs."""
import io, json, re
import requests
import pandas as pd

BASE = "https://data.orthodb.org/v12"
REFERENCES = {"SMG5": "Q9UPR3", "SMG6": "Q86US8", "SMG7": "Q92540"}
EUKARYOTA = "2759"


def get_json(endpoint, params):
    r = requests.get(f"{BASE}/{endpoint}", params=params, timeout=180)
    r.raise_for_status()
    return r.json()


def get_text(endpoint, params):
    r = requests.get(f"{BASE}/{endpoint}", params=params, timeout=180)
    r.raise_for_status()
    return r.text


def map_accession(accession):
    search = get_json("genesearch", {"query": accession})
    gene_ids = sorted(set(re.findall(r"\b9606_\d+:[0-9a-fA-F]+\b", json.dumps(search))))
    if not gene_ids:
        raise RuntimeError(f"No human OrthoDB gene found for {accession}")
    odb_gene = gene_ids[0]

    details = get_json("ogdetails", {"id": odb_gene})
    ncbi_ids = re.findall(r"\b\d+\b", json.dumps(details.get("data", details).get("entrezgene", "")))
    if not ncbi_ids:
        raise RuntimeError(f"No Entrez GeneID found for {accession}")

    result = get_json("search", {"gid": ncbi_ids[0], "take": 10000})
    ogs = result.get("data", [])
    if isinstance(ogs, str):
        ogs = [ogs]

    euk_ogs = []
    for og in sorted(set(ogs)):
        if not og.endswith(f"at{EUKARYOTA}"):
            continue
        table = pd.read_csv(io.StringIO(get_text("tab", {"id": og})), sep="\t", dtype=str)
        euk_ogs.append(og)
        table.to_csv(f"{accession}_{og}_members.tsv", sep="\t", index=False)
    return euk_ogs


rows = []
for protein, accession in REFERENCES.items():
    ogs = map_accession(accession)
    rows.append({"protein": protein, "UniProtKB": accession, "Eukaryota_OG": ";".join(ogs)})

pd.DataFrame(rows).to_csv("human_reference_Eukaryota_OGs.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
