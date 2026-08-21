# MSc pseudoenzymes — selected analysis code

This folder contains a  small selection of scripts supporting the computational analyses described in the dissertation Methods. It is not intended to reproduce every exploratory or manual curation step.

Included scripts:

1. `01_orthodb_mapping.py` — maps the reviewed human SMG5/SMG6/SMG7 references to Eukaryota-level OrthoDB v12 orthogroups.
2. `02_parnas_sampling.py` — evaluates PARNAS diversity for k = 2–35 and selects the final 25 representatives (5 fixed + 20 selected).
3. `03_build_PIN32_dataset.py` — constructs the 32-sequence PIN-domain dataset from the PARNAS subset plus NMD4 and SWT1.
4. `04_phylogenetic_workflows.sh` — records the final MAFFT, ClipKIT and IQ-TREE 3 commands for the 90-protein full-length and 32-protein PIN-domain analyses; the exploratory 14-3-3-like analysis is included as an optional final block because it is reported briefly in Results.
5. `05_map_functional_residues.py` — maps the experimentally characterised SMG5/SMG6 residues onto an alignment and records exact residue conservation.

Manual/database steps described in the dissertation (PANTHER inspection, NCBI BLAST/PSI-BLAST, EukProt searches, InterPro web screening, Domain Architecture Search and Jalview inspection) are not represented here. Plotting-only code and superseded exploratory analyses are also omitted.

Software versions reported in the dissertation include MAFFT 7.490, PARNAS 0.1.7, ClipKIT 2.14.0, IQ-TREE 3.1.3 and Biopython 1.88.
