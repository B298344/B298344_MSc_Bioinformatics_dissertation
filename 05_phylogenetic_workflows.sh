#!/usr/bin/env bash
set -euo pipefail

# Final phylogenetic workflows reported in the dissertation.
# MAFFT 7.490; ClipKIT 2.14.0; IQ-TREE 3.1.3.

THREADS=${THREADS:-8}

# 1) Full-length dataset: 90 PIN-containing proteins
mafft --genafpair --maxiterate 1000 --ep 0 --thread "$THREADS" \
  Final_full_length_90.fasta > Final_full_length_90_EINSI.fasta
clipkit Final_full_length_90_EINSI.fasta -m gappy -g 0.90 -s aa \
  -o Final_full_length_90_EINSI_gappy090.fasta
# ModelFinder selected Q.YEAST+F+R7 in the final analysis.
iqtree3 -s Final_full_length_90_EINSI_gappy090.fasta -st AA -m MFP \
  -alrt 1000 -B 1000 -bnni -T AUTO --threads-max "$THREADS" \
  --prefix Final_full_length_90

# 2) PIN-domain dataset: 32 proteins
mafft --localpair --maxiterate 1000 --thread "$THREADS" \
  Final_PIN_32_unaligned.fasta > Final_PIN_32_LINSI.fasta
clipkit Final_PIN_32_LINSI.fasta -m gappy -g 0.90 -s aa \
  -o Final_PIN_32_LINSI_gappy090.fasta
# ModelFinder selected Q.PFAM+R4 in the final analysis.
iqtree3 -s Final_PIN_32_LINSI_gappy090.fasta -st AA -m MFP \
  -alrt 1000 -B 1000 -bnni -T AUTO --threads-max "$THREADS" \
  --prefix Final_PIN_32

# 3) Exploratory 14-3-3-like analysis reported briefly in Results
# (not used for the principal evolutionary inferences).
# mafft --localpair --maxiterate 1000 --thread "$THREADS" Final_36_14-3-3-like.fasta > Final_36_14-3-3-like_LINSI.fasta
# clipkit Final_36_14-3-3-like_LINSI.fasta -m gappy -g 0.90 -s aa -o Final_36_14-3-3-like_gappy090.fasta
# iqtree3 -s Final_36_14-3-3-like_gappy090.fasta -st AA -m MFP -alrt 1000 -B 1000 -bnni -T AUTO --threads-max "$THREADS" --prefix Final_36_14-3-3-like
