#!/bin/bash

featureCounts \
  -a "../reference/Mus_musculus.GRCm38.99.gtf" \
  -o "counts.txt" \
  -F "GTF" \
  -t "exon" \
  -g "gene_id" \
  -p \
  -s 0 \
  -T 1 \
  ../3_mapping/*.bam
