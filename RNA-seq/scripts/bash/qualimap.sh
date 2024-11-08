#!/bin/bash

# get output filename prefix
prefix=$( basename "$1" .bam)

unset DISPLAY

qualimap rnaseq -pe \
  -bam $1 \
  -gtf "../reference/Mus_musculus.GRCm38.99-19.gtf" \
  -outdir "../4_qualimap/${prefix}/" \
  -outfile "$prefix" \
  -outformat "HTML" \
  --java-mem-size=6G >& "${prefix}-qualimap.log"
