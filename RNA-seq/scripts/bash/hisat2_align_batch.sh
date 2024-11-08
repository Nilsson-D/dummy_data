#!/bin/bash

# Loop through each read1 file
for r1 in ../1_raw/*1.fq.gz; do
  # Define the corresponding read2 file by replacing "r1" with "r2" in the filename
  r2="${r1/r1/r2}"


  # Extract the base name of the sample
  sample_name=$(basename "$r1" | sed 's/_1.fastq.gz//')

  # Define output files
  out_sam="./${sample_name}.sam"

  ## use the sample name and add suffix (_1.fq.gz or _2.fq.gz)
  echo "Mapping ${r1} and ${r2} ..."
  bash ../new_scripts/rnaseq-hisat2_align-pe.sh ${r1} ${r2}
done
