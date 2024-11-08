#!/bin/bash
## 2020 Roy Francis

# $1 paired-end read1
# $2 paired-end read2

if [ -z "$1" ]; then
    echo "read1 not provided."
    exit 1
fi

if [ -z "$2" ]; then
    echo "read2 not provided."
    exit 1
fi

# create output file name
prefix="${1##*/}"
prefix="${prefix/_*/}"

hisat2 \
-p 1 \
-x ../reference/mouse_chr19_hisat2/mouse_chr19_hisat2 \
--summary-file "${prefix}.summary" \
-1 $1 \
-2 $2 | samtools sort -O BAM > "${prefix}.bam"
