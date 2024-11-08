#!/bin/bash

# run from directory /2_fastqc/
# $1 is a fastq filename
# $2 is the number of cores

fastqc -o . ../1_raw/*.fq.gz
