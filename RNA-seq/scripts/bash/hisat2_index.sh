#!/bin/bash
## 2020 Roy Francis

# run from directory /reference directly

hisat2-build \
-p 1 \
Mus_musculus.GRCm38.dna.chromosome.19.fa \
mouse_chr19_hisat2/mouse_chr19_hisat2
