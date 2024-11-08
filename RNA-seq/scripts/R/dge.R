#!/usr/bin/env Rscript
#Script from NBIS rnaseq-intro
## R script to run differential gene expression


library(edgeR)
library(DESeq2)
#library(methods)

# input
cat("Reading counts and annotation files ...\n")

# count data
file_counts <- "counts_full.txt"

# metadata
table_meta <- data.frame(accession=c("SRR3222409","SRR3222410","SRR3222411","SRR3222412","SRR3222413","SRR3222414"),condition=c(rep(c("KO","Wt"),each=3)),replicate=rep(1:3,2))
table_meta$condition <- factor(as.character(table_meta$condition),levels=c("Wt","KO"))
#table_meta$condition <- relevel(table_meta$condition,"Wt")
rownames(table_meta) <- table_meta$accession

# reading in data
data_counts <- read.delim(file_counts, skip=1)
table_counts <- as.matrix(data_counts[,-c(1:6)])
row.names(table_counts) <- data_counts$Geneid

# clean-up column labels
# base r solution
mth <- regexpr("SRR[0-9]+",colnames(table_counts))
colnames(table_counts) <- substr(colnames(table_counts),start=mth,stop=(mth+attr(mth,"match.length"))-1)
# easier tidyverse solution
# colnames(table_counts) <- stringr::str_extract(colnames(table_counts),"SRR[0-9]+")

# match order of counts and metadata
mth <- match(colnames(table_counts),rownames(table_meta))
table_counts <- table_counts[,mth]
all.equal(rownames(table_meta),colnames(table_counts))

# remove genes with low counts
# keep genes that have minimum 1 CPM across 3 samples (since group has three replicates)
keepgenes <- rowSums(edgeR::cpm(table_counts)>1) >= 3
table_counts <- table_counts[keepgenes,]

cat("\nPreview of count table:\n")
print(head(table_counts))
cat("\nView of metadata table:\n")
print(table_meta)

## deseq2
cat("\nRunning DESeq2 ...\n")

# run deseq2
d <- DESeqDataSetFromMatrix(countData=table_counts,colData=table_meta,design=~condition)
d <- DESeq2::estimateSizeFactors(d,type="ratio")
d <- DESeq2::estimateDispersions(d)

# export vst counts
cv <- as.data.frame(assay(varianceStabilizingTransformation(d,blind=T)),check.names=F)
write.table(cv,"counts_vst_full.txt",sep="\t",dec=".",quote=FALSE)
saveRDS(cv,"counts_vst_full.Rds")

dg <- nbinomWaldTest(d)
print(resultsNames(dg))
res <- results(dg,contrast=c("condition","KO","Wt"),alpha=0.05)
cat("\nSummary of DEGs:\n")
summary(res)

# lfc shrink
res1 <- lfcShrink(dg, contrast=c("condition","KO","Wt"), res=res, type="normal")

# convert table to data.frame
table_res <- as.data.frame(res1)
table_res$ensembl_gene_id <- rownames(table_res)

cat("Preview of DEG table:\n")
print(head(table_res))

cat("\nExporting results ...\n")
write.table(table_res, "dge_results_full.txt", sep="\t", quote=F, row.names=F)
saveRDS(table_res,"dge_results_full.Rds", version = 2)

cat("Completed ...\n")

