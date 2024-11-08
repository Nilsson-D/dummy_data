#!/usr/bin/env Rscript
## Script for annotating differential expression results for RNA-seq data 
## Author Daniel Nilsson
## This should be run in the funannot/ directory


suppressMessages({
library(goseq)
library(GO.db)
library(reactome.db)
library(org.Mm.eg.db)


# Read the annotation file
anno_data_path <- "../reference/mm-biomart99-genes.txt"
annotations <- read.delim(anno_data_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

deseq2_results_path <- "../4_dge/dge_results_full.txt"
deseq2_results <- read.delim(deseq2_results_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Merge the Count Table with the Annotations
annotated_results <- deseq2_results %>%
  left_join(annotations, by = join_by("ensembl_gene_id") )

# Prepare Gene Lists for Up- and Down-Regulated Genes
lfc_threshold <- 0.5
padj_threshold <- 0.05

# Filter for up- and down-regulated genes
up_genes <- annotated_results %>%
  dplyr::filter(log2FoldChange > lfc_threshold, padj < padj_threshold) %>%
  pull(ensembl_gene_id)

down_genes <- annotated_results %>%
  dplyr::filter(log2FoldChange < -lfc_threshold, padj < padj_threshold) %>%
  pull(ensembl_gene_id)

# Run GO Enrichment Analysis ---------------------
#Up regulated-----------------
# Prepare gene vector with 1 for selected genes and 0 for all others
all_genes <- as.integer(annotations$ensembl_gene_id %in% up_genes)
names(all_genes) <- annotations$ensembl_gene_id

# Remove duplicates to prevent the row name error 
all_genes <- all_genes[!duplicated(names(all_genes))]

# GO analysis for up-regulated genes
pwf <- nullp(all_genes, "mm9", "ensGene",plot.fit=FALSE) 
go_up <- goseq(pwf = pwf, genome = "mm9", id = "ensGene", test.cats = c("GO:CC", "GO:BP", "GO:MF"))

#Down regulated----------------
# GO analysis for down-regulated genes
all_genes <- as.integer(annotations$ensembl_gene_id %in% down_genes)
names(all_genes) <- annotations$ensembl_gene_id
# Remove duplicates to prevent the row name error
all_genes <- all_genes[!duplicated(names(all_genes))]

# GO analysis for up-regulated genes
pwf <- nullp(all_genes, "mm9", "ensGene",plot.fit=FALSE)

go_down <- goseq(pwf = pwf, genome = "mm9", id = "ensGene", test.cats = c("GO:CC", "GO:BP", "GO:MF"))


# Reactome---------------------------------------------------
# Map genes to reactome pathways 
all_entrez_ids <- keys(org.Mm.eg.db, keytype = "ENTREZID")

# Map Entrez IDs to reactome pathway IDs
# Retrieve Entrez ID, reactome ID, and Pathway Name
entrez_to_reactome <- AnnotationDbi::select(reactome.db, keys = all_entrez_ids, keytype = "ENTREZID", columns = c("REACTOMEID", "PATHNAME"))
entrez_to_reactome <- na.omit(entrez_to_reactome) # Remove any rows with NA values


#Make a similar connection for mapping between ensemble and reactome as 
#in the  https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf "grepKEGG"
en2eg <- as.list(org.Mm.egENSEMBL2EG)
eg2reactome <- split(entrez_to_reactome$REACTOMEID, entrez_to_reactome$ENTREZID)

grepReactome <- function(id, mapkeys) {
  unique(unlist(mapkeys[id], use.names = FALSE))
}
reactome_ensemble_database <- lapply(en2eg, grepReactome, eg2reactome)

# Test up regulated with reactome ---------------------------------------------------------------
all_genes <- as.integer(annotations$ensembl_gene_id %in% up_genes)
names(all_genes) <- annotations$ensembl_gene_id
all_genes <- all_genes[!duplicated(names(all_genes))]

pwf <- nullp(all_genes, "mm9", "ensGene", plot.fit = FALSE)
reactome_up <- goseq(pwf, gene2cat = reactome_ensemble_database)

# Test down regulated with reactome ---------------------------------------------------------------
all_genes <- as.integer(annotations$ensembl_gene_id %in% down_genes)
names(all_genes) <- annotations$ensembl_gene_id
all_genes <- all_genes[!duplicated(names(all_genes))]

pwf <- nullp(all_genes, "mm9", "ensGene", plot.fit = FALSE)
reactome_down <- goseq(pwf, gene2cat = reactome_ensemble_database)


# Just add the path name to the reactome output as the category do not say much
reactome_pathname <- entrez_to_reactome %>% dplyr::select(REACTOMEID, PATHNAME)

reactome_up <- reactome_up %>%
  left_join(reactome_pathname, by = join_by("category" == "REACTOMEID")) %>% distinct()
reactome_down <- reactome_down %>%
  left_join(reactome_pathname, by = join_by("category" == "REACTOMEID")) %>% distinct()

# Save the annotated tables to a file
write.table(go_up, "go_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(go_down, "go_down.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(reactome_up, "reactome_up.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(reactome_down, "reactome_down.txt", sep = "\t", quote = FALSE, row.names = FALSE)

})
