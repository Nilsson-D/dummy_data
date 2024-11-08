#!/usr/bin/env Rscript
## Plotting functions for RNA-seq data 
## Author Daniel Nilsson
## This should be run in the plots/ directory


suppressMessages(library(tidyverse))
suppressMessages(library(pheatmap))
suppressMessages(library(EnsDb.Mmusculus.v79))

#taken from the NBIS rnaseq
table_meta <- data.frame(accession=c("SRR3222409","SRR3222410","SRR3222411","SRR3222412","SRR3222413","SRR3222414"),condition=c(rep(c("KO","Wt"),each=3)),replicate=rep(1:3,2))
table_meta$condition <- factor(as.character(table_meta$condition),levels=c("Wt","KO"))

#deseq2 results (our DEGs)
deseq2_results_path <- "../5_dge/dge_results_full.txt"
deseq2_results <- read.delim(deseq2_results_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#variance Stabilizing Transformation assay for plotting
count_matrix_path <- "../5_dge/counts_vst_full.txt" 
count_table <- read.delim(count_matrix_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


## PCA--------------------------------------
# Perform PCA on transposed count matrix
pca <- prcomp(t(count_table), scale. = TRUE)
pca_data <- as.data.frame(pca$x)

# Add metadata to PCA data
pca_data$Accession <- rownames(pca_data)
pca_data <- merge(pca_data %>% rownames_to_column("accession"), table_meta, by = "accession")

# Define colors for conditions 
colors <- c("Wt" = "red", "KO" = "skyblue")


png("pca_plot.png", width = 480, height = 480)
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = accession)) +
  geom_point(size = 3) + 
  geom_text(vjust = -0.5, size = 3) + 
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(title = "", x = "PC1", y = "PC2") +
  theme(legend.position = "top")
dev.off()


## MA plot--------------------------------------
padj_threshold <- 0.05
# Add a significance column based on an adjusted p-value threshold
deseq2_results$Significance <- ifelse(deseq2_results$padj < padj_threshold, "Sig", "NotSig")

#remove NAs
deseq2_results <- deseq2_results %>%
  dplyr::filter(!is.na(Significance))

png("ma_plot.png", width = 480, height = 480)

#  MA plot
ggplot(deseq2_results, aes(x = log10(baseMean), y = log2FoldChange, color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("NotSig" = "gray", "Sig" = "skyblue")) +
  labs(title = "MA Plot",
       x = "Log10 Mean expression",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "top")
dev.off()

## Volcano Plot --------------------------------------
# Add a significance column based on an adjusted p-value threshold and fold change
deseq2_results$Significance <- ifelse(deseq2_results$padj < 0.05 & abs(deseq2_results$log2FoldChange) > 1, "Significant", "Not Significant")

# Define color scheme
colors <- c("Not Significant" = "gray", "Significant" = "dodgerblue")


png("volcano_plot.png", width = 480, height = 480)
# Create the volcano plot
ggplot(deseq2_results, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = colors) +
  
  # Add significance threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red", linewidth = 0.5) +
  
  # axis labels and title
  labs(
    title = "Volcano Plot of Differential Expression",
    x = expression("Log"[2]*" Fold Change"),
    y = expression("-Log"[10]*" P-value"),
    color = "Gene Significance"
  ) +
  theme_classic() +
  scale_y_continuous(trans = 'log1p') 
dev.off()


## Pheatmap -------------------------------------------
# Subset to include only the most variable genes (e.g., top 50)
# Get the top 50 genes with the lowest p-values
top_genes <- deseq2_results %>%
  arrange(pvalue) %>%
  head(50) %>%
  pull(ensembl_gene_id)

heatmap_data <- count_table[top_genes, ] #subset the top 50

#change to gene symbols
gene_info <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                                   keys = top_genes,
                                   keytype = "GENEID",
                                   columns = c("GENEID", "SYMBOL"))

#change the rownames to SYMBOLs
heatmap_data <- heatmap_data %>% 
  rownames_to_column("GENEID") %>%
  left_join(gene_info, by = join_by("GENEID")) %>%
  column_to_rownames("SYMBOL") %>%
  dplyr::select(-GENEID)

# Set up annotation for columns (samples)
annotation_col <- table_meta[, c("replicate", "condition")]
rownames(annotation_col) <- table_meta$accession  

#convert columns to factors for vizualtion
annotation_col$replicate <- as.factor(annotation_col$replicate)
annotation_col$condition <- as.factor(annotation_col$condition)

# Define colors for the  annotations
replicate_colors <- c("1" = "plum3", "2" = "#92C5DE", "3" = "#0571B0")
condition_colors <- c("Wt" = "coral1", "KO" = "aquamarine2")

# Create annotation colors
annotation_colors <- list(
  condition = condition_colors,
  replicate = replicate_colors
)

png("heatmap_top_50.png", width = 480, height = 800)

# Generate the heatmap
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,      
  cluster_cols = TRUE,      
  scale = "row",            # Scale each gene (row) for better visualization?
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 10,        
  fontsize_col = 10,         
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,        # Remove borders 
  width = 4,                
  height = 8               
)
dev.off()