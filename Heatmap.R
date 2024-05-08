# Load necessary libraries
library(DESeq2)
library(WGCNA)
library(data.table)
library(boot)
library(pheatmap)

# Set the working directory to where the featureCounts files are located
setwd("C:/Users/theoa/OneDrive/Documents/SRR lists/Featurecounts")

# Metadata for DESeq2 analysis
metadata <- data.frame(
  sampleName = c("SRR26296713", "SRR26296714", "SRR26296715", "SRR26296716", "SRR26296717", "SRR26296718", 
                 "SRR26296719", "SRR26296720", "SRR26296721", "SRR26296722", "SRR26296723", "SRR26296724"),
  genotype = c("deltapstS", "deltapstS", "deltaphoB", "deltaphoB", "CB15N", "CB15N", 
               "deltapstS", "deltapstS", "deltaphoB", "deltaphoB", "CB15N", "CB15N"),
  treatment = c("conditional_expression", "conditional_expression", "conditional_expression", "conditional_expression",
                "conditional_expression", "conditional_expression", "not", "not", "not", "not",
                "not", "not"),
  stringsAsFactors = FALSE
)

# Reading and preprocessing featureCounts files
files <- list.files(pattern="*.tabular", full.names=TRUE)
data_list <- lapply(files, function(file) {
  dt <- fread(file, header = TRUE)
  setnames(dt, old = colnames(dt)[2], new = gsub(".tabular", "", basename(file)))
  dt
})

# Combine all data using Geneid
combined_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), data_list)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = combined_counts[, -1, with = FALSE], colData = metadata, design = ~ genotype + treatment)
rownames(dds) <- combined_counts$Geneid

# Run DESeq analysis
dds <- DESeq(dds)

# Extract and write normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(normalized_counts), file = "Normalized_DESeq_counts.csv")

# Creating annotations for samples
sample_annotations <- data.frame(
  Genotype = metadata$genotype,
  Treatment = metadata$treatment
)
rownames(sample_annotations) <- metadata$sampleName
col_ann_colors <- list(
  Genotype = c(deltapstS = "blue", deltaphoB = "red", CB15N = "green"),
  Treatment = c(conditional_expression = "orange", not = "purple")
)

# Extracting top 30 significant genes for each comparison and creating heatmaps
create_top_genes_heatmap <- function(dds, comparison_name, annotations, ann_colors) {
  res <- results(dds)
  # Ordering by p-value and filtering top 30
  top_genes <- head(res[order(res$pvalue), ], 30)
  top_gene_data <- normalized_counts[rownames(normalized_counts) %in% rownames(top_genes), ]
  
  # Create heatmap
  pheatmap(log2(top_gene_data + 1),
           clustering_distance_rows = "euclidean", 
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           color = colorRampPalette(c("navy", "white", "firebrick"))(50),
           show_rownames = FALSE, 
           show_colnames = TRUE,
           annotation_col = annotations,
           annotation_colors = ann_colors,
           main = paste("Top 30 Genes Heatmap:", comparison_name))
}

# Analyzing and generating heatmap for each comparison
# Conditional Expression vs Not
dds_cond_exp_vs_not <- dds[dds$treatment %in% c("conditional_expression", "not"), ]
dds_cond_exp_vs_not <- DESeq(dds_cond_exp_vs_not)
create_top_genes_heatmap(dds_cond_exp_vs_not, "Conditional Expression vs Not", sample_annotations, col_ann_colors)

# deltapstS vs CB15N
dds_deltapstS_vs_CB15N <- dds[dds$genotype %in% c("deltapstS", "CB15N"), ]
dds_deltapstS_vs_CB15N <- DESeq(dds_deltapstS_vs_CB15N)
create_top_genes_heatmap(dds_deltapstS_vs_CB15N, "deltapstS vs CB15N", sample_annotations, col_ann_colors)

# deltaphoB vs CB15N
dds_deltaphoB_vs_CB15N <- dds[dds$genotype %in% c("deltaphoB", "CB15N"), ]
dds_deltaphoB_vs_CB15N <- DESeq(dds_deltaphoB_vs_CB15N)
create_top_genes_heatmap(dds_deltaphoB_vs_CB15N, "deltaphoB vs CB15N", sample_annotations, col_ann_colors)
