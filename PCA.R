library(DESeq2)
library(ggplot2)
library(RColorBrewer)

# Metadata for DESeq2 analysis
# Metadata for DESeq2 analysis
metadata <- data.frame(
  sampleName = c("SRR26296713", "SRR26296714", "SRR26296715", "SRR26296716", "SRR26296717", "SRR26296718", 
                 "SRR26296719", "SRR26296720", "SRR26296721", "SRR26296722", "SRR26296723", "SRR26296724"),
  genotype = c("deltapstS", "deltapstS", "deltaphoB", "deltaphoB", "CB15N", "CB15N", 
               "deltapstS", "deltapstS", "deltaphoB", "deltaphoB", "CB15N", "CB15N"),
  treatment = c("conditional_expression", "conditional_expression", "conditional_expression", "conditional_expression",
                "conditional_expression", "conditional_expression", "not", "not", "not", "not",
                "not", "not"),
  stringsAsFactors = TRUE
)

# Load feature counts
# Adjust the path to match the directory on your system
# Assuming the files have a ".tabular" extension and correcting the path as per your setup
count_files <- list.files(path = "C:/Users/theoa/OneDrive/Documents/SRR lists/FeatureCounts", pattern = "*.tabular", full.names = TRUE)

# Check if files were loaded successfully
if(length(count_files) == 0) {
  stop("No files found. Check the directory path and file pattern.")
}

# Assign names to the files based on the metadata sample names
names(count_files) <- metadata$sampleName

# Read counts
# Read counts, specifying tab delimiter and using the first column as row names
count_list <- lapply(count_files, read.table, header = TRUE, sep="\t", row.names=1)
counts <- do.call(cbind, count_list)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~ genotype + treatment)

# Run the DESeq2 analysis
dds <- DESeq(dds)

# Normalizing counts
norm_counts <- counts(dds, normalized = TRUE)

write.csv(norm_counts, file = "C:/Users/theoa/OneDrive/Documents/Normalized_DESeq2_Results.csv", row.names = TRUE)

# PCA on normalized counts
pca_data <- prcomp(t(norm_counts))

# Calculate percentage of variance explained by PCs
var_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100

# Data frame for ggplot
pca_df <- data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2],
                     Genotype = metadata$genotype, Treatment = metadata$treatment,
                     PC1_Var = sprintf("%.1f%%", var_explained[1]),
                     PC2_Var = sprintf("%.1f%%", var_explained[2]))

# Plotting PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Genotype, shape = Treatment)) +
  geom_point(size = 4) +
  xlab(paste("PC1 - ", pca_df$PC1_Var[1])) +
  ylab(paste("PC2 - ", pca_df$PC2_Var[1])) +
  labs(title = "PCA of DESeq2 Normalized Data") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank(), plot.title = element_text(hjust = 0.5))

# Extracting loadings
loadings <- pca_data$rotation

# Define a cutoff for significant contribution, e.g., top 5% of the absolute loading values
cutoff_pc1 <- quantile(abs(loadings[,1]), 0.95)
cutoff_pc2 <- quantile(abs(loadings[,2]), 0.95)

# Get genes exceeding the cutoff for PC1 and PC2
significant_genes_pc1 <- names(loadings[abs(loadings[,1]) >= cutoff_pc1, 1])
significant_genes_pc2 <- names(loadings[abs(loadings[,2]) >= cutoff_pc2, 2])

# Create a data frame of the significant genes for each component
significant_genes_df <- data.frame(
  Gene_PC1 = significant_genes_pc1,
  Contribution_PC1 = loadings[significant_genes_pc1, 1],
  Gene_PC2 = significant_genes_pc2,
  Contribution_PC2 = loadings[significant_genes_pc2, 2]
)

# Save to CSV file
write.csv(significant_genes_df, "Significant_Contributing_Genes_PC1_PC2.csv", row.names = FALSE)