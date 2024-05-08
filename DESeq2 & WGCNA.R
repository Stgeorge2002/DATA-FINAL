# Load necessary libraries
library(DESeq2)
library(WGCNA)
library(data.table)
library(boot)

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
  stringsAsFactors = TRUE
)

# Reading and preprocessing featureCounts files
files <- list.files(pattern="*.tabular", full.names=TRUE)
names <- gsub(".tabular", "", basename(files))

data_list <- lapply(files, function(file) {
  dt <- fread(file, header = TRUE)
  setnames(dt, old = colnames(dt)[2], new = gsub(".tabular", "", basename(file)))
  dt
})

combined_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), data_list)
rownames(combined_counts) <- combined_counts$Geneid
combined_counts$Geneid <- NULL  # Remove the Geneid column

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = combined_counts, colData = metadata, design = ~ genotype + treatment)
dds <- DESeq(dds)
res <- results(dds)

# Save DESeq2 results
write.csv(as.data.frame(res), file="DESeq2_results.csv")

# Prepare data for WGCNA
filtered_counts <- combined_counts[rowSums(combined_counts) >= 10, ]
datExpr <- t(as.matrix(filtered_counts))  # Transpose for WGCNA
datExpr <- apply(datExpr, 2, as.numeric)  # Convert each column to numeric

# Check if all values in datExpr are numeric
if (!all(sapply(datExpr, is.numeric))) {
  stop("Non-numeric values found in datExpr. Please check the data.")
}

# Define the vector of powers explicitly and simply
powers <- as.numeric(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20))

# Selecting soft-threshold power
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate

# Plot the scale-free topology fit and mean connectivity
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.90, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = "red")

# Network construction ensuring all inputs are numeric
net <- blockwiseModules(datExpr, power = softPower, maxBlockSize = 5000,
                        TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE)

# Plot the dendrogram and module colors
plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Bootstrap analysis for module stability
nBootstrap <- 10
bootResults <- boot(datExpr, statistic = function(data, indices) {
  sampleData <- data[indices,]
  net <- blockwiseModules(sampleData, power = softPower, maxBlockSize = 5000,
                          TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE)
  return(net$colors)
}, R = nBootstrap)

# Save Bootstrap results for further analysis
save(bootResults, file = "WGCNA_bootstrap_results.RData")
