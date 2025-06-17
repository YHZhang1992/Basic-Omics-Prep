# Install DESeq2 and supporting packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("org.Hs.eg.db")      # For gene annotation (change for other species)
BiocManager::install("AnnotationDbi")
install.packages("pheatmap")
install.packages("ggplot2")

# Load libraries
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)

# Load count data (genes as rows, samples as columns)
counts <- read.csv("counts.csv", row.names = 1)

# Load metadata (must include a column for condition or group)
coldata <- read.csv("metadata.csv", row.names = 1)

# Check alignment of sample names
all(colnames(counts) == rownames(coldata))  # Should be TRUE

# Construct DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)  # Replace 'condition' with your variable

# Pre-filter low count genes (optional but recommended)
dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)  # Main DESeq2 function: normalization, dispersion estimation, model fitting

# Transform counts for visualization (VST or rlog)
vsd <- vst(dds, blind = FALSE)

# PCA plot
plotPCA(vsd, intgroup = "condition")

# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists), annotation_col = coldata)

# Get DE results: change condition2 and condition1 as needed
res <- results(dds, contrast = c("condition", "treated", "control"))

# Shrink log fold change (optional but helps visualization)
resLFC <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm")

# View summary
summary(resLFC)

# Order by p-value
resOrdered <- resLFC[order(resLFC$pvalue), ]

# Filter significant genes
resSig <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)

# Assume rownames of resSig are Ensembl IDs, remove version numbers if present
resSig$ensembl <- gsub("\\..*", "", rownames(resSig))

# Map Ensembl to gene symbols
resSig$symbol <- mapIds(org.Hs.eg.db,
                        keys = resSig$ensembl,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# MA plot
plotMA(resLFC, ylim = c(-5, 5))

# Volcano plot
resSig$logP <- -log10(resSig$padj)
ggplot(resSig, aes(x = log2FoldChange, y = logP, label = symbol)) +
  geom_point(alpha = 0.5) +
  geom_text(aes(label=ifelse(logP > 5 & abs(log2FoldChange) > 2, symbol, "")), hjust=1.1, vjust=1.1) +
  theme_minimal() +
  xlab("Log2 Fold Change") + ylab("-log10 adjusted p-value")

# Heatmap of top genes
topGenes <- head(order(resSig$padj), 30)
pheatmap(assay(vsd)[topGenes, ], cluster_rows = TRUE, show_rownames = TRUE,
         cluster_cols = TRUE, annotation_col = coldata)

# Export to CSV
write.csv(as.data.frame(resSig), file = "DEGs_annotated.csv")
