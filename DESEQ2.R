# Load required packages
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(tximport)
library(apeglm)
library(RColorBrewer)

# -------------------------------
# 1. Load Gene Count Data
# -------------------------------
# Read the gene count matrix (rows: genes, columns: samples)
counts <- read.csv("gene_counts_matrix.csv", row.names=1)

# Read sample metadata (with sample names, conditions, and batch info if needed)
metadata <- read.csv("metadata.csv", row.names=1)

# Ensure the column names of counts match row names of metadata
all(colnames(counts) == rownames(metadata))

# -------------------------------
# 2. Create DESeq2 Dataset
# -------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)  # Change as needed (e.g., ~ batch + condition)

# -------------------------------
# 3. Preprocessing
# -------------------------------

# Remove low-expressed genes (filter genes with very low counts)
dds <- dds[rowSums(counts(dds)) > 10, ]

# Normalize counts using variance stabilization
vsd <- vst(dds, blind=FALSE)

# Log-transformed counts for visualization
rld <- rlog(dds, blind=FALSE)

# -------------------------------
# 4. Run DESeq2 Analysis
# -------------------------------
dds <- DESeq(dds)

# Get results with log2 fold change (LFC) shrinkage
res <- results(dds, contrast=c("condition", "treated", "control"))
resLFC <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm")

# Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Save results as CSV
write.csv(as.data.frame(resOrdered), "DESeq2_gene_results.csv")

# -------------------------------
# 5. PCA Plot for Sample Clustering
# -------------------------------
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

# -------------------------------
# 6. Heatmap of Top 50 Genes
# -------------------------------
top50 <- head(order(res$padj, decreasing=FALSE), 50)
pheatmap(assay(vsd)[top50,], cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=TRUE, annotation_col=metadata,
         color=colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(100))

# -------------------------------
# 7. Volcano Plot of Differentially Expressed Genes
# -------------------------------
EnhancedVolcano(res,
                lab=rownames(res),
                x="log2FoldChange",
                y="pvalue",
                title="Differential Gene Expression",
                pCutoff=0.05,
                FCcutoff=1,
                colAlpha=0.75,
                legendPosition="right")

# -------------------------------
# 8. Isoform-Based Analysis (Using tximport)
# -------------------------------
# Read transcript-level quantifications (Salmon output)
files <- list.files("results/salmon", pattern="quant.sf", full.names=TRUE)
names(files) <- gsub("results/salmon/", "", dirname(files))  # Extract sample names

# Load transcript-to-gene mapping file
tx2gene <- read.csv("transcript_to_gene_mapping.csv")

# Import transcript counts using tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)

# Create DESeq2 dataset
dds_iso <- DESeqDataSetFromTximport(txi, colData=metadata, design=~ condition)

# Run DESeq2 on isoform data
dds_iso <- DESeq(dds_iso)
res_iso <- results(dds_iso, contrast=c("condition", "treated", "control"))
write.csv(as.data.frame(res_iso), "DESeq2_isoform_results.csv")

# -------------------------------
# 9. Save Normalized Counts
# -------------------------------
write.csv(assay(vsd), "normalized_counts.csv")
