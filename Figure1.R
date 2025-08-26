### RNA-seq Analysis Pipeline ###
# Author: Zhen Li

# Description: RNA-seq analysis including PCA, clustering, DEG analysis, and enrichment

# ----------------------------
# Section 1: Load Libraries
# ----------------------------
rm(list = ls())
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(apeglm)
library(dplyr)
library(rio)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(DOSE)
library(AnnotationHub)
library(AnnotationDbi)
library(rtracklayer)
library(ChIPseeker)
library(GenomicFeatures)
library(txdbmaker)
library(ggprism)
library(Mfuzz)

# ----------------------------
# Section 2: PCA and Heatmap
# ----------------------------

# Set working directory
setwd("~/project/RNA-seq")

# Read count data and sample information
dat <- read.table("GC_count.txt", header = TRUE, row.names = 1)
group <- read.table("follicle_rna_group.txt", sep = "\t")

# Prepare count matrix
sample_names <- group$V1
countdata <- dplyr::select(dat, all_of(sample_names))

# Filter genes with zero counts
countdata <- countdata[rowSums(countdata) > 0, ]

# Create condition vector
condition <- factor(group$V2)
coldata <- data.frame(row.names = colnames(countdata), condition)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countdata, 
                              colData = coldata, 
                              design = ~ condition)

# Variance stabilizing transformation
rld <- vst(dds, blind = FALSE)

# Generate PCA plot
pdf("PCA_follicle.pdf")
plotPCA(rld, intgroup = "condition")
dev.off()

# Custom PCA plot
pcadata <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
pcadata$group <- factor(pcadata$group, levels = c('SWF', 'LWF', 'SYF', 'F5', 'F3', "F1", "POF"))

p <- ggplot(pcadata, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group), size = 3) +
  scale_color_manual(values = c("#3951A2", "#72AACF", "#CAE8F2", "#FEFBBA", "#FDB96B", "#EC5D3B", "#A80326")) +
  theme_bw() +
  labs(title = "RNA-seq PCA", x = "PC1", y = "PC2") +
  theme(panel.grid = element_blank())
print(p)

# Correlation heatmap
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

ann_colors <- list(condition = c(SWF = "#3951A2", LWF = "#72AACF", SYF = "#CAE8F2",
                                 F5 = "#FEFBBA", F3 = "#FDB96B", F1 = "#EC5D3B", POF = "#A80326"))

pheatmap(rld_cor, annotation = coldata,
         annotation_colors = ann_colors,
         show_rownames = FALSE, show_colnames = FALSE)

# ----------------------------
# Section 3: DEG Analysis
# ----------------------------

# Set working directory
setwd("~/project/RNA-seq/DEG")

# Get DEG files and merge
file_list <- list.files(pattern = "*.csv")
all_genes <- data.frame()

for (file in file_list) {
  data <- read.csv(file, stringsAsFactors = FALSE)
  filtered_data <- data[abs(data$log2FoldChange) > 1 & data$padj < 0.05, ]
  genes <- filtered_data[[1]]
  all_genes <- rbind(all_genes, data.frame(gene = genes))
}

unique_genes <- unique(all_genes$gene)
write.csv(unique_genes, "filtered_DEG.csv", row.names = FALSE)

# Merge with TPM data
df_wide <- read.table('GC_tpm_no_syf.txt', header = TRUE, sep = '\t')
df_gene <- read.csv('filtered_DEG.csv', header = TRUE)
merged_df <- merge(df_gene, df_wide, by = "gene_id", all.x = TRUE)
write.csv(merged_df, "DEG_tpm.csv", row.names = FALSE)

# ----------------------------
# Section 4: Z-score Normalization and Clustering
# ----------------------------
setwd("~/project/RNA-seq")

# Z-score normalization
df_wide <- read.table('GC_tpm_no_syf.txt', header = TRUE, row.names = 1, sep = '\t')
rowmean <- apply(df_wide, 1, mean)
rowsd <- apply(df_wide, 1, sd)
data1 <- sweep(df_wide, 1, rowmean)
data2 <- sweep(data1, 1, rowsd, '/')
write.table(data2, "follicle_DEG_zscore_tpm.txt", sep = "\t")

# K-means clustering
z_df <- read.table('follicle_DEG_zscore_tpm.txt', header = TRUE, sep = '\t')
rownames(z_df) <- z_df[, 1]
data_scale <- z_df[, -1]
data_scale1 <- as.matrix(na.omit(data_scale))

pdf("kmeans_8cluster.pdf", width = 4, height = 6)
p <- Heatmap(data_scale1, name = "z-score", 
             km = 8,
             column_names_side = "bottom",
             col = colorRamp2(c(-2, 0, 2), c("#808080", "#ffffff", "red")),
             cluster_columns = FALSE,
             show_row_names = FALSE)
dev.off()

# Extract cluster information
ht_list <- draw(p)
row_order_result <- row_order(ht_list)
clu_df <- lapply(1:length(row_order_result), function(i) {
  data.frame(GeneID = rownames(data_scale1[row_order_result[[i]], ]),
             Cluster = paste0("cluster", i))
})
cluster_assignments <- do.call(rbind, clu_df)
write.table(cluster_assignments, "DEG_kmeans-8_order.txt", sep = "\t", quote = FALSE)

# ----------------------------
# Section 5: Enrichment Analysis
# ----------------------------

# KEGG enrichment
data_scale <- read.table("DEG_kmeans-8.txt", header = TRUE, sep = "\t")
duck_symbols <- read.table("duck_symbols.txt", header = TRUE, sep = "\t")
df <- merge(data_scale, duck_symbols, by.x = "Gene", by.y = "gene_id", all.x = TRUE)

xx <- compareCluster(GeneID ~ cluster, 
                     data = df,
                     fun = "enrichKEGG",
                     organism = "apla",
                     pvalueCutoff = 0.05)
p1 <- dotplot(xx, showCategory = 10, label_format = 50)
p1 + scale_color_continuous(low = "#FF3333", high = "#77FFEE")

# GO enrichment
txdb <- loadDb("Anas_platyrhynchos.orgDb")
data_GO <- compareCluster(GeneID ~ cluster, 
                          data = df, 
                          fun = "enrichGO", 
                          OrgDb = txdb,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)

data_GO_sim <- simplify(data_GO, cutoff = 0.7)
p1 <- dotplot(data_GO_sim, showCategory = 5)
p1 + scale_color_continuous(low = "#FF3333", high = "#77FFEE")


# ----------------------------
# Section 6: Statistical Tests
# ----------------------------

# Fisher's exact test for DEG ratios
data <- read.table("DEG_ratio.txt", header = TRUE)
pairs <- combn(data$GROUP, 2, simplify = FALSE)
p_values <- sapply(pairs, function(pair) {
  subset <- data[data$GROUP %in% pair, ]
  fisher.test(as.matrix(subset[, c("DEG", "Non_DEG")]))$p.value
})
adjusted_p <- p.adjust(p_values, method = "BH")
significant_results <- data.frame(
  Comparison = sapply(pairs, paste, collapse = " vs "),
  Adjusted_P = adjusted_p
)
write.csv(significant_results, "DEG_significance.csv", row.names = FALSE)