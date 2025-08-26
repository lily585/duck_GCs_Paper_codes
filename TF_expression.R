#===========================================================================
#TF, Enhancer and Silencer Target Gene Expression Analysis
# Author: Zhen Li
#===========================================================================

#-------------------------
# 1. Load Required Libraries
#-------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(ggprism)
library(reshape2)

#-------------------------
# 2. Core TF Gene Expression Analysis
#-------------------------

# Set working directory and read data
setwd("C:/Users/35471/Desktop/follcile/相邻阶段转录网络/差异增强子表达量")
target_genes <- read.table("SYF_F5_5核心TF_gene.txt", header = FALSE)
expression_data <- read.table(
  "C:/Users/35471/Desktop/follcile/rna-seq/follicle_gene_zscore_tpm.txt",
  header = TRUE, 
  sep = "\t"
)

# Merge data
merged_data <- merge(
  target_genes, 
  expression_data, 
  by.x = "V1", 
  by.y = "Gene", 
  all.x = TRUE
)

# Select relevant columns
plot_data <- subset(merged_data, select = c(V1, F3, F1, POF))

# Calculate mean expression
mean_F3 <- mean(plot_data$F3, na.rm = TRUE)
mean_F1 <- mean(plot_data$F1, na.rm = TRUE)
mean_POF <- mean(plot_data$POF, na.rm = TRUE)

# Prepare data for plotting
melted_data <- melt(
  plot_data, 
  id.vars = "V1",
  variable.name = "tissue",
  value.name = "zscore_tpm"
)

# Define colors
tissue_colors <- c("#ea8379", "#8fc9e2", "#b395bd")

# Create violin plot
tf_plot <- ggplot(melted_data, aes(x = tissue, y = zscore_tpm)) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_violin(aes(fill = tissue), outlier.colour = "white", width = 0.8) +
  geom_boxplot(width = 0.15) +
  scale_fill_manual(values = tissue_colors) +
  labs(
    title = "Expression of Core TF Target Genes",
    x = "Developmental Stage", 
    y = "Z-score Normalized TPM"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  stat_compare_means(
    comparisons = list(c("F3", "F1"), c("F1", "POF")),
    label = "p.signif",
    method = "t.test",
    size = 5
  )

# Save plot
ggsave("Core_TF_Target_Expression.png", tf_plot, width = 8, height = 6, dpi = 300)

#-------------------------
# 3. Esrrb Target Gene Expression Analysis
#-------------------------

# Load required library
library(pheatmap)
library(RColorBrewer)

# Set working directory and read data
setwd("C:/Users/35471/Desktop")
overlap_genes <- read.table("F1_F3_ATAC_H3K27AC_overlap.txt", header = TRUE)
expression_data <- read.table(
  "C:/Users/35471/Desktop/follcile/rna-seq/follicle_gene_zscore_tpm.txt",
  header = TRUE, 
  sep = "\t"
)

# Merge data
esrrb_data <- merge(
  overlap_genes, 
  expression_data, 
  by.x = "gene", 
  by.y = "Gene", 
  all.x = TRUE
)

# Select relevant columns
esrrb_plot_data <- subset(esrrb_data, select = c(gene, F3, F1))

# Prepare data for plotting
esrrb_melted <- melt(
  esrrb_plot_data, 
  id.vars = "gene",
  variable.name = "tissue",
  value.name = "zscore_tpm"
)

# Define colors
esrrb_colors <- c("#FDB96B", "#EC5D3B")

# Create violin plot for Esrrb targets
esrrb_plot <- ggplot(esrrb_melted, aes(x = tissue, y = zscore_tpm)) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_violin(aes(fill = tissue), outlier.colour = "white", width = 0.8) +
  geom_boxplot(width = 0.15) +
  scale_fill_manual(values = esrrb_colors) +
  labs(
    title = "Expression of ATAC/H3K27ac Overlap Genes",
    x = "Developmental Stage", 
    y = "Z-score Normalized TPM"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  geom_signif(
    comparisons = list(c("F3", "F1")),
    map_signif_level = TRUE,
    test = t.test,
    tip_length = 0.05,
    size = 0.8,
    color = "black"
  )

# Save plot
ggsave("ATAC_H3K27ac_Overlap_Expression.png", esrrb_plot, width = 7, height = 6, dpi = 300)

