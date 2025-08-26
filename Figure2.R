#Chromatin State Analysis
# Author: Zhen Li

# Load required libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(RColorBrewer)

#-------------------------
# 1. PCA Plot
#-------------------------

# Read data
pca_data <- read.table("atac_PCA_total.txt", header = TRUE, sep = "\t", row.names = 1)
group_info <- read.table("atac_list.txt", header = TRUE, sep = "\t", row.names = 1)

# Prepare data
pca_df <- as.data.frame(t(pca_data[1:2, ]))
colnames(pca_df) <- c("PC1", "PC2")

# Define color palette and factor levels
tissue_colors <- c("#3951A2", "#72AACF", "#CAE8F2", "#FEFBBA", "#FDB96B", "#EC5D3B", "#A80326")
group_info$Tissue <- factor(group_info$Tissue, 
                            levels = c('SWF', 'LWF', 'SYF', 'F5', 'F3', "F1", "POF"))

# Create PCA plot
pca_plot <- ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group_info$Tissue), size = 3) +
  scale_color_manual(name = "Tissue", values = tissue_colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = "black"), 
    axis.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    panel.grid = element_blank()
  ) +
  labs(
    title = "ATAC-Seq Principal Component Analysis",
  )

# Save plot
ggsave("ATAC_PCA_plot.png", pca_plot, width = 10, height = 8, dpi = 300)

#-------------------------
# 2. Stacked Bar Plot (Genome Coverage)
#-------------------------

# Read data
coverage_data <- read.table("genome_coverage.txt", header = TRUE, sep = "\t")

# Prepare data
melted_df <- melt(coverage_data, 
                  id.vars = "State", 
                  variable.name = "tissue",
                  value.name = "value")

# Define chromatin state colors and levels
chromatin_colors <- c(
  "TssA" = "#FF0000", "TssAHet" = "#E54C5E",
  "TxFlnk" = "#006400", "TxFlinkWk" = "#588E31",
  "TxWk" = "#ADD88D", "TxFLInkHet1" = "#C9E4B4",
  "TxFLInkHet2" = "#C9E4B4", "EnhA" = "#FFFF00",
  "EnhAmMe" = "#EE822F", "EnhAWk" = "#F4B382",
  "EnhPois" = "#F9CBAA", "TssBiv" = "#CD5C5C",
  "Repr" = "#808080", "reprwk" = "#C0C0C0", "Qui" = "#EAEAEA"
)

chromatin_levels <- c(
  'TssA', 'TssAHet', 'TxFlnk', 'TxFlinkWk', 'TxWk', 
  "TxFLInkHet1", "TxFLInkHet2", "EnhA", "EnhAmMe", 
  "EnhAWk", "EnhPois", "TssBiv", "Repr", "reprwk", "Qui"
)

melted_df$State <- factor(melted_df$State, levels = chromatin_levels)

# Create stacked bar plot
stacked_plot <- ggplot(data = melted_df, aes(x = tissue, y = value, fill = State)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = chromatin_colors) +
  theme(
    panel.background = element_blank(),        
    panel.grid.major = element_blank(),        
    panel.grid.minor = element_blank(),        
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),        
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),        
    axis.text.y = element_text(size = 12, colour = "black"),        
    axis.title.y = element_text(size = 14)
  ) +
  labs(x = NULL, y = "Genome Coverage (%)") +
  guides(fill = guide_legend(title = "Chromatin State", nrow = 3))

# Save plot
ggsave("Genome_Coverage_Stacked.png", stacked_plot, width = 12, height = 8, dpi = 300)

#-------------------------
# 3. Enhancer Target Genes Expression
#-------------------------

# Read data
enhancer_genes <- read.table('EnhA_gene.txt', sep = '\t', header = FALSE)
expression_data <- read.table('average_tpm.txt', header = TRUE, sep = '\t')

# Merge datasets
merged_data <- merge(enhancer_genes, expression_data, 
                     by.x = "V1", by.y = "gene_id", all.x = TRUE)

# Create expression categories
merged_data$category <- cut(merged_data$V2, 
                            breaks = c(-Inf, 3, 10, 21, 50, Inf), 
                            labels = c("Low", "NA", "Medium", "NA", "High"),
                            right = FALSE)

# Filter data
filtered_data <- filter(merged_data, category != "NA")

# Define colors
expr_colors <- c("#f08f92", "#9cbedb", "#a9d5a5")

# Create boxplot
expr_plot <- ggplot(filtered_data, aes(x = category, y = log)) +
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.8) +
  geom_boxplot(aes(fill = category), outlier.colour = "white", size = 0.8) +
  scale_fill_manual(values = expr_colors) +
  theme(
    panel.background = element_blank(),        
    panel.grid.major = element_blank(),        
    panel.grid.minor = element_blank(),        
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    legend.position = "none",        
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),        
    axis.text.x = element_text(size = 14, colour = "black"),        
    axis.text.y = element_text(size = 14, colour = "black"),        
    axis.title.y = element_text(size = 16)
  ) +
  labs(
    title = "Expression of EnhAmMe Target Genes",
    x = "Expression Level", 
    y = "Log(TPM)"
  ) +
  geom_signif(
    comparisons = list(c("Low", "Medium"), c("Medium", "High"), c("Low", "High")),
    map_signif_level = TRUE,
    test = t.test,
    tip_length = 0.05,
    size = 0.8,
    color = "black"
  )

# Save plot
ggsave("Enhancer_Target_Expression.png", expr_plot, width = 10, height = 8, dpi = 300)

#-------------------------
# 4. FRIP Analysis
#-------------------------

# Read data
frip_data <- read.table("follicle_FRIP.txt", header = TRUE, sep = "\t", row.names = 1)

# Set factor levels
frip_data$Tissue <- factor(frip_data$Tissue, 
                           levels = c('SWF', 'LWF', 'SYF', 'F5', 'F3', "F1", "POF"))

# Create bar plot
frip_plot <- ggplot(frip_data, aes(x = Factor, y = FRiP, fill = Tissue)) + 
  geom_col(position = position_dodge(0.8), width = 0.5) +
  scale_fill_manual(values = tissue_colors) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  labs(x = "Factor", y = "FRiP (Fraction of Reads in Peaks)")

# Save plot
ggsave("FRIP_Analysis.png", frip_plot, width = 10, height = 6, dpi = 300)

#-------------------------
# 5. Mapped Ratio Analysis
#-------------------------

# Read data
mapped_data <- read.table("mapped_ratio.txt", header = TRUE, sep = "\t")

# Define color palette and factor levels
data_type_colors <- c("#C47273", "#4387B5", "#A9D179", "#F5AE6B", "#BCB8D3", "#84CAC0")
mapped_data$Type <- factor(mapped_data$Type, 
                           levels = c('ATAC', 'H3K27AC', 'H3K4ME3', 'H3K4ME1', 'H3K27ME3', "RNA"))

# Create violin plot
mapped_plot <- ggplot(mapped_data, aes(x = Type, y = Mapped_ratio)) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_violin(aes(fill = Type), outlier.colour = "white", cex = 1.2, width = 0.8) +
  geom_boxplot(width = 0.15, cex = 1.2) +
  scale_fill_manual(values = data_type_colors) +
  theme(
    panel.background = element_blank(),        
    panel.grid.major = element_blank(),        
    panel.grid.minor = element_blank(),        
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    legend.position = "none",        
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),        
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),        
    axis.text.y = element_text(size = 12, colour = "black"),        
    axis.title.y = element_text(size = 14)
  ) +
  labs(x = "Data Type", y = "Mapped Ratio (%)")

# Save plot
ggsave("Mapped_Ratio_Violin.png", mapped_plot, width = 10, height = 8, dpi = 300)