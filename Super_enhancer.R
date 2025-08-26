#===========================================================================
# Super Enhancer

#===========================================================================

#-------------------------
# 1. Load Required Libraries
#-------------------------
library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(DOSE)
library(ComplexHeatmap)
library(circlize)
library(aplot)
library(ggpubr)
library(ggsignif)
library(ggprism)

#-------------------------
# 2. Signal Scatter Plot
#-------------------------

# Set working directory and read data
setwd("C:/Users/35471/Desktop/")
df_gene <- read.table('signal_SS.txt', header = TRUE, sep = '\t')

# Prepare data for plotting
melted_data <- melt(df_gene, 
                    id.vars = "RANK",
                    variable.name = "group",
                    value.name = "signal")

# Define color palette
tissue_colors <- c("#3951A2","#72AACF","#CAE8F2","#FDB96B","#EC5D3B","#A80326")

# Create scatter plot
signal_plot <- ggplot(melted_data, aes(x = RANK, y = signal, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = tissue_colors) +
  theme_classic() +
  labs(title = "Signal Distribution by Rank",
       x = "Rank", 
       y = "Signal Intensity") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Save plot
ggsave("SS_signal.pdf", signal_plot, device = cairo_pdf, width = 7, height = 6)

#-------------------------
# 3. SE and SS Count Bar Plot
#-------------------------

# Read data
count_data <- read.table("SE_SS.number.txt", header = TRUE, sep = "\t")

# Set factor levels
count_data$State <- factor(count_data$State, 
                           levels = c("SWF", "LWF", "SYF", "F5", "F3", "F1", "POF"))

# Define colors
count_colors <- c("#E6846d", "#8dcdd5")

# Create stacked bar plot
count_plot <- ggplot(count_data, aes(x = State, y = Number, fill = Type)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = count_colors) +
  labs(title = "SE and SS Counts Across Tissues",
       x = "Tissue", 
       y = "Count") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "bottom"
  )

# Save plot
ggsave("SE_SS_counts.pdf", count_plot, width = 8, height = 6)

#-------------------------
# 4. Enhancer/Silencer Target Gene Expression
#-------------------------

# Read and combine data
setwd("C:/Users/35471/Desktop/")
df1 <- read.table('merged_SS_SE.txt', sep = '\t')
df2 <- read.table('total.silencer_gene.txt', sep = '\t')
df3 <- read.table('total_enhancer_gene.txt', sep = '\t')
combined_df <- rbind(df1, df2, df3)

# Read expression data
expression_data <- read.table(
  'C:/Users/35471/Desktop/follcile/rna-seq/follicle_gene_zscore_tpm.txt',
  header = TRUE, 
  sep = '\t'
)

# Merge data
merged_data <- merge(combined_df, expression_data, 
                     by.x = "V2", by.y = "Gene", all.x = TRUE)
merged_data <- na.omit(merged_data)

# Prepare for plotting
melted_data <- melt(merged_data, 
                    id.vars = c("V2", "V1", "V3", "V4"),
                    variable.name = "tissue",
                    value.name = "zscore_tpm")

# Filter matching tissues
filtered_data <- melted_data %>% filter(V3 == tissue)

# Define colors
enhancer_colors <- c("#f08f92", "#9cbedb", "#a9d5a5", "#FDB96B")

# Create expression comparison plot
expression_plot <- ggplot(filtered_data, aes(x = V4, y = zscore_tpm, fill = V4)) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_violin(alpha = 0.8, width = 0.8, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = enhancer_colors) +
  labs(title = "Expression of Enhancer/Silencer Target Genes",
       x = "Enhancer Type", 
       y = "Z-score Normalized TPM") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  geom_signif(
    comparisons = list(
      c("TE", "SE"),
      c("TS", "SS"),
      c("TE", "TS"),
      c("SE", "SS")
    ),
    map_signif_level = TRUE,
    test = "t.test",
    tip_length = 0.01,
    size = 0.8,
    textsize = 4,
    y_position = c(5, 6, 7, 8) # Adjust based on your data range
  )

# Save plot
ggsave("Enhancer_Silencer_Expression.pdf", expression_plot, width = 10, height = 8)


#-------------------------
# 5. Tissue Specificity Index (TSI) Analysis
#-------------------------

# Read expression data
expression_data <- read.table(
  "C:/Users/35471/Desktop/follcile/rna-seq/GC_tpm_no_syf.txt",
  header = TRUE, 
  row.names = 1
)

# Calculate TSI
TSI <- apply(expression_data, 1, function(row) {
  max_val <- max(row)
  sum_val <- sum(row)
  max_val / sum_val
})

# Create TSI data frame
tsi_df <- data.frame(Gene = names(TSI), TSI = TSI)

# Save TSI results
write.csv(tsi_df, "gene_TSI_results.csv", row.names = FALSE)

# Plot TSI distribution
tsi_plot <- ggplot(tsi_df, aes(x = TSI)) + 
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_bw() +
  labs(title = "Distribution of Tissue Specificity Index (TSI)",
       x = "TSI", 
       y = "Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

# Save TSI plot
ggsave("TSI_Distribution.pdf", tsi_plot, width = 8, height = 6)

#-------------------------
# 6. FOXO1 Target Gene Expression
#-------------------------

# Read FOXO1 target data
foxo1_data <- read.table('FOXO1.txt', sep = '\t', header = FALSE)

# Read expression data
expression_data <- read.table(
  'C:/Users/35471/Desktop/follcile/rna-seq/follicle_gene_zscore_tpm.txt',
  header = TRUE, 
  sep = '\t'
)

# Merge data
merged_foxo1 <- merge(foxo1_data, expression_data, 
                      by.x = "V1", by.y = "Gene", all.x = TRUE)

# Prepare for heatmap
rownames(merged_foxo1) <- merged_foxo1$V2
heatmap_matrix <- as.matrix(merged_foxo1[, -c(1:2)])  # Remove ID columns

# Create heatmap
foxo1_heatmap <- pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#64A2CE", "white", "#F0282D"))(100),
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  border_color = NA,
  main = "FOXO1 Target Gene Expression"
)

# Save heatmap
ggsave("FOXO1_Target_Expression.png", foxo1_heatmap, width = 8, height = 10, dpi = 300)