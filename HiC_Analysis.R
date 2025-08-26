#===========================================================================
# Hi-C Data Analysis: Compartments, TADs, and Chromatin Interactions

#===========================================================================

#-------------------------
# 1. Load Required Libraries
#-------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(reshape2)
library(clusterProfiler)
library(DOSE)
library(AnnotationHub)
library(rtracklayer)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

#-------------------------
# 2. A/B Compartment Analysis
#-------------------------

# Set working directory and read data
setwd("D:/work/hic/duck-F1/AB")

# Read PC1 data for each tissue
read_pc1 <- function(file) {
  data <- read.table(file)
  colnames(data) <- c("peakID", "chr", "start", "end", "strand", "PC1")
  return(data)
}

f1_pc1 <- read_pc1("F1.40k.PC1.txt")
f5_pc1 <- read_pc1("F5_40k.PC1.txt")
syf_pc1 <- read_pc1("SYF_40k.PC1.txt")

# Merge PC1 data
merge_pc1 <- function(df1, df2, name1, name2) {
  merged <- merge(df1, df2, by = c("peakID", "chr", "start", "end", "strand"))
  colnames(merged)[6:7] <- c(name1, name2)
  return(merged)
}

pc1_f1_f5 <- merge_pc1(f1_pc1, f5_pc1, "F1", "F5")
pc1_syf_f5 <- merge_pc1(syf_pc1, f5_pc1, "SYF", "F5")
pc1_data <- merge(pc1_f1_f5, pc1_syf_f5, 
                  by = c("peakID", "chr", "start", "end", "strand", "F5"))

# Assign compartment labels
assign_compartment <- function(pc_value) {
  ifelse(pc_value < 0, "B", "A")
}

pc1_data <- pc1_data %>%
  mutate(
    SYF_AB = assign_compartment(SYF),
    F5_AB = assign_compartment(F5),
    F1_AB = assign_compartment(F1)
  )

# Save results
write.table(pc1_data, "ab_compartments.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Calculate compartment transitions
pc1_trans <- pc1_data %>%
  mutate(ABtrans = paste0(SYF_AB, F5_AB, F1_AB)) %>%
  group_by(ABtrans) %>%
  summarise(n = n()) %>%
  mutate(percentage = round(n / sum(n), 2),
         label = paste0(ABtrans, " (", percentage * 100, "%)"))

# Save transitions
write.table(pc1_trans, "ab_transitions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot compartment transitions
mycol <- c("#BD6263", "#4387B5", "#8EA325", "#A9D179", 
           "#84CAC0", "#F5AE6B", "#BCB8D3", "#91D1C2")

pc1_trans$ABtrans <- factor(pc1_trans$ABtrans,
                            levels = c("AAA", "AAB", "ABA", "ABB", 
                                       "BAA", "BAB", "BBA", "BBB"))

trans_plot <- ggplot(pc1_trans, aes(x = "", y = percentage, fill = ABtrans)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "A/B Compartment Transitions",
       x = "Transition Type",
       y = "Percentage") +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16))

# Save plot
ggsave("AB_Compartment_Transitions.pdf", trans_plot, width = 10, height = 6)

#-------------------------
# 3. Compartment Transition Gene Expression
#-------------------------

# Read expression data
expr_data <- read.table("AB_TPM_long.txt", header = TRUE, sep = "\t")

# Set factor levels
expr_data$type <- factor(expr_data$type, 
                         levels = c("Stable", "AB", "BA"))

# Create boxplot
expr_plot <- ggplot(expr_data, aes(x = type, y = fold, fill = type)) +
  facet_wrap("~tissue") +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("Stable", "AB"), c("Stable", "BA")),
    map_signif_level = TRUE, 
    test = t.test,
    y_position = c(10, 11, 12),
    tip_length = 0
  ) +
  scale_fill_manual(values = c("#8582BD", "#4F99C9", "#A8D3A0")) +
  labs(title = "Gene Expression in Compartment Transitions",
       x = "Transition Type",
       y = "log2(TPM fold change)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, face = "bold")
  )

# Save plot
ggsave("Compartment_Transition_Expression.pdf", expr_plot, width = 12, height = 8)

#-------------------------
# 4. TAD Size and Interaction Analysis
#-------------------------

# Set working directory
setwd("D:/work/hic/tad/tad_size")

# Read TAD size data
tad_size <- read.table("TAD_size.txt", header = TRUE)
melted_size <- melt(tad_size, variable.name = "tissue", value.name = "size")

# Plot TAD size
tad_size_plot <- ggplot(melted_size, aes(x = tissue, y = size, fill = tissue)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_violin(alpha = 0.7) +
  stat_compare_means(comparisons = list(c("SYF", "F5"), c("F5", "F1"), c("SYF", "F1"))) +
  scale_fill_manual(values = c("#8582BD", "#4F99C9", "#A8D3A0")) +
  labs(title = "TAD Size Distribution",
       x = "Tissue",
       y = "TAD Size (bp)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

# Save plot
ggsave("TAD_Size_Distribution.pdf", tad_size_plot, width = 8, height = 6)

# Read TAD interaction data
tad_interaction <- read.table("TAD_value.txt", header = TRUE)
melted_interaction <- melt(tad_interaction, variable.name = "tissue", value.name = "value")

# Plot TAD interactions
tad_interaction_plot <- ggplot(melted_interaction, aes(x = tissue, y = value, fill = tissue)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_violin(alpha = 0.7) +
  stat_compare_means(comparisons = list(c("SYF", "F5"), c("F5", "F1"), c("SYF", "F1"))) +
  scale_fill_manual(values = c("#8582BD", "#4F99C9", "#A8D3A0")) +
  labs(title = "Inter-TAD Interaction Strength",
       x = "Tissue",
       y = "Interaction Value") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

# Save plot
ggsave("TAD_Interaction_Strength.pdf", tad_interaction_plot, width = 8, height = 6)

#-------------------------
# 5. TAD-associated Gene Expression
#-------------------------

# Set working directory
setwd("C:/Users/35471/Desktop/")

# Read TAD differential genes
tad_syf_f5 <- read.table("D:/work/hic/tad/SYF_F5_diff_gene.txt", header = TRUE)
tad_f1_f5 <- read.table("D:/work/hic/tad/F1_F5_DIFF_gene.txt", header = TRUE)

# Read expression data
expr_data <- read.table("follicle_gene_zscore_tpm.txt", header = TRUE, sep = "\t")

# Merge and prepare data
prepare_tad_expr <- function(tad_data, expr_data, cols) {
  merged <- merge(tad_data, expr_data, by.x = "gene", by.y = "Gene", all.x = TRUE)
  merged <- merged %>% select(all_of(cols))
  melted <- melt(merged, id.vars = "gene", variable.name = "state", value.name = "tpm")
  return(melted)
}

f1_f5_expr <- prepare_tad_expr(tad_f1_f5, expr_data, c("gene", "F5", "F1"))
syf_f5_expr <- prepare_tad_expr(tad_syf_f5, expr_data, c("GENE", "SYF", "F5"))

# Plot expression
plot_tad_expr <- function(data, title) {
  ggplot(data, aes(x = state, y = tpm, fill = state)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot() +
    stat_compare_means(comparisons = list(c("F5", "F1")), label = "p.signif") +
    scale_fill_manual(values = c("#E6846d", "#8dcdd5")) +
    labs(title = title,
         x = "State",
         y = "Z-score Normalized TPM") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
}

f1_f5_plot <- plot_tad_expr(f1_f5_expr, "F5 vs F1 TAD-associated Genes")
syf_f5_plot <- plot_tad_expr(syf_f5_expr, "SYF vs F5 TAD-associated Genes")

# Save plots
ggsave("F5_F1_TAD_Gene_Expression.pdf", f1_f5_plot, width = 6, height = 6)
ggsave("SYF_F5_TAD_Gene_Expression.pdf", syf_f5_plot, width = 6, height = 6)

#-------------------------
# 6. Enhancer-Gene Distance Analysis
#-------------------------

# Set working directory
setwd("C:/Users/35471/Desktop/follcile/HIC/loop")

# Read distance data
dist_data <- read.table("enhancer_distance.txt", header = TRUE, sep = "\t")

# Plot distance distribution
dist_plot <- ggplot(dist_data, aes(x = distance, fill = type)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("#E6846d", "#8dcdd5")) +
  labs(title = "Enhancer-Gene Distance Distribution",
       x = "Distance (Kb)",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

# Save plot
ggsave("Enhancer_Gene_Distance.pdf", dist_plot, width = 8, height = 6)

