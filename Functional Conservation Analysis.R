#===========================================================================
# Functional Conservation Analysis of Chicken-Duck Enhancers
#===========================================================================

#-------------------------
# 1. Load Required Libraries
#-------------------------
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(pheatmap)
library(dplyr)
library(ggpubr)

#-------------------------
# 1. Homologous Gene Expression Analysis
#-------------------------

# Set working directory and read data
setwd("C:/Users/35471/Desktop")
homolog_data <- read.table("chicken_homelogus_gene.txt", header = TRUE)
expression_data <- read.table("chicken_GC_zscore_tpm.txt", header = TRUE, sep = "\t")

# Merge data
merged_data <- merge(homolog_data, expression_data, 
                     by.x = "gene", by.y = "gene_id", all.x = TRUE)

# Preserve original order
merged_data <- merged_data[match(homolog_data$gene, merged_data$gene), ]
rownames(merged_data) <- merged_data$gene
heatmap_matrix <- as.matrix(merged_data[, -1])  # Remove gene column

# Create heatmap using ComplexHeatmap
pdf("Homologous_Gene_Expression_ComplexHeatmap.pdf", width = 6, height = 8)
heatmap_plot <- Heatmap(
  heatmap_matrix,
  name = "Expression",
  km = 4,
  column_names_side = "bottom",
  col = colorRamp2(c(-2, 0, 2), c("#aecdd7", "#ffffff", "#c85d4d")),
  cluster_columns = FALSE,
  row_dend_side = "left",
  show_row_names = FALSE,
  row_title = "Clustered Genes",
  column_title = "Homologous Gene Expression in Chicken"
)
draw(heatmap_plot)
dev.off()

# Extract clusters
ht_list <- draw(heatmap_plot)
row_order_result <- row_order(ht_list)

# Create cluster assignments
cluster_df <- lapply(seq_along(row_order_result), function(i) {
  data.frame(
    GeneID = rownames(heatmap_matrix[row_order_result[[i]], ]),
    Cluster = paste0("cluster", i),
    stringsAsFactors = FALSE
  )
})

cluster_assignments <- do.call(rbind, cluster_df)
write.table(cluster_assignments, "duck_homelogene_order.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------------
#2. KEGG Enrichment Analysis
#-------------------------

# Read gene symbols
symbol_data <- read.table(
  "D:/work/RNA-seq/final_symbols_name.txt",
  header = TRUE, 
  sep = "\t", 
  quote = ""
)

# Merge with cluster data
merged_symbols <- merge(cluster_assignments, symbol_data, 
                        by.x = "GeneID", by.y = "gene_id", all.x = TRUE)

# Perform KEGG enrichment
kegg_enrichment <- compareCluster(
  Gene_id ~ Cluster, 
  data = merged_symbols,
  fun = "enrichKEGG",
  organism = "apla",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Plot KEGG results
kegg_plot <- dotplot(kegg_enrichment, 
                     showCategory = 15,
                     font.size = 12,
                     label_format = 50) +
  scale_color_continuous(low = "#FF3333", high = "#77FFEE") +
  labs(title = "KEGG Pathway Enrichment of Homologous Genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

# Save KEGG plot
ggsave("KEGG_Enrichment_Homologous_Genes.tiff", kegg_plot, 
       width = 8, height = 10, dpi = 900)

# Save KEGG results
kegg_results <- kegg_enrichment@result
write.table(kegg_results, "kegg_results_0.05.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------------
# 3. Conserved Enhancer Analysis
#-------------------------

# Read conserved enhancer data
conserved_data <- read.table("conserved_enh.txt", header = TRUE)

# Set factor levels
conserved_data$type <- factor(
  conserved_data$type,
  levels = c("sequence_not_conserved", 
             "sequence_but_not_usage_conserved", 
             "usage_conserved")
)

# Define colors
conserved_colors <- c("#cadeec", "#a4c3dd", "#6697c3")

# Create stacked bar plot
conserved_plot <- ggplot(conserved_data, aes(x = species, y = number, fill = type)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_fill_manual(values = conserved_colors) +
  labs(title = "Conservation of Enhancers Across Species",
       x = "Species", 
       y = "Proportion") +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Save plot
ggsave("Conserved_Enhancer_Proportions.pdf", conserved_plot, width = 8, height = 6)

#-------------------------
# 4. Transcription Factor Density Analysis
#-------------------------

# Read TF density data (adjust reading method as needed)
# tf_data <- read.table('clipboard', header = T)  # For clipboard data
tf_data <- read.table("tf_density_data.txt", header = TRUE)  # For file data

# Define colors
tf_colors <- c("#C5B0D5", "#f08f92", "#9cbedb", "#a9d5a5", "#FDB96B")

# Create bar plot
tf_plot <- ggplot(tf_data, aes(x = species, y = Density_of_motif, fill = species)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  scale_fill_manual(values = tf_colors) +
  labs(title = "Transcription Factor Motif Density",
       x = "Species", 
       y = "Density of motif/1000kb Enhancers") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  # Add significance comparisons if needed
  geom_signif(
    comparisons = list(c("Chicken", "Duck"), 
                       map_signif_level = TRUE,
                       y_position = max(tf_data$Density_of_motif) * 1.1
    )
    
    # Save plot
    ggsave("TF_Motif_Density.pdf", tf_plot, width = 8, height = 6)