#===========================================================================
# Differential Chromatin Accessibility Analysis
#===========================================================================

#-------------------------
# 1. DiffBind Analysis
#-------------------------

# Load required libraries
library(DiffBind)
library(edgeR)
library(BiocParallel)

# Register parallel processing (using serial for stability)
register(SerialParam())

# Set working directory
setwd("/storage-05/poultrylab5/workspace/epi_duck/follicle_result")

# Create DBA object
dbObj <- dba(sampleSheet = "atac_diffibind.csv", minOverlap = 2)

# Perform read counting
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE, minOverlap = 2)

# Normalize using library size factors
dbObj <- dba.normalize(dbObj, normalize = DBA_NORM_LIB)

# Set up contrasts
dbObj <- dba.contrast(dbObj, minMembers = 2)

# Run differential analysis
dbObj <- dba.analyze(dbObj, method = DBA_ALL_METHODS)

# Summary of results
dba.show(dbObj, bContrasts = TRUE)

#-------------------------
# 2. Visualization
#-------------------------

# PCA plot
pdf("plotPCA.pdf")
dba.plotPCA(dbObj, label = DBA_TREATMENT)
dev.off()


#-------------------------
# 3. Export Results
#-------------------------

# Export all contrast results
for (i in 1:length(dbObj$contrasts)) {
  report <- dba.report(dbObj, method = DBA_DESEQ2, contrast = i, th = 1)
  output_file <- paste0("diffbind_results_contrast_", i, ".csv")
  write.csv(report, file = output_file, row.names = FALSE)
  message("Exported contrast ", i, " to ", output_file)
}

#-------------------------
# 4. Expression Analysis of Target Genes
#-------------------------

# Load required libraries
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(reshape2)

# Read expression data
setwd("C:/Users/35471/Desktop/DIFFBIND新/SYF_F5")
target_genes <- read.table('rna_up_gene.txt', sep = '\t', header = FALSE)
expression_data <- read.table(
  'C:/Users/35471/Desktop/follcile/rna-seq/GC_tpm_no_syf.txt',
  header = TRUE, 
  sep = '\t'
)

# Merge data
merged_data <- merge(
  target_genes, 
  expression_data, 
  by.x = "V1", 
  by.y = "gene_id", 
  all.x = TRUE
)

# Prepare for plotting
plot_data <- merged_data %>% select("V1", "SYF", "F5")
melted_data <- melt(
  plot_data, 
  id.vars = "V1",
  variable.name = "tissue",
  value.name = "tpm"
)

# Define colors
tissue_colors <- c("#F29093", "#9CBEDB")

# Create expression plot
expression_plot <- ggplot(melted_data, aes(x = tissue, y = tpm)) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_violin(aes(fill = tissue), outlier.colour = "white", width = 0.8) +
  geom_boxplot(width = 0.15) +
  scale_fill_manual(values = tissue_colors) +
  labs(title = "Expression of Target Genes", 
       x = "Tissue", 
       y = "TPM") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14)
  ) +
  geom_signif(
    comparisons = list(c("SYF", "F5")),
    map_signif_level = TRUE,
    test = t.test,
    tip_length = 0.05,
    size = 0.8,
    color = "black"
  )

# Save plot
ggsave("Target_Gene_Expression_SYF_vs_F5.png", expression_plot, 
       width = 8, height = 6, dpi = 300)

#-------------------------
# 5. Chromatin State Dynamics Analysis
#-------------------------

# Load required libraries
library(tidyverse)
library(data.table)
library(GenomicRanges)

# Set up directories
setwd("C:/Users/35471/Desktop/follcile/染色质变化/state")
output_dir <- "chromatin_dynamics"
dir.create(output_dir, showWarnings = FALSE)

# Read chromosome sizes
chrom_sizes <- fread("cau_ASM874695v1_112.size", col.names = c("chr", "size"))
genome_size <- sum(chrom_sizes$size)

# Define developmental stages
stages <- c("SWF", "LWF", "SYF", "F5", "F3", "F1", "POF")

# State mapping function
state_mapping <- function(state) {
  case_when(
    state %in% c("F9", "F8", "F13") ~ "TSS",
    state %in% c("F15", "F14", "F5", "F4", "F10") ~ "TX",
    state %in% c("F11", "F12", "F7", "F2") ~ "Enh",
    state %in% c("F3", "F6") ~ "RePr",
    state == "F1" ~ "Qui",
    TRUE ~ "Other"
  )
}

# Compute state proportions
compute_state_proportions <- function(stage) {
  state_bed <- fread(
    paste0(stage, "_15_segments.bed"), 
    col.names = c("chr", "start", "end", "state")
  )
  
  state_bed %>%
    mutate(
      state_merged = state_mapping(state),
      width = end - start
    ) %>%
    group_by(state_merged) %>%
    summarise(total_length = sum(width)) %>%
    mutate(
      proportion = total_length / genome_size,
      stage = stage
    ) %>%
    select(stage, state_merged, proportion)
}

# Compute chromatin state changes
compute_transition_changes <- function(prev_stage, next_stage) {
  prev_bed <- fread(
    paste0(prev_stage, "_15_segments.bed"), 
    col.names = c("chr", "start", "end", "state")
  )
  next_bed <- fread(
    paste0(next_stage, "_15_segments.bed"), 
    col.names = c("chr", "start", "end", "state")
  )
  
  prev_bed <- prev_bed %>% mutate(state_merged = state_mapping(state))
  next_bed <- next_bed %>% mutate(state_merged = state_mapping(state))
  
  prev_gr <- makeGRangesFromDataFrame(prev_bed, keep.extra.columns = TRUE)
  next_gr <- makeGRangesFromDataFrame(next_bed, keep.extra.columns = TRUE)
  
  overlaps <- findOverlaps(prev_gr, next_gr)
  
  changes_df <- data.frame(
    prev_chr = seqnames(prev_gr)[queryHits(overlaps)],
    prev_start = start(prev_gr)[queryHits(overlaps)],
    prev_end = end(prev_gr)[queryHits(overlaps)],
    prev_state = prev_gr$state_merged[queryHits(overlaps)],
    next_state = next_gr$state_merged[subjectHits(overlaps)]
  ) %>%
    mutate(
      overlap_start = pmax(prev_start, start(next_gr)[subjectHits(overlaps)]),
      overlap_end = pmin(prev_end, end(next_gr)[subjectHits(overlaps)]),
      width = overlap_end - overlap_start
    )
  
  changes_df %>%
    filter(prev_state != next_state) %>%
    group_by(prev_state) %>%
    summarise(changed_length = sum(width)) %>%
    mutate(
      transition = paste0(prev_stage, "→", next_stage),
      proportion = changed_length / genome_size
    ) %>%
    select(transition, prev_state, proportion)
}

# Compute all proportions
all_proportions <- map_dfr(stages, compute_state_proportions)

# Compute all transitions
stage_pairs <- list(
  c("SWF", "LWF"), c("LWF", "SYF"), c("SYF", "F5"),
  c("F5", "F3"), c("F3", "F1"), c("F1", "POF")
)
all_changes <- map_dfr(stage_pairs, ~ compute_transition_changes(.[1], .[2]))

# Define state colors
state_colors <- c(
  "TSS" = "#E41A1C", "TX" = "#377EB8", 
  "Enh" = "#4DAF4A", "RePr" = "#984EA3", 
  "Qui" = "#FF7F00", "Other" = "#A9A9A9"
)

# Plot state proportions
prop_plot <- all_proportions %>%
  mutate(state_merged = factor(state_merged, levels = names(state_colors))) %>%
  ggplot(aes(x = stage, y = proportion, fill = state_merged)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = state_colors) +
  labs(title = "Chromatin State Proportions Across Development",
       y = "Genomic Proportion", 
       x = "Developmental Stage",
       fill = "Chromatin State") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "chromatin_state_proportions.pdf"), 
       prop_plot, width = 10, height = 6)

# Plot state changes
change_plot <- all_changes %>%
  filter(prev_state != "Qui") %>%
  mutate(transition = factor(
    transition, 
    levels = sapply(stage_pairs, function(x) paste0(x[1], "→", x[2]))
  )) %>%
  ggplot(aes(x = transition, y = proportion, fill = prev_state)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = state_colors) +
  labs(title = "Chromatin State Changes Between Stages",
       y = "Proportion of Genome with Changed State", 
       x = "Developmental Transition",
       fill = "Previous State") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "chromatin_state_changes.pdf"), 
       change_plot, width = 10, height = 6)