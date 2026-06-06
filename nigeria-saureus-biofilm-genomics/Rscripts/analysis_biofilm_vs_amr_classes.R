#!/usr/bin/env Rscript
# Correlation between biofilm score and AMR gene classes

library(dplyr)
library(ggplot2)
library(corrplot)

# ----------------------------------------------------------------------
# 1. Read data
# ----------------------------------------------------------------------
df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)
amr <- read.csv('amr_genes.csv', stringsAsFactors = FALSE)

# ----------------------------------------------------------------------
# 2. Define AMR classes based on gene name patterns
# ----------------------------------------------------------------------
classify_amr <- function(gene) {
  gene_lower <- tolower(gene)
  if (grepl('mec', gene_lower)) return('Beta_lactam')
  if (grepl('bla', gene_lower)) return('Beta_lactam')
  if (grepl('erm', gene_lower)) return('Macrolide')
  if (grepl('msr', gene_lower)) return('Macrolide')
  if (grepl('mph', gene_lower)) return('Macrolide')
  if (grepl('aac|aph|ant', gene_lower)) return('Aminoglycoside')
  if (grepl('str', gene_lower)) return('Aminoglycoside')
  if (grepl('tet', gene_lower)) return('Tetracycline')
  if (grepl('dfr', gene_lower)) return('Trimethoprim')
  if (grepl('cat', gene_lower)) return('Chloramphenicol')
  if (grepl('lnu|lin', gene_lower)) return('Lincosamide')
  if (grepl('qac', gene_lower)) return('Biocide')
  if (grepl('fos', gene_lower)) return('Fosfomycin')
  if (grepl('van', gene_lower)) return('Vancomycin')
  if (grepl('mup', gene_lower)) return('Mupirocin')
  return('Other')
}

# ----------------------------------------------------------------------
# 3. Build per‑isolate AMR class counts
# ----------------------------------------------------------------------
amr_classes <- amr %>%
  mutate(Class = sapply(Gene, classify_amr)) %>%
  group_by(Class) %>%
  do({
    # For each class, extract all sample IDs
    samples <- unlist(strsplit(.$Genomes, ';'))
    data.frame(Sample = samples, stringsAsFactors = FALSE)
  }) %>%
  group_by(Class, Sample) %>%
  summarise(count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = Class, values_from = count, values_fill = 0)

# Merge with biofilm data
df_biofilm <- df[, c('staphscope_id', 'biofilm_score')]
colnames(amr_classes)[2] <- 'staphscope_id'  # adjust column name
df_merged <- merge(df_biofilm, amr_classes, by = 'staphscope_id', all.x = TRUE)
df_merged[is.na(df_merged)] <- 0

# ----------------------------------------------------------------------
# 4. Correlation matrix (Spearman) between biofilm score and AMR classes
# ----------------------------------------------------------------------
# Exclude the first column (staphscope_id)
cor_matrix <- cor(df_merged[, -1], method = "spearman")
# Extract only correlations with biofilm_score
cor_biofilm <- cor_matrix['biofilm_score', , drop = FALSE]
print(cor_biofilm)

# Convert to data frame for saving
cor_df <- data.frame(
  AMR_Class = colnames(cor_biofilm)[-1],  # exclude biofilm_score itself
  Spearman_rho = as.numeric(cor_biofilm[1, -1])
)
cor_df <- cor_df[order(-abs(cor_df$Spearman_rho)), ]
print(cor_df)

write.csv(cor_df, "biofilm_vs_amr_class_correlations.csv", row.names = FALSE)

# ----------------------------------------------------------------------
# 5. Barplot of top correlations
# ----------------------------------------------------------------------
p <- ggplot(cor_df[1:min(10, nrow(cor_df)), ], aes(x = reorder(AMR_Class, Spearman_rho), y = Spearman_rho, fill = Spearman_rho > 0)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(title = "Spearman correlation: biofilm score vs AMR class count",
       x = "AMR class", y = "Spearman's rho") +
  theme_minimal() +
  theme(legend.position = 'none')

ggsave("fig8_biofilm_vs_amr_class_corr.pdf", plot = p, width = 6, height = 5, dpi = 300)
cat("Correlation plot saved as fig8_biofilm_vs_amr_class_corr.pdf\n")
