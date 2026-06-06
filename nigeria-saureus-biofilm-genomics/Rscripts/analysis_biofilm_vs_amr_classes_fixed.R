#!/usr/bin/env Rscript
# Correlation between biofilm score and AMR gene classes (fixed)

library(dplyr)
library(ggplot2)

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
    samples <- unlist(strsplit(.$Genomes, ';'))
    data.frame(Sample = samples, stringsAsFactors = FALSE)
  }) %>%
  group_by(Class, Sample) %>%
  summarise(count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = Class, values_from = count, values_fill = 0)

# Ensure all columns except Sample are numeric
amr_classes <- as.data.frame(amr_classes)
for (col in setdiff(names(amr_classes), 'Sample')) {
  amr_classes[[col]] <- as.numeric(amr_classes[[col]])
}

# Merge with biofilm data
df_biofilm <- df[, c('staphscope_id', 'biofilm_score')]
colnames(amr_classes)[colnames(amr_classes) == 'Sample'] <- 'staphscope_id'
df_merged <- merge(df_biofilm, amr_classes, by = 'staphscope_id', all.x = TRUE)
df_merged[is.na(df_merged)] <- 0

# ----------------------------------------------------------------------
# 4. Correlation matrix (Spearman) between biofilm score and AMR classes
# ----------------------------------------------------------------------
# Identify numeric columns (exclude staphscope_id and biofilm_score for the matrix)
numeric_cols <- sapply(df_merged, is.numeric)
# Keep only AMR class columns (exclude biofilm_score itself)
amr_class_cols <- setdiff(names(df_merged)[numeric_cols], 'biofilm_score')
if (length(amr_class_cols) == 0) {
  cat("No AMR class columns found. Check classification.\n")
  quit()
}
cor_matrix <- cor(df_merged[, amr_class_cols], df_merged$biofilm_score, method = "spearman")
cor_df <- data.frame(
  AMR_Class = rownames(cor_matrix),
  Spearman_rho = as.numeric(cor_matrix[,1])
)
cor_df <- cor_df[order(-abs(cor_df$Spearman_rho)), ]
print(cor_df)

write.csv(cor_df, "biofilm_vs_amr_class_correlations.csv", row.names = FALSE)

# ----------------------------------------------------------------------
# 5. Barplot of top correlations (if any)
# ----------------------------------------------------------------------
if (nrow(cor_df) > 0) {
  p <- ggplot(cor_df[1:min(10, nrow(cor_df)), ], aes(x = reorder(AMR_Class, Spearman_rho), y = Spearman_rho, fill = Spearman_rho > 0)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    labs(title = "Spearman correlation: biofilm score vs AMR class count",
         x = "AMR class", y = "Spearman's rho") +
    theme_minimal() +
    theme(legend.position = 'none')
  
  ggsave("fig8_biofilm_vs_amr_class_corr.pdf", plot = p, width = 6, height = 5, dpi = 300)
  cat("Correlation plot saved as fig8_biofilm_vs_amr_class_corr.pdf\n")
} else {
  cat("No AMR classes to plot.\n")
}
