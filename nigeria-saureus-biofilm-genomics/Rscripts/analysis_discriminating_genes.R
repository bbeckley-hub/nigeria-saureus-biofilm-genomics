#!/usr/bin/env Rscript
# Fisher's exact test for each biofilm gene vs MRSA status
# Handles genes with no variation (present in all isolates)

library(dplyr)

df <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)

# All biofilm gene columns (2:25)
gene_cols <- names(df)[2:25]

results <- list()
for (gene in gene_cols) {
  tab <- table(df[[gene]], df$mrsa_binary)
  # Check if table has at least 2 rows and 2 columns
  if (nrow(tab) >= 2 && ncol(tab) >= 2) {
    test <- fisher.test(tab, simulate.p.value = TRUE)
    p_val <- test$p.value
  } else {
    # No variation: p-value = 1
    p_val <- 1
  }
  results[[gene]] <- data.frame(
    Gene = gene,
    p_value = p_val,
    MRSA_pos = sum(df[[gene]] == 1 & df$mrsa_binary == 1),
    MRSA_neg = sum(df[[gene]] == 1 & df$mrsa_binary == 0),
    MSSA_pos = sum(df[[gene]] == 0 & df$mrsa_binary == 1),
    MSSA_neg = sum(df[[gene]] == 0 & df$mrsa_binary == 0)
  )
}
df_genes <- bind_rows(results)
df_genes <- df_genes[order(df_genes$p_value), ]
print(df_genes)

write.csv(df_genes, "discriminating_biofilm_genes_expanded.csv", row.names = FALSE)
cat("Results saved to discriminating_biofilm_genes_expanded.csv\n")
