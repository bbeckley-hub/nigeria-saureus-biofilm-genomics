#!/usr/bin/env Rscript
# Chi-square tests for each biofilm gene vs MRSA status

library(dplyr)
library(broom)

df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)

# List of biofilm gene columns
gene_cols <- c('icaA', 'icaB', 'icaC', 'icaD', 'icaR',
               'fnbA', 'fnbB', 'clfA', 'clfB', 'sdrC', 'sdrD', 'sdrE',
               'ebp', 'cna')

# Function to run chi-square test for each gene
results <- list()
for (gene in gene_cols) {
  tab <- table(df[[gene]], df$mrsa_binary)
  test <- chisq.test(tab, simulate.p.value = TRUE)  # simulate for small counts
  results[[gene]] <- data.frame(
    Gene = gene,
    p_value = test$p.value,
    MRSA_pos = sum(df[[gene]] == 1 & df$mrsa_binary == 1),
    MRSA_neg = sum(df[[gene]] == 1 & df$mrsa_binary == 0),
    MSSA_pos = sum(df[[gene]] == 0 & df$mrsa_binary == 1),
    MSSA_neg = sum(df[[gene]] == 0 & df$mrsa_binary == 0)
  )
}

df_genes <- bind_rows(results)
df_genes <- df_genes[order(df_genes$p_value), ]
print(df_genes)

# Save to CSV
write.csv(df_genes, "discriminating_biofilm_genes.csv", row.names = FALSE)
cat("Results saved to discriminating_biofilm_genes.csv\n")
