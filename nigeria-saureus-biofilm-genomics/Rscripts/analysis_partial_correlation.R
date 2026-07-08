df <- read.csv("biofilm_analysis_dataset_expanded.csv")
lm_fit <- lm(biofilm_score ~ as.factor(MLST), data = df)
residuals_biofilm <- residuals(lm_fit)
cor_partial <- cor.test(residuals_biofilm, df$amr_gene_count, method = "spearman")
cor_raw <- cor.test(df$biofilm_score, df$amr_gene_count, method = "spearman")
results <- data.frame(
  method = c("Raw Spearman", "Partial (controlling for MLST)"),
  rho = c(cor_raw$estimate, cor_partial$estimate),
  p = c(cor_raw$p.value, cor_partial$p.value)
)
write.csv(results, "partial_correlation_results.csv", row.names = FALSE)
print(results)
