#!/usr/bin/env Rscript
# Correlation between biofilm score and total AMR gene count

library(ggplot2)
library(dplyr)

df <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)

# Overall Spearman correlation
cor_test <- cor.test(df$biofilm_score, df$amr_gene_count, method = "spearman")
print(cor_test)

# Scatter plot
p_all <- ggplot(df, aes(x = amr_gene_count, y = biofilm_score)) +
  geom_point(aes(color = MRSA_Status), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", fill = "gray80") +
  scale_color_manual(values = c("MSSA" = "#2E86AB", "MRSA" = "#D64933")) +
  labs(x = "Total AMR gene count", y = "Biofilm score (0-24)", color = "MRSA status") +
  theme_minimal(base_size = 14)

# Add correlation annotation
rho <- round(cor_test$estimate, 3)
p_val <- cor_test$p.value
p_label <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 4)))
anno_text <- paste0("Spearman rho = ", rho, "\n", p_label)
p_all <- p_all + annotate("text", x = max(df$amr_gene_count) * 0.8, y = min(df$biofilm_score) + 0.5, 
                          label = anno_text, size = 4, hjust = 0)

print(p_all)

# Stratified correlations
cor_mssa <- cor.test(df$biofilm_score[df$MRSA_Status == "MSSA"], 
                     df$amr_gene_count[df$MRSA_Status == "MSSA"], 
                     method = "spearman")
cor_mrsa <- cor.test(df$biofilm_score[df$MRSA_Status == "MRSA"], 
                     df$amr_gene_count[df$MRSA_Status == "MRSA"], 
                     method = "spearman")
cat("\nStratified correlations:\n")
cat("MSSA (n = 49): rho =", round(cor_mssa$estimate, 3), ", p =", format(cor_mssa$p.value, scientific = TRUE), "\n")
cat("MRSA (n = 58): rho =", round(cor_mrsa$estimate, 3), ", p =", format(cor_mrsa$p.value, scientific = TRUE), "\n")

ggsave("fig4_biofilm_vs_AMR_count.pdf", plot = p_all, width = 7, height = 6, dpi = 300)
cat("Plot saved as fig4_biofilm_vs_AMR_count.pdf\n")
