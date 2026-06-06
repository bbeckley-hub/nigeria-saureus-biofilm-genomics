#!/usr/bin/env Rscript
# Analysis 1: Compare biofilm scores between MRSA and MSSA

library(ggplot2)
library(dplyr)

# ----------------------------------------------------------------------
# 1. Read data
# ----------------------------------------------------------------------
df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)

# Convert MRSA_Status to factor with appropriate levels
df$MRSA_Status <- factor(df$MRSA_Status, levels = c('MSSA', 'MRSA'))

# ----------------------------------------------------------------------
# 2. Summary statistics
# ----------------------------------------------------------------------
summary_stats <- df %>%
  group_by(MRSA_Status) %>%
  summarise(
    n = n(),
    mean_score = mean(biofilm_score, na.rm = TRUE),
    median_score = median(biofilm_score, na.rm = TRUE),
    sd_score = sd(biofilm_score, na.rm = TRUE),
    min_score = min(biofilm_score, na.rm = TRUE),
    max_score = max(biofilm_score, na.rm = TRUE)
  )
print(summary_stats)

# ----------------------------------------------------------------------
# 3. Statistical test (Wilcoxon / Mann-Whitney)
# ----------------------------------------------------------------------
test_result <- wilcox.test(biofilm_score ~ MRSA_Status, data = df)
print(test_result)

# ----------------------------------------------------------------------
# 4. Boxplot with jittered points
# ----------------------------------------------------------------------
p <- ggplot(df, aes(x = MRSA_Status, y = biofilm_score, fill = MRSA_Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("MSSA" = "#2E86AB", "MRSA" = "#D64933")) +
  labs(
    title = "Biofilm genetic score by MRSA status",
    x = NULL,
    y = "Biofilm score (sum of 14 genes)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Add p-value annotation
p_value <- test_result$p.value
p_label <- ifelse(p_value < 0.001, "p < 0.001", paste("p =", round(p_value, 4)))
p <- p + annotate("text", x = 1.5, y = max(df$biofilm_score) + 0.3, label = p_label, size = 4.5)

# Display plot
print(p)

# ----------------------------------------------------------------------
# 5. Save plot as PDF
# ----------------------------------------------------------------------
ggsave("fig1_biofilm_score_MRSA_vs_MSSA.pdf", plot = p, width = 5, height = 6, dpi = 300)
cat("Plot saved as fig1_biofilm_score_MRSA_vs_MSSA.pdf\n")
