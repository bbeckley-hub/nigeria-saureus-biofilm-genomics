#!/usr/bin/env Rscript
# Compare biofilm scores across MLST sequence types

library(ggplot2)
library(dplyr)
library(rstatix)

df <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)

# Keep only STs with at least 5 isolates
st_counts <- df %>%
  filter(MLST != 'Not typed') %>%
  group_by(MLST) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

print("ST counts (>=5 isolates):")
print(st_counts %>% filter(n >= 5))

top_sts <- st_counts %>% filter(n >= 5) %>% pull(MLST)
df_st <- df %>% filter(MLST %in% top_sts)
df_st$MLST <- factor(df_st$MLST, levels = top_sts)

# Summary statistics by ST
summary_st <- df_st %>%
  group_by(MLST) %>%
  summarise(
    n = n(),
    mean_score = mean(biofilm_score, na.rm = TRUE),
    median_score = median(biofilm_score, na.rm = TRUE),
    sd_score = sd(biofilm_score, na.rm = TRUE)
  )
print(summary_st)

# Kruskal-Wallis test
kw <- kruskal.test(biofilm_score ~ MLST, data = df_st)
print(kw)

# Post-hoc Dunn test with Bonferroni correction
if(kw$p.value < 0.05) {
  dunn <- df_st %>%
    dunn_test(biofilm_score ~ MLST, p.adjust.method = "bonferroni")
  print(dunn)
  # Compact letter display
  library(multcompView)
  pairs <- paste(dunn$group1, dunn$group2, sep = "-")
  p_vals <- dunn$p.adj
  names(p_vals) <- pairs
  letters <- multcompLetters(p_vals)$Letters
  print("Compact letter display:")
  print(letters)
}

# Boxplot
p <- ggplot(df_st, aes(x = MLST, y = biofilm_score, fill = MLST)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "MLST", y = "Biofilm score (0-24)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_label <- paste("Kruskal-Wallis p =", format(kw$p.value, scientific = TRUE, digits = 3))
p <- p + annotate("text", x = 1.5, y = max(df_st$biofilm_score) + 0.3, label = p_label, size = 4)

print(p)
ggsave("fig2_biofilm_score_by_MLST.pdf", plot = p, width = 8, height = 6, dpi = 300)
cat("Plot saved as fig2_biofilm_score_by_MLST.pdf\n")
