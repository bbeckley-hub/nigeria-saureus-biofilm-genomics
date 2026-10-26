#!/usr/bin/env Rscript
# Compare biofilm scores across sample sources

library(ggplot2)
library(dplyr)
library(rstatix)

df <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)

# Show original categories
print("Original SAMPLE categories:")
print(unique(df$SAMPLE))

# Create simplified source categories
df <- df %>%
  mutate(
    Source = case_when(
      grepl('wound', SAMPLE, ignore.case = TRUE) ~ 'Wound',
      grepl('blood', SAMPLE, ignore.case = TRUE) ~ 'Blood',
      grepl('urine|catheter', SAMPLE, ignore.case = TRUE) ~ 'Urine/Catheter',
      grepl('aspirate|ulcer|abscess', SAMPLE, ignore.case = TRUE) ~ 'Aspirate/Abscess',
      grepl('ear|eye|trachea|sputum|HVS|ECS|genitals|discharge', SAMPLE, ignore.case = TRUE) ~ 'Other',
      TRUE ~ 'Other'
    )
  )
df$Source <- factor(df$Source, levels = c('Wound', 'Blood', 'Urine/Catheter', 'Aspirate/Abscess', 'Other'))

print("Source counts:")
print(table(df$Source))

# Summary per source
summary_source <- df %>%
  group_by(Source) %>%
  summarise(
    n = n(),
    mean_score = mean(biofilm_score, na.rm = TRUE),
    median_score = median(biofilm_score, na.rm = TRUE),
    sd_score = sd(biofilm_score, na.rm = TRUE)
  )
print(summary_source)

# Kruskal-Wallis test
kw <- kruskal.test(biofilm_score ~ Source, data = df)
print(kw)

# Post-hoc if significant
if(kw$p.value < 0.05) {
  dunn <- df %>% dunn_test(biofilm_score ~ Source, p.adjust.method = "bonferroni")
  print(dunn)
}

# Boxplot
p <- ggplot(df, aes(x = Source, y = biofilm_score, fill = Source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Source", y = "Biofilm score (0-24)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_label <- paste("Kruskal-Wallis p =", format(kw$p.value, scientific = TRUE, digits = 3))
p <- p + annotate("text", x = 2, y = max(df$biofilm_score) + 0.3, label = p_label, size = 4)

print(p)
ggsave("fig3_biofilm_score_by_source.pdf", plot = p, width = 7, height = 6, dpi = 300)
cat("Plot saved as fig3_biofilm_score_by_source.pdf\n")
