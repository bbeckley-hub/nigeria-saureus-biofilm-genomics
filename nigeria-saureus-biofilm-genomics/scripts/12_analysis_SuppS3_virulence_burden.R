#!/usr/bin/env Rscript
# Supplementary Figure S3: Virulence burden by biofilm-associated gene carriage score group
# Computes total VFDB virulence gene count per isolate from virulence_genes.csv

library(dplyr)
library(tidyr)
library(ggplot2)

df          <- read.csv("biofilm_analysis_dataset_expanded.csv", stringsAsFactors = FALSE)
vir_genes   <- read.csv("virulence_genes.csv", stringsAsFactors = FALSE)

# --- Compute per-isolate virulence gene count from VFDB table ---
if (!"Genomes" %in% names(vir_genes)) {
  stop("virulence_genes.csv must contain a 'Genomes' column with semicolon-separated isolate IDs")
}

vir_per_isolate <- vir_genes %>%
  separate_rows(Genomes, sep = ";") %>%
  rename(staphscope_id = Genomes) %>%
  distinct(staphscope_id, Gene) %>%
  count(staphscope_id, name = "virulence_count")

cat("Isolates in virulence_genes.csv:", nrow(vir_per_isolate), "\n")

# --- Merge with biofilm dataset ---
merged <- df %>%
  dplyr::select(staphscope_id, biofilm_score, MRSA_Status) %>%
  left_join(vir_per_isolate, by = "staphscope_id") %>%
  mutate(virulence_count = replace_na(virulence_count, 0))

cat("Isolates merged:", nrow(merged), "\n")
cat("Biofilm score range:", range(merged$biofilm_score), "\n")
cat("Virulence count range:", range(merged$virulence_count), "\n\n")

# --- Spearman correlation ---
ct <- cor.test(merged$biofilm_score, merged$virulence_count,
               method = "spearman", exact = FALSE)
cat("Spearman (biofilm score vs virulence count):\n")
cat(sprintf("  rho = %.3f, p = %.3f, n = %d\n\n",
            ct$estimate, ct$p.value, nrow(merged)))

# --- High vs low groups ---
merged$group <- factor(
  ifelse(merged$biofilm_score >= 23, "High (\u226523)", "Low (\u226422)"),
  levels = c("Low (\u226422)", "High (\u226523)")
)

wt <- wilcox.test(virulence_count ~ group, data = merged)
cat("Wilcoxon rank-sum test:\n")
cat(sprintf("  W = %.1f, p = %.3f\n\n", wt$statistic, wt$p.value))

cat("Summary by group:\n")
print(merged %>% group_by(group) %>%
  summarise(n = n(), median = median(virulence_count),
            mean = round(mean(virulence_count), 1), .groups = "drop"))

# --- Figure ---
p <- ggplot(merged, aes(x = group, y = virulence_count, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.55) +
  scale_fill_manual(values = c("Low (\u226422)" = "#2E86AB",
                                "High (\u226523)" = "#D64933"),
                    guide = "none") +
  labs(x = "Biofilm-associated gene carriage score group",
       y = "Total VFDB virulence gene count per isolate",
       subtitle = sprintf("Wilcoxon p = %.2f; Spearman rho = %.2f",
                          wt$p.value, ct$estimate)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("figS3_virulence_burden.pdf", p, width = 6, height = 5, dpi = 300)
ggsave("figS3_virulence_burden.tif", p, width = 6, height = 5,
       dpi = 300, compression = "lzw")
cat("\nSaved: figS3_virulence_burden.{pdf,tif}\n")

# Bootstrap CI for Spearman rho
library(boot)
set.seed(42)
b <- boot(merged, function(d, i) cor(d$biofilm_score[i], d$virulence_count[i],
                                     method = "spearman"), R = 2000)
ci <- boot.ci(b, type = "bca")$bca[4:5]
cat(sprintf("Bootstrap 95%% CI: %.3f to %.3f\n", ci[1], ci[2]))
