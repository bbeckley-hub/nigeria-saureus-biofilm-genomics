#!/usr/bin/env Rscript
# Figure 6: Biofilm-associated gene carriage score vs AMR functional classes
# Uses cleaned AMR data from the three-step pipeline
# Adds Intrinsic_efflux class; filters figure to classes with >=10 positive isolates

library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)

select <- dplyr::select

amr_clean <- read.csv('amr_clean_pass2_long.csv', stringsAsFactors = FALSE)
df        <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)

classify_amr <- function(gene) {
  g <- tolower(trimws(gene))
  if (grepl('^(mepa|mepb|lmrs|nora|norb|norc)$', g)) return('Intrinsic_efflux')
  if (grepl('^mec', g) || grepl('^bla', g) || grepl('^pc1', g)) return('Beta_lactam')
  if (grepl('^erm', g) || grepl('^msr', g) || grepl('^mph', g)) return('Macrolide')
  if (grepl('^aac', g) || grepl('^aph', g) || grepl('^ant', g) ||
      grepl('^aad', g) || grepl('^str$', g) || grepl('^sat', g)) return('Aminoglycoside')
  if (grepl('^tet', g)) return('Tetracycline')
  if (grepl('^dfr', g)) return('Trimethoprim')
  if (grepl('^cat', g)) return('Chloramphenicol')
  if (grepl('^lnu', g) || grepl('^lin', g)) return('Lincosamide')
  if (grepl('^qac', g)) return('Biocide')
  if (grepl('^fos', g)) return('Fosfomycin')
  if (grepl('^mup', g) || grepl('^iles', g)) return('Mupirocin')
  return('Other')
}

amr_class_counts <- amr_clean %>%
  mutate(Class = sapply(Gene, classify_amr)) %>%
  count(Isolate, Class) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  rename(staphscope_id = Isolate)

merged <- df %>%
  dplyr::select(staphscope_id, biofilm_score) %>%
  left_join(amr_class_counts, by = 'staphscope_id') %>%
  mutate(across(-c(staphscope_id, biofilm_score), ~ replace_na(.x, 0)))

class_cols <- setdiff(names(amr_class_counts), 'staphscope_id')
cat("Classes detected:", length(class_cols), "\n")
cat(paste(" -", class_cols, collapse = "\n"), "\n\n")

spearman_boot <- function(x, y, R = 2000) {
  d <- data.frame(x = x, y = y); d <- d[complete.cases(d), ]
  if (nrow(d) < 5 || sd(d$x) == 0 || sd(d$y) == 0) {
    return(list(rho = NA_real_, ci = c(NA, NA), p = NA_real_, n = nrow(d)))
  }
  stat <- function(d, i) cor(d$x[i], d$y[i], method = "spearman")
  set.seed(42)
  b  <- boot(d, stat, R = R)
  ci <- tryCatch(boot.ci(b, type = "bca")$bca[4:5],
                 error = function(e) boot.ci(b, type = "perc")$perc[4:5])
  ct <- cor.test(d$x, d$y, method = "spearman", exact = FALSE)
  list(rho = unname(ct$estimate), ci = ci, p = ct$p.value, n = nrow(d))
}

results <- lapply(class_cols, function(cls) {
  s <- spearman_boot(merged$biofilm_score, merged[[cls]])
  data.frame(
    AMR_Class    = cls,
    Spearman_rho = round(s$rho, 3),
    CI_low       = round(s$ci[1], 3),
    CI_high      = round(s$ci[2], 3),
    p_value      = signif(s$p, 3),
    n            = s$n,
    n_positive   = sum(merged[[cls]] > 0),
    mean_count   = round(mean(merged[[cls]]), 2),
    stringsAsFactors = FALSE
  )
})
cor_df <- do.call(rbind, results) %>% arrange(desc(abs(Spearman_rho)))

cat("=== Per-class correlations (ALL classes) ===\n")
print(cor_df)
cat("\n")

write.csv(cor_df, "Supplementary_Table_S3_correlations.csv", row.names = FALSE)
cat("Saved: Supplementary_Table_S3_correlations.csv\n\n")

# Filter to n_positive >= 10 for the manuscript figure
MIN_POSITIVE <- 10
plot_df <- cor_df %>%
  filter(!is.na(Spearman_rho), n_positive >= MIN_POSITIVE)
dropped <- setdiff(cor_df$AMR_Class, plot_df$AMR_Class)

cat("=== Figure filter (n_positive >= ", MIN_POSITIVE, ") ===\n", sep = "")
cat("Classes retained for figure:", nrow(plot_df), "of", nrow(cor_df), "\n")
if (length(dropped) > 0) {
  cat("Classes dropped:\n")
  for (d in dropped) {
    cat(sprintf("  - %s (n_positive = %d)\n", d,
                cor_df$n_positive[cor_df$AMR_Class == d]))
  }
}
cat("\n")

if (nrow(plot_df) == 0) stop("No reliable classes to plot.")

# ============ FIGURE 6 (manuscript; n_positive >= 10) ============
p_fig6 <- ggplot(plot_df, aes(x = reorder(AMR_Class, Spearman_rho),
                                y = Spearman_rho,
                                fill = Spearman_rho > 0)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.25) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#D64933"),
                    guide = "none") +
  labs(x = "AMR functional class",
       y = "Spearman's rho (95% bootstrap CI)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("fig6_biofilm_vs_amr_class_corr.pdf", p_fig6, width = 7, height = 5, dpi = 300)
ggsave("fig6_biofilm_vs_amr_class_corr.tif", p_fig6, width = 7, height = 5,
       dpi = 300, compression = "lzw")
ggsave("fig6_biofilm_vs_amr_class_corr.png", p_fig6, width = 7, height = 5, dpi = 300)
cat("Saved: fig6_biofilm_vs_amr_class_corr.{pdf,tif,png}\n\n")

# ============ SUPPLEMENTARY FIGURE S4 (all 12 classes) ============
p_figS4 <- ggplot(cor_df, aes(x = reorder(AMR_Class, Spearman_rho),
                                y = Spearman_rho,
                                fill = Spearman_rho > 0)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.25) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#D64933"),
                    guide = "none") +
  labs(x = "AMR functional class",
       y = "Spearman's rho (95% bootstrap CI)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("figS4_biofilm_vs_amr_class_all.pdf", p_figS4, width = 7, height = 5.5, dpi = 300)
ggsave("figS4_biofilm_vs_amr_class_all.tif", p_figS4, width = 7, height = 5.5,
       dpi = 300, compression = "lzw")
cat("Saved: figS4_biofilm_vs_amr_class_all.{pdf,tif}\n\n")

# ============ Sensitivity analysis ============
merged_noeff <- merged %>%
  dplyr::select(-any_of("Intrinsic_efflux")) %>%
  mutate(AMR_count_noeff = rowSums(across(-c(staphscope_id, biofilm_score))))

s_total <- spearman_boot(merged$biofilm_score,
                         rowSums(merged[, class_cols, drop = FALSE]))
s_noeff <- spearman_boot(merged_noeff$biofilm_score,
                         merged_noeff$AMR_count_noeff)

cat("=== Sensitivity: with vs without intrinsic efflux ===\n")
cat(sprintf("With intrinsic efflux:    rho = %.3f, p = %.3g, n = %d\n",
            s_total$rho, s_total$p, s_total$n))
cat(sprintf("Without intrinsic efflux: rho = %.3f, p = %.3g, n = %d\n",
            s_noeff$rho, s_noeff$p, s_noeff$n))
