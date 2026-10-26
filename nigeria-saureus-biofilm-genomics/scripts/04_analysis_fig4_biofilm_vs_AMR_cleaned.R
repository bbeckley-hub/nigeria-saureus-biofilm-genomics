#!/usr/bin/env Rscript
# Figure 4 + printed correlations: Biofilm-associated gene carriage score
# vs cleaned AMR determinant count (deduplicated, no biofilm/virulence genes)

library(dplyr)
library(tidyr)
library(boot)
library(ppcor)
library(ggplot2)

# ---------- Load cleaned data ----------
amr <- read.csv('amr_clean_pass2_per_isolate.csv', stringsAsFactors = FALSE)
df  <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)

# ---------- Column detection ----------
find_col <- function(df, candidates) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  return(NA)
}
id_col   <- find_col(df, c("staphscope_id","Staphscope_id","isolate","Isolate"))
mrsa_col <- find_col(df, c("MRSA_Status","mrsa_status","MRSA_status"))
mlst_col <- find_col(df, c("MLST","mlst","ST"))

stopifnot(!is.na(id_col), !is.na(mrsa_col), !is.na(mlst_col))
df[[mrsa_col]] <- as.character(df[[mrsa_col]])
df[[mlst_col]] <- as.character(df[[mlst_col]])

# ---------- Merge ----------
merged <- df %>%
  left_join(amr, by = setNames("Isolate", id_col))

cat("Merged:", nrow(merged), "| Missing AMR count:",
    sum(is.na(merged$AMR_gene_count)), "\n")
cat("MRSA statuses:\n"); print(table(merged[[mrsa_col]]))
cat("\n")

# ---------- Bootstrap Spearman helper ----------
spearman_boot <- function(x, y, R = 2000) {
  d <- data.frame(x = x, y = y); d <- d[complete.cases(d), ]
  stat <- function(d, i) cor(d$x[i], d$y[i], method = "spearman")
  set.seed(42)
  b  <- boot(d, stat, R = R)
  ci <- tryCatch(boot.ci(b, type = "bca")$bca[4:5],
                 error = function(e) boot.ci(b, type = "perc")$perc[4:5])
  ct <- cor.test(d$x, d$y, method = "spearman", exact = FALSE)
  list(rho = unname(ct$estimate), ci = ci, p = ct$p.value, n = nrow(d))
}

# ---------- Total ----------
sb <- spearman_boot(merged$biofilm_score, merged$AMR_gene_count)
cat("=== Total AMR correlation ===\n")
cat(sprintf("rho = %.3f | 95%% CI: %.3f to %.3f | p = %.3g | n = %d\n\n",
            sb$rho, sb$ci[1], sb$ci[2], sb$p, sb$n))

# ---------- Stratified ----------
cat("=== Stratified by MRSA status ===\n")
for (status in sort(unique(merged[[mrsa_col]]))) {
  sub <- merged[merged[[mrsa_col]] == status, ]
  if (nrow(sub) < 3) next
  s <- spearman_boot(sub$biofilm_score, sub$AMR_gene_count)
  cat(sprintf("%-6s: rho = %.3f | 95%% CI: %.3f to %.3f | p = %.3g | n = %d\n",
              status, s$rho, s$ci[1], s$ci[2], s$p, s$n))
}

# ---------- Partial (MLST-adjusted) ----------
cat("\n=== Partial correlation controlling for MLST ===\n")
merged$mlst_numeric <- as.numeric(factor(merged[[mlst_col]]))
pc <- pcor.test(merged$biofilm_score, merged$AMR_gene_count,
                merged$mlst_numeric, method = "spearman")
cat(sprintf("partial rho = %.3f | p = %.3g\n\n", pc$estimate, pc$p.value))

# ---------- Save merged dataset ----------
write.csv(merged, "biofilm_analysis_dataset_cleaned_AMR.csv", row.names = FALSE)
cat("Saved: biofilm_analysis_dataset_cleaned_AMR.csv\n")

# ---------- Figure 4 ----------
p_caption <- sprintf(
  "Spearman rho = %.2f (95%% CI %.2f\u2013%.2f), p = %.1e; n = %d",
  sb$rho, sb$ci[1], sb$ci[2], sb$p, sb$n
)

p_fig4 <- ggplot(merged,
                 aes(x = AMR_gene_count, y = biofilm_score,
                     colour = .data[[mrsa_col]])) +
  geom_jitter(width = 0.4, height = 0.12, size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, colour = "black",
              linetype = "dashed", linewidth = 0.7,
              aes(group = 1)) +
  scale_colour_manual(values = c("MRSA" = "#D64933", "MSSA" = "#2E86AB")) +
  labs(x = "AMR determinant count (unique, deduplicated)",
       y = "Biofilm-associated gene carriage score (0\u201324)",
       colour = "MRSA status",
       caption = p_caption) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("fig4_biofilm_vs_AMR_count.pdf", p_fig4, width = 7, height = 5, dpi = 300)
ggsave("fig4_biofilm_vs_AMR_count.tif", p_fig4, width = 7, height = 5,
       dpi = 300, compression = "lzw")
ggsave("fig4_biofilm_vs_AMR_count.png", p_fig4, width = 7, height = 5, dpi = 300)
cat("Saved: fig4_biofilm_vs_AMR_count.{pdf,tif,png}\n")
