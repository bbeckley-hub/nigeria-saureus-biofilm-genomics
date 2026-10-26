#!/usr/bin/env Rscript
# Post-canonicalisation filter for housekeeping genes that survived Layer 1
# (e.g., dha1_1 → dha1 because the Layer 1 anchor ^dha1$ didn't match)

library(dplyr)

amr <- read.csv('amr_clean_pass2_long.csv', stringsAsFactors = FALSE)

housekeeping_canonical <- c("dha1", "dhap", "kdpd", "rlmh", "sepa")

before <- nrow(amr)
amr_filtered <- amr %>% filter(!(Gene %in% housekeeping_canonical))
removed <- before - nrow(amr_filtered)

cat("Rows before:", before, "| after:", nrow(amr_filtered),
    "| removed:", removed, "\n")
cat("Housekeeping genes still present after pass2:\n")
print(intersect(unique(amr$Gene), housekeeping_canonical))

if (removed == 0) {
  cat("\nNo surviving housekeeping genes. Nothing to change.\n")
} else {
  per_isolate <- amr_filtered %>%
    count(Isolate, name = "AMR_gene_count") %>%
    arrange(desc(AMR_gene_count))
  cat("\nNew per-isolate count:\n")
  cat("  Range:", range(per_isolate$AMR_gene_count), "\n")
  cat("  Median:", median(per_isolate$AMR_gene_count), "\n")
  cat("  Mean:", round(mean(per_isolate$AMR_gene_count), 1), "\n")
  write.csv(amr_filtered, "amr_clean_pass2_long.csv", row.names = FALSE)
  write.csv(per_isolate, "amr_clean_pass2_per_isolate.csv", row.names = FALSE)
  cat("\nOverwrote: amr_clean_pass2_long.csv, amr_clean_pass2_per_isolate.csv\n")
}
