#!/usr/bin/env Rscript
# Supplementary Table S6: Partial correlation between biofilm-associated gene
# carriage score and AMR determinant count, controlling for MLST
# Also computes partial correlations for each AMR functional class

library(dplyr)
library(tidyr)
library(boot)
library(ppcor)

# Prevent MASS (loaded by ppcor) from masking dplyr::select
select <- dplyr::select

# ---------- Load data ----------
amr_long <- read.csv('amr_clean_pass2_long.csv', stringsAsFactors = FALSE)
df       <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)

# ---------- Classifier (matches script 07) ----------
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

# ---------- Per-isolate counts ----------
amr_total <- amr_long %>%
  count(Isolate, name = "AMR_gene_count") %>%
  rename(staphscope_id = Isolate)

amr_classes <- amr_long %>%
  mutate(Class = sapply(Gene, classify_amr)) %>%
  count(Isolate, Class) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  rename(staphscope_id = Isolate)

# ---------- Merge with biofilm data ----------
merged <- df %>%
  dplyr::select(staphscope_id, biofilm_score, MLST) %>%
  left_join(amr_total,   by = 'staphscope_id') %>%
  left_join(amr_classes, by = 'staphscope_id') %>%
  mutate(across(where(is.numeric) & -c(biofilm_score),
                ~ replace_na(.x, 0))) %>%
  mutate(MLST = as.character(MLST))

cat("Merged:", nrow(merged), "isolates\n")
cat("Missing AMR count:", sum(is.na(merged$AMR_gene_count)), "\n\n")

# ---------- Bootstrap partial Spearman helper ----------
partial_spearman <- function(x, y, z, R = 2000) {
  d <- data.frame(x = x, y = y, z = z)
  d <- d[complete.cases(d), ]
  if (nrow(d) < 10) return(list(rho = NA, ci = c(NA, NA), p = NA, n = nrow(d)))
  # Skip if x or y is constant
  if (sd(d$x, na.rm = TRUE) == 0 || sd(d$y, na.rm = TRUE) == 0) {
    return(list(rho = NA, ci = c(NA, NA), p = NA, n = nrow(d)))
  }
  pc <- pcor.test(d$x, d$y, d$z, method = "spearman")
  stat <- function(d, i) {
    di <- d[i, ]
    if (sd(di$x) == 0 || sd(di$y) == 0) return(NA_real_)
    pcor.test(di$x, di$y, di$z, method = "spearman")$estimate
  }
  set.seed(42)
  b  <- boot(d, stat, R = R)
  ci <- tryCatch(boot.ci(b, type = "bca")$bca[4:5],
                 error = function(e) boot.ci(b, type = "perc")$perc[4:5])
  list(rho = pc$estimate, ci = ci, p = pc$p.value, n = nrow(d))
}

# ---------- Encode MLST as numeric ----------
merged$mlst_numeric <- as.numeric(factor(merged$MLST))

# ---------- Raw (unadjusted) correlation for comparison ----------
raw_ct <- cor.test(merged$biofilm_score, merged$AMR_gene_count,
                   method = "spearman", exact = FALSE)

# ---------- Total AMR count (partial) ----------
res_total <- partial_spearman(merged$biofilm_score,
                              merged$AMR_gene_count,
                              merged$mlst_numeric)

# ---------- Per-class partial correlations ----------
class_cols <- intersect(setdiff(names(amr_classes), 'staphscope_id'),
                        names(merged))

cat("Testing partial correlation for", length(class_cols), "classes\n\n")

class_results <- lapply(class_cols, function(cls) {
  raw_rho <- tryCatch(
    suppressWarnings(cor(merged$biofilm_score, merged[[cls]],
                         method = "spearman")),
    error = function(e) NA_real_
  )
  r <- partial_spearman(merged$biofilm_score, merged[[cls]],
                        merged$mlst_numeric)
  data.frame(
    Comparison  = paste0("Biofilm score vs ", cls),
    Raw_rho     = round(raw_rho, 3),
    Partial_rho = round(r$rho, 3),
    CI_low      = round(r$ci[1], 3),
    CI_high     = round(r$ci[2], 3),
    p_value     = signif(r$p, 3),
    n           = r$n,
    Control     = "MLST (as factor)",
    stringsAsFactors = FALSE
  )
})

result_df <- rbind(
  data.frame(
    Comparison  = "Biofilm score vs total AMR determinant count",
    Raw_rho     = round(unname(raw_ct$estimate), 3),
    Partial_rho = round(res_total$rho, 3),
    CI_low      = round(res_total$ci[1], 3),
    CI_high     = round(res_total$ci[2], 3),
    p_value     = signif(res_total$p, 3),
    n           = res_total$n,
    Control     = "MLST (as factor)",
    stringsAsFactors = FALSE
  ),
  do.call(rbind, class_results)
)

cat("=== Partial correlations controlling for MLST ===\n")
print(result_df)
cat("\n")

# ---------- Save ----------
write.csv(result_df, "Supplementary_Table_S6_partial_correlation.csv",
          row.names = FALSE)
cat("Saved: Supplementary_Table_S6_partial_correlation.csv\n")

# ---------- Interpretation summary ----------
cat("\n=== Interpretation summary ===\n")
cat(sprintf("Total AMR count:\n"))
cat(sprintf("  Raw Spearman rho       = %.3f\n", unname(raw_ct$estimate)))
cat(sprintf("  Partial rho (MLST-adj) = %.3f (95%% CI %.3f-%.3f, p = %.3g)\n",
            res_total$rho, res_total$ci[1], res_total$ci[2], res_total$p))
cat(sprintf("  Attenuation from raw to partial: %.3f\n",
            unname(raw_ct$estimate) - res_total$rho))
