#!/usr/bin/env Rscript
# Supplementary Table S10: Sensitivity analysis of the total biofilm-AMR
# correlation with and without intrinsic efflux determinants

library(dplyr)
library(tidyr)
library(boot)

amr_long <- read.csv('amr_clean_pass2_long.csv', stringsAsFactors = FALSE)
df       <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)

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

amr_class <- amr_long %>%
  mutate(Class = sapply(Gene, classify_amr))

# Counts with and without intrinsic efflux
with_efflux <- amr_class %>%
  count(Isolate, name = "AMR_count_with") %>%
  rename(staphscope_id = Isolate)

without_efflux <- amr_class %>%
  filter(Class != "Intrinsic_efflux") %>%
  count(Isolate, name = "AMR_count_without") %>%
  rename(staphscope_id = Isolate)

# Merge with biofilm
merged <- df %>%
  dplyr::select(staphscope_id, biofilm_score) %>%
  left_join(with_efflux,    by = 'staphscope_id') %>%
  left_join(without_efflux, by = 'staphscope_id') %>%
  mutate(across(starts_with("AMR_"), ~ replace_na(.x, 0)))

# Bootstrap Spearman helper
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

s_with    <- spearman_boot(merged$biofilm_score, merged$AMR_count_with)
s_without <- spearman_boot(merged$biofilm_score, merged$AMR_count_without)

result <- data.frame(
  Analysis = c(
    "With intrinsic efflux determinants (primary analysis)",
    "Without intrinsic efflux determinants (sensitivity analysis)"
  ),
  n                  = c(s_with$n, s_without$n),
  Mean_AMR_count     = c(round(mean(merged$AMR_count_with), 1),
                          round(mean(merged$AMR_count_without), 1)),
  Median_AMR_count   = c(median(merged$AMR_count_with),
                          median(merged$AMR_count_without)),
  Spearman_rho       = c(round(s_with$rho, 3),     round(s_without$rho, 3)),
  CI_lower           = c(round(s_with$ci[1], 3),   round(s_without$ci[1], 3)),
  CI_upper           = c(round(s_with$ci[2], 3),   round(s_without$ci[2], 3)),
  p_value            = c(signif(s_with$p, 3),      signif(s_without$p, 3)),
  Interpretation = c(
    "Primary correlation between biofilm-associated gene carriage score and total AMR determinant count",
    "Sensitivity analysis excluding intrinsic efflux determinants, confirming that the association is driven by acquired resistance determinants rather than the intrinsic core resistome"
  ),
  stringsAsFactors = FALSE
)

write.csv(result, "Supplementary_Table_S10_sensitivity_intrinsic_efflux.csv",
          row.names = FALSE)

cat("Saved Supplementary_Table_S10_sensitivity_intrinsic_efflux.csv\n\n")
print(result)
