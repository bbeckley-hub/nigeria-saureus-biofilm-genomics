#!/usr/bin/env Rscript
# Figure 7: ROC curves for MRSA status
# Model 1: biofilm-associated gene carriage score only
# Model 2: biofilm score + age + sex + sample source
# Both models evaluated on the same 99 complete cases for fair comparison

library(pROC)
library(dplyr)

df <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)

# Prepare covariates
df$age_years <- df$AGE.in.years
df$age_years[df$age_years < 1] <- NA
df$sex <- factor(df$SEX, levels = c('Male', 'Female'))
df$source_group <- case_when(
  grepl('wound', df$SAMPLE, ignore.case = TRUE) ~ 'Wound',
  grepl('Blood', df$SAMPLE, ignore.case = TRUE) ~ 'Blood',
  grepl('urine|catheter', df$SAMPLE, ignore.case = TRUE) ~ 'Urine/Catheter',
  grepl('aspirate|abscess', df$SAMPLE, ignore.case = TRUE) ~ 'Aspirate/Abscess',
  TRUE ~ 'Other'
)
df$source_group <- factor(df$source_group,
                          levels = c('Wound','Blood','Urine/Catheter',
                                     'Aspirate/Abscess','Other'))

# Complete cases
df_complete <- df %>% filter(!is.na(age_years), !is.na(sex),
                              !is.na(source_group), !is.na(biofilm_score))
cat("Complete cases:", nrow(df_complete), "of", nrow(df), "\n\n")

# ================= FAIR COMPARISON (both models on df_complete) =================
m1_fair <- glm(mrsa_binary ~ biofilm_score, data = df_complete, family = binomial)
m2      <- glm(mrsa_binary ~ biofilm_score + age_years + sex + source_group,
               data = df_complete, family = binomial)

r1_fair <- roc(df_complete$mrsa_binary,
               predict(m1_fair, type = 'response'), quiet = TRUE)
r2      <- roc(df_complete$mrsa_binary,
               predict(m2, type = 'response'), quiet = TRUE)
a1_fair <- ci.auc(r1_fair); a2 <- ci.auc(r2)

cat("=== Fair comparison (both models on n =", nrow(df_complete), ") ===\n")
cat(sprintf("Model 1 (biofilm only):          AUC = %.3f (95%% CI %.3f-%.3f)\n",
            a1_fair[2], a1_fair[1], a1_fair[3]))
cat(sprintf("Model 2 (biofilm + covariates):  AUC = %.3f (95%% CI %.3f-%.3f)\n",
            a2[2], a2[1], a2[3]))
cat(sprintf("Delta AUC:                        %.3f\n\n", a2[2] - a1_fair[2]))

# DeLong test on the fair comparison
dl <- roc.test(r1_fair, r2, method = "delong")
cat(sprintf("DeLong test (Model 1 vs Model 2): p = %.4f\n\n", dl$p.value))

# ================= ORIGINAL (Model 1 on full n = 107) =================
m1_full <- glm(mrsa_binary ~ biofilm_score, data = df, family = binomial)
r1_full <- roc(df$mrsa_binary, predict(m1_full, type = 'response'), quiet = TRUE)
a1_full <- ci.auc(r1_full)

cat("=== Model 1 on full dataset (n = 107), for reference ===\n")
cat(sprintf("Model 1 (biofilm only):          AUC = %.3f (95%% CI %.3f-%.3f)\n\n",
            a1_full[2], a1_full[1], a1_full[3]))

# ================= PLOT =================
pdf("fig7_roc_curves.pdf", width = 7, height = 6)
plot(r1_full, col = "gray40", lwd = 1.5, lty = 2, main = "",
     xlab = "1 - Specificity", ylab = "Sensitivity")
lines(r1_fair, col = "#2E86AB", lwd = 2)
lines(r2,      col = "#D64933", lwd = 2)
abline(a = 0, b = 1, lty = 3, col = "gray70")
legend("bottomright",
       legend = c(
         sprintf("Biofilm only (n = 107, AUC = %.3f)", a1_full[2]),
         sprintf("Biofilm only (n = 99, AUC = %.3f)",  a1_fair[2]),
         sprintf("Biofilm + age + sex + source (n = 99, AUC = %.3f)", a2[2])
       ),
       col = c("gray40","#2E86AB","#D64933"),
       lwd = c(1.5, 2, 2), lty = c(2, 1, 1), bty = "n", cex = 0.85)
dev.off()

# ================= SAVE =================
write.csv(data.frame(
  Model    = c("Biofilm only (n=107)",
               "Biofilm only (n=99)",
               "Biofilm + age + sex + source (n=99)"),
  n        = c(nrow(df), nrow(df_complete), nrow(df_complete)),
  AUC      = c(a1_full[2], a1_fair[2], a2[2]),
  CI_lower = c(a1_full[1], a1_fair[1], a2[1]),
  CI_upper = c(a1_full[3], a1_fair[3], a2[3]),
  stringsAsFactors = FALSE
), "Supplementary_Table_S11_ROC_covariate_comparison.csv", row.names = FALSE)

# Save TIFF version for BMC
tiff("fig7_roc_curves.tif", width = 7, height = 6, units = "in",
     res = 300, compression = "lzw")
plot(r1_full, col = "gray40", lwd = 1.5, lty = 2, main = "",
     xlab = "1 - Specificity", ylab = "Sensitivity")
lines(r1_fair, col = "#2E86AB", lwd = 2)
lines(r2,      col = "#D64933", lwd = 2)
abline(a = 0, b = 1, lty = 3, col = "gray70")
legend("bottomright",
       legend = c(
         sprintf("Biofilm only (n = 107, AUC = %.3f)", a1_full[2]),
         sprintf("Biofilm only (n = 99, AUC = %.3f)",  a1_fair[2]),
         sprintf("Biofilm + age + sex + source (n = 99, AUC = %.3f)", a2[2])
       ),
       col = c("gray40","#2E86AB","#D64933"),
       lwd = c(1.5, 2, 2), lty = c(2, 1, 1), bty = "n", cex = 0.85)
dev.off()
cat("Saved: fig7_roc_curves.tif\n")

cat("Saved: fig7_roc_curves.pdf, Supplementary_Table_S11_ROC_covariate_comparison.csv\n")
