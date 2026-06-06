#!/usr/bin/env Rscript
# Multivariate logistic regression with complete cases only

library(dplyr)
library(pROC)

# Read data
df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)

# Clean age: exclude <1
df$age_years <- df$AGE.in.years
df$age_years[df$age_years < 1] <- NA

# Sex
df$sex <- factor(df$SEX, levels = c('Male', 'Female'))

# Source grouping
df$source_group <- case_when(
  grepl('wound', df$SAMPLE, ignore.case = TRUE) ~ 'Wound',
  grepl('Blood', df$SAMPLE, ignore.case = TRUE) ~ 'Blood',
  grepl('urine|catheter', df$SAMPLE, ignore.case = TRUE) ~ 'Urine/Catheter',
  grepl('aspirate|abscess', df$SAMPLE, ignore.case = TRUE) ~ 'Aspirate/Abscess',
  TRUE ~ 'Other'
)
df$source_group <- factor(df$source_group, levels = c('Wound', 'Blood', 'Urine/Catheter', 'Aspirate/Abscess', 'Other'))

# Keep only complete cases for all variables
df_complete <- df %>%
  filter(!is.na(age_years), !is.na(sex), !is.na(source_group), !is.na(biofilm_score))

cat("Complete cases for multivariate model:", nrow(df_complete), "out of", nrow(df), "\n")

# Multivariate logistic regression
model_mv <- glm(mrsa_binary ~ biofilm_score + age_years + sex + source_group,
                data = df_complete, family = binomial)

summary(model_mv)

# Odds ratios
or_mv <- exp(coef(model_mv))
ci_mv <- exp(confint(model_mv))

results_mv <- data.frame(
  Predictor = names(coef(model_mv)),
  OR = round(or_mv, 3),
  CI_lower = round(ci_mv[,1], 3),
  CI_upper = round(ci_mv[,2], 3),
  p_value = round(summary(model_mv)$coefficients[,4], 4)
)
print(results_mv)
write.csv(results_mv, "model_multivariate_odds_ratios.csv", row.names = FALSE)

# ROC curve
pred_mv <- predict(model_mv, type = 'response')
roc_mv <- roc(df_complete$mrsa_binary, pred_mv)
auc_mv <- ci.auc(roc_mv)
cat("\nMultivariate model AUC =", round(auc_mv[2], 3), "(95% CI:", round(auc_mv[1], 3), "-", round(auc_mv[3], 3), ")\n")

pdf("fig7_multivariate_roc.pdf", width = 6, height = 6)
plot(roc_mv, col = "#2E86AB", lwd = 2, main = "ROC curve: multivariate model")
abline(a = 0, b = 1, lty = 2, col = "gray50")
legend("bottomright", legend = paste("AUC =", round(auc_mv[2], 3)), col = "#2E86AB", lwd = 2, bty = "n")
dev.off()
cat("Multivariate ROC plot saved as fig7_multivariate_roc.pdf\n")
