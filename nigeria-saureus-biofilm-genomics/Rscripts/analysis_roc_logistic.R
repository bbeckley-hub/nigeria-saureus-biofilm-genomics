#!/usr/bin/env Rscript
# Additional analysis: ROC curve and logistic regression for MRSA prediction

library(pROC)
library(ggplot2)
library(dplyr)

# ----------------------------------------------------------------------
# 1. Read data
# ----------------------------------------------------------------------
df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)

# Convert MRSA_Status to binary (1 = MRSA, 0 = MSSA)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)

# Remove 'Not typed' MLST for logistic regression (optional)
df_log <- df %>% filter(MLST != 'Not typed')

# ----------------------------------------------------------------------
# 2. Logistic regression models
# ----------------------------------------------------------------------
# Model 1: biofilm_score only
model1 <- glm(mrsa_binary ~ biofilm_score, data = df, family = binomial)
summary(model1)

# Model 2: biofilm_score + MLST
model2 <- glm(mrsa_binary ~ biofilm_score + MLST, data = df_log, family = binomial)
summary(model2)

# ----------------------------------------------------------------------
# 3. ROC curves and AUC
# ----------------------------------------------------------------------
# Predict probabilities
pred1 <- predict(model1, type = 'response')
pred2 <- predict(model2, type = 'response')

# ROC objects
roc1 <- roc(df$mrsa_binary, pred1)
roc2 <- roc(df_log$mrsa_binary, pred2)

# AUC with confidence intervals
auc1 <- ci.auc(roc1)
auc2 <- ci.auc(roc2)

cat("\n=== ROC Analysis ===\n")
cat("Model 1 (biofilm_score only): AUC =", round(auc1[2], 3), 
    " (95% CI:", round(auc1[1], 3), "-", round(auc1[3], 3), ")\n")
cat("Model 2 (biofilm_score + MLST): AUC =", round(auc2[2], 3), 
    " (95% CI:", round(auc2[1], 3), "-", round(auc2[3], 3), ")\n")

# ----------------------------------------------------------------------
# 4. Plot ROC curves using pROC's built-in plotting
# ----------------------------------------------------------------------
# Open PDF device
pdf("fig6_roc_curves.pdf", width = 7, height = 6)

# Plot ROC curve for model1
plot(roc1, col = "#2E86AB", lwd = 2, main = "ROC curves for predicting MRSA status",
     xlab = "1 - Specificity (False positive rate)",
     ylab = "Sensitivity (True positive rate)")

# Add ROC curve for model2
lines(roc2, col = "#D64933", lwd = 2)

# Add diagonal reference line
abline(a = 0, b = 1, lty = 2, col = "gray50")

# Add legend
legend("bottomright", 
       legend = c(paste("Biofilm score only (AUC =", round(auc1[2], 3), ")"),
                  paste("Biofilm score + MLST (AUC =", round(auc2[2], 3), ")")),
       col = c("#2E86AB", "#D64933"), lwd = 2, bty = "n")

dev.off()
cat("ROC plot saved as fig6_roc_curves.pdf\n")

# ----------------------------------------------------------------------
# 5. Logistic regression table (odds ratios)
# ----------------------------------------------------------------------
or_table <- function(model) {
  coefs <- coef(model)
  ors <- exp(coefs)
  ci <- exp(confint(model))
  result <- data.frame(
    Predictor = names(coefs),
    OR = round(ors, 3),
    CI_lower = round(ci[,1], 3),
    CI_upper = round(ci[,2], 3),
    p_value = round(summary(model)$coefficients[,4], 4)
  )
  return(result)
}

cat("\n=== Odds ratios for Model 1 ===\n")
print(or_table(model1))

cat("\n=== Odds ratios for Model 2 ===\n")
print(or_table(model2))

# Save table as CSV
write.csv(or_table(model1), "model1_odds_ratios.csv", row.names = FALSE)
write.csv(or_table(model2), "model2_odds_ratios.csv", row.names = FALSE)
cat("Odds ratio tables saved as CSV files.\n")
