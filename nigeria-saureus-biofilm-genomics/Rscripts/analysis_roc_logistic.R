#!/usr/bin/env Rscript
# ROC and logistic regression for MRSA prediction

library(pROC)
library(dplyr)

df <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)
df_log <- df %>% filter(MLST != 'Not typed')

# Model 1: biofilm score only
model1 <- glm(mrsa_binary ~ biofilm_score, data = df, family = binomial)
summary(model1)

# Model 2: biofilm score + MLST
model2 <- glm(mrsa_binary ~ biofilm_score + MLST, data = df_log, family = binomial)
summary(model2)

# Predictions and ROC
pred1 <- predict(model1, type = 'response')
pred2 <- predict(model2, type = 'response')
roc1 <- roc(df$mrsa_binary, pred1)
roc2 <- roc(df_log$mrsa_binary, pred2)
auc1 <- ci.auc(roc1)
auc2 <- ci.auc(roc2)

cat("\n=== ROC Analysis ===\n")
cat("Model 1 (biofilm_score only): AUC =", round(auc1[2], 3), 
    " (95% CI:", round(auc1[1], 3), "-", round(auc1[3], 3), ")\n")
cat("Model 2 (biofilm_score + MLST): AUC =", round(auc2[2], 3), 
    " (95% CI:", round(auc2[1], 3), "-", round(auc2[3], 3), ")\n")

# ROC plot
pdf("fig6_roc_curves.pdf", width = 7, height = 6)
plot(roc1, col = "#2E86AB", lwd = 2, main = "",
     xlab = "1 - Specificity", ylab = "Sensitivity")
lines(roc2, col = "#D64933", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray50")
legend("bottomright", 
       legend = c(paste("Biofilm score only (AUC =", round(auc1[2], 3), ")"),
                  paste("Biofilm score + MLST (AUC =", round(auc2[2], 3), ")")),
       col = c("#2E86AB", "#D64933"), lwd = 2, bty = "n")
dev.off()
cat("ROC plot saved as fig6_roc_curves.pdf\n")

# Odds ratio tables
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

write.csv(or_table(model1), "model1_odds_ratios.csv", row.names = FALSE)
write.csv(or_table(model2), "model2_odds_ratios.csv", row.names = FALSE)
cat("Odds ratio tables saved.\n")
