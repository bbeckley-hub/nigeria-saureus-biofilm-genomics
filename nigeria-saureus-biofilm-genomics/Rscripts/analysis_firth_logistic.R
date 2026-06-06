#!/usr/bin/env Rscript
# Firth's penalised logistic regression for MRSA prediction

library(logistf)
library(dplyr)

# Read data
df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)

# Remove 'Not typed' MLST
df_log <- df %>% filter(MLST != 'Not typed')
df_log$MLST <- factor(df_log$MLST)

# Firth's penalised logistic regression
fit_firth <- logistf(mrsa_binary ~ biofilm_score + MLST, data = df_log)

# Summary
summary(fit_firth)

# Odds ratios
or_firth <- exp(coef(fit_firth))
ci_firth <- exp(confint(fit_firth))

results_firth <- data.frame(
  Predictor = names(coef(fit_firth)),
  OR = round(or_firth, 3),
  CI_lower = round(ci_firth[,1], 3),
  CI_upper = round(ci_firth[,2], 3),
  p_value = round(fit_firth$prob, 4)
)
print(results_firth)

write.csv(results_firth, "model_firth_odds_ratios.csv", row.names = FALSE)
cat("Firth's penalised logistic regression results saved.\n")
