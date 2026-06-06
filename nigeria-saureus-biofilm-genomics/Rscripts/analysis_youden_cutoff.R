#!/usr/bin/env Rscript
# Youden index for biofilm score threshold to predict MRSA

library(pROC)
library(dplyr)

df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)

# Use the simple model (biofilm_score only) from earlier
model1 <- glm(mrsa_binary ~ biofilm_score, data = df, family = binomial)
pred <- predict(model1, type = 'response')

roc_obj <- roc(df$mrsa_binary, pred)
youden <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
print(youden)

# Also compute confusion matrix at that threshold
threshold <- youden$threshold
pred_class <- ifelse(pred >= threshold, 1, 0)
cm <- table(Predicted = pred_class, Actual = df$mrsa_binary)
print(cm)

# Save results
write.csv(youden, "youden_cutoff.csv", row.names = FALSE)
cat("Optimal cut-off saved to youden_cutoff.csv\n")
