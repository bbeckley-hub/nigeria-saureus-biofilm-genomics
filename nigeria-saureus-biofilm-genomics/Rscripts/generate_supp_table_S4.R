# generate_supp_table_S4.R
library(pROC)
df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == 'MRSA', 1, 0)
model <- glm(mrsa_binary ~ biofilm_score, data = df, family = binomial)
pred <- predict(model, type = 'response')
roc_obj <- roc(df$mrsa_binary, pred)
coords <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy"))
pred_class <- ifelse(pred >= coords$threshold, 1, 0)
cm <- table(Predicted = pred_class, Actual = df$mrsa_binary)

write.csv(data.frame(metric = names(coords), value = as.numeric(coords)), 
          "Supplementary_Table_S4_Youden.csv", row.names = FALSE)
write.csv(as.data.frame.matrix(cm), 
          "Supplementary_Table_S4_confusion_matrix.csv", row.names = FALSE)
print("Supplementary_Table_S4_Youden.csv and Supplementary_Table_S4_confusion_matrix.csv created")
