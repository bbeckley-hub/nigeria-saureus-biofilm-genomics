library(randomForest)
df <- read.csv("biofilm_analysis_dataset_expanded.csv")
gene_cols <- names(df)[2:25]  # columns 2‑25 contain the 27 biofilm genes
if (length(gene_cols) == 0) stop("No biofilm gene columns found.")
genes <- df[, gene_cols]
genes$MRSA <- as.factor(df$MRSA_Status)
set.seed(123)
rf <- randomForest(MRSA ~ ., data = genes, importance = TRUE, ntree = 500)
imp <- importance(rf)
write.csv(imp, "random_forest_importance.csv")
pdf("rf_importance.pdf")
varImpPlot(rf)
dev.off()
print(imp)
