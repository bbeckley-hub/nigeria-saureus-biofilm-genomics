#!/usr/bin/env Rscript
# Supplementary Figure S2 and Supplementary Table S8:
# Random forest variable importance for 24 biofilm-associated genes

library(randomForest)

df <- read.csv("biofilm_analysis_dataset_expanded.csv", stringsAsFactors = FALSE)

# Columns 2 to 25 are the 24 biofilm-associated genes
gene_cols <- names(df)[2:25]
if (length(gene_cols) == 0) stop("No biofilm-associated gene columns found.")
cat("Using", length(gene_cols), "biofilm-associated genes\n")

genes <- df[, gene_cols]
genes$MRSA <- as.factor(df$MRSA_Status)

set.seed(123)
rf <- randomForest(MRSA ~ ., data = genes, importance = TRUE, ntree = 500)

imp <- importance(rf)
write.csv(imp, "Supplementary_Table_S8_RF_importance.csv", row.names = TRUE)
cat("Saved: Supplementary_Table_S8_RF_importance.csv\n")

# Top five by overall MeanDecreaseAccuracy
overall <- imp[, "MeanDecreaseAccuracy"]
top5 <- sort(overall, decreasing = TRUE)[1:5]
cat("\nTop five genes by MeanDecreaseAccuracy:\n")
print(top5)

# PDF
pdf("figS2_random_forest_importance.pdf", width = 7, height = 6)
varImpPlot(rf, main = "", type = 1)
dev.off()

# TIFF (BMC-compatible)
tiff("figS2_random_forest_importance.tif", width = 7, height = 6, units = "in",
     res = 300, compression = "lzw")
varImpPlot(rf, main = "", type = 1)
dev.off()

cat("Saved: figS2_random_forest_importance.{pdf,tif}\n")
