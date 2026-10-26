#!/usr/bin/env Rscript
# Figure 8: ROC comparison of the 24-gene biofilm-associated gene carriage score,
# sdrD presence/absence, and MLST for predicting MRSA status.

library(pROC)

df <- read.csv("biofilm_analysis_dataset_expanded.csv", stringsAsFactors = FALSE)
df$mrsa_binary <- ifelse(df$MRSA_Status == "MRSA", 1, 0)

m1 <- glm(mrsa_binary ~ biofilm_score, data = df, family = binomial)
m2 <- glm(mrsa_binary ~ sdrD,           data = df, family = binomial)
m3 <- glm(mrsa_binary ~ MLST,           data = df, family = binomial)

r1 <- roc(df$mrsa_binary, predict(m1, type = "response"), quiet = TRUE)
r2 <- roc(df$mrsa_binary, predict(m2, type = "response"), quiet = TRUE)
r3 <- roc(df$mrsa_binary, predict(m3, type = "response"), quiet = TRUE)

comp12 <- roc.test(r1, r2, method = "delong")
comp13 <- roc.test(r1, r3, method = "delong")
comp23 <- roc.test(r2, r3, method = "delong")

cat(sprintf("AUC 24-gene score: %.3f\n", auc(r1)))
cat(sprintf("AUC sdrD:          %.3f\n", auc(r2)))
cat(sprintf("AUC MLST:          %.3f\n", auc(r3)))
cat(sprintf("Score vs sdrD:     p = %.3g\n", comp12$p.value))
cat(sprintf("Score vs MLST:     p = %.3g\n", comp13$p.value))
cat(sprintf("sdrD vs MLST:      p = %.3g\n", comp23$p.value))

# Save PDF
pdf("fig8_roc_comparison.pdf", width = 7, height = 6)
plot(r1, col = "#2E86AB", lwd = 2, main = "",
     xlab = "1 - Specificity", ylab = "Sensitivity")
lines(r2, col = "#D64933", lwd = 2)
lines(r3, col = "#4A934A", lwd = 2)
abline(a = 0, b = 1, lty = 3, col = "gray70")
legend("bottomright",
       legend = c(sprintf("24-gene score (AUC = %.3f)", auc(r1)),
                  sprintf("sdrD only (AUC = %.3f)",      auc(r2)),
                  sprintf("MLST (AUC = %.3f)",          auc(r3))),
       col = c("#2E86AB", "#D64933", "#4A934A"),
       lwd = 2, bty = "n", cex = 0.9)
dev.off()

# Save TIFF (BMC-compatible)
tiff("fig8_roc_comparison.tif", width = 7, height = 6, units = "in",
     res = 300, compression = "lzw")
plot(r1, col = "#2E86AB", lwd = 2, main = "",
     xlab = "1 - Specificity", ylab = "Sensitivity")
lines(r2, col = "#D64933", lwd = 2)
lines(r3, col = "#4A934A", lwd = 2)
abline(a = 0, b = 1, lty = 3, col = "gray70")
legend("bottomright",
       legend = c(sprintf("24-gene score (AUC = %.3f)", auc(r1)),
                  sprintf("sdrD only (AUC = %.3f)",      auc(r2)),
                  sprintf("MLST (AUC = %.3f)",          auc(r3))),
       col = c("#2E86AB", "#D64933", "#4A934A"),
       lwd = 2, bty = "n", cex = 0.9)
dev.off()

# Save results table
write.csv(
  data.frame(
    Model = c("24-gene score", "sdrD", "MLST"),
    AUC   = round(c(auc(r1), auc(r2), auc(r3)), 3),
    stringsAsFactors = FALSE
  ),
  "fig8_roc_comparison_results.csv", row.names = FALSE
)

cat("Saved: fig8_roc_comparison.pdf, fig8_roc_comparison.tif, fig8_roc_comparison_results.csv\n")
