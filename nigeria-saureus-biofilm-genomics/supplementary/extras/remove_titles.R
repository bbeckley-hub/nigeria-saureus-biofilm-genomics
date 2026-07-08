# -----------------------------------------------------------------------------
# Remove titles from all main figures
# This recreates each figure with NO main title using the existing data
# -----------------------------------------------------------------------------

library(ggplot2)
library(pheatmap)
library(pROC)
library(dplyr)

df <- read.csv("biofilm_analysis_dataset_expanded.csv")
df$MRSA_Status <- factor(df$MRSA_Status, levels = c("MSSA", "MRSA"))

# -----------------------------------------------------------------------------
# Figure 1: MRSA vs MSSA
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = MRSA_Status, y = biofilm_score, fill = MRSA_Status)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "MRSA Status", y = "Biofilm Score (0–27)") +
  theme_bw() +
  theme(legend.position = "none")
ggsave("fig1_biofilm_score_MRSA_vs_MSSA.pdf", p, width = 5, height = 6)

# -----------------------------------------------------------------------------
# Figure 2: Biofilm score by MLST
# -----------------------------------------------------------------------------
st_keep <- names(table(df$MLST)[table(df$MLST) >= 5])
df_st <- df[df$MLST %in% st_keep, ]
df_st$MLST <- factor(df_st$MLST, levels = st_keep[order(st_keep)])

p <- ggplot(df_st, aes(x = MLST, y = biofilm_score, fill = MLST)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "Sequence Type (MLST)", y = "Biofilm Score (0–27)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fig2_biofilm_score_by_MLST.pdf", p, width = 8, height = 6)

# -----------------------------------------------------------------------------
# Figure 3: Biofilm score by source
# -----------------------------------------------------------------------------
df$Source <- factor(df$Source, levels = c("Wound", "Blood", "Urine/Catheter", "Aspirate/Abscess", "Other"))

p <- ggplot(df, aes(x = Source, y = biofilm_score, fill = Source)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "Sample Source", y = "Biofilm Score (0–27)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fig3_biofilm_score_by_source.pdf", p, width = 6, height = 6)

# -----------------------------------------------------------------------------
# Figure 4: AMR correlation
# -----------------------------------------------------------------------------
p <- ggplot(df, aes(x = amr_gene_count, y = biofilm_score, color = MRSA_Status)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black", linetype = "dashed") +
  labs(x = "Total AMR Gene Count", y = "Biofilm Score (0–27)") +
  theme_bw()
ggsave("fig4_biofilm_vs_AMR_count.pdf", p, width = 6, height = 5)

# -----------------------------------------------------------------------------
# Figure 5: Heatmap (all isolates)
# -----------------------------------------------------------------------------
gene_cols <- names(df)[2:25]
mat <- df[, gene_cols]
rownames(mat) <- df$staphscope_id
annotation <- df[, c("MRSA_Status", "MLST")]
annotation$MLST <- ifelse(annotation$MLST %in% st_keep, annotation$MLST, "Other")

pdf("fig5_biofilm_gene_heatmap.pdf", width = 10, height = 8)
pheatmap(mat,
         annotation_row = annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_col = 8,
         main = NA)
dev.off()

# -----------------------------------------------------------------------------
# Figure 6: ROC curves (logistic models)
# -----------------------------------------------------------------------------
df$mrsa_binary <- ifelse(df$MRSA_Status == "MRSA", 1, 0)
m1 <- glm(mrsa_binary ~ biofilm_score, data = df, family = binomial)
m2 <- glm(mrsa_binary ~ biofilm_score + Age + Sex + Source, data = df, family = binomial)

r1 <- roc(df$mrsa_binary, predict(m1, type = "response"))
r2 <- roc(df$mrsa_binary, predict(m2, type = "response"))

pdf("fig6_roc_curves.pdf", width = 6, height = 6)
plot(r1, col = "blue", main = NA, legacy.axes = TRUE)
lines(r2, col = "red")
legend("bottomright",
       legend = c("Biofilm score only (AUC = 0.791)",
                  "Multivariate (AUC = 0.797)"),
       col = c("blue", "red"), lwd = 2, bty = "n")
dev.off()

# -----------------------------------------------------------------------------
# Figure 7: AMR class correlation (use existing data)
# -----------------------------------------------------------------------------
cor_data <- read.csv("biofilm_vs_amr_class_correlations.csv")
cor_data <- cor_data[order(-abs(cor_data$Spearman_rho)), ]
cor_data$AMR_Class <- factor(cor_data$AMR_Class, levels = cor_data$AMR_Class)

p <- ggplot(cor_data, aes(x = AMR_Class, y = Spearman_rho, fill = Spearman_rho > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightcoral"), guide = "none") +
  labs(x = "AMR Functional Class", y = "Spearman's ρ") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fig8_biofilm_vs_amr_class_corr.pdf", p, width = 7, height = 5)

# -----------------------------------------------------------------------------
# Figure 8: ROC comparison (score vs sdrD vs MLST)
# -----------------------------------------------------------------------------
m_sdrD <- glm(mrsa_binary ~ sdrD, data = df, family = binomial)
m_mlst <- glm(mrsa_binary ~ MLST, data = df, family = binomial)
r_sdrD <- roc(df$mrsa_binary, predict(m_sdrD, type = "response"))
r_mlst <- roc(df$mrsa_binary, predict(m_mlst, type = "response"))

pdf("roc_comparison.pdf", width = 6, height = 6)
plot(r1, col = "blue", main = NA, legacy.axes = TRUE)
lines(r_sdrD, col = "red")
lines(r_mlst, col = "green")
legend("bottomright",
       legend = c("27-gene score (AUC = 0.791)",
                  "sdrD alone (AUC = 0.644)",
                  "MLST (AUC = 0.955)"),
       col = c("blue", "red", "green"), lwd = 2, bty = "n")
dev.off()

cat("All figures regenerated without titles.\n")
