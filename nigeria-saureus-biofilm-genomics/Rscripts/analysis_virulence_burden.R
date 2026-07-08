main <- read.csv("biofilm_analysis_dataset_expanded.csv")
so <- read.csv("sample_overview.csv")
df <- merge(main, so, by.x = "staphscope_id", by.y = "Sample")

count_vir <- function(x) {
  if (is.na(x) || x == "") return(0)
  genes <- unlist(strsplit(as.character(x), "[;,]"))
  length(genes[genes != ""])
}
df$virulence_count <- sapply(df$VFDB_Genes, count_vir)

cor_test <- cor.test(df$biofilm_score, df$virulence_count, method = "spearman")
print(cor_test)

df$group <- ifelse(df$biofilm_score >= median(df$biofilm_score), "High", "Low")
wt <- wilcox.test(df$virulence_count ~ df$group)
print(wt)

write.csv(data.frame(rho = cor_test$estimate, p = cor_test$p.value),
          "virulence_correlation.csv", row.names = FALSE)

pdf("virulence_boxplot.pdf")
boxplot(virulence_count ~ group, data = df,
        main = "Virulence burden by biofilm score",
        xlab = "Biofilm score group (High = 23, Low = 20-22)",
        ylab = "Number of virulence genes")
dev.off()
