#!/usr/bin/env Rscript
# Heatmap of biofilm gene presence/absence with visible gene names

library(dplyr)
library(pheatmap)

df <- read.csv('biofilm_analysis_dataset_expanded.csv', stringsAsFactors = FALSE)

# Select all biofilm gene columns (2:25)
gene_cols <- names(df)[2:25]
cat("Gene columns found:\n")
print(gene_cols)

mat <- as.matrix(df[, gene_cols])
rownames(mat) <- df$staphscope_id

# Annotation rows
annotation_row <- data.frame(
  MRSA_Status = df$MRSA_Status,
  MLST = df$MLST,
  Source = df$SAMPLE
)
rownames(annotation_row) <- df$staphscope_id

annotation_row$MLST_group <- ifelse(annotation_row$MLST %in% c('ST772', 'ST152', 'ST8', 'ST30', 'ST508', 'ST1', 'ST15', 'ST121', 'ST789'),
                                    annotation_row$MLST, 'Other')
annotation_row$Source_group <- case_when(
  grepl('wound', annotation_row$Source, ignore.case = TRUE) ~ 'Wound',
  grepl('Blood', annotation_row$Source, ignore.case = TRUE) ~ 'Blood',
  grepl('urine|catheter', annotation_row$Source, ignore.case = TRUE) ~ 'Urine/Catheter',
  grepl('aspirate|abscess', annotation_row$Source, ignore.case = TRUE) ~ 'Aspirate/Abscess',
  TRUE ~ 'Other'
)

annotation_for_plot <- annotation_row[, c('MRSA_Status', 'MLST_group', 'Source_group')]
colnames(annotation_for_plot) <- c('MRSA', 'MLST', 'Source')

ann_colors <- list(
  MRSA = c('MSSA' = '#2E86AB', 'MRSA' = '#D64933'),
  MLST = c('ST772' = '#1B9E77', 'ST152' = '#D95F02', 'ST8' = '#7570B3',
           'ST30' = '#E7298A', 'ST508' = '#66A61E', 'ST1' = '#E6AB02',
           'ST15' = '#A6761D', 'ST121' = '#666666', 'ST789' = '#E41A1C',
           'Other' = 'gray80'),
  Source = c('Wound' = '#F8766D', 'Blood' = '#00BA38', 'Urine/Catheter' = '#619CFF',
             'Aspirate/Abscess' = '#F564E3', 'Other' = '#B79F00')
)

# Full heatmap with visible column names
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_row = annotation_for_plot,
  annotation_colors = ann_colors,
  color = c("white", "#2c7bb6"),
  legend = TRUE,
  fontsize_col = 7,
  fontsize_row = 6,
  angle_col = 45,                     # Rotate column labels
  filename = "fig5_biofilm_gene_heatmap.pdf",
  width = 12,                         # Increased width
  height = 14,                        # Increased height
  cellwidth = 18,
  cellheight = 7,
  main = NA
)
cat("Heatmap saved as fig5_biofilm_gene_heatmap.pdf\n")

# Heatmap for major STs only
major_sts <- c('ST772', 'ST152', 'ST30', 'ST508', 'ST8', 'ST121', 'ST789', 'ST1', 'ST15')
df_major <- df %>% filter(MLST %in% major_sts)
mat_major <- as.matrix(df_major[, gene_cols])
rownames(mat_major) <- df_major$staphscope_id

annotation_row_major <- annotation_row[rownames(mat_major), c('MRSA_Status', 'MLST_group')]
colnames(annotation_row_major) <- c('MRSA', 'MLST')

pheatmap(
  mat_major,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_row = annotation_row_major,
  annotation_colors = ann_colors[c('MRSA', 'MLST')],
  color = c("white", "#2c7bb6"),
  fontsize_col = 8,
  angle_col = 45,
  filename = "fig5_biofilm_gene_heatmap_major_STs.pdf",
  width = 10,
  height = 12,
  cellwidth = 20,
  cellheight = 8,
  main = NA
)
cat("Heatmap for major STs saved as fig5_biofilm_gene_heatmap_major_STs.pdf\n")
