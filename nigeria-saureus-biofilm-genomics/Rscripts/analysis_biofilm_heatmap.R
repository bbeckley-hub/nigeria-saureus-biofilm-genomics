#!/usr/bin/env Rscript
# Analysis 5: Heatmap of biofilm gene presence/absence

library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)

# ----------------------------------------------------------------------
# 1. Read data
# ----------------------------------------------------------------------
df <- read.csv('biofilm_analysis_dataset.csv', stringsAsFactors = FALSE)

# ----------------------------------------------------------------------
# 2. Prepare biofilm gene matrix (rows = isolates, columns = genes)
# ----------------------------------------------------------------------
gene_cols <- c('icaA', 'icaB', 'icaC', 'icaD', 'icaR',
               'fnbA', 'fnbB', 'clfA', 'clfB', 'sdrC', 'sdrD', 'sdrE',
               'ebp', 'cna')

# Convert to matrix (values 0/1)
mat <- as.matrix(df[, gene_cols])
rownames(mat) <- df$staphscope_id

# ----------------------------------------------------------------------
# 3. Create annotation data frame for columns (optional: gene categories)
# ----------------------------------------------------------------------
# For simplicity, we'll annotate rows with metadata
annotation_row <- data.frame(
  MRSA_Status = df$MRSA_Status,
  MLST = df$MLST,
  Source = df$SAMPLE
)
rownames(annotation_row) <- df$staphscope_id

# Simplify MLST for colour display (highlight major STs)
annotation_row$MLST_group <- ifelse(annotation_row$MLST %in% c('ST772', 'ST152', 'ST8', 'ST30', 'ST508', 'ST1', 'ST15', 'ST121', 'ST789'),
                                    annotation_row$MLST, 'Other')

# Simplify source categories (as previously defined)
annotation_row$Source_group <- case_when(
  grepl('wound', annotation_row$Source, ignore.case = TRUE) ~ 'Wound',
  grepl('Blood', annotation_row$Source, ignore.case = TRUE) ~ 'Blood',
  grepl('urine|catheter', annotation_row$Source, ignore.case = TRUE) ~ 'Urine/Catheter',
  grepl('aspirate|abscess', annotation_row$Source, ignore.case = TRUE) ~ 'Aspirate/Abscess',
  TRUE ~ 'Other'
)

# Select annotation columns for heatmap
annotation_for_plot <- annotation_row[, c('MRSA_Status', 'MLST_group', 'Source_group')]
colnames(annotation_for_plot) <- c('MRSA', 'MLST', 'Source')

# Define annotation colours
ann_colors <- list(
  MRSA = c('MSSA' = '#2E86AB', 'MRSA' = '#D64933'),
  MLST = c('ST772' = '#1B9E77', 'ST152' = '#D95F02', 'ST8' = '#7570B3',
           'ST30' = '#E7298A', 'ST508' = '#66A61E', 'ST1' = '#E6AB02',
           'ST15' = '#A6761D', 'ST121' = '#666666', 'ST789' = '#E41A1C',
           'Other' = 'gray80'),
  Source = c('Wound' = '#F8766D', 'Blood' = '#00BA38', 'Urine/Catheter' = '#619CFF',
             'Aspirate/Abscess' = '#F564E3', 'Other' = '#B79F00')
)

# ----------------------------------------------------------------------
# 4. Generate heatmap
# ----------------------------------------------------------------------
pheatmap(
  mat,
  main = "Biofilm gene presence/absence across isolates",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,          # Too many isolates to show names
  show_colnames = TRUE,
  annotation_row = annotation_for_plot,
  annotation_colors = ann_colors,
  color = c("white", "#2c7bb6"),
  legend = TRUE,
  fontsize_col = 10,
  filename = "fig5_biofilm_gene_heatmap.pdf",
  width = 10,
  height = 12,
  cellwidth = 20,
  cellheight = 8
)

cat("Heatmap saved as fig5_biofilm_gene_heatmap.pdf\n")

# ----------------------------------------------------------------------
# 5. Generate a simplified heatmap with only key isolates (e.g., top STs)
# ----------------------------------------------------------------------
# Filter to isolates from major STs (≥5 isolates)
major_sts <- c('ST772', 'ST152', 'ST30', 'ST508', 'ST8', 'ST121', 'ST789', 'ST1', 'ST15')
df_major <- df %>% filter(MLST %in% major_sts)
mat_major <- as.matrix(df_major[, gene_cols])
rownames(mat_major) <- df_major$staphscope_id

annotation_row_major <- annotation_row[rownames(mat_major), c('MRSA_Status', 'MLST_group')]
colnames(annotation_row_major) <- c('MRSA', 'MLST')

pheatmap(
  mat_major,
  main = "Biofilm gene presence/absence (major STs only)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_row = annotation_row_major,
  annotation_colors = ann_colors[c('MRSA', 'MLST')],
  color = c("white", "#2c7bb6"),
  filename = "fig5_biofilm_gene_heatmap_major_STs.pdf",
  width = 8,
  height = 10,
  cellwidth = 20,
  cellheight = 8
)
cat("Heatmap for major STs saved as fig5_biofilm_gene_heatmap_major_STs.pdf\n")
