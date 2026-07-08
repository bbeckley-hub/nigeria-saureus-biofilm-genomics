#!/bin/bash
# Rename figures to logical order

echo "Renaming main figures..."

# Main figures
mv fig1_biofilm_score_MRSA_vs_MSSA.pdf Figure_1_MRSA_vs_MSSA.pdf
mv fig2_biofilm_score_by_MLST.pdf Figure_2_MLST_boxplot.pdf
mv fig3_biofilm_score_by_source.pdf Figure_3_source_boxplot.pdf
mv fig4_biofilm_vs_AMR_count.pdf Figure_4_AMR_correlation.pdf
mv fig5_biofilm_gene_heatmap.pdf Figure_5_heatmap_all.pdf
mv fig6_roc_curves.pdf Figure_6_ROC_logistic.pdf
mv fig8_biofilm_vs_amr_class_corr.pdf Figure_7_AMR_class_corr.pdf
mv roc_comparison.pdf Figure_8_ROC_comparison.pdf

echo "Renaming supplementary figures..."

# Supplementary figures
mv fig5_biofilm_gene_heatmap_major_STs.pdf Supplementary_Figure_S1_heatmap_major_STs.pdf
mv rf_importance.pdf Supplementary_Figure_S2_RF_importance.pdf
mv virulence_boxplot.pdf Supplementary_Figure_S3_virulence_burden.pdf

# Optional: if Supplementary_Figure_S5.pdf exists and is redundant, remove or rename
if [ -f Supplementary_Figure_S5.pdf ]; then
    echo "Supplementary_Figure_S5.pdf already exists – check if it should be kept."
fi

echo "Done. List of renamed files:"
ls -l Figure_*.pdf Supplementary_Figure_S*.pdf 2>/dev/null
