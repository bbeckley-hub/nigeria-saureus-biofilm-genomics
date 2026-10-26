#!/usr/bin/env bash
# Final repository organization for BMC submission
# Creates figures/ (TIFF only) and supplementary/ folders
# Archives superseded files, deletes duplicates

set -euo pipefail

echo "=== Creating target directories ==="
mkdir -p figures
mkdir -p supplementary
mkdir -p _archive/figures_preview
mkdir -p _archive/stale

echo ""
echo "=== Deleting -1.png duplicates (12 files) ==="
rm -f fig1_biofilm_score_MRSA_vs_MSSA-1.png
rm -f fig2_biofilm_score_by_MLST-1.png
rm -f fig3_biofilm_score_by_source-1.png
rm -f fig4_biofilm_vs_AMR_count-1.png
rm -f fig5_biofilm_gene_heatmap-1.png
rm -f fig6_biofilm_vs_amr_class_corr-1.png
rm -f fig7_roc_curves-1.png
rm -f fig8_roc_comparison-1.png
rm -f figS1_biofilm_gene_heatmap_major_STs-1.png
rm -f figS2_random_forest_importance-1.png
rm -f figS3_virulence_burden-1.png
rm -f figS4_biofilm_vs_amr_class_all-1.png

echo ""
echo "=== Deleting PNG previews (2 files) ==="
rm -f fig4_biofilm_vs_AMR_count.png
rm -f fig6_biofilm_vs_amr_class_corr.png

echo ""
echo "=== Archiving stale outputs (4 files) ==="
[ -f random_forest_importance.csv ] && mv random_forest_importance.csv _archive/stale/
[ -f rf_importance.pdf ]           && mv rf_importance.pdf           _archive/stale/
[ -f virulence_boxplot.pdf ]       && mv virulence_boxplot.pdf       _archive/stale/
[ -f virulence_correlation.csv ]   && mv virulence_correlation.csv   _archive/stale/

echo ""
echo "=== Archiving source metadata (2 files) ==="
[ -f "Meta data_AREF RDF GENOMIC STUDIES.xlsx" ] && mv "Meta data_AREF RDF GENOMIC STUDIES.xlsx" _archive/
[ -f sample_overview.csv ]                     && mv sample_overview.csv                     _archive/stale/

echo ""
echo "=== Moving main figure TIFFs to figures/ (8 files) ==="
mv fig1_biofilm_score_MRSA_vs_MSSA.tif figures/
mv fig2_biofilm_score_by_MLST.tif       figures/
mv fig3_biofilm_score_by_source.tif     figures/
mv fig4_biofilm_vs_AMR_count.tif        figures/
mv fig5_biofilm_gene_heatmap.tif        figures/
mv fig6_biofilm_vs_amr_class_corr.tif   figures/
mv fig7_roc_curves.tif                  figures/
mv fig8_roc_comparison.tif              figures/

echo ""
echo "=== Moving main figure PDFs to _archive/figures_preview/ (8 files) ==="
mv fig1_biofilm_score_MRSA_vs_MSSA.pdf _archive/figures_preview/
mv fig2_biofilm_score_by_MLST.pdf      _archive/figures_preview/
mv fig3_biofilm_score_by_source.pdf    _archive/figures_preview/
mv fig4_biofilm_vs_AMR_count.pdf       _archive/figures_preview/
mv fig5_biofilm_gene_heatmap.pdf       _archive/figures_preview/
mv fig6_biofilm_vs_amr_class_corr.pdf  _archive/figures_preview/
mv fig7_roc_curves.pdf                 _archive/figures_preview/
mv fig8_roc_comparison.pdf             _archive/figures_preview/

echo ""
echo "=== Moving supplementary tables to supplementary/ (14 files) ==="
mv Supplementary_Table_S1_metadata.csv                        supplementary/
mv Supplementary_Table_S1_metadata.xlsx                       supplementary/
mv Supplementary_Table_S2_Dunn_test.csv                       supplementary/
mv Supplementary_Table_S3_correlations.csv                    supplementary/
mv Supplementary_Table_S4_confusion_matrix.csv                supplementary/
mv Supplementary_Table_S4_Youden_cutoff.csv                   supplementary/
mv Supplementary_Table_S5_Firth.csv                           supplementary/
mv Supplementary_Table_S6_partial_correlation.csv             supplementary/
mv Supplementary_Table_S7_ROC_comparison.csv                  supplementary/
mv Supplementary_Table_S8_RF_importance.csv                   supplementary/
mv Supplementary_Table_S9_AMR_other_class_definition.csv      supplementary/
mv Supplementary_Table_S10_sensitivity_intrinsic_efflux.csv   supplementary/
mv Supplementary_Table_S11_ROC_covariate_comparison.csv       supplementary/
mv Supplementary_Table_S12_ROC_comparison_results.csv         supplementary/

echo ""
echo "=== Moving supplementary figures to supplementary/ (8 files) ==="
mv figS1_biofilm_gene_heatmap_major_STs.pdf supplementary/
mv figS1_biofilm_gene_heatmap_major_STs.tif supplementary/
mv figS2_random_forest_importance.pdf       supplementary/
mv figS2_random_forest_importance.tif       supplementary/
mv figS3_virulence_burden.pdf               supplementary/
mv figS3_virulence_burden.tif               supplementary/
mv figS4_biofilm_vs_amr_class_all.pdf       supplementary/
mv figS4_biofilm_vs_amr_class_all.tif       supplementary/

echo ""
echo "=================================================================="
echo "=== FINAL STRUCTURE ==="
echo "=================================================================="
echo ""
echo "--- figures/ (main figures, TIFF only: expect 8) ---"
ls -1 figures/
echo ""
echo "--- supplementary/ (tables + supplementary figures: expect 22) ---"
ls -1 supplementary/
echo ""
echo "--- Root directory (scripts + core data) ---"
ls -1 *.R *.sh *.py *.csv *.txt *.html 2>/dev/null
echo ""
echo "--- _archive/ structure ---"
ls -1 _archive/
echo "  _archive/figures_preview/: $(ls _archive/figures_preview 2>/dev/null | wc -l) files"
echo "  _archive/stale/:           $(ls _archive/stale 2>/dev/null | wc -l) files"
echo ""
echo "=== Done. ==="
