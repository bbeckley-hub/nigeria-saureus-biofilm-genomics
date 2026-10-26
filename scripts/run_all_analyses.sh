#!/usr/bin/env bash
# Master runner for the complete analysis pipeline.
# Runs all R scripts from data/ (where inputs live), then moves outputs
# to figures/ and supplementary/.
#
# Usage (from repository root):
#   bash scripts/run_all_analyses.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT/data"

echo "Working directory: $(pwd)"
echo ""

run() {
  echo "--- $1 ---"
  Rscript "../scripts/$1"
}

# AMR cleaning pipeline
run 00a_clean_AMR_layer1_3.R
run 00b_clean_AMR_pass2.R
run 19_postfilter_housekeeping.R

# Figures 1-8
run 01_analysis_fig1_mrsa_vs_mssa.R
run 02_analysis_fig2_biofilm_by_MLST.R
run 03_analysis_fig3_biofilm_by_source.R
run 04_analysis_fig4_biofilm_vs_AMR_cleaned.R
run 05_analysis_fig5_biofilm_heatmap.R
run 06_analysis_fig6_AMR_class_correlations.R
run 07_analysis_fig7_ROC_covariates.R
run 08_analysis_fig8_ROC_comparison.R

# Supplementary tables and figures
run 09_analysis_table1_gene_prevalence.R
run 10_analysis_SuppS6_partial_correlation.R
run 11_analysis_SuppS2_random_forest.R
run 12_analysis_SuppS3_virulence_burden.R
run 13_analysis_SuppS5_firth.R
run 14_analysis_SuppS4_youden.R
run 15_supp_table_S2_dunn_pairs.R
run 16_supp_table_S4_confusion.R
run 17_supp_table_S9_other_class.R
run 18_supp_table_S10_sensitivity.R

# Move outputs to canonical folders
echo ""
echo "=== Distributing outputs ==="
mv fig[1-8]_*.tif              "$REPO_ROOT/figures/" 2>/dev/null       || true
mv figS*.pdf figS*.tif         "$REPO_ROOT/supplementary/" 2>/dev/null || true
mv Supplementary_Table_*.csv   "$REPO_ROOT/supplementary/" 2>/dev/null || true

# Clean up intermediate PNGs
rm -f fig*-1.png
rm -f fig4_biofilm_vs_AMR_count.png
rm -f fig6_biofilm_vs_amr_class_corr.png

echo ""
echo "=== Complete ==="
echo "Main figures:        $REPO_ROOT/figures/"
echo "Supplementary files: $REPO_ROOT/supplementary/"
echo "Working data:        $REPO_ROOT/data/"
