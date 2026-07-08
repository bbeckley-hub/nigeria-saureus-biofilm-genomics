#!/bin/bash
# Complete analysis pipeline: Python build -> all R scripts
# Usage: bash run_all_analyses.sh

set -e  # exit immediately if any command fails

echo "========================================="
echo "Starting complete analysis pipeline"
echo "========================================="

# ------------------------------------------------------------
# Step 1: Python – build the expanded biofilm dataset
# ------------------------------------------------------------
echo ""
echo "[1/3] Running Python: build_biofilm_dataset_expanded.py"
python3 build_biofilm_dataset_expanded.py

# ------------------------------------------------------------
# Step 2: Baseline R analyses (figures & main tables)
# ------------------------------------------------------------
echo ""
echo "[2/3] Running baseline R analyses..."
echo "  > analysis_mrsa_vs_mssa.R"
Rscript analysis_mrsa_vs_mssa.R

echo "  > analysis_biofilm_by_st.R"
Rscript analysis_biofilm_by_st.R

echo "  > analysis_biofilm_by_sample_source.R"
Rscript analysis_biofilm_by_sample_source.R

echo "  > analysis_biofilm_vs_amr.R"
Rscript analysis_biofilm_vs_amr.R

echo "  > analysis_biofilm_heatmap.R"
Rscript analysis_biofilm_heatmap.R

echo "  > analysis_biofilm_vs_amr_classes_fixed.R"
Rscript analysis_biofilm_vs_amr_classes_fixed.R

echo "  > analysis_roc_logistic.R"
Rscript analysis_roc_logistic.R

echo "  > analysis_discriminating_genes.R"
Rscript analysis_discriminating_genes.R

# ------------------------------------------------------------
# Step 3: New supplementary analyses (post-revision)
# ------------------------------------------------------------
echo ""
echo "[3/3] Running new supplementary analyses..."
echo "  > analysis_partial_correlation.R"
Rscript analysis_partial_correlation.R

echo "  > analysis_roc_comparison.R"
Rscript analysis_roc_comparison.R

echo "  > analysis_random_forest.R"
Rscript analysis_random_forest.R

echo "  > analysis_virulence_burden.R"
Rscript analysis_virulence_burden.R

# ------------------------------------------------------------
# Done
# ------------------------------------------------------------
echo ""
echo "========================================="
echo "All analyses completed successfully!"
echo "========================================="

# Show the most important output files
echo ""
echo "Key outputs generated:"
ls -1 fig*.pdf roc_comparison.pdf rf_importance.pdf virulence_boxplot.pdf \
      *correlation*.csv roc_comparison_results.csv random_forest_importance.csv \
      2>/dev/null | sed 's/^/  - /'
