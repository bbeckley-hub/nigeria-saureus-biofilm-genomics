#!/bin/bash
# Remove all plot titles from original R scripts

echo "Removing titles from R scripts..."

# List of scripts to fix
scripts=(
    "analysis_biofilm_by_sample_source.R"
    "analysis_biofilm_by_st.R"
    "analysis_biofilm_vs_amr.R"
    "analysis_biofilm_vs_amr_classes_fixed.R"
    "analysis_mrsa_vs_mssa.R"
    "analysis_roc_logistic.R"
    "analysis_roc_comparison.R"
)

for script in "${scripts[@]}"; do
    if [ -f "$script" ]; then
        # Remove title lines containing 'title = "..."' and replace with no title
        # For ggplots: replace title = "..." with title = NA or remove it
        sed -i 's/title = "[^"]*"/title = NA/g' "$script"
        sed -i 's/labs(title = "[^"]*", /labs(/g' "$script"
        sed -i 's/main = "[^"]*"/main = NA/g' "$script"
        sed -i '/plot.title = element_text/d' "$script"
        echo "Fixed: $script"
    fi
done

echo "All scripts updated. Re-run your analyses to generate clean figures."
