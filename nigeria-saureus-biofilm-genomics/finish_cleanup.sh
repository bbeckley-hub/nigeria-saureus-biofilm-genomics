#!/usr/bin/env bash
set -euo pipefail

echo "=== Archiving old biofilm_analysis_dataset.csv ==="
if [ -f data/biofilm_analysis_dataset.csv ]; then
  mv data/biofilm_analysis_dataset.csv _archive/superseded_folders/
  echo "  Moved: data/biofilm_analysis_dataset.csv -> _archive/superseded_folders/"
fi

echo ""
echo "=== Removing duplicate archived files (identical to data/) ==="
for f in amr_genes.csv bacmet_genes.csv staphscope_ultimate_report.html virulence_genes.csv; do
  if [ -f "_archive/superseded_folders/$f" ]; then
    rm "_archive/superseded_folders/$f"
    echo "  Deleted duplicate: $f"
  fi
done

echo ""
echo "=== Adding .gitignore ==="
cat > .gitignore <<'EOF'
# Intermediate R outputs
*-1.png
Rplots.pdf
.DS_Store

# Editor files
*.swp
*~

# R session files
.RData
.Rhistory
EOF
echo "  Created: .gitignore"

echo ""
echo "=== Verification ==="
echo "--- data/ ($(ls data | wc -l) files, expected 11) ---"
ls -1 data/
echo ""
echo "--- _archive/superseded_folders/ ($(ls _archive/superseded_folders | wc -l) files) ---"
ls -1 _archive/superseded_folders/
