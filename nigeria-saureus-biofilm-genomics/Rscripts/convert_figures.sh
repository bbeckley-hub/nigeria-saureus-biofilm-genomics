#!/bin/bash
# Convert all Fig*.pdf to TIFF, PNG, SVG at 300 dpi

for f in fig*.pdf; do
    base="${f%.pdf}"
    echo "Converting $f ..."
    # TIFF (300 dpi)
    pdftoppm -tiff -r 300 "$f" "${base}"
    # PNG (300 dpi)
    pdftoppm -png -r 300 "$f" "${base}"
    # SVG (using pdf2svg)
    pdf2svg "$f" "${base}.svg"
done

# Rename the converted files (pdftoppm adds -1 suffix)
for f in fig*.tif; do
    mv "$f" "${f/-1/}"
done
for f in Fig*.png; do
    mv "$f" "${f/-1/}"
done
