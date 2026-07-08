# Supplementary Data for: Genomic profiling of biofilm‑associated genes in clinical *Staphylococcus aureus* from Nigeria

This repository contains all supplementary data, analysis scripts, and results for the manuscript:

**Genomic profiling of biofilm‑associated genes in clinical *Staphylococcus aureus* from Nigeria: association with multidrug resistance and clonal lineages**

**Authors:** Brown Beckley^1,2*†^, Christiana Jesumirhewe^3^, Solomon Babajide Bamiro^4^, Vincent Amarh^1^, Abdurrahman Hassan Jibril^5^, Innocent Afeke^6^, Adesola Olalekan^7*†^

**Affiliations:**
^1^ Department of Medical Biochemistry, University of Ghana Medical School, Accra, Ghana  
^2^ Department of Biochemistry and Biotechnology, Kwame Nkrumah University of Science and Technology (KNUST), Kumasi, Ghana  
^3^ Department of Pharmaceutical Microbiology and Biotechnology, Faculty of Pharmacy, University of Benin, Benin City, Edo State, Nigeria  
^4^ Department of Obstetrics and Gynaecology, College of Medicine, University of Lagos, Nigeria  
^5^ Department of Veterinary Public Health and Preventive Medicine, Faculty of Veterinary Medicine, Usmanu Danfodiyo University, Sokoto, Nigeria  
^6^ Department of Medical Laboratory Sciences, School of Allied Health Sciences, University of Health and Allied Sciences, Ho, Ghana  
^7^ Department of Medical Laboratory Science, College of Medicine of the University of Lagos  

*Corresponding authors: brownbeckley94@gmail.com, adesolaolalekan@unilag.edu.ng  
†These authors contributed equally to this work

---

## Repository Structure

```
.
├── Rscripts/                          # All analysis scripts
│   ├── analysis_*.R                   # R scripts for statistical analyses
│   ├── build_biofilm_dataset_expanded.py  # Python script to build dataset
│   ├── convert_figures.sh             # Convert figures to multiple formats
│   ├── R_packages.txt                 # Required R packages
│   └── run_all_analyses.sh            # Master script to run all analyses
│
├── supplementary/                     # Supplementary data and outputs
│   ├── Supplementary files/           # Organised supplementary tables & figures
│   │   ├── Supplementary Data/        # Raw StaphScope output files
│   │   ├── Supplementary_Table_S*.csv # Supplementary tables (S1–S9)
│   │   └── Supplementary_Figure_S*.pdf# Supplementary figures (S1–S3)
│   │
│   ├── figures/                       # Main manuscript figures (1–8)
│   │   ├── Figure_1_*.pdf/png/tif
│   │   ├── Figure_2_*.pdf/png/tif
│   │   ├── ...
│   │   └── Figure_8_*.pdf/png/tif
│   │
│   └── extras/                        # Intermediate files and additional outputs
│       ├── amr_genes.csv
│       ├── bacmet_genes.csv
│       ├── biofilm_analysis_dataset_expanded.csv
│       ├── sample_overview.csv
│       └── virulence_genes.csv
│
└── README.md                          # This file
```

---

## Supplementary Tables (CSV)

| Table | Description |
|-------|-------------|
| **S1** | Clinical metadata for all 107 isolates (age, sex, sample source, hospital, diagnosis) |
| **S2** | Pairwise Dunn's post‑hoc comparisons for biofilm score by MLST |
| **S3** | Spearman correlations between biofilm score and AMR functional classes |
| **S4** | Confusion matrix and performance metrics at the Youden cut‑off |
| **S5** | Full output of Firth's penalised logistic regression |
| **S6** | Partial Spearman correlation results between biofilm score and AMR gene count, controlling for MLST |
| **S7** | ROC comparison: 24‑gene biofilm score vs. *sdrD* alone vs. MLST |
| **S8** | Random forest variable importance for 24 biofilm‑associated genes |
| **S9** | Definition of AMR genes classified as the "Other" functional class |

---

## Supplementary Figures

| Figure | Description |
|--------|-------------|
| **S1** | Heatmap of biofilm gene presence/absence for major STs (≥5 isolates) |
| **S2** | Random forest variable importance plot for 24 biofilm‑associated genes |
| **S3** | Virulence burden by biofilm score group (high vs low) |

---

## Supplementary Data (StaphScope Output Files)

| File | Description |
|------|-------------|
| `amr_genes.csv` | Antimicrobial resistance genes detected |
| `bacmet_genes.csv` | Biocide and heavy metal resistance genes detected |
| `biofilm_analysis_dataset.csv` | Merged analysis dataset |
| `sample_overview.csv` | Sample-level summary of typing, AMR, and virulence genes |
| `staphscope_ultimate_report.html` | Full StaphScope HTML report |
| `virulence_genes.csv` | Virulence genes detected |

---

## How to Reproduce the Analysis

### 1. Install Required R Packages

```r
install.packages(c("ggplot2", "dplyr", "pROC", "logistf", "rstatix", "pheatmap", "corrplot", "ppcor", "randomForest"))
```

Or use the provided `R_packages.txt`:

```bash
Rscript -e "install.packages(readLines('R_packages.txt'), repos='https://cloud.r-project.org')"
```

### 2. Run the Complete Analysis Pipeline

```bash
cd Rscripts
bash run_all_analyses.sh
```

This will:
- Build the expanded biofilm dataset (Python)
- Run all baseline R analyses (Figures 1–7, Table 1)
- Run supplementary analyses (partial correlation, ROC comparison, random forest, virulence burden)

### 3. Convert Figures to Multiple Formats (Optional)

```bash
bash convert_figures.sh
```

This converts PDF figures to PNG, TIFF, and SVG formats using ImageMagick and pdf2svg.

---

## Main Figures

| Figure | Description | Generated By |
|--------|-------------|--------------|
| **1** | Biofilm genetic score in MRSA and MSSA | `analysis_mrsa_vs_mssa.R` |
| **2** | Biofilm genetic score by MLST sequence type | `analysis_biofilm_by_st.R` |
| **3** | Biofilm genetic score by sample source | `analysis_biofilm_by_sample_source.R` |
| **4** | Correlation between biofilm score and total AMR gene count | `analysis_biofilm_vs_amr.R` |
| **5** | Heatmap of biofilm gene presence/absence (all isolates) | `analysis_biofilm_heatmap.R` |
| **6** | ROC curves for predicting MRSA status | `analysis_roc_logistic.R` |
| **7** | Correlation between biofilm score and AMR functional classes | `analysis_biofilm_vs_amr_classes.R` |
| **8** | ROC comparison: 24‑gene score vs. *sdrD* vs. MLST | `analysis_roc_comparison.R` |

---

## Citation

If you use these data or scripts, please cite the original manuscript (DOI to be added) and this repository:

> Beckley B, Jesumirhewe C, Bamiro SB, Amarh V, Jibril AH, Afeke I, Olalekan A. Supplementary data for: Genomic profiling of biofilm‑associated genes in clinical *Staphylococcus aureus* from Nigeria: association with multidrug resistance and clonal lineages. GitHub; 2026. Available from: https://github.com/bbeckley-hub/nigeria-saureus-biofilm-genomics

### Data Availability

Raw sequencing reads are available at NCBI under BioProject **PRJNA1235457**.

All processed data and code are available at:  
https://github.com/bbeckley-hub/nigeria-saureus-biofilm-genomics

---

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

---

## Contact

For questions about the data or analysis scripts, please contact:

- Brown Beckley: brownbeckley94@gmail.com
- Adesola Olalekan: adesolaolalekan@unilag.edu.ng

---

**Last updated:** July 2026

