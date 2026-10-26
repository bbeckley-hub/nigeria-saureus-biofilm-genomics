# Supplementary Data and Analysis Repository

**Manuscript:** Genomic profiling of biofilm-associated genes in clinical 
*Staphylococcus aureus* from Nigeria: association with multidrug resistance 
and clonal lineages

**Authors:** Brown Beckley^1,2,*,dagger^, Christiana Jesumirhewe^3^, 
Solomon Babajide Bamiro^4^, Vincent Amarh^1^, Abdurrahman Hassan Jibril^5^, 
Innocent Afeke^6^, Adesola Olalekan^7,*,dagger^

**Affiliations:**
1. Department of Medical Biochemistry, University of Ghana Medical School, Accra, Ghana
2. Department of Biochemistry and Biotechnology, Kwame Nkrumah University of 
   Science and Technology (KNUST), Kumasi, Ghana
3. Department of Pharmaceutical Microbiology and Biotechnology, Faculty of Pharmacy, 
   University of Benin, Benin City, Edo State, Nigeria
4. Department of Obstetrics and Gynaecology, College of Medicine, University of Lagos, Nigeria
5. Department of Veterinary Public Health and Preventive Medicine, Faculty of 
   Veterinary Medicine, Usmanu Danfodiyo University, Sokoto, Nigeria
6. Department of Medical Laboratory Sciences, School of Allied Health Sciences, 
   University of Health and Allied Sciences, Ho, Ghana
7. Department of Medical Laboratory Science, College of Medicine of the University 
   of Lagos, Lagos, Nigeria

\* Corresponding authors: brownbeckley94@gmail.com, adesolaolalekan@unilag.edu.ng
\dagger These authors contributed equally to this work

**Repository:** https://github.com/bbeckley-hub/nigeria-saureus-biofilm-genomics
**Raw sequencing reads:** NCBI BioProject PRJNA1235457 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1235457)
**Licence:** Creative Commons Attribution 4.0 International (CC BY 4.0)

---

## Revision note (this version)

This repository accompanies the revised submission of the manuscript to 
BMC Microbiology. Substantive changes from the original submission:

1. **AMR data cleaning.** The original multi-database AMR output contained
   cross-database duplication, allele-name variants, virulence genes, 
   chromosomal point mutations and five housekeeping genes (dha1, dhap, kdpd,
   rlmh, sepa). The cleaned pipeline (scripts 00a, 00b, 19) applies category
   filtering, name canonicalisation, cross-database deduplication and a final 
   housekeeping-gene removal. The resulting per-isolate AMR determinant count 
   ranges from 9 to 27 (median 21, mean 18.8), down from the previously 
   reported 49 to 133.

2. **coa correction.** *coa* is significantly associated with MRSA status 
   (54/58, 93.1% in MRSA vs 33/49, 67.3% in MSSA; p = 0.001), not *sdrD* alone.
   All sections updated.

3. **AMR-biofilm association reframed.** The unadjusted correlation is 
   rho = 0.42 (95% CI 0.22-0.57, p = 6.0e-06). Within-group correlations were
   not significant (MRSA rho = 0.22, p = 0.10; MSSA rho = 0.15, p = 0.30),
   so the association is driven by the between-group MRSA-MSSA contrast.
   After MLST adjustment, partial rho = 0.25 (95% CI 0.04-0.41, p = 0.009).

4. **Discrimination analyses relabelled exploratory.** No cross-validation 
   or independent validation was performed. All AUCs are apparent in-sample
   performance.

5. **Figure renumbering.** Manuscript Figure 6 is now AMR functional class
   correlations; Figure 7 is now ROC curves (biofilm-only vs biofilm + 
   covariates); Figure 8 is ROC comparison of the 24-gene score vs sdrD vs MLST.

6. **New supplementary items.** Tables S10 (sensitivity analysis), S11 
   (ROC covariate comparison) and S12 (ROC comparison results) are new.
   Figure S4 (all 12 AMR class correlations, for transparency) is new.

---

## Repository structure

```
.
|-- scripts/                          Analysis scripts (numbered 00-19)
|   |-- 00a_clean_AMR_layer1_3.R      AMR cleaning: category filter + canonicalisation
|   |-- 00b_clean_AMR_pass2.R         Second-pass name consolidation
|   |-- 00c_build_biofilm_dataset.py  Build merged biofilm + metadata dataset
|   |-- 01_analysis_fig1_mrsa_vs_mssa.R
|   |-- 02_analysis_fig2_biofilm_by_MLST.R
|   |-- 03_analysis_fig3_biofilm_by_source.R
|   |-- 04_analysis_fig4_biofilm_vs_AMR_cleaned.R
|   |-- 05_analysis_fig5_biofilm_heatmap.R
|   |-- 06_analysis_fig6_AMR_class_correlations.R
|   |-- 07_analysis_fig7_ROC_covariates.R
|   |-- 08_analysis_fig8_ROC_comparison.R
|   |-- 09_analysis_table1_gene_prevalence.R
|   |-- 10_analysis_SuppS6_partial_correlation.R
|   |-- 11_analysis_SuppS2_random_forest.R
|   |-- 12_analysis_SuppS3_virulence_burden.R
|   |-- 13_analysis_SuppS5_firth.R
|   |-- 14_analysis_SuppS4_youden.R
|   |-- 15_supp_table_S2_dunn_pairs.R
|   |-- 16_supp_table_S4_confusion.R
|   |-- 17_supp_table_S9_other_class.R
|   |-- 18_supp_table_S10_sensitivity.R
|   |-- 19_postfilter_housekeeping.R  Final pass: remove dha1, dhap, kdpd, rlmh, sepa
|   |-- run_all_analyses.sh           Master runner
|   |-- convert_figures.sh            PDF -> TIFF conversion helper
|   |-- R_packages.txt                Required R packages with verification block
|
|-- figures/                          Main manuscript figures (TIFF only)
|   |-- fig1_biofilm_score_MRSA_vs_MSSA.tif
|   |-- fig2_biofilm_score_by_MLST.tif
|   |-- fig3_biofilm_score_by_source.tif
|   |-- fig4_biofilm_vs_AMR_count.tif
|   |-- fig5_biofilm_gene_heatmap.tif
|   |-- fig6_biofilm_vs_amr_class_corr.tif
|   |-- fig7_roc_curves.tif
|   |-- fig8_roc_comparison.tif
|
|-- supplementary/                    Tables and supplementary figures
|   |-- Supplementary_Table_S1_metadata.csv / .xlsx
|   |-- Supplementary_Table_S2_Dunn_test.csv
|   |-- Supplementary_Table_S3_correlations.csv
|   |-- Supplementary_Table_S4_confusion_matrix.csv
|   |-- Supplementary_Table_S4_Youden_cutoff.csv
|   |-- Supplementary_Table_S5_Firth.csv
|   |-- Supplementary_Table_S6_partial_correlation.csv
|   |-- Supplementary_Table_S7_ROC_comparison.csv
|   |-- Supplementary_Table_S8_RF_importance.csv
|   |-- Supplementary_Table_S9_AMR_other_class_definition.csv
|   |-- Supplementary_Table_S10_sensitivity_intrinsic_efflux.csv
|   |-- Supplementary_Table_S11_ROC_covariate_comparison.csv
|   |-- Supplementary_Table_S12_ROC_comparison_results.csv
|   |-- figS1_biofilm_gene_heatmap_major_STs.pdf / .tif
|   |-- figS2_random_forest_importance.pdf / .tif
|   |-- figS3_virulence_burden.pdf / .tif
|   |-- figS4_biofilm_vs_amr_class_all.pdf / .tif
|
|-- data/                             Input and cleaned datasets
|   |-- amr_genes.csv                       Raw multi-database AMR output
|   |-- amr_clean_pass2_long.csv            Cleaned AMR determinants (long format)
|   |-- amr_clean_pass2_per_isolate.csv     Cleaned AMR determinant counts per isolate
|   |-- biofilm_analysis_dataset_expanded.csv       Merged biofilm + metadata
|   |-- biofilm_analysis_dataset_cleaned_AMR.csv    Merged biofilm + cleaned AMR counts
|   |-- virulence_genes.csv                 VFDB virulence gene table
|   |-- bacmet_genes.csv                    BacMet biocide/metal gene table
|   |-- discriminating_biofilm_genes_expanded.csv   Fisher's test per gene (Table 1)
|   |-- staphscope_ultimate_report.html     Full StaphScope HTML report
|   |-- R_packages.txt                      Required R packages
|-- README.md                         This file
```

---

## Main figures and generating scripts

| Figure | Description | Script |
|--------|-------------|--------|
| 1 | Biofilm-associated gene carriage score in MRSA vs MSSA | `01_analysis_fig1_mrsa_vs_mssa.R` |
| 2 | Biofilm-associated gene carriage score by MLST sequence type | `02_analysis_fig2_biofilm_by_MLST.R` |
| 3 | Biofilm-associated gene carriage score by sample source | `03_analysis_fig3_biofilm_by_source.R` |
| 4 | Correlation between score and total AMR determinant count | `04_analysis_fig4_biofilm_vs_AMR_cleaned.R` |
| 5 | Heatmap of 24 biofilm-associated genes (all isolates) | `05_analysis_fig5_biofilm_heatmap.R` |
| 6 | Correlation between score and AMR functional classes | `06_analysis_fig6_AMR_class_correlations.R` |
| 7 | ROC curves: biofilm-only vs biofilm + age + sex + source | `07_analysis_fig7_ROC_covariates.R` |
| 8 | ROC comparison: 24-gene score vs sdrD vs MLST | `08_analysis_fig8_ROC_comparison.R` |

---

## Supplementary tables

| Table | Description | Generating script |
|-------|-------------|-------------------|
| S1 | Anonymised clinical metadata | (compiled from source records; provided as CSV and XLSX) |
| S2 | Pairwise Dunn's post-hoc comparisons by MLST | `15_supp_table_S2_dunn_pairs.R` |
| S3 | Spearman correlations with 12 AMR functional classes | `06_analysis_fig6_AMR_class_correlations.R` |
| S4 | Confusion matrix and Youden cut-off metrics | `16_supp_table_S4_confusion.R` |
| S5 | Firth's penalised logistic regression output | `13_analysis_SuppS5_firth.R` |
| S6 | Partial Spearman correlations controlling for MLST | `10_analysis_SuppS6_partial_correlation.R` |
| S7 | ROC comparison: score vs sdrD vs MLST | `08_analysis_fig8_ROC_comparison.R` |
| S8 | Random forest variable importance | `11_analysis_SuppS2_random_forest.R` |
| S9 | AMR determinants classified as "Other" | `17_supp_table_S9_other_class.R` |
| S10 | Sensitivity analysis (with and without intrinsic efflux) | `18_supp_table_S10_sensitivity.R` |
| S11 | ROC covariate comparison on n = 99 complete cases | `07_analysis_fig7_ROC_covariates.R` |
| S12 | ROC comparison results (AUC values and DeLong tests) | `08_analysis_fig8_ROC_comparison.R` |

---

## Supplementary figures

| Figure | Description | Generating script |
|--------|-------------|-------------------|
| S1 | Heatmap of biofilm-associated genes for major STs (>=5 isolates) | `05_analysis_fig5_biofilm_heatmap.R` |
| S2 | Random forest variable importance plot | `11_analysis_SuppS2_random_forest.R` |
| S3 | Virulence burden by score group (high vs low) | `12_analysis_SuppS3_virulence_burden.R` |
| S4 | Biofilm vs AMR functional class correlations (all 12 classes) | `06_analysis_fig6_AMR_class_correlations.R` |

---

## Key results

| Quantity | Value | Source |
|----------|-------|--------|
| Isolates | 107 (58 MRSA, 49 MSSA) | S1 |
| Sequence types | 16 (ST772 dominant at n = 35, 32.7%) | S1 |
| Biofilm score range | 20-23 (median 23) | Fig 1 |
| MRSA vs MSSA biofilm score | 23 vs 22, p = 3.4e-09 | Fig 1 |
| sdrD prevalence | 98.3% MRSA vs 69.4% MSSA, p < 0.001 | Table 1 |
| coa prevalence | 93.1% MRSA vs 67.3% MSSA, p = 0.001 | Table 1 |
| AMR determinants per isolate | 9-27 (median 21, mean 18.8) | amr_clean_pass2_per_isolate.csv |
| Total biofilm-AMR correlation | rho = 0.42 (95% CI 0.22-0.57, p = 6.0e-06) | Fig 4 |
| Within-MRSA correlation | rho = 0.22 (95% CI -0.15 to 0.49, p = 0.10) | Fig 4 |
| Within-MSSA correlation | rho = 0.15 (95% CI -0.07 to 0.37, p = 0.30) | Fig 4 |
| Partial correlation (MLST-adjusted) | partial rho = 0.25 (95% CI 0.04-0.41, p = 0.009) | S6 |
| Biofilm score apparent AUC | 0.791 (95% CI 0.711-0.872) | Fig 7 |
| Covariate model apparent AUC (n=99) | 0.851 (95% CI 0.775-0.926) | Fig 7 |
| MLST apparent AUC | 0.955 | Fig 8 |
| Virulence burden correlation | rho = -0.01 (95% CI -0.24 to 0.19, p = 0.89) | S3 |

---

## How to reproduce

### Requirements

- R version 4.3.3 (or newer 4.x)
- Python 3.13.5 with pandas, numpy, openpyxl
- For figure conversion (optional): ImageMagick, Ghostscript

### 1. Install R packages

The full dependency list (12 packages) and a verification block are in
`data/R_packages.txt`. To install everything:

```r
install.packages(c(
  "pROC", "logistf", "ppcor", "boot", "MASS",
  "dplyr", "tidyr",
  "ggplot2", "pheatmap", "corrplot", "randomForest",
  "rstatix"
))
```

Verify installation:

```r
source("data/R_packages.txt")
```

This prints installed versions and, if any package is missing, the exact
install command to run.

### 2. Install Python packages

```bash
pip install pandas numpy openpyxl
```

### 3. Run the analysis

The scripts must run in numerical order. The master script handles this:

```bash
cd scripts
bash run_all_analyses.sh
```

Or run each script individually:

```bash
cd scripts
Rscript 00a_clean_AMR_layer1_3.R
Rscript 00b_clean_AMR_pass2.R
Rscript 19_postfilter_housekeeping.R
Rscript 04_analysis_fig4_biofilm_vs_AMR_cleaned.R
Rscript 06_analysis_fig6_AMR_class_correlations.R
Rscript 12_analysis_SuppS3_virulence_burden.R
# ... etc (scripts 01-18)
```

Script 00c (Python) is a one-off merge; rerun only if the raw inputs change.

### 4. Convert figures (optional)

```bash
bash scripts/convert_figures.sh
```

Produces PNG previews of the TIFF figures for quick viewing. TIFFs in
`figures/` are the submission-ready versions.

---

## Analysis pipeline overview

```
amr_genes.csv (260 rows, 230 unique gene names)
       |
       |  00a_clean_AMR_layer1_3.R   (category filter + canonicalisation)
       |  00b_clean_AMR_pass2.R      (second-pass name consolidation)
       |  19_postfilter_housekeeping.R (final pass: remove 5 housekeeping genes)
       v
amr_clean_pass2_long.csv (2,014 unique gene-isolate pairs, 53 canonical genes)
amr_clean_pass2_per_isolate.csv (range 9-27, median 21, mean 18.8)
       |
       |  04_analysis_fig4_biofilm_vs_AMR_cleaned.R (merge + correlations)
       v
biofilm_analysis_dataset_cleaned_AMR.csv
       |
       |  All downstream analyses (Figures 1-8, S1-S4, Tables S2-S12)
       v
figures/  supplementary/  _archive/
```

---

## Software versions

| Tool | Version |
|------|---------|
| R | 4.3.3 |
| Python | 3.13.5 |
| StaphScope | 1.2.2 |
| pROC | 1.18.x |
| logistf | 1.26.x |
| ppcor | 1.1 |
| boot | 1.3-30 |
| dplyr | 1.1.x |
| tidyr | 1.3.x |
| ggplot2 | 3.5.x |
| pheatmap | 1.0.12 |
| randomForest | 4.7-1.2 |
| rstatix | 0.7.x |

Exact versions are printed when running `data/R_packages.txt`.

---

## Data availability

- **Raw sequencing reads:** NCBI BioProject PRJNA1235457
- **Processed data and code:** This repository
- **Correspondence:** brownbeckley94@gmail.com, adesolaolalekan@unilag.edu.ng

---

## Licence

This work is licensed under a Creative Commons Attribution 4.0 International
License (CC BY 4.0). See https://creativecommons.org/licenses/by/4.0/

---

## Citation

If you use these data or scripts, please cite:

Beckley B, Jesumirhewe C, Bamiro SB, Amarh V, Jibril AH, Afeke I, Olalekan A.
Supplementary data for: Genomic profiling of biofilm-associated genes in
clinical Staphylococcus aureus from Nigeria: association with multidrug
resistance and clonal lineages. GitHub; 2026.
https://github.com/bbeckley-hub/nigeria-saureus-biofilm-genomics

---

**Last updated:** October 2026
**Repository version:** v2.0 (revised submission)

