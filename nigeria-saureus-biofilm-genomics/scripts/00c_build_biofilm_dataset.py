#!/usr/bin/env python3
"""
Build expanded biofilm genetic score dataset from StaphScope virulence_genes.csv
and clinical metadata for S. aureus Nigeria isolates.
Expanded panel includes 27 genes (ica operon, adhesins, regulators, enzymes, immune evasion).
"""

import pandas as pd
import re

# ----------------------------------------------------------------------
# 1. Read metadata (Excel)
# ----------------------------------------------------------------------
metadata = pd.read_excel('Meta data_AREF RDF GENOMIC STUDIES.xlsx')
metadata = metadata.rename(columns={'isolate number': 'isolate_num'})
metadata['isolate_num'] = metadata['isolate_num'].astype(int)

# ----------------------------------------------------------------------
# 2. Helper: extract numeric suffix from StaphScope sample ID
# ----------------------------------------------------------------------
def extract_numeric_suffix(sample_id: str) -> int:
    parts = sample_id.rsplit('_', 1)
    if len(parts) < 2:
        raise ValueError(f"Invalid format: {sample_id}")
    suffix = parts[1]
    match = re.search(r'\d+', suffix)
    if match is None:
        raise ValueError(f"No numeric suffix in: {sample_id}")
    return int(match.group())

# ----------------------------------------------------------------------
# 3. Load sample overview (MLST, spa, SCCmec, MRSA status)
# ----------------------------------------------------------------------
overview = pd.read_csv('sample_overview.csv')
overview['isolate_num'] = overview['Sample'].apply(extract_numeric_suffix)
overview_subset = overview[['isolate_num', 'MLST', 'spa_Type', 'SCCmec_Type', 'MRSA_Status']]

# ----------------------------------------------------------------------
# 4. Build biofilm gene presence matrix from virulence_genes.csv
# ----------------------------------------------------------------------
vir = pd.read_csv('virulence_genes.csv')

# Expanded biofilm gene list (27 genes)
biofilm_genes = [
    # Matrix synthesis (5)
    'icaA', 'icaB', 'icaC', 'icaD', 'icaR',
    # Adhesins (10)
    'fnbA', 'fnbB', 'clfA', 'clfB', 'sdrC', 'sdrD', 'sdrE', 'ebp', 'cna', 'coa',
    # Regulation (2)
    'sarA', 'hld',
    # Enzymes (5)
    'geh', 'lip', 'sspA', 'sspB', 'aur',
    # Immune evasion (3)
    'scn', 'sbi', 'adsA'
]

existing_genes = [g for g in biofilm_genes if g in vir['Gene'].values]
print(f"Found {len(existing_genes)} biofilm genes: {existing_genes}")

gene_to_samples = {}
for _, row in vir.iterrows():
    gene = row['Gene']
    if gene not in existing_genes:
        continue
    samples = set(row['Genomes'].split(';'))
    gene_to_samples[gene] = samples

all_samples = set()
for s_set in gene_to_samples.values():
    all_samples.update(s_set)
all_samples = sorted(all_samples)

biofilm_binary = pd.DataFrame(0, index=all_samples, columns=existing_genes)
for gene, samples in gene_to_samples.items():
    for s in samples:
        biofilm_binary.loc[s, gene] = 1

biofilm_binary['biofilm_score'] = biofilm_binary.sum(axis=1)
biofilm_binary['isolate_num'] = biofilm_binary.index.map(extract_numeric_suffix)

# ----------------------------------------------------------------------
# 5. Load AMR gene counts from amr_genes.csv
# ----------------------------------------------------------------------
amr = pd.read_csv('amr_genes.csv')
amr_count = {}
for _, row in amr.iterrows():
    samples = row['Genomes'].split(';')
    for s in samples:
        amr_count[s] = amr_count.get(s, 0) + 1
amr_df = pd.DataFrame(list(amr_count.items()), columns=['Sample', 'amr_gene_count'])
amr_df['isolate_num'] = amr_df['Sample'].apply(extract_numeric_suffix)
amr_df = amr_df[['isolate_num', 'amr_gene_count']]

# ----------------------------------------------------------------------
# 6. Merge all data frames
# ----------------------------------------------------------------------
merged = biofilm_binary.reset_index().rename(columns={'index': 'staphscope_id'})
merged = merged.merge(metadata, on='isolate_num', how='left')
merged = merged.merge(overview_subset, on='isolate_num', how='left')
merged = merged.merge(amr_df, on='isolate_num', how='left')

# ----------------------------------------------------------------------
# 7. Save final dataset
# ----------------------------------------------------------------------
output_file = 'biofilm_analysis_dataset_expanded.csv'
merged.to_csv(output_file, index=False)
print(f"Saved {merged.shape[0]} rows and {merged.shape[1]} columns to {output_file}")
print("Expanded biofilm genes included:", existing_genes)
