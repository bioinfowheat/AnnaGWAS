#!/usr/bin/env python3
"""
Prepare phenotype file for PLINK2 GWAS analysis.
Converts user phenotype file (0/1 coding) to PLINK format (1/2 coding).
"""

import pandas as pd

# Read sample IDs from FAM file
fam = pd.read_csv(snakemake.input.fam, sep='\s+', header=None,
                  names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])

# Read user phenotype file
# Expected format: sample_id<tab>phenotype (0 or 1)
pheno = pd.read_csv(snakemake.input.pheno, sep='\t', header=0)
pheno.columns = ['IID', 'PHENO1']

# Convert 0/1 to 1/2 for PLINK (1=control, 2=case)
pheno['PHENO1'] = pheno['PHENO1'] + 1

# Merge with FAM file to ensure correct sample order
merged = fam[['FID', 'IID']].merge(pheno, on='IID', how='left')

# Fill missing phenotypes with -9 (PLINK missing code)
merged['PHENO1'] = merged['PHENO1'].fillna(-9).astype(int)

# Save in PLINK phenotype format
merged[['FID', 'IID', 'PHENO1']].to_csv(
    snakemake.output[0], sep='\t', index=False, header=True
)

print(f"Phenotype file created: {snakemake.output[0]}")
print(f"  Total samples: {len(merged)}")
print(f"  Controls (1): {sum(merged['PHENO1'] == 1)}")
print(f"  Cases (2): {sum(merged['PHENO1'] == 2)}")
print(f"  Missing (-9): {sum(merged['PHENO1'] == -9)}")
