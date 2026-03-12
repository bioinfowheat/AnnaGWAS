# GWAS Pipeline for *Pararge aegeria* Diapause Phenotype

Snakemake-based genome-wide association study pipeline for identifying genetic variants associated with diapause (developmental dormancy) in the speckled wood butterfly (*Pararge aegeria*).

## Study Design

- **Organism**: *Pararge aegeria* (speckled wood butterfly)
- **Reference genome**: GCF_905163445.1 (ilParAegt1.1) with standardised chromosome names
- **Samples**: 94 individuals (44 non-diapause NN controls, 50 diapause DD cases)
- **Sequencing**: Paired-end whole-genome resequencing

## Pipeline Overview

The pipeline runs 10 phases through Snakemake:

| Phase | Step | Tools |
|-------|------|-------|
| 1 | Read QC & trimming | FastQC, Fastp |
| 2 | Alignment | BWA-mem2, Samtools |
| 3 | Post-alignment QC | Picard MarkDuplicates, Mosdepth |
| 4 | Variant calling | Bcftools mpileup/call (parallelised by chromosome) |
| 5 | Variant filtering | Bcftools (QUAL>20, DP>5, biallelic SNPs) |
| 6 | Genotype QC | PLINK2 (MAF, genotyping rate, HWE, missingness) |
| 7 | Population structure | PLINK2 LD pruning + PCA (10 components) |
| 8 | Association testing | PLINK2 logistic regression with Firth fallback |
| 9 | Visualisation | R (Manhattan plots, Q-Q plots, PCA) |
| 10 | Performance reporting | Pipeline metrics |

## Key Analyses

Three GWAS filtering strategies are run in parallel to test robustness:

1. **Standard** &mdash; MAF &ge; 0.05, genotyping rate &ge; 0.9, HWE p > 1e-6
2. **No HWE filter** &mdash; same as standard but without Hardy-Weinberg filtering (important for diapause loci under strong selection)
3. **Relaxed MAF/geno** &mdash; MAF &ge; 0.10, genotyping rate &ge; 0.8

Additional downstream analyses include:

- **QTL-seq** &mdash; allele frequency differentiation between DD and NN groups
- **BayPass GWAS** &mdash; alternative association framework
- **EHH / iHS / XP-EHH** &mdash; extended haplotype homozygosity tests for selection on chromosome 3
- **LD decay plots** &mdash; linkage disequilibrium structure around candidate regions

## Repository Structure

```
├── Snakefile                 # Main pipeline definition
├── config.yaml               # Pipeline configuration
├── samples.tsv               # Sample metadata
├── data/phenotypes.txt        # Diapause phenotype coding
├── envs/                     # Conda environment definitions
│   ├── align.yaml
│   ├── gwas.yaml
│   ├── qc.yaml
│   └── variant.yaml
├── scripts/                  # Analysis scripts
│   ├── plot_gwas_results.R
│   ├── EHH_analysis_chr3.Rmd
│   ├── chr3_detailed_analysis.py
│   ├── chr3_region_analysis.py
│   ├── qtl_seq_analysis.py
│   └── performance_report.py
├── genome/                   # Reference genome metadata & QUAST comparison
├── reports/                  # MultiQC HTML reports (reads, alignment, variants, GWAS)
├── results/                  # Analysis outputs (HTML reports, plots, summary tables)
└── old_html/                 # Archived QC reports
```

> **Note:** Large files (FASTQ, BAM, VCF, PLINK binary genotypes, genome FASTA) are excluded from this repository via `.gitignore`. Only HTML reports, plots, scripts, and small summary files are tracked.

## Running the Pipeline

```bash
# Create conda environments and run all phases
snakemake --use-conda --cores 22 all
```

Requires Snakemake, Conda/Mamba, and the tools specified in `envs/*.yaml`.

## HTML Reports

Key HTML reports included in this repository:

- `GWAS_workflow_v2.html` &mdash; Pipeline overview and workflow documentation
- `EHH_analysis_chr3.html` &mdash; Extended haplotype homozygosity analysis on chromosome 3
- `Need_for_Ref_polarization.html` &mdash; Reference allele polarisation analysis
- `results/gwas_maf0.1geno0.2/summary_report.html` &mdash; GWAS summary with relaxed filters
- `results/BayPass.gwas_maf0.1geno0.2/summary_report.html` &mdash; BayPass GWAS summary
- `results/QTLseq_chrom3/summary_report.html` &mdash; QTL-seq chromosome 3 analysis
- `results/LDplot/ld_report.html` &mdash; LD decay analysis
- `reports/qc*/multiqc_*.html` &mdash; MultiQC quality control reports
