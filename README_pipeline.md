# GWAS Snakemake Pipeline

A complete workflow from FASTQ files to GWAS association results.

## Requirements

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- [Snakemake](https://snakemake.readthedocs.io/) (version 7.0+)

## Installation

```bash
# Install Snakemake (if not already installed)
conda install -c bioconda -c conda-forge snakemake

# Or with mamba (faster)
mamba install -c bioconda -c conda-forge snakemake
```

## Directory Structure

```
Annas_GWAS/
├── Snakefile              # Main workflow
├── config.yaml            # Configuration (edit this!)
├── README_pipeline.md     # This file
├── envs/                  # Conda environments
│   ├── qc.yaml
│   ├── align.yaml
│   ├── variant.yaml
│   ├── gwas.yaml
│   └── plotting.yaml
├── scripts/               # Helper scripts
│   ├── prepare_phenotypes.py
│   └── plot_gwas_results.R
└── data/                  # Your input data (create this)
    ├── fastq/             # FASTQ files
    ├── reference/         # Reference genome
    └── phenotypes.txt     # Phenotype file
```

## Setup

### 1. Prepare Input Data

Create the data directory structure:
```bash
mkdir -p data/fastq data/reference
```

**FASTQ files** should be named:
- `{sample}_R1.fastq.gz` (forward reads)
- `{sample}_R2.fastq.gz` (reverse reads)

**Reference genome**:
```bash
# Place your reference FASTA in data/reference/
cp your_genome.fasta data/reference/genome.fasta

# Index it (pipeline will also do this, but doing it first saves time)
samtools faidx data/reference/genome.fasta
```

**Phenotype file** (`data/phenotypes.txt`):
```
sample_id	phenotype
sample1	0
sample2	1
sample3	0
...
```
- Tab-separated
- Phenotype: 0 = control, 1 = case

### 2. Edit Configuration

Edit `config.yaml` to set:
- Path to FASTQ directory
- Path to reference genome
- Path to phenotype file
- QC thresholds (MAF, missingness, HWE)

## Running the Pipeline

### Dry Run (see what will be executed)
```bash
snakemake -n
```

### Visualize the workflow
```bash
snakemake --dag | dot -Tpdf > workflow_dag.pdf
```

### Run the pipeline
```bash
# Run with 24 cores and conda environments
snakemake --cores 24 --use-conda

# Or specify fewer cores
snakemake --cores 12 --use-conda
```

### Run specific steps only
```bash
# Only run up to alignment
snakemake --cores 24 --use-conda results/aligned/{sample}.dedup.bam

# Only run variant calling
snakemake --cores 24 --use-conda results/variants/filtered_variants.vcf.gz
```

### Resume after failure
```bash
# Snakemake automatically resumes from where it left off
snakemake --cores 24 --use-conda
```

### Generate report
```bash
snakemake --report report.html
```

## Output Files

```
results/
├── qc/
│   ├── fastqc/           # FastQC reports
│   ├── fastp/            # Trimming reports
│   ├── alignment/        # Alignment statistics
│   ├── dedup/            # Duplicate metrics
│   ├── variants/         # Variant statistics
│   └── samples/          # Sample QC metrics
├── trimmed/              # Trimmed FASTQ files
├── aligned/              # BAM files
├── variants/
│   ├── raw_variants.vcf
│   └── filtered_variants.vcf.gz
├── plink/                # PLINK binary files
├── pca/                  # PCA results
├── gwas/
│   └── association_results.PHENO1.glm.logistic
├── plots/
│   ├── manhattan_plot.pdf
│   ├── qq_plot.pdf
│   └── pca_plot.pdf
└── multiqc/
    └── multiqc_report.html
```

## Cleaning Up

```bash
# Remove intermediate files (keeps final results)
snakemake clean_intermediates

# Remove all results (start fresh)
snakemake clean
```

## Troubleshooting

### Conda environment issues on M3 Mac
Some bioinformatics tools may need Rosetta 2:
```bash
# Install Rosetta 2 if not already installed
softwareupdate --install-rosetta

# Run Snakemake with architecture specification if needed
arch -x86_64 snakemake --cores 24 --use-conda
```

### Out of memory
Reduce parallelization:
```bash
snakemake --cores 8 --use-conda
```

### Check what went wrong
```bash
# See the log for a specific rule
cat .snakemake/log/*.snakemake.log
```

## Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `maf_threshold` | 0.05 | Minor allele frequency filter |
| `geno_threshold` | 0.1 | Variant missingness threshold |
| `hwe_threshold` | 1e-6 | Hardy-Weinberg p-value threshold |
| `mind_threshold` | 0.1 | Sample missingness threshold |
