# How to Update GitHub for This Project

**Remote:** https://github.com/bioinfowheat/AnnaGWAS.git
**Branch:** `main`
**Local path:** `/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS`

## Quick Steps

```bash
cd /Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS
git add <specific files>
git commit -m "descriptive message"
git push origin main
```

## What Gets Included vs Excluded

The `.gitignore` is carefully configured. The key principle is: **scripts, configs, HTML reports, plots, and small data files are tracked; large binary data files are excluded.**

### INCLUDED (track these)
- **Scripts:** `scripts/*.R`, `scripts/*.py`, `scripts/*.Rmd`, `*.sh`
- **Pipeline files:** `Snakefile`, `config.yaml`, `samples.tsv`, `envs/*.yaml`
- **Phenotype data:** `data/phenotypes.txt`, `data/Annas_phenotypes.tsv`, `data/*.xlsx`
- **HTML reports:** All `.html` files everywhere (GWAS summaries, QTL-seq, EHH, LD, MultiQC reports)
- **Plots:** `.pdf`, `.png`, `.svg` plot outputs in `results/` and `reports/`
- **Small result files:** PCA eigenvectors/eigenvalues, QTL-seq windows, log files, sample lists
- **Genome metadata:** `genome/genome_source.md`, `genome/quast_comparison/`, BWA `.amb`/`.ann` files

### EXCLUDED by .gitignore (never commit these)
- **Genome FASTA:** `genome/*.fna`, `*.fa`, `*.fasta`, `*.gff`, `*.gtf`, `*.dict`, `*.fai`
- **BWA indexes:** `*.0123`, `*.bwt.2bit.64`, `*.pac`
- **Sequencing data:** `*.fastq.gz`, `*.fq.gz`
- **Alignments:** `*.bam`, `*.bam.bai`, `*.cram`, `*.crai`
- **Variants:** `*.vcf`, `*.vcf.gz`, `*.vcf.gz.tbi`, `*.bcf`, `*.csi`
- **PLINK binary:** `*.bed`, `*.bim`, `*.fam`
- **Association results:** `*.glm.logistic.hybrid`, `*.glm.logistic`, `*.glm.linear`
- **PLINK intermediate:** `*.prune.out`, `*.prune.in`, `*.afreq`, `*.vmiss`, `*.nosex`
- **LD data:** `*.vcor2`
- **Large tables:** `results/BayPass.gwas_maf0.1geno0.2/qtlseq/snp_index.tsv`, `results/qtl_seq/qtl_seq_per_snp.tsv.gz`
- **Third-party tools:** `results/SparsePainter/`, `results/*/plink2`
- **Intermediate outputs:** `results/variants/per_region/`, `results/qc/fastqc/`
- **Package caches:** `results/*/micromamba/`, `.snakemake/`
- **MultiQC raw data:** `reports/*/multiqc_*_data/`
- **Other:** `*.parquet`, `.DS_Store`, `__pycache__/`, `*.pyc`

## Before Committing: Checklist

1. **Check what's new:** `git status` — review untracked and modified files
2. **Check for oversized files:** No single file should exceed 100 MB (GitHub hard limit); files >50 MB get warnings
3. **Stage specific files** rather than `git add .` if you want to be selective
4. **If adding a new file type** that's large (new analysis output format), add a pattern to `.gitignore` first
5. **Total repo size** should stay under ~1 GB ideally

## If .gitignore Needs Updating

If a new analysis produces large output files in a new format, add the pattern to `.gitignore` before committing:

```bash
# Example: adding a new large format
echo "*.newformat" >> .gitignore
git add .gitignore
git commit -m "Exclude *.newformat files from tracking"
```

## Common Scenarios

### Pushing new analysis results
```bash
git add results/new_analysis/summary_report.html results/new_analysis/*.png scripts/new_analysis.py
git commit -m "Add new analysis results and script"
git push origin main
```

### Updating phenotype data
```bash
git add data/Annas_phenotypes.tsv data/phenotype_summary.html
git commit -m "Update phenotype file with new traits"
git push origin main
```

### Pushing everything new at once
```bash
# First verify nothing huge slipped through
git add --dry-run . 2>&1 | sed "s/^add '//;s/'$//" | while read f; do
  [ -f "$f" ] && size=$(stat -f "%z" "$f"); [ "$size" -gt 52428800 ] && echo "WARNING >50MB: $f ($((size/1048576))MB)"
done

# If no warnings, stage and push
git add .
git commit -m "Add latest analysis outputs"
git push origin main
```
