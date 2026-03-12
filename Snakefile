# =============================================================================
# GWAS Pipeline v2: FASTQ → Variant Calling → GWAS Association
# Project: YF-4363 / XI-4103 (snpseq01427) — Pararge aegeria
# NN (non-diapause) vs DD (diapause) — 10 samples (5 NN + 5 DD)
# =============================================================================
#
# Usage:
#   snakemake -n                                              # dry run
#   snakemake --cores 22 --use-conda --conda-frontend micromamba  # execute
#   snakemake --report report.html                            # report
#
# =============================================================================

configfile: "config.yaml"

import pandas as pd

# Read sample table
samples_df = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples_df["sample"].tolist()

# Build sample→file lookup
SAMPLE_R1 = dict(zip(samples_df["sample"], samples_df["r1"]))
SAMPLE_R2 = dict(zip(samples_df["sample"], samples_df["r2"]))

# Reference genome path
REF = config["reference"]

# Get chromosome/scaffold names from .fai index for per-region variant calling
# The .fai is created by samtools faidx in the index rule
import os
FAI_FILE = REF + ".fai"
if os.path.exists(FAI_FILE):
    with open(FAI_FILE) as f:
        REGIONS = [line.split("\t")[0] for line in f]
else:
    # Will be populated after indexing; for dry-run use a placeholder
    REGIONS = ["placeholder"]

# =============================================================================
# Target rule
# =============================================================================
rule all:
    input:
        # QC Checkpoint 1: Read quality
        "reports/qc1_reads/multiqc_reads.html",
        # QC Checkpoint 2: Alignment quality
        "reports/qc2_alignment/multiqc_alignment.html",
        # QC Checkpoint 3: Variant quality
        "reports/qc3_variants/multiqc_variants.html",
        # QC Checkpoint 4: GWAS-ready data
        "reports/qc4_gwas/multiqc_gwas.html",
        # Filtered variants
        "results/variants/filtered_variants.vcf.gz",
        # GWAS results
        "results/gwas/association_results.PHENO1.glm.logistic.hybrid",
        # Plots
        "results/plots/manhattan_plot.pdf",
        "results/plots/qq_plot.pdf",
        "results/plots/pca_plot.pdf",
        # Performance report
        "reports/performance_report.html",
        # GWAS without HWE filter
        "results/gwas_noHWEfilt/association_results.PHENO1.glm.logistic.hybrid",
        "results/plots_noHWEfilt/manhattan_plot.pdf",
        "results/plots_noHWEfilt/qq_plot.pdf",
        # GWAS with MAF 0.10, geno 0.2, no HWE
        "results/gwas_maf0.1geno0.2/association_results.PHENO1.glm.logistic.hybrid",
        "results/plots_maf0.1geno0.2/manhattan_plot.pdf",
        "results/plots_maf0.1geno0.2/qq_plot.pdf",

# =============================================================================
# Phase 1: Read Quality Control
# =============================================================================
rule fastqc_raw:
    input:
        r1 = lambda wc: SAMPLE_R1[wc.sample],
        r2 = lambda wc: SAMPLE_R2[wc.sample]
    output:
        zip1 = "results/qc/fastqc/{sample}_R1_fastqc.zip",
        zip2 = "results/qc/fastqc/{sample}_R2_fastqc.zip",
        html1 = "results/qc/fastqc/{sample}_R1_fastqc.html",
        html2 = "results/qc/fastqc/{sample}_R2_fastqc.html"
    params:
        outdir = "results/qc/fastqc"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        # FastQC names outputs based on input filename, so we need symlinks
        # with our desired sample names
        ln -sf $(readlink -f {input.r1}) {params.outdir}/{wildcards.sample}_R1.fastq.gz
        ln -sf $(readlink -f {input.r2}) {params.outdir}/{wildcards.sample}_R2.fastq.gz
        fastqc -t {threads} -o {params.outdir} \
            {params.outdir}/{wildcards.sample}_R1.fastq.gz \
            {params.outdir}/{wildcards.sample}_R2.fastq.gz
        rm -f {params.outdir}/{wildcards.sample}_R1.fastq.gz \
              {params.outdir}/{wildcards.sample}_R2.fastq.gz
        """

rule fastp:
    input:
        r1 = lambda wc: SAMPLE_R1[wc.sample],
        r2 = lambda wc: SAMPLE_R2[wc.sample]
    output:
        r1 = "results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.trimmed.fastq.gz",
        html = "results/qc/fastp/{sample}_fastp.html",
        json = "results/qc/fastp/{sample}_fastp.json"
    threads: 4
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --html {output.html} --json {output.json} \
            --thread {threads} \
            --qualified_quality_phred 20 \
            --length_required 50 \
            --detect_adapter_for_pe
        """

# =============================================================================
# Phase 2: Alignment
# =============================================================================
rule bwa_mem2_index:
    input:
        REF
    output:
        REF + ".0123",
        REF + ".amb",
        REF + ".ann",
        REF + ".bwt.2bit.64",
        REF + ".pac"
    conda:
        "envs/align.yaml"
    shell:
        """
        bwa-mem2 index {input}
        """

rule samtools_faidx:
    input:
        REF
    output:
        REF + ".fai"
    conda:
        "envs/align.yaml"
    shell:
        """
        samtools faidx {input}
        """

rule bwa_mem2_align:
    input:
        r1 = "results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.trimmed.fastq.gz",
        ref = REF,
        idx = REF + ".0123",
        fai = REF + ".fai"
    output:
        temp("results/aligned/{sample}.sorted.bam")
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\tLB:{sample}"
    threads: 7
    conda:
        "envs/align.yaml"
    shell:
        """
        bwa-mem2 mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} | \
        samtools sort -@ 2 -m 10G -o {output} -
        """

# =============================================================================
# Phase 3: Post-Alignment Processing
# =============================================================================

rule mark_duplicates:
    input:
        "results/aligned/{sample}.sorted.bam"
    output:
        bam = "results/aligned/{sample}.dedup.bam",
        metrics = "results/qc/dedup/{sample}_dedup_metrics.txt"
    threads: 4
    conda:
        "envs/align.yaml"
    shell:
        """
        picard MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            --REMOVE_DUPLICATES false \
            --VALIDATION_STRINGENCY LENIENT
        """

rule index_bam:
    input:
        "results/aligned/{sample}.dedup.bam"
    output:
        "results/aligned/{sample}.dedup.bam.bai"
    threads: 2
    conda:
        "envs/align.yaml"
    shell:
        """
        samtools index -@ {threads} {input}
        """

rule alignment_stats:
    input:
        bam = "results/aligned/{sample}.dedup.bam",
        bai = "results/aligned/{sample}.dedup.bam.bai",
        ref = REF
    output:
        stats = "results/qc/alignment/{sample}_stats.txt",
        flagstat = "results/qc/alignment/{sample}_flagstat.txt",
        idxstats = "results/qc/alignment/{sample}_idxstats.txt"
    conda:
        "envs/align.yaml"
    shell:
        """
        samtools stats --ref-seq {input.ref} {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        samtools idxstats {input.bam} > {output.idxstats}
        """

rule mosdepth:
    input:
        bam = "results/aligned/{sample}.dedup.bam",
        bai = "results/aligned/{sample}.dedup.bam.bai"
    output:
        summary = "results/qc/mosdepth/{sample}.mosdepth.summary.txt",
        dist = "results/qc/mosdepth/{sample}.mosdepth.global.dist.txt"
    params:
        prefix = "results/qc/mosdepth/{sample}"
    threads: 4
    conda:
        "envs/align.yaml"
    shell:
        """
        mosdepth --threads {threads} --no-per-base {params.prefix} {input.bam}
        """

# =============================================================================
# Phase 4: Variant Calling (bcftools mpileup/call — per region)
# =============================================================================
rule create_bam_list:
    input:
        bams = expand("results/aligned/{sample}.dedup.bam", sample=SAMPLES),
        bais = expand("results/aligned/{sample}.dedup.bam.bai", sample=SAMPLES)
    output:
        "results/variants/bam_list.txt"
    run:
        with open(output[0], 'w') as f:
            for bam in input.bams:
                f.write(bam + '\n')

rule bcftools_call:
    input:
        bam_list = "results/variants/bam_list.txt",
        ref = REF,
        fai = REF + ".fai",
        bams = expand("results/aligned/{sample}.dedup.bam", sample=SAMPLES),
        bais = expand("results/aligned/{sample}.dedup.bam.bai", sample=SAMPLES)
    output:
        temp("results/variants/per_region/{region}.vcf.gz")
    conda:
        "envs/variant.yaml"
    shell:
        """
        bcftools mpileup \
            --fasta-ref {input.ref} \
            --bam-list {input.bam_list} \
            --regions {wildcards.region} \
            --min-MQ 20 \
            --min-BQ 20 \
            --annotate FORMAT/AD,FORMAT/DP \
            --max-depth 500 \
        | bcftools call \
            --multiallelic-caller \
            --variants-only \
            -Oz -o {output}

        bcftools index {output}
        """

rule concat_vcfs:
    input:
        vcfs = expand("results/variants/per_region/{region}.vcf.gz", region=REGIONS),
    output:
        "results/variants/raw_variants.vcf.gz"
    conda:
        "envs/variant.yaml"
    shell:
        """
        bcftools concat {input.vcfs} -a -Oz -o {output}
        bcftools index {output}
        """

# =============================================================================
# Phase 5: Variant Filtering
# =============================================================================
rule filter_variants:
    input:
        vcf = "results/variants/raw_variants.vcf.gz",
        ref = REF
    output:
        "results/variants/filtered_variants.vcf.gz"
    conda:
        "envs/variant.yaml"
    shell:
        """
        bcftools view -i 'QUAL>20 && INFO/DP>5' {input.vcf} | \
        bcftools norm -f {input.ref} -m -any | \
        bcftools view -m2 -M2 -v snps -Oz -o {output}

        bcftools index {output}
        """

# =============================================================================
# Phase 6: Convert to PLINK format and Sample QC
# =============================================================================
rule vcf_to_plink:
    input:
        "results/variants/filtered_variants.vcf.gz"
    output:
        bed = "results/plink/variants.bed",
        bim = "results/plink/variants.bim",
        fam = "results/plink/variants.fam"
    params:
        prefix = "results/plink/variants"
    # plink2 installed natively at ~/bin/plink2 (ARM64, no conda)
    shell:
        """
        plink2 --vcf {input} \
            --make-bed \
            --out {params.prefix} \
            --chr-set 27 --allow-extra-chr \
            --set-missing-var-ids '@:#:\$r:\$a'
        """

rule variant_qc:
    input:
        bed = "results/plink/variants.bed",
        bim = "results/plink/variants.bim",
        fam = "results/plink/variants.fam"
    output:
        bed = "results/plink/variants_qc.bed",
        bim = "results/plink/variants_qc.bim",
        fam = "results/plink/variants_qc.fam"
    params:
        input_prefix = "results/plink/variants",
        output_prefix = "results/plink/variants_qc",
        maf = config["maf_threshold"],
        geno = config["geno_threshold"],
        hwe = config["hwe_threshold"]
    # plink2 installed natively at ~/bin/plink2 (ARM64, no conda)
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --maf {params.maf} \
            --geno {params.geno} \
            --hwe {params.hwe} \
            --make-bed \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

rule sample_qc:
    input:
        bed = "results/plink/variants_qc.bed",
        bim = "results/plink/variants_qc.bim",
        fam = "results/plink/variants_qc.fam"
    output:
        bed = "results/plink/variants_qc_samples.bed",
        bim = "results/plink/variants_qc_samples.bim",
        fam = "results/plink/variants_qc_samples.fam",
        het = "results/qc/samples/heterozygosity.het",
        missing = "results/qc/samples/missing.smiss"
    params:
        input_prefix = "results/plink/variants_qc",
        output_prefix = "results/plink/variants_qc_samples",
        mind = config["mind_threshold"]
    # plink2 installed natively at ~/bin/plink2 (ARM64, no conda)
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --freq \
            --out results/qc/samples/allele_freq \
            --chr-set 27 --allow-extra-chr

        plink2 --bfile {params.input_prefix} \
            --het \
            --out results/qc/samples/heterozygosity \
            --chr-set 27 --allow-extra-chr

        plink2 --bfile {params.input_prefix} \
            --missing \
            --out results/qc/samples/missing \
            --chr-set 27 --allow-extra-chr

        plink2 --bfile {params.input_prefix} \
            --mind {params.mind} \
            --make-bed \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

# =============================================================================
# Phase 7: Population Structure (PCA)
# =============================================================================
rule ld_prune:
    input:
        bed = "results/plink/variants_qc_samples.bed",
        bim = "results/plink/variants_qc_samples.bim",
        fam = "results/plink/variants_qc_samples.fam"
    output:
        prune_in = "results/plink/pruned.prune.in",
        prune_out = "results/plink/pruned.prune.out"
    params:
        input_prefix = "results/plink/variants_qc_samples",
        output_prefix = "results/plink/pruned"
    # plink2 installed natively at ~/bin/plink2 (ARM64, no conda)
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --indep-pairwise 50 5 0.2 \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

rule pca:
    input:
        bed = "results/plink/variants_qc_samples.bed",
        bim = "results/plink/variants_qc_samples.bim",
        fam = "results/plink/variants_qc_samples.fam",
        prune_in = "results/plink/pruned.prune.in"
    output:
        eigenvec = "results/pca/pca.eigenvec",
        eigenval = "results/pca/pca.eigenval"
    params:
        input_prefix = "results/plink/variants_qc_samples",
        output_prefix = "results/pca/pca"
    # plink2 installed natively at ~/bin/plink2 (ARM64, no conda)
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --extract {input.prune_in} \
            --pca 10 \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

# =============================================================================
# Phase 8: Association Analysis
# =============================================================================
rule prepare_phenotypes:
    input:
        fam = "results/plink/variants_qc_samples.fam",
        pheno = config["phenotype_file"]
    output:
        "results/gwas/phenotypes.txt"
    # plink2 installed natively at ~/bin/plink2 (ARM64, no conda)
    script:
        "scripts/prepare_phenotypes.py"

rule gwas_association:
    input:
        bed = "results/plink/variants_qc_samples.bed",
        bim = "results/plink/variants_qc_samples.bim",
        fam = "results/plink/variants_qc_samples.fam",
        pheno = "results/gwas/phenotypes.txt",
        covar = "results/pca/pca.eigenvec"
    output:
        "results/gwas/association_results.PHENO1.glm.logistic.hybrid"
    params:
        input_prefix = "results/plink/variants_qc_samples",
        output_prefix = "results/gwas/association_results"
    # plink2 installed natively at ~/bin/plink2 (ARM64, no conda)
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --pheno {input.pheno} \
            --covar {input.covar} \
            --covar-col-nums 3-12 \
            --glm firth-fallback hide-covar \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

# =============================================================================
# Phase 8b: Association Analysis without HWE filter
# The focal locus (diapause vs non-diapause) is expected to be homozygous
# ref vs. alt without heterozygotes, violating HWE. This branch keeps those
# variants in the analysis. PCA covariates from the HWE-filtered set are
# reused (standard for population structure estimation).
# =============================================================================
rule variant_qc_noHWE:
    input:
        bed = "results/plink/variants.bed",
        bim = "results/plink/variants.bim",
        fam = "results/plink/variants.fam"
    output:
        bed = "results/plink/variants_qc_noHWE.bed",
        bim = "results/plink/variants_qc_noHWE.bim",
        fam = "results/plink/variants_qc_noHWE.fam"
    params:
        input_prefix = "results/plink/variants",
        output_prefix = "results/plink/variants_qc_noHWE",
        maf = config["maf_threshold"],
        geno = config["geno_threshold"]
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --maf {params.maf} \
            --geno {params.geno} \
            --make-bed \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

rule sample_qc_noHWE:
    input:
        bed = "results/plink/variants_qc_noHWE.bed",
        bim = "results/plink/variants_qc_noHWE.bim",
        fam = "results/plink/variants_qc_noHWE.fam"
    output:
        bed = "results/plink/variants_qc_samples_noHWE.bed",
        bim = "results/plink/variants_qc_samples_noHWE.bim",
        fam = "results/plink/variants_qc_samples_noHWE.fam"
    params:
        input_prefix = "results/plink/variants_qc_noHWE",
        output_prefix = "results/plink/variants_qc_samples_noHWE",
        mind = config["mind_threshold"]
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --mind {params.mind} \
            --make-bed \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

rule prepare_phenotypes_noHWE:
    input:
        fam = "results/plink/variants_qc_samples_noHWE.fam",
        pheno = config["phenotype_file"]
    output:
        "results/gwas_noHWEfilt/phenotypes.txt"
    script:
        "scripts/prepare_phenotypes.py"

rule gwas_association_noHWE:
    input:
        bed = "results/plink/variants_qc_samples_noHWE.bed",
        bim = "results/plink/variants_qc_samples_noHWE.bim",
        fam = "results/plink/variants_qc_samples_noHWE.fam",
        pheno = "results/gwas_noHWEfilt/phenotypes.txt",
        covar = "results/pca/pca.eigenvec"
    output:
        "results/gwas_noHWEfilt/association_results.PHENO1.glm.logistic.hybrid"
    params:
        input_prefix = "results/plink/variants_qc_samples_noHWE",
        output_prefix = "results/gwas_noHWEfilt/association_results"
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --pheno {input.pheno} \
            --covar {input.covar} \
            --covar-col-nums 3-12 \
            --glm firth-fallback hide-covar \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

rule plot_results_noHWE:
    input:
        gwas = "results/gwas_noHWEfilt/association_results.PHENO1.glm.logistic.hybrid",
        pca_vec = "results/pca/pca.eigenvec",
        pca_val = "results/pca/pca.eigenval"
    output:
        manhattan = "results/plots_noHWEfilt/manhattan_plot.pdf",
        qq = "results/plots_noHWEfilt/qq_plot.pdf",
        pca = "results/plots_noHWEfilt/pca_plot.pdf"
    conda:
        "envs/plotting.yaml"
    shell:
        """
        Rscript scripts/plot_gwas_results.R \
            {input.gwas} {input.pca_vec} {input.pca_val} \
            {output.manhattan} {output.qq} {output.pca}
        """

# =============================================================================
# Phase 8c: GWAS with MAF 0.10, geno 0.2, no HWE filter
# =============================================================================
rule variant_qc_maf01geno02:
    input:
        bed = "results/plink/variants.bed",
        bim = "results/plink/variants.bim",
        fam = "results/plink/variants.fam"
    output:
        bed = "results/plink/variants_qc_maf01geno02.bed",
        bim = "results/plink/variants_qc_maf01geno02.bim",
        fam = "results/plink/variants_qc_maf01geno02.fam"
    params:
        input_prefix = "results/plink/variants",
        output_prefix = "results/plink/variants_qc_maf01geno02"
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --maf 0.10 \
            --geno 0.2 \
            --make-bed \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

rule sample_qc_maf01geno02:
    input:
        bed = "results/plink/variants_qc_maf01geno02.bed",
        bim = "results/plink/variants_qc_maf01geno02.bim",
        fam = "results/plink/variants_qc_maf01geno02.fam"
    output:
        bed = "results/plink/variants_qc_samples_maf01geno02.bed",
        bim = "results/plink/variants_qc_samples_maf01geno02.bim",
        fam = "results/plink/variants_qc_samples_maf01geno02.fam"
    params:
        input_prefix = "results/plink/variants_qc_maf01geno02",
        output_prefix = "results/plink/variants_qc_samples_maf01geno02",
        mind = config["mind_threshold"]
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --mind {params.mind} \
            --make-bed \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

rule prepare_phenotypes_maf01geno02:
    input:
        fam = "results/plink/variants_qc_samples_maf01geno02.fam",
        pheno = config["phenotype_file"]
    output:
        "results/gwas_maf0.1geno0.2/phenotypes.txt"
    script:
        "scripts/prepare_phenotypes.py"

rule gwas_association_maf01geno02:
    input:
        bed = "results/plink/variants_qc_samples_maf01geno02.bed",
        bim = "results/plink/variants_qc_samples_maf01geno02.bim",
        fam = "results/plink/variants_qc_samples_maf01geno02.fam",
        pheno = "results/gwas_maf0.1geno0.2/phenotypes.txt",
        covar = "results/pca/pca.eigenvec"
    output:
        "results/gwas_maf0.1geno0.2/association_results.PHENO1.glm.logistic.hybrid"
    params:
        input_prefix = "results/plink/variants_qc_samples_maf01geno02",
        output_prefix = "results/gwas_maf0.1geno0.2/association_results"
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --pheno {input.pheno} \
            --covar {input.covar} \
            --covar-col-nums 3-12 \
            --glm firth-fallback hide-covar \
            --out {params.output_prefix} \
            --chr-set 27 --allow-extra-chr
        """

rule plot_results_maf01geno02:
    input:
        gwas = "results/gwas_maf0.1geno0.2/association_results.PHENO1.glm.logistic.hybrid",
        pca_vec = "results/pca/pca.eigenvec",
        pca_val = "results/pca/pca.eigenval"
    output:
        manhattan = "results/plots_maf0.1geno0.2/manhattan_plot.pdf",
        qq = "results/plots_maf0.1geno0.2/qq_plot.pdf",
        pca = "results/plots_maf0.1geno0.2/pca_plot.pdf"
    conda:
        "envs/plotting.yaml"
    shell:
        """
        Rscript scripts/plot_gwas_results.R \
            {input.gwas} {input.pca_vec} {input.pca_val} \
            {output.manhattan} {output.qq} {output.pca}
        """

# =============================================================================
# Phase 9: Visualization
# =============================================================================
rule plot_results:
    input:
        gwas = "results/gwas/association_results.PHENO1.glm.logistic.hybrid",
        pca_vec = "results/pca/pca.eigenvec",
        pca_val = "results/pca/pca.eigenval"
    output:
        manhattan = "results/plots/manhattan_plot.pdf",
        qq = "results/plots/qq_plot.pdf",
        pca = "results/plots/pca_plot.pdf"
    conda:
        "envs/plotting.yaml"
    shell:
        """
        Rscript scripts/plot_gwas_results.R \
            {input.gwas} {input.pca_vec} {input.pca_val} \
            {output.manhattan} {output.qq} {output.pca}
        """

# =============================================================================
# Phase 10: Performance Report
# =============================================================================
rule performance_report:
    input:
        # Depend on final BAM and VCF outputs to ensure this runs after them
        expand("results/aligned/{sample}.dedup.bam", sample=SAMPLES),
        "results/variants/filtered_variants.vcf.gz"
    output:
        "reports/performance_report.html"
    params:
        log = ".snakemake/log/2026-03-05T100003.181364.snakemake.log",
        results = "results"
    run:
        # Find the main pipeline log (the largest one, which is the full run)
        import glob as gl
        logs = gl.glob(".snakemake/log/*.snakemake.log")
        # Use the specified log if it exists, otherwise find the largest
        main_log = params.log
        if not os.path.exists(main_log):
            main_log = max(logs, key=os.path.getsize) if logs else params.log
        shell("python scripts/performance_report.py {main_log} {params.results} {output}")

# =============================================================================
# MultiQC Checkpoints
# =============================================================================

# QC Checkpoint 1: Read Quality
rule multiqc_reads:
    input:
        expand("results/qc/fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("results/qc/fastqc/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand("results/qc/fastp/{sample}_fastp.json", sample=SAMPLES)
    output:
        "reports/qc1_reads/multiqc_reads.html"
    params:
        outdir = "reports/qc1_reads",
        title = "QC1: Read Quality Assessment"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc results/qc/fastqc results/qc/fastp \
            -o {params.outdir} \
            -n multiqc_reads.html \
            --title "{params.title}" \
            --force
        """

# QC Checkpoint 2: Alignment Quality
rule multiqc_alignment:
    input:
        expand("results/qc/alignment/{sample}_stats.txt", sample=SAMPLES),
        expand("results/qc/alignment/{sample}_flagstat.txt", sample=SAMPLES),
        expand("results/qc/dedup/{sample}_dedup_metrics.txt", sample=SAMPLES),
        expand("results/qc/mosdepth/{sample}.mosdepth.summary.txt", sample=SAMPLES),
        expand("results/qc/mosdepth/{sample}.mosdepth.global.dist.txt", sample=SAMPLES)
    output:
        "reports/qc2_alignment/multiqc_alignment.html"
    params:
        outdir = "reports/qc2_alignment",
        title = "QC2: Alignment Quality Assessment"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc results/qc/alignment results/qc/dedup results/qc/mosdepth \
            -o {params.outdir} \
            -n multiqc_alignment.html \
            --title "{params.title}" \
            --force
        """

# QC Checkpoint 3: Variant Quality
rule variant_stats_raw:
    input:
        "results/variants/raw_variants.vcf.gz"
    output:
        "results/qc/variants/raw_variants.bcftools_stats.txt"
    conda:
        "envs/variant.yaml"
    shell:
        """
        bcftools stats {input} > {output}
        """

rule variant_stats_filtered:
    input:
        "results/variants/filtered_variants.vcf.gz"
    output:
        "results/qc/variants/filtered_variants.bcftools_stats.txt"
    conda:
        "envs/variant.yaml"
    shell:
        """
        bcftools stats -s - {input} > {output}
        """

rule multiqc_variants:
    input:
        "results/qc/variants/raw_variants.bcftools_stats.txt",
        "results/qc/variants/filtered_variants.bcftools_stats.txt"
    output:
        "reports/qc3_variants/multiqc_variants.html"
    params:
        outdir = "reports/qc3_variants",
        title = "QC3: Variant Quality Assessment"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc results/qc/variants \
            -o {params.outdir} \
            -n multiqc_variants.html \
            --title "{params.title}" \
            --force
        """

# QC Checkpoint 4: GWAS-Ready Data
rule gwas_qc_stats:
    input:
        bed = "results/plink/variants_qc_samples.bed",
        bim = "results/plink/variants_qc_samples.bim",
        fam = "results/plink/variants_qc_samples.fam",
        het = "results/qc/samples/heterozygosity.het",
        missing = "results/qc/samples/missing.smiss"
    output:
        summary = "results/qc/gwas/gwas_summary_mqc.yaml",
    run:
        import statistics

        with open(input.bim) as f:
            n_variants = sum(1 for _ in f)
        with open(input.fam) as f:
            n_samples = sum(1 for _ in f)

        het_vals = []
        with open(input.het) as f:
            next(f)
            for line in f:
                parts = line.split()
                if len(parts) >= 6:
                    het_vals.append(float(parts[5]))

        if het_vals and len(het_vals) > 1:
            mean_f = statistics.mean(het_vals)
            sd_f = statistics.stdev(het_vals)
            fail_het = sum(1 for v in het_vals if v < mean_f - 3*sd_f or v > mean_f + 3*sd_f)
        else:
            fail_het = 0

        fail_miss = 0
        with open(input.missing) as f:
            next(f)
            for line in f:
                parts = line.split()
                if len(parts) >= 5 and float(parts[4]) > 0.1:
                    fail_miss += 1

        with open(output.summary, 'w') as f:
            lines = [
                "id: 'gwas_summary'",
                "section_name: 'GWAS QC Summary'",
                "description: 'Summary statistics for GWAS-ready data'",
                "plot_type: 'table'",
                "pconfig:",
                "  id: 'gwas_summary_table'",
                "  title: 'GWAS Quality Control Summary'",
                "data:",
                "  GWAS_Ready_Data:",
                f"    Total_Variants: {n_variants}",
                f"    Total_Samples: {n_samples}",
                f"    Samples_Fail_Missingness: {fail_miss}",
                f"    Samples_Fail_Heterozygosity: {fail_het}",
            ]
            f.write('\n'.join(lines) + '\n')

rule multiqc_gwas:
    input:
        "results/qc/gwas/gwas_summary_mqc.yaml",
        "results/qc/samples/heterozygosity.het",
        "results/qc/samples/missing.smiss"
    output:
        "reports/qc4_gwas/multiqc_gwas.html"
    params:
        outdir = "reports/qc4_gwas",
        title = "QC4: GWAS-Ready Data Assessment"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc results/qc/gwas results/qc/samples \
            -o {params.outdir} \
            -n multiqc_gwas.html \
            --title "{params.title}" \
            --force
        """
