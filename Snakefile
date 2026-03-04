# =============================================================================
# GWAS Pipeline: From FASTQ to Association Analysis
# Snakemake workflow for 96 samples, 350 Mbp genome, 15X coverage
# =============================================================================
#
# Usage:
#   1. Edit config.yaml with your paths and parameters
#   2. Dry run: snakemake -n
#   3. Run: snakemake --cores 24 --use-conda
#   4. Generate report: snakemake --report report.html
#
# =============================================================================

configfile: "config.yaml"

# Get sample names from FASTQ directory
import os
SAMPLES = [f.replace("_R1.fastq.gz", "")
           for f in os.listdir(config["fastq_dir"])
           if f.endswith("_R1.fastq.gz")]

# =============================================================================
# Target rule - defines final outputs
# =============================================================================
rule all:
    input:
        # QC Checkpoint Reports (4 stages)
        "reports/qc1_reads/multiqc_reads.html",
        "reports/qc2_alignment/multiqc_alignment.html",
        "reports/qc3_variants/multiqc_variants.html",
        "reports/qc4_gwas/multiqc_gwas.html",
        # Association results
        "results/gwas/association_results.PHENO1.glm.logistic",
        # Visualizations
        "results/plots/manhattan_plot.pdf",
        "results/plots/qq_plot.pdf"

# =============================================================================
# Phase 1: Read Quality Control
# =============================================================================
rule fastqc_raw:
    input:
        r1 = config["fastq_dir"] + "/{sample}_R1.fastq.gz",
        r2 = config["fastq_dir"] + "/{sample}_R2.fastq.gz"
    output:
        html1 = "results/qc/fastqc/{sample}_R1_fastqc.html",
        html2 = "results/qc/fastqc/{sample}_R2_fastqc.html",
        zip1 = "results/qc/fastqc/{sample}_R1_fastqc.zip",
        zip2 = "results/qc/fastqc/{sample}_R2_fastqc.zip"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o results/qc/fastqc {input.r1} {input.r2}
        """

rule fastp:
    input:
        r1 = config["fastq_dir"] + "/{sample}_R1.fastq.gz",
        r2 = config["fastq_dir"] + "/{sample}_R2.fastq.gz"
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
# Phase 2: Read Alignment
# =============================================================================
rule bwa_index:
    input:
        config["reference"]
    output:
        config["reference"] + ".bwt"
    conda:
        "envs/align.yaml"
    shell:
        """
        bwa index {input}
        """

rule bwa_mem:
    input:
        r1 = "results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.trimmed.fastq.gz",
        ref = config["reference"],
        idx = config["reference"] + ".bwt"
    output:
        temp("results/aligned/{sample}.unsorted.bam")
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    threads: 8
    conda:
        "envs/align.yaml"
    shell:
        """
        bwa mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} | \
        samtools view -bS - > {output}
        """

# =============================================================================
# Phase 3: Post-Alignment Processing
# =============================================================================
rule sort_bam:
    input:
        "results/aligned/{sample}.unsorted.bam"
    output:
        temp("results/aligned/{sample}.sorted.bam")
    threads: 4
    conda:
        "envs/align.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

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
        samtools markdup -@ {threads} --write-index \
            -f {output.metrics} \
            {input} {output.bam}
        """

rule index_bam:
    input:
        "results/aligned/{sample}.dedup.bam"
    output:
        "results/aligned/{sample}.dedup.bam.bai"
    conda:
        "envs/align.yaml"
    shell:
        """
        samtools index {input}
        """

rule alignment_stats:
    input:
        bam = "results/aligned/{sample}.dedup.bam",
        bai = "results/aligned/{sample}.dedup.bam.bai"
    output:
        stats = "results/qc/alignment/{sample}_stats.txt",
        flagstat = "results/qc/alignment/{sample}_flagstat.txt"
    conda:
        "envs/align.yaml"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

# =============================================================================
# Phase 4: Variant Calling (FreeBayes - joint calling)
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

rule freebayes_call:
    input:
        bam_list = "results/variants/bam_list.txt",
        ref = config["reference"],
        bams = expand("results/aligned/{sample}.dedup.bam", sample=SAMPLES),
        bais = expand("results/aligned/{sample}.dedup.bam.bai", sample=SAMPLES)
    output:
        "results/variants/raw_variants.vcf"
    params:
        min_alt_fraction = 0.2,
        min_coverage = 5
    threads: 24
    conda:
        "envs/variant.yaml"
    shell:
        """
        freebayes \
            -f {input.ref} \
            -L {input.bam_list} \
            --min-alternate-fraction {params.min_alt_fraction} \
            --min-coverage {params.min_coverage} \
            --pooled-continuous \
            > {output}
        """

# =============================================================================
# Phase 5: Variant Quality Control
# =============================================================================
rule filter_variants:
    input:
        vcf = "results/variants/raw_variants.vcf",
        ref = config["reference"]
    output:
        "results/variants/filtered_variants.vcf.gz"
    params:
        qual = 20,
        depth = 5
    conda:
        "envs/variant.yaml"
    shell:
        """
        bcftools view -i 'QUAL>{params.qual} && INFO/DP>{params.depth}' {input.vcf} | \
        bcftools norm -f {input.ref} -m -any | \
        bcftools view -m2 -M2 -v snps | \
        bgzip -c > {output}

        tabix -p vcf {output}
        """

# Note: variant_stats rules moved to MultiQC Checkpoints section

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
    conda:
        "envs/gwas.yaml"
    shell:
        """
        plink2 --vcf {input} \
            --make-bed \
            --out {params.prefix} \
            --allow-extra-chr \
            --set-missing-var-ids @:#
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
    conda:
        "envs/gwas.yaml"
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --maf {params.maf} \
            --geno {params.geno} \
            --hwe {params.hwe} \
            --make-bed \
            --out {params.output_prefix} \
            --allow-extra-chr
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
    conda:
        "envs/gwas.yaml"
    shell:
        """
        # Calculate heterozygosity
        plink2 --bfile {params.input_prefix} \
            --het \
            --out results/qc/samples/heterozygosity \
            --allow-extra-chr

        # Calculate missingness
        plink2 --bfile {params.input_prefix} \
            --missing \
            --out results/qc/samples/missing \
            --allow-extra-chr

        # Filter samples by missingness
        plink2 --bfile {params.input_prefix} \
            --mind {params.mind} \
            --make-bed \
            --out {params.output_prefix} \
            --allow-extra-chr
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
    conda:
        "envs/gwas.yaml"
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --indep-pairwise 50 5 0.2 \
            --out {params.output_prefix} \
            --allow-extra-chr
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
    conda:
        "envs/gwas.yaml"
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --extract {input.prune_in} \
            --pca 10 \
            --out {params.output_prefix} \
            --allow-extra-chr
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
    conda:
        "envs/gwas.yaml"
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
        "results/gwas/association_results.PHENO1.glm.logistic"
    params:
        input_prefix = "results/plink/variants_qc_samples",
        output_prefix = "results/gwas/association_results"
    conda:
        "envs/gwas.yaml"
    shell:
        """
        plink2 --bfile {params.input_prefix} \
            --pheno {input.pheno} \
            --covar {input.covar} \
            --covar-col-nums 3-7 \
            --glm firth-fallback hide-covar \
            --out {params.output_prefix} \
            --allow-extra-chr
        """

# =============================================================================
# Phase 9: Visualization
# =============================================================================
rule plot_results:
    input:
        gwas = "results/gwas/association_results.PHENO1.glm.logistic",
        pca_vec = "results/pca/pca.eigenvec",
        pca_val = "results/pca/pca.eigenval"
    output:
        manhattan = "results/plots/manhattan_plot.pdf",
        qq = "results/plots/qq_plot.pdf",
        pca = "results/plots/pca_plot.pdf"
    conda:
        "envs/plotting.yaml"
    script:
        "scripts/plot_gwas_results.R"

# =============================================================================
# MultiQC Checkpoints - QC reports at each stage
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
        expand("results/qc/dedup/{sample}_dedup_metrics.txt", sample=SAMPLES)
    output:
        "reports/qc2_alignment/multiqc_alignment.html"
    params:
        outdir = "reports/qc2_alignment",
        title = "QC2: Alignment Quality Assessment"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc results/qc/alignment results/qc/dedup \
            -o {params.outdir} \
            -n multiqc_alignment.html \
            --title "{params.title}" \
            --force
        """

# QC Checkpoint 3: Variant Quality
rule variant_stats_raw:
    input:
        "results/variants/raw_variants.vcf"
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
        bcftools stats {input} > {output}
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
        het_plot = "results/qc/gwas/het_miss_plot.txt"
    params:
        input_prefix = "results/plink/variants_qc_samples"
    run:
        import subprocess

        # Count variants and samples
        with open(input.bim) as f:
            n_variants = sum(1 for _ in f)
        with open(input.fam) as f:
            n_samples = sum(1 for _ in f)

        # Count samples failing heterozygosity (beyond 3 SD)
        het_vals = []
        with open(input.het) as f:
            next(f)  # skip header
            for line in f:
                parts = line.split()
                if len(parts) >= 6:
                    het_vals.append(float(parts[5]))

        if het_vals:
            import statistics
            mean_f = statistics.mean(het_vals)
            sd_f = statistics.stdev(het_vals) if len(het_vals) > 1 else 0
            fail_het = sum(1 for f in het_vals if f < mean_f - 3*sd_f or f > mean_f + 3*sd_f)
        else:
            fail_het = 0

        # Count samples failing missingness
        fail_miss = 0
        with open(input.missing) as f:
            next(f)  # skip header
            for line in f:
                parts = line.split()
                if len(parts) >= 5 and float(parts[4]) > 0.1:
                    fail_miss += 1

        # Write YAML for MultiQC custom content
        with open(output.summary, 'w') as f:
            f.write(f"""id: 'gwas_summary'
section_name: 'GWAS QC Summary'
description: 'Summary statistics for GWAS-ready data'
plot_type: 'table'
pconfig:
    id: 'gwas_summary_table'
    title: 'GWAS Quality Control Summary'
data:
    GWAS_Ready_Data:
        Total_Variants: {n_variants}
        Total_Samples: {n_samples}
        Samples_Fail_Missingness: {fail_miss}
        Samples_Fail_Heterozygosity: {fail_het}
""")

        # Write placeholder for het plot data
        with open(output.het_plot, 'w') as f:
            f.write("Het vs Missingness data prepared\n")

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

# =============================================================================
# Utility rules
# =============================================================================
rule clean:
    shell:
        """
        rm -rf results/
        """

rule clean_intermediates:
    shell:
        """
        rm -f results/aligned/*.unsorted.bam
        rm -f results/aligned/*.sorted.bam
        rm -f results/trimmed/*.fastq.gz
        """
