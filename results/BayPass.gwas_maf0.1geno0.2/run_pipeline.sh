#!/bin/bash
set -euo pipefail

# =============================================================================
# FreeBayes GWAS Pipeline
# Calls SNPs+indels with FreeBayes, runs GWAS + QTLseq, focuses on chr3
# =============================================================================

PROJ="/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/BayPass.gwas_maf0.1geno0.2"
BASE="/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS"
FB="${PROJ}/micromamba/envs/freebayes/bin"
REF="${BASE}/genome/GCF_905163445.1_ilParAegt1.1_genomic.edited.fna"
BAMLIST="${PROJ}/bam_list.txt"
PLINK2="${PROJ}/plink2"

FREEBAYES="${FB}/freebayes"
BCFTOOLS="${FB}/bcftools"
SAMTOOLS="${FB}/samtools"
BGZIP="${FB}/bgzip"
TABIX="${FB}/tabix"

# Chromosomes to process (autosomes 1-27 + Z + W)
CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 Z W"

# Max parallel FreeBayes jobs (conservative for memory)
MAX_JOBS=3

echo "============================================"
echo "FreeBayes GWAS Pipeline"
echo "Started: $(date)"
echo "Samples: $(wc -l < $BAMLIST)"
echo "Reference: $REF"
echo "============================================"

# =========================================================================
# STEP 1: Run FreeBayes per chromosome
# =========================================================================
echo ""
echo "=== STEP 1: FreeBayes variant calling ==="

run_freebayes_chr() {
    local chr=$1
    local outvcf="${PROJ}/vcf/per_chr/${chr}.vcf.gz"

    if [ -f "$outvcf" ] && [ -f "${outvcf}.tbi" ]; then
        echo "  [SKIP] chr${chr} — already done"
        return 0
    fi

    echo "  [START] chr${chr} — $(date '+%H:%M:%S')"

    # FreeBayes with haplotype-aware calling
    # --min-mapping-quality 20: skip poorly mapped reads
    # --min-base-quality 20: skip low quality bases
    # --min-alternate-count 2: require at least 2 reads supporting alt
    # --min-alternate-fraction 0.05: minimum alt allele fraction
    # --use-best-n-alleles 4: limit alleles for speed
    # --report-monomorphic off (default): skip invariant sites
    ${FREEBAYES} \
        --fasta-reference "${REF}" \
        --bam-list "${BAMLIST}" \
        --region "${chr}" \
        --min-mapping-quality 20 \
        --min-base-quality 20 \
        --min-alternate-count 2 \
        --min-alternate-fraction 0.05 \
        --use-best-n-alleles 4 \
        2>/dev/null \
    | ${BCFTOOLS} view -i 'QUAL>20' \
    | ${BGZIP} -c > "${outvcf}" \
    && ${TABIX} -p vcf "${outvcf}"

    local n_var=$(${BCFTOOLS} view -H "${outvcf}" 2>/dev/null | wc -l | tr -d ' ')
    echo "  [DONE] chr${chr} — ${n_var} variants — $(date '+%H:%M:%S')"
}

# Run chromosome 3 FIRST (priority), then others
echo "  Priority: chromosome 3"
run_freebayes_chr 3

echo "  Running remaining chromosomes (max ${MAX_JOBS} parallel)..."
job_count=0
for chr in $CHROMOSOMES; do
    if [ "$chr" = "3" ]; then continue; fi

    run_freebayes_chr "$chr" &
    job_count=$((job_count + 1))

    if [ $job_count -ge $MAX_JOBS ]; then
        wait -n 2>/dev/null || wait
        job_count=$((job_count - 1))
    fi
done
wait
echo "  All chromosomes complete."

# =========================================================================
# STEP 2: Merge per-chromosome VCFs
# =========================================================================
echo ""
echo "=== STEP 2: Merging VCFs ==="

MERGED="${PROJ}/vcf/freebayes_all.vcf.gz"
if [ -f "$MERGED" ] && [ -f "${MERGED}.tbi" ]; then
    echo "  [SKIP] Already merged"
else
    # Build list of per-chr VCFs in chromosome order
    vcf_list=""
    for chr in $CHROMOSOMES; do
        f="${PROJ}/vcf/per_chr/${chr}.vcf.gz"
        if [ -f "$f" ]; then
            vcf_list="${vcf_list} ${f}"
        fi
    done

    ${BCFTOOLS} concat $vcf_list -a -Oz -o "${MERGED}"
    ${TABIX} -p vcf "${MERGED}"
    n_total=$(${BCFTOOLS} view -H "${MERGED}" | wc -l | tr -d ' ')
    echo "  Merged VCF: ${n_total} total variants"
fi

# =========================================================================
# STEP 3: Filter variants
# =========================================================================
echo ""
echo "=== STEP 3: Filtering variants ==="

FILTERED="${PROJ}/vcf/freebayes_filtered.vcf.gz"
if [ -f "$FILTERED" ] && [ -f "${FILTERED}.tbi" ]; then
    echo "  [SKIP] Already filtered"
else
    # Filter: biallelic, QUAL>30, at least 5 total reads, normalize
    ${BCFTOOLS} view "${MERGED}" \
        -m2 -M2 \
        -i 'QUAL>30 && INFO/DP>50' \
    | ${BCFTOOLS} norm -m -any -f "${REF}" \
    | ${BGZIP} -c > "${FILTERED}"
    ${TABIX} -p vcf "${FILTERED}"

    n_snps=$(${BCFTOOLS} view -H -v snps "${FILTERED}" | wc -l | tr -d ' ')
    n_indels=$(${BCFTOOLS} view -H -v indels "${FILTERED}" | wc -l | tr -d ' ')
    echo "  Filtered: ${n_snps} SNPs, ${n_indels} indels"
fi

# =========================================================================
# STEP 4: Convert to PLINK format
# =========================================================================
echo ""
echo "=== STEP 4: PLINK conversion + QC ==="

# Convert VCF to PLINK binary
${PLINK2} \
    --vcf "${FILTERED}" \
    --make-bed \
    --set-all-var-ids '@:#:\$r:\$a' \
    --chr-set 27 --allow-extra-chr \
    --out "${PROJ}/plink/variants_raw" \
    2>&1 | tail -5

# Variant QC: MAF 0.10, geno 0.2 (matching the gwas_maf0.1geno0.2 settings)
${PLINK2} \
    --bfile "${PROJ}/plink/variants_raw" \
    --maf 0.10 \
    --geno 0.2 \
    --make-bed \
    --chr-set 27 --allow-extra-chr \
    --out "${PROJ}/plink/variants_qc" \
    2>&1 | tail -5

# Sample QC (mind 0.1)
${PLINK2} \
    --bfile "${PROJ}/plink/variants_qc" \
    --mind 0.1 \
    --make-bed \
    --chr-set 27 --allow-extra-chr \
    --out "${PROJ}/plink/variants_qc_samples" \
    2>&1 | tail -5

echo "  Final dataset:"
wc -l "${PROJ}/plink/variants_qc_samples.bim" "${PROJ}/plink/variants_qc_samples.fam"

# =========================================================================
# STEP 5: PCA for population structure
# =========================================================================
echo ""
echo "=== STEP 5: PCA ==="

# LD pruning
${PLINK2} \
    --bfile "${PROJ}/plink/variants_qc_samples" \
    --indep-pairwise 50 5 0.2 \
    --chr-set 27 --allow-extra-chr \
    --out "${PROJ}/plink/ld_prune" \
    2>&1 | tail -3

# PCA on pruned set
${PLINK2} \
    --bfile "${PROJ}/plink/variants_qc_samples" \
    --extract "${PROJ}/plink/ld_prune.prune.in" \
    --pca 10 \
    --chr-set 27 --allow-extra-chr \
    --out "${PROJ}/plink/pca" \
    2>&1 | tail -3

# =========================================================================
# STEP 6: Create phenotype file
# =========================================================================
echo ""
echo "=== STEP 6: Phenotype file ==="

# Create phenotype file matching PLINK sample IDs
# Read samples.tsv: NN=1 (control), DD=2 (case) for PLINK
python3 -c "
import sys
samples_tsv = '${BASE}/samples.tsv'
fam_file = '${PROJ}/plink/variants_qc_samples.fam'
drop = set('nr89 nr22 nr36 nr64 nr67 nr81 nr65 nr88 nr76 nr50'.split())

# Read group assignments
groups = {}
with open(samples_tsv) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        groups[parts[0]] = parts[2]  # NN or DD

# Read FAM to get actual sample IDs in the dataset
with open(fam_file) as f:
    fam_samples = [line.strip().split()[1] for line in f]

# Write phenotype file
with open('${PROJ}/phenotypes.txt', 'w') as out:
    out.write('FID\tIID\tPHENO1\n')
    nn, dd = 0, 0
    for s in fam_samples:
        g = groups.get(s, '')
        pheno = 1 if g == 'NN' else 2 if g == 'DD' else -9
        if pheno == 1: nn += 1
        elif pheno == 2: dd += 1
        out.write(f'0\t{s}\t{pheno}\n')

print(f'  Phenotype file: {nn} NN (control=1), {dd} DD (case=2)')
"

# =========================================================================
# STEP 7: GWAS association (logistic regression with Firth fallback)
# =========================================================================
echo ""
echo "=== STEP 7: GWAS ==="

${PLINK2} \
    --bfile "${PROJ}/plink/variants_qc_samples" \
    --pheno "${PROJ}/phenotypes.txt" \
    --covar "${PROJ}/plink/pca.eigenvec" \
    --covar-col-nums 3-12 \
    --glm firth-fallback hide-covar \
    --out "${PROJ}/gwas/association" \
    --chr-set 27 --allow-extra-chr \
    2>&1 | tail -5

echo "  GWAS results:"
gwas_file="${PROJ}/gwas/association.PHENO1.glm.logistic.hybrid"
if [ -f "$gwas_file" ]; then
    total=$(wc -l < "$gwas_file")
    sig=$(awk '$NF != "NA" && $NF != "P" && $NF+0 < 5e-8' "$gwas_file" | wc -l | tr -d ' ')
    sug=$(awk '$NF != "NA" && $NF != "P" && $NF+0 < 1e-5' "$gwas_file" | wc -l | tr -d ' ')
    echo "  Total variants: $((total - 1))"
    echo "  Genome-wide significant (p<5e-8): $sig"
    echo "  Suggestive (p<1e-5): $sug"
fi

# =========================================================================
# STEP 8: QTLseq analysis
# =========================================================================
echo ""
echo "=== STEP 8: QTLseq analysis ==="

python3 "${PROJ}/qtlseq_analysis.py"

echo ""
echo "============================================"
echo "Pipeline complete: $(date)"
echo "============================================"
