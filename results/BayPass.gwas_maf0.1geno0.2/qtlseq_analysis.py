#!/usr/bin/env python3
"""
QTLseq Analysis — computes SNP-index and delta SNP-index between NN and DD pools.
Uses the FreeBayes VCF with individual genotype data to compute per-pool allele frequencies.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import base64
import io
import gzip
import sys
import os

PROJ = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/BayPass.gwas_maf0.1geno0.2"
BASE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS"
VCF = f"{PROJ}/vcf/freebayes_filtered.vcf.gz"
SAMPLES_TSV = f"{BASE}/samples.tsv"
PHENO_FILE = f"{PROJ}/phenotypes.txt"

WINDOW_SIZE = 1_000_000  # 1 Mb sliding window
STEP_SIZE = 100_000      # 100 kb step

print("=== QTLseq Analysis ===")

# Read sample groups
groups = {}
with open(SAMPLES_TSV) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        groups[parts[0]] = parts[2]  # NN or DD

# Read phenotype file to get actual samples in dataset
pheno_samples = {}
with open(PHENO_FILE) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        pheno_samples[parts[1]] = int(parts[2])  # 1=NN, 2=DD

# Parse VCF header to get sample columns
print("Reading VCF...")
vcf_samples = []
nn_idx = []
dd_idx = []

opener = gzip.open if VCF.endswith('.gz') else open

with opener(VCF, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            vcf_samples = fields[9:]
            for i, s in enumerate(vcf_samples):
                g = groups.get(s, '')
                if g == 'NN':
                    nn_idx.append(i)
                elif g == 'DD':
                    dd_idx.append(i)
            break

print(f"  VCF samples: {len(vcf_samples)} (NN={len(nn_idx)}, DD={len(dd_idx)})")

# Parse genotypes and compute per-pool allele counts
# Store: chr, pos, nn_ref, nn_alt, dd_ref, dd_alt
data = []  # list of (chr, pos, nn_ref, nn_alt, dd_ref, dd_alt)
chr_order = [str(i) for i in range(1, 28)] + ['Z', 'W']

with opener(VCF, 'rt') as f:
    count = 0
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])

        # Skip non-chromosome scaffolds
        if chrom not in chr_order:
            continue

        # Get genotype field index (GT is always first in FORMAT)
        fmt = fields[8].split(':')
        gt_idx = fmt.index('GT') if 'GT' in fmt else 0

        # Count alleles per pool
        nn_ref_count = 0
        nn_alt_count = 0
        dd_ref_count = 0
        dd_alt_count = 0

        genos = fields[9:]
        for i in nn_idx:
            if i >= len(genos):
                continue
            gt = genos[i].split(':')[gt_idx]
            if '.' in gt:
                continue
            alleles = gt.replace('|', '/').split('/')
            for a in alleles:
                if a == '0':
                    nn_ref_count += 1
                else:
                    nn_alt_count += 1

        for i in dd_idx:
            if i >= len(genos):
                continue
            gt = genos[i].split(':')[gt_idx]
            if '.' in gt:
                continue
            alleles = gt.replace('|', '/').split('/')
            for a in alleles:
                if a == '0':
                    dd_ref_count += 1
                else:
                    dd_alt_count += 1

        # Only keep sites with data in both pools
        nn_total = nn_ref_count + nn_alt_count
        dd_total = dd_ref_count + dd_alt_count
        if nn_total >= 10 and dd_total >= 10:  # minimum depth filter
            data.append((chrom, pos, nn_ref_count, nn_alt_count, dd_ref_count, dd_alt_count))

        count += 1
        if count % 500000 == 0:
            print(f"  Processed {count} variants, kept {len(data)}...")

print(f"  Total variants for QTLseq: {len(data)}")

if len(data) == 0:
    print("ERROR: No variants passed filters!")
    sys.exit(1)

# Convert to arrays
chroms = [d[0] for d in data]
positions = np.array([d[1] for d in data])
nn_ref = np.array([d[2] for d in data], dtype=float)
nn_alt = np.array([d[3] for d in data], dtype=float)
dd_ref = np.array([d[4] for d in data], dtype=float)
dd_alt = np.array([d[5] for d in data], dtype=float)

# Compute SNP-index per pool
# SNP-index = alt_count / total_count (frequency of alt allele)
nn_snp_index = nn_alt / (nn_ref + nn_alt)
dd_snp_index = dd_alt / (dd_ref + dd_alt)

# Delta SNP-index = DD_index - NN_index
delta_snp_index = dd_snp_index - nn_snp_index

# =========================================================================
# Sliding window analysis
# =========================================================================
print("Computing sliding windows...")

window_data = {}  # chr -> list of (center, mean_delta, n_snps)

for chrom in chr_order:
    mask = np.array([c == chrom for c in chroms])
    if mask.sum() == 0:
        continue

    chr_pos = positions[mask]
    chr_delta = delta_snp_index[mask]

    min_pos = chr_pos.min()
    max_pos = chr_pos.max()

    windows = []
    for start in range(int(min_pos), int(max_pos), STEP_SIZE):
        end = start + WINDOW_SIZE
        w_mask = (chr_pos >= start) & (chr_pos < end)
        n = w_mask.sum()
        if n >= 5:  # minimum SNPs per window
            center = start + WINDOW_SIZE // 2
            mean_delta = chr_delta[w_mask].mean()
            windows.append((center, mean_delta, n))

    if windows:
        window_data[chrom] = windows

# =========================================================================
# Save results
# =========================================================================
print("Saving results...")

# Per-SNP results
with open(f"{PROJ}/qtlseq/snp_index.tsv", 'w') as f:
    f.write("CHR\tPOS\tNN_REF\tNN_ALT\tDD_REF\tDD_ALT\tNN_INDEX\tDD_INDEX\tDELTA\n")
    for i in range(len(data)):
        f.write(f"{chroms[i]}\t{positions[i]}\t{int(nn_ref[i])}\t{int(nn_alt[i])}\t"
                f"{int(dd_ref[i])}\t{int(dd_alt[i])}\t{nn_snp_index[i]:.4f}\t"
                f"{dd_snp_index[i]:.4f}\t{delta_snp_index[i]:.4f}\n")

# Window results
with open(f"{PROJ}/qtlseq/window_delta.tsv", 'w') as f:
    f.write("CHR\tCENTER\tMEAN_DELTA\tN_SNPS\n")
    for chrom in chr_order:
        if chrom in window_data:
            for center, mean_d, n in window_data[chrom]:
                f.write(f"{chrom}\t{center}\t{mean_d:.6f}\t{n}\n")

# =========================================================================
# Generate plots
# =========================================================================
print("Generating QTLseq plots...")

# Build cumulative positions for genome-wide plot
chr_cum_offset = {}
cum = 0
for chrom in chr_order:
    mask = np.array([c == chrom for c in chroms])
    if mask.sum() == 0:
        continue
    chr_cum_offset[chrom] = cum
    cum += positions[mask].max() + 2_000_000

# --- Genome-wide delta SNP-index Manhattan ---
fig, ax = plt.subplots(figsize=(16, 5))

chr_centers = {}
for chrom in chr_order:
    mask = np.array([c == chrom for c in chroms])
    if mask.sum() == 0 or chrom not in chr_cum_offset:
        continue

    x = positions[mask] + chr_cum_offset[chrom]
    y = delta_snp_index[mask]
    chr_num = chr_order.index(chrom)
    color = '#2166ac' if chr_num % 2 == 0 else '#67a9cf'
    ax.scatter(x / 1e6, y, s=0.5, alpha=0.3, c=color, edgecolors='none', rasterized=True)
    chr_centers[chrom] = x.mean() / 1e6

    # Plot window means
    if chrom in window_data:
        wx = np.array([w[0] + chr_cum_offset[chrom] for w in window_data[chrom]]) / 1e6
        wy = np.array([w[1] for w in window_data[chrom]])
        ax.plot(wx, wy, color='red', linewidth=0.8, alpha=0.8)

ax.axhline(0, color='black', linewidth=0.3)
ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('Δ SNP-index (DD − NN)', fontsize=12)
ax.set_title('QTLseq — Genome-wide Δ SNP-index', fontsize=14, fontweight='bold')
ax.set_xticks(list(chr_centers.values()))
ax.set_xticklabels(list(chr_centers.keys()), fontsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
fig.savefig(f"{PROJ}/plots/qtlseq_genomewide.png", dpi=200, bbox_inches='tight')

buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
buf.seek(0)
b64_gw = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# --- Chromosome 3 zoom ---
mask3 = np.array([c == '3' for c in chroms])
if mask3.sum() > 0:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True,
                                     gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.1})

    pos3 = positions[mask3] / 1e6
    delta3 = delta_snp_index[mask3]
    nn3 = nn_snp_index[mask3]
    dd3 = dd_snp_index[mask3]

    # Panel 1: SNP-index per pool
    ax1.scatter(pos3, nn3, s=1, alpha=0.3, c='#3498db', label='NN', edgecolors='none', rasterized=True)
    ax1.scatter(pos3, dd3, s=1, alpha=0.3, c='#e74c3c', label='DD', edgecolors='none', rasterized=True)
    ax1.set_ylabel('SNP-index', fontsize=11)
    ax1.set_title('QTLseq — Chromosome 3', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=9, markerscale=5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_ylim(-0.05, 1.05)

    # Panel 2: Delta SNP-index
    ax2.scatter(pos3, delta3, s=1, alpha=0.3, c='#2c3e50', edgecolors='none', rasterized=True)
    ax2.axhline(0, color='black', linewidth=0.3)

    if '3' in window_data:
        wx = np.array([w[0] for w in window_data['3']]) / 1e6
        wy = np.array([w[1] for w in window_data['3']])
        ax2.plot(wx, wy, color='red', linewidth=1.5, alpha=0.9, label='1 Mb window mean')
        ax2.legend(fontsize=9)

    ax2.set_xlabel('Position (Mb)', fontsize=11)
    ax2.set_ylabel('Δ SNP-index (DD − NN)', fontsize=11)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(f"{PROJ}/plots/qtlseq_chr3.png", dpi=200, bbox_inches='tight')

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
    buf.seek(0)
    b64_chr3 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
else:
    b64_chr3 = ""

# =========================================================================
# Summary statistics
# =========================================================================
print("\n=== QTLseq Summary ===")
print(f"Total variants analyzed: {len(data)}")
print(f"Mean delta SNP-index: {delta_snp_index.mean():.4f}")
print(f"Std delta SNP-index: {delta_snp_index.std():.4f}")

# Find top windows
print("\nTop 10 windows by |delta SNP-index|:")
all_windows = []
for chrom in chr_order:
    if chrom in window_data:
        for center, mean_d, n in window_data[chrom]:
            all_windows.append((chrom, center, mean_d, n))

all_windows.sort(key=lambda x: abs(x[2]), reverse=True)
for chrom, center, mean_d, n in all_windows[:10]:
    print(f"  chr{chrom}:{center/1e6:.1f}Mb  delta={mean_d:+.4f}  n_snps={n}")

# Save base64 images for HTML report generation
with open(f"{PROJ}/qtlseq/plot_b64_gw.txt", 'w') as f:
    f.write(b64_gw)
if b64_chr3:
    with open(f"{PROJ}/qtlseq/plot_b64_chr3.txt", 'w') as f:
        f.write(b64_chr3)

print("\nQTLseq analysis complete!")
