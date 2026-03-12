#!/usr/bin/env python3
"""
Generate comprehensive HTML report for FreeBayes GWAS + QTLseq analysis.
Includes: Manhattan plot, QQ plot, chr3 zoom, QTLseq plots, summary stats.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import base64
import io
import os

PROJ = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/BayPass.gwas_maf0.1geno0.2"
BASE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS"
GWAS = f"{PROJ}/gwas/association_3pc.PHENO1.glm.logistic.hybrid"
GWAS_NOPC = f"{PROJ}/gwas_no_covar/association_noPC.PHENO1.glm.logistic.hybrid"
GWAS_V2 = f"{PROJ}/gwas_v2/association_v2.PHENO1.glm.logistic.hybrid"

def read_gwas(path):
    """Read GWAS results, return chroms, positions, pvals arrays."""
    cl, pl, pv = [], [], []
    with open(path) as f:
        header = f.readline()
        cols = header.strip().split('\t')
        p_idx = cols.index('P') if 'P' in cols else len(cols) - 1
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= p_idx:
                continue
            p_str = parts[p_idx]
            if p_str == 'NA':
                continue
            try:
                p = float(p_str)
                if p > 0:
                    cl.append(parts[0])
                    pl.append(int(parts[1]))
                    pv.append(p)
            except ValueError:
                continue
    return np.array(cl), np.array(pl), np.array(pv)

# ===== Read GWAS results =====
print("Reading GWAS results (3 PCs)...")
chroms, positions, pvals = read_gwas(GWAS)
logp = -np.log10(pvals)
# Keep old-style lists for top-hits table
chrom_list = list(chroms)
pos_list = list(positions)
p_list = list(pvals)
print(f"  {len(pvals)} variants with valid p-values")

print("Reading GWAS results (no covariates)...")
chroms_nc, positions_nc, pvals_nc = read_gwas(GWAS_NOPC)
logp_nc = -np.log10(pvals_nc)
print(f"  {len(pvals_nc)} variants with valid p-values")

print("Reading GWAS results (v2: no LD prune PCA)...")
chroms_v2, positions_v2, pvals_v2 = read_gwas(GWAS_V2)
logp_v2 = -np.log10(pvals_v2)
print(f"  {len(pvals_v2)} variants with valid p-values")

chr_order = [str(i) for i in range(1, 28)] + ['Z', 'W']

# ===== Cumulative positions =====
chr_cum_offset = {}
cum = 0
for c in chr_order:
    mask = chroms == c
    if mask.sum() == 0:
        continue
    chr_cum_offset[c] = cum
    cum += positions[mask].max() + 2_000_000

# ===== Manhattan Plot =====
print("Generating Manhattan plot...")
fig, ax = plt.subplots(figsize=(16, 5))
chr_centers = {}
colors = ['#1f77b4', '#aec7e8']

for c in chr_order:
    mask = chroms == c
    if mask.sum() == 0 or c not in chr_cum_offset:
        continue
    x = (positions[mask] + chr_cum_offset[c]) / 1e6
    y = logp[mask]
    ci = chr_order.index(c) % 2
    ax.scatter(x, y, s=1, alpha=0.4, c=colors[ci], edgecolors='none', rasterized=True)
    chr_centers[c] = x.mean()

ax.axhline(-np.log10(5e-8), color='red', linewidth=0.8, linestyle='--', label='Genome-wide (5e-8)')
ax.axhline(-np.log10(1e-5), color='blue', linewidth=0.5, linestyle='--', alpha=0.5, label='Suggestive (1e-5)')
ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('-log10(p)', fontsize=12)
ax.set_title('GWAS Manhattan Plot — FreeBayes + PLINK2 (3 PCs)', fontsize=14, fontweight='bold')
ax.set_xticks(list(chr_centers.values()))
ax.set_xticklabels(list(chr_centers.keys()), fontsize=7)
ax.legend(fontsize=8, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
buf.seek(0)
manhattan_b64 = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# ===== QQ Plot =====
print("Generating QQ plot...")
n = len(pvals)
expected = -np.log10(np.arange(1, n + 1) / (n + 1))
observed = np.sort(logp)[::-1]

# Genomic inflation factor
chi2 = stats.chi2.isf(pvals, 1)
lambda_gc = np.nanmedian(chi2) / stats.chi2.ppf(0.5, 1)

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(expected, observed, s=1, alpha=0.3, c='#2c3e50', edgecolors='none', rasterized=True)
max_val = max(expected.max(), observed.max()) + 0.5
ax.plot([0, max_val], [0, max_val], 'r-', linewidth=0.8, alpha=0.7)
ax.set_xlabel('Expected -log10(p)', fontsize=12)
ax.set_ylabel('Observed -log10(p)', fontsize=12)
ax.set_title(f'QQ Plot (λ_GC = {lambda_gc:.3f})', fontsize=14, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
buf.seek(0)
qq_b64 = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# ===== No-covariate GWAS: compute cumulative positions =====
chr_cum_offset_nc = {}
cum_nc = 0
for c in chr_order:
    mask = chroms_nc == c
    if mask.sum() == 0:
        continue
    chr_cum_offset_nc[c] = cum_nc
    cum_nc += positions_nc[mask].max() + 2_000_000

# Lambda for no-covariate
chi2_nc = stats.chi2.isf(pvals_nc, 1)
lambda_gc_nc = np.nanmedian(chi2_nc) / stats.chi2.ppf(0.5, 1)

n_sig_nc = (pvals_nc < 5e-8).sum()
n_sug_nc = (pvals_nc < 1e-5).sum()

# v2: compute cumulative positions, lambda, counts
chr_cum_offset_v2 = {}
cum_v2 = 0
for c in chr_order:
    mask = chroms_v2 == c
    if mask.sum() == 0:
        continue
    chr_cum_offset_v2[c] = cum_v2
    cum_v2 += positions_v2[mask].max() + 2_000_000

chi2_v2 = stats.chi2.isf(pvals_v2, 1)
lambda_gc_v2 = np.nanmedian(chi2_v2) / stats.chi2.ppf(0.5, 1)
n_sig_v2 = (pvals_v2 < 5e-8).sum()
n_sug_v2 = (pvals_v2 < 1e-5).sum()

# ===== V2 COMPARISON: 3-panel Manhattan (original vs v2 vs no-covar) =====
print("Generating v2 comparison plots...")
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 12), sharex=False)

# Shared y-limit
ymax_all = max(logp.max(), logp_nc.max(), logp_v2.max()) + 1

# Panel 1: Original (LD-pruned PCA + 3 PCs)
for c in chr_order:
    mask = chroms == c
    if mask.sum() == 0 or c not in chr_cum_offset:
        continue
    x = (positions[mask] + chr_cum_offset[c]) / 1e6
    y = logp[mask]
    ci = chr_order.index(c) % 2
    ax1.scatter(x, y, s=0.8, alpha=0.4, c=colors[ci], edgecolors='none', rasterized=True)
ax1.axhline(-np.log10(5e-8), color='red', linewidth=0.8, linestyle='--')
ax1.axhline(-np.log10(1e-5), color='blue', linewidth=0.5, linestyle='--', alpha=0.5)
ax1.set_ylabel('-log10(p)', fontsize=10)
ax1.set_title(f'Original: 3 PCs from LD-pruned PCA  (λ={lambda_gc:.3f}, {(pvals < 5e-8).sum()} sig.)', fontsize=11, fontweight='bold')
ax1.set_xticks(list(chr_centers.values()))
ax1.set_xticklabels(list(chr_centers.keys()), fontsize=6)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_ylim(0, ymax_all)

# Panel 2: v2 (no LD prune PCA + 3 PCs)
chr_centers_v2 = {}
for c in chr_order:
    mask = chroms_v2 == c
    if mask.sum() == 0 or c not in chr_cum_offset_v2:
        continue
    x = (positions_v2[mask] + chr_cum_offset_v2[c]) / 1e6
    y = logp_v2[mask]
    ci = chr_order.index(c) % 2
    ax2.scatter(x, y, s=0.8, alpha=0.4, c=colors[ci], edgecolors='none', rasterized=True)
    chr_centers_v2[c] = x.mean()
ax2.axhline(-np.log10(5e-8), color='red', linewidth=0.8, linestyle='--')
ax2.axhline(-np.log10(1e-5), color='blue', linewidth=0.5, linestyle='--', alpha=0.5)
ax2.set_ylabel('-log10(p)', fontsize=10)
ax2.set_title(f'v2: 3 PCs from non-LD-pruned PCA  (λ={lambda_gc_v2:.3f}, {n_sig_v2} sig.)', fontsize=11, fontweight='bold')
ax2.set_xticks(list(chr_centers_v2.values()))
ax2.set_xticklabels(list(chr_centers_v2.keys()), fontsize=6)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_ylim(0, ymax_all)

# Panel 3: No covariates
chr_centers_nc = {}
for c in chr_order:
    mask = chroms_nc == c
    if mask.sum() == 0 or c not in chr_cum_offset_nc:
        continue
    x = (positions_nc[mask] + chr_cum_offset_nc[c]) / 1e6
    y = logp_nc[mask]
    ci = chr_order.index(c) % 2
    ax3.scatter(x, y, s=0.8, alpha=0.4, c=colors[ci], edgecolors='none', rasterized=True)
    chr_centers_nc[c] = x.mean()
ax3.axhline(-np.log10(5e-8), color='red', linewidth=0.8, linestyle='--')
ax3.axhline(-np.log10(1e-5), color='blue', linewidth=0.5, linestyle='--', alpha=0.5)
ax3.set_xlabel('Chromosome', fontsize=10)
ax3.set_ylabel('-log10(p)', fontsize=10)
ax3.set_title(f'No covariates  (λ={lambda_gc_nc:.3f}, {n_sig_nc} sig.)', fontsize=11, fontweight='bold')
ax3.set_xticks(list(chr_centers_nc.values()))
ax3.set_xticklabels(list(chr_centers_nc.keys()), fontsize=6)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.set_ylim(0, ymax_all)

plt.tight_layout()
buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
buf.seek(0)
v2_manhattan_b64 = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# ===== V2: 3-panel QQ =====
print("Generating v2 QQ comparison...")
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

qmax_all = max(np.sort(logp)[-1], np.sort(logp_v2)[-1], np.sort(logp_nc)[-1]) + 1
exp_max = max(
    -np.log10(1 / (len(pvals) + 1)),
    -np.log10(1 / (len(pvals_v2) + 1)),
    -np.log10(1 / (len(pvals_nc) + 1))
) + 0.5

for ax_i, (pv, lp, lam, title, col) in enumerate([
    (pvals, logp, lambda_gc, 'Original (LD-pruned PCA)', '#2c3e50'),
    (pvals_v2, logp_v2, lambda_gc_v2, 'v2 (no LD prune PCA)', '#27ae60'),
    (pvals_nc, logp_nc, lambda_gc_nc, 'No covariates', '#e74c3c'),
]):
    ax = [ax1, ax2, ax3][ax_i]
    n_i = len(pv)
    exp_i = -np.log10(np.arange(1, n_i + 1) / (n_i + 1))
    obs_i = np.sort(lp)[::-1]
    ax.scatter(exp_i, obs_i, s=1, alpha=0.3, c=col, edgecolors='none', rasterized=True)
    ax.plot([0, exp_max], [0, exp_max], 'r-', linewidth=0.8, alpha=0.7)
    ax.set_xlabel('Expected -log10(p)', fontsize=10)
    ax.set_ylabel('Observed -log10(p)', fontsize=10)
    ax.set_title(f'{title}\n(λ = {lam:.3f})', fontsize=10, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, exp_max)
    ax.set_ylim(0, qmax_all)

plt.tight_layout()
buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
buf.seek(0)
v2_qq_b64 = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# ===== V2: Chr3 3-panel comparison =====
print("Generating v2 chr3 comparison...")
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

mask3_orig = chroms == '3'
mask3_v2 = chroms_v2 == '3'
mask3_nopc = chroms_nc == '3'
ymax_chr3 = max(logp[mask3_orig].max(), logp_v2[mask3_v2].max(), logp_nc[mask3_nopc].max()) + 1

for ax, m, pos, lp, title, col in [
    (ax1, mask3_orig, positions, logp, f'Original: LD-pruned PCA ({(pvals[mask3_orig] < 5e-8).sum()} sig.)', '#2c3e50'),
    (ax2, mask3_v2, positions_v2, logp_v2, f'v2: no LD prune PCA ({(pvals_v2[mask3_v2] < 5e-8).sum()} sig.)', '#27ae60'),
    (ax3, mask3_nopc, positions_nc, logp_nc, f'No covariates ({(pvals_nc[mask3_nopc] < 5e-8).sum()} sig.)', '#e74c3c'),
]:
    ax.scatter(pos[m] / 1e6, lp[m], s=2, alpha=0.5, c=col, edgecolors='none', rasterized=True)
    ax.axhline(-np.log10(5e-8), color='red', linewidth=0.8, linestyle='--')
    ax.axhline(-np.log10(1e-5), color='blue', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.set_ylabel('-log10(p)', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, ymax_chr3)

ax3.set_xlabel('Position (Mb)', fontsize=10)
plt.tight_layout()
buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
buf.seek(0)
v2_chr3_b64 = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# v2 per-chromosome comparison table (3 columns)
v2_compare_table = ""
for ch in chr_order:
    m_orig = chroms == ch
    m_v2 = chroms_v2 == ch
    m_nc = chroms_nc == ch
    h_orig = (pvals[m_orig] < 5e-8).sum() if m_orig.sum() > 0 else 0
    h_v2 = (pvals_v2[m_v2] < 5e-8).sum() if m_v2.sum() > 0 else 0
    h_nc = (pvals_nc[m_nc] < 5e-8).sum() if m_nc.sum() > 0 else 0
    if h_orig > 0 or h_v2 > 0 or h_nc > 0:
        v2_compare_table += f"<tr><td>{ch}</td><td>{h_orig}</td><td>{h_v2}</td><td>{h_nc}</td></tr>\n"

# ===== Chromosome 3 Manhattan =====
print("Generating chr3 Manhattan...")
mask3 = chroms == '3'
fig, ax = plt.subplots(figsize=(14, 4))
x3 = positions[mask3] / 1e6
y3 = logp[mask3]
ax.scatter(x3, y3, s=2, alpha=0.5, c='#2c3e50', edgecolors='none', rasterized=True)
ax.axhline(-np.log10(5e-8), color='red', linewidth=0.8, linestyle='--', label='5e-8')
ax.axhline(-np.log10(1e-5), color='blue', linewidth=0.5, linestyle='--', alpha=0.5, label='1e-5')
ax.set_xlabel('Position (Mb)', fontsize=12)
ax.set_ylabel('-log10(p)', fontsize=12)
ax.set_title('Chromosome 3 — GWAS Manhattan', fontsize=14, fontweight='bold')
ax.legend(fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
buf.seek(0)
chr3_b64 = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# ===== PCA Plot =====
print("Generating PCA plot...")
eigenvec = {}
with open(f"{PROJ}/plink/pca.eigenvec") as f:
    header_pca = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        eigenvec[parts[1]] = [float(x) for x in parts[2:]]

# Read groups
groups = {}
with open("/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/samples.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        groups[parts[0]] = parts[2]

fig, ax = plt.subplots(figsize=(9, 7))
nn_x, nn_y, nn_names = [], [], []
dd_x, dd_y, dd_names = [], [], []
outlier_samples = []  # PC2 > 0.1 and PC1 < -0.1

for sample, pcs in eigenvec.items():
    g = groups.get(sample, '?')
    if g == 'NN':
        nn_x.append(pcs[0])
        nn_y.append(pcs[1])
        nn_names.append(sample)
    elif g == 'DD':
        dd_x.append(pcs[0])
        dd_y.append(pcs[1])
        dd_names.append(sample)
    # Check outlier criteria
    if pcs[1] > 0.1 and pcs[0] < -0.1:
        outlier_samples.append((sample, g, pcs[0], pcs[1]))

ax.scatter(nn_x, nn_y, s=30, alpha=0.7, c='#3498db', label='NN (non-diapause)', edgecolors='white', linewidths=0.3)
ax.scatter(dd_x, dd_y, s=30, alpha=0.7, c='#e74c3c', label='DD (diapause)', edgecolors='white', linewidths=0.3)

# Label outlier samples
from matplotlib.offsetbox import AnnotationBbox
for sample, g, pc1, pc2 in outlier_samples:
    ax.annotate(sample, (pc1, pc2), fontsize=6, ha='left', va='bottom',
                xytext=(3, 3), textcoords='offset points', color='#2c3e50',
                fontweight='bold')

# Draw box around outlier region
import matplotlib.patches as mpatches
rect = mpatches.FancyBboxPatch((-0.20, 0.10), 0.12, 0.09,
                                boxstyle="round,pad=0.005",
                                facecolor='none', edgecolor='#e67e22',
                                linewidth=1.5, linestyle='--')
ax.add_patch(rect)
ax.annotate('Outlier cluster\n(PC1<-0.1, PC2>0.1)', xy=(-0.14, 0.19),
            fontsize=8, color='#e67e22', fontweight='bold', ha='center', va='bottom')

ax.set_xlabel('PC1', fontsize=12)
ax.set_ylabel('PC2', fontsize=12)
ax.set_title('PCA (PC1 vs PC2)', fontsize=14, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=200, bbox_inches='tight')
buf.seek(0)
pca_b64 = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# Build outlier table HTML with QC and family info
outlier_samples.sort(key=lambda x: x[0])

# Get family/cage info from QC file names
import re
sample_meta = {}
qc1_path = f"{BASE}/reports/qc1_reads/multiqc_reads_data/multiqc_general_stats.txt"
try:
    with open(qc1_path) as qf:
        qf.readline()
        for line in qf:
            m = re.search(r'(nr\d+)-(C\d+)-(NN|DD)-(\d+)', line)
            if m:
                sample_meta[m.group(1)] = {'cage': m.group(2), 'family': m.group(4)}
except FileNotFoundError:
    pass

# Get coverage from QC2
qc2_path = f"{BASE}/reports/qc2_alignment/multiqc_alignment_data/multiqc_general_stats.txt"
sample_cov = {}
try:
    with open(qc2_path) as qf:
        hdr = qf.readline().strip().split('\t')
        cov_idx = hdr.index('mosdepth-mean_coverage') if 'mosdepth-mean_coverage' in hdr else None
        dup_idx = hdr.index('picard_duplicates_mark_duplicates-PERCENT_DUPLICATION') if 'picard_duplicates_mark_duplicates-PERCENT_DUPLICATION' in hdr else None
        for line in qf:
            parts = line.strip().split('\t')
            m = re.match(r'^(nr\d+)$', parts[0])
            if m:
                sid = m.group(1)
                sample_cov[sid] = {
                    'mean_cov': float(parts[cov_idx]) if cov_idx and cov_idx < len(parts) and parts[cov_idx] else None,
                    'pct_dup': float(parts[dup_idx]) if dup_idx and dup_idx < len(parts) and parts[dup_idx] else None,
                }
except FileNotFoundError:
    pass

pca_outlier_table = ""
for sample, g, pc1, pc2 in outlier_samples:
    meta = sample_meta.get(sample, {})
    cov_info = sample_cov.get(sample, {})
    cage = meta.get('cage', '?')
    family = meta.get('family', '?')
    mcov = cov_info.get('mean_cov')
    mcov_str = f"{mcov:.1f}x" if mcov else '?'
    pca_outlier_table += f"<tr><td>{sample}</td><td>{g}</td><td>{cage}</td><td>{family}</td><td>{mcov_str}</td><td>{pc1:.4f}</td><td>{pc2:.4f}</td></tr>\n"
pca_outlier_count = len(outlier_samples)

# Family summary for outliers vs non-outliers
outlier_ids = set(s for s, _, _, _ in outlier_samples)
outlier_families = sorted(set(sample_meta.get(s, {}).get('family', '?') for s in outlier_ids))
non_outlier_nn_families = sorted(set(
    sample_meta.get(s, {}).get('family', '?')
    for s, info in sample_meta.items()
    if info.get('cage') and s not in outlier_ids and 'NN' in str(groups.get(s, ''))
))
outlier_family_str = ', '.join(outlier_families)
non_outlier_family_str = ', '.join(non_outlier_nn_families)

# ===== Read QTLseq plots =====
print("Reading QTLseq plots...")
qtl_gw_b64 = ""
qtl_chr3_b64 = ""
try:
    with open(f"{PROJ}/qtlseq/plot_b64_gw.txt") as f:
        qtl_gw_b64 = f.read().strip()
except FileNotFoundError:
    pass
try:
    with open(f"{PROJ}/qtlseq/plot_b64_chr3.txt") as f:
        qtl_chr3_b64 = f.read().strip()
except FileNotFoundError:
    pass

# ===== Compute summary statistics =====
print("Computing summary stats...")
n_sig = (pvals < 5e-8).sum()
n_sug = (pvals < 1e-5).sum()

# Top hits table
top_idx = np.argsort(pvals)[:30]
top_table = ""
for i in top_idx:
    top_table += f"<tr><td>{chrom_list[i]}</td><td>{pos_list[i]:,}</td><td>{logp[i]:.2f}</td><td>{pvals[i]:.2e}</td></tr>\n"

# ===== Generate HTML =====
print("Writing HTML report...")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>FreeBayes GWAS + QTLseq Report</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; background: #fafafa; }}
h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
h2 {{ color: #34495e; margin-top: 40px; border-bottom: 1px solid #bdc3c7; padding-bottom: 8px; }}
h3 {{ color: #5a6c7d; margin-top: 25px; }}
.summary-box {{ background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 20px; margin: 15px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
.stat {{ display: inline-block; margin: 10px 20px 10px 0; padding: 10px 20px; background: #ecf0f1; border-radius: 5px; }}
.stat .value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
.stat .label {{ font-size: 12px; color: #7f8c8d; }}
img {{ max-width: 100%; border: 1px solid #ddd; border-radius: 4px; margin: 10px 0; }}
table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
th {{ background: #3498db; color: white; }}
tr:nth-child(even) {{ background: #f2f2f2; }}
.highlight {{ background: #fff3cd; border: 1px solid #ffc107; padding: 15px; border-radius: 5px; margin: 15px 0; }}
.method {{ background: #e8f4fd; border-left: 4px solid #3498db; padding: 15px; margin: 15px 0; }}
pre.code-block {{ background: #1e1e1e; color: #d4d4d4; padding: 16px; border-radius: 6px; overflow-x: auto; font-size: 13px; line-height: 1.5; margin: 10px 0; }}
pre.code-block .comment {{ color: #6a9955; }}
pre.code-block .result {{ color: #9cdcfe; font-style: italic; }}
.code-section {{ background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 20px; margin: 15px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
.step-label {{ display: inline-block; background: #3498db; color: white; padding: 2px 10px; border-radius: 3px; font-weight: bold; font-size: 13px; margin-bottom: 8px; }}
</style>
</head>
<body>

<h1>FreeBayes GWAS + QTLseq Report</h1>
<p><em>Pararge aegeria — NN (non-diapause) vs DD (diapause)</em></p>

<div class="summary-box">
<div class="stat"><div class="value">84</div><div class="label">Samples (42 NN + 42 DD)</div></div>
<div class="stat"><div class="value">{len(pvals):,}</div><div class="label">Variants tested</div></div>
<div class="stat"><div class="value">{n_sig}</div><div class="label">Genome-wide significant</div></div>
<div class="stat"><div class="value">{n_sug}</div><div class="label">Suggestive (p&lt;1e-5)</div></div>
<div class="stat"><div class="value">{lambda_gc:.3f}</div><div class="label">Genomic inflation (λ_GC)</div></div>
</div>

<div class="method">
<strong>Methods summary:</strong> Variants called with <strong>FreeBayes v1.3.10</strong> (haplotype-based caller) from 84 BAM files.
Filtered: biallelic SNPs+indels, QUAL&gt;30, DP&gt;50, normalized with bcftools.
QC: MAF&gt;0.10, genotyping rate&gt;80%.
GWAS: PLINK2 logistic regression with Firth fallback, 3 PCA covariates.
QTLseq: per-pool (NN vs DD) allele frequency analysis with 1 Mb sliding windows.
</div>

<h2>Analysis Code</h2>

<div class="code-section">
<h3>1. Variant Calling (FreeBayes)</h3>
<p>FreeBayes was run per chromosome on all 84 BAM files (29 chromosomes: 1&ndash;27, Z, W):</p>
<pre class="code-block"><span class="comment"># FreeBayes per-chromosome variant calling</span>
freebayes \\
    --fasta-reference GCF_905163445.1_ilParAegt1.1_genomic.edited.fna \\
    --bam-list bam_list.txt \\
    --region ${{chr}} \\
    --min-mapping-quality 20 \\
    --min-base-quality 20 \\
    --min-alternate-count 2 \\
    --min-alternate-fraction 0.05 \\
    --use-best-n-alleles 4 \\
  | bcftools view -i 'QUAL&gt;20' \\
  | bgzip -c &gt; per_chr/${{chr}}.vcf.gz

<span class="comment"># Merge all per-chromosome VCFs</span>
bcftools concat per_chr/{{1..27}}.vcf.gz per_chr/Z.vcf.gz per_chr/W.vcf.gz \\
    -a -Oz -o freebayes_all.vcf.gz
<span class="result"># Result: 5,798,154 variants (4.3M SNPs, 495K MNPs, 833K indels)</span>
</pre>

<h3>2. Variant Filtering (bcftools)</h3>
<pre class="code-block"><span class="comment"># Biallelic only, QUAL&gt;30, total DP&gt;50, normalize</span>
bcftools view freebayes_all.vcf.gz \\
    -m2 -M2 \\
    -i 'QUAL&gt;30 &amp;&amp; INFO/DP&gt;50' \\
  | bcftools norm -m -any -f reference.fna \\
  | bgzip -c &gt; freebayes_filtered.vcf.gz
<span class="result"># Result: 4,748,634 variants</span>
</pre>

<h3>3. PLINK2 Conversion &amp; QC</h3>
<pre class="code-block"><span class="comment"># Step 1: VCF &rarr; PLINK binary format</span>
plink2 --vcf freebayes_filtered.vcf.gz \\
    --make-bed \\
    --set-all-var-ids '@:#:$r:$a' \\
    --new-id-max-allele-len 80 \\
    --chr-set 27 --allow-extra-chr \\
    --out variants_raw
<span class="result"># 4,748,634 variants, 84 samples loaded</span>

<span class="comment"># Step 2: Variant QC &mdash; MAF and genotyping rate</span>
plink2 --bfile variants_raw \\
    --maf 0.10 \\
    --geno 0.2 \\
    --make-bed \\
    --chr-set 27 --allow-extra-chr \\
    --out variants_qc
<span class="result"># 95,447 variants removed due to missing genotype data (--geno 0.2)
# 1,532,987 variants removed due to allele frequency threshold (--maf 0.10)
# 3,120,200 variants remaining</span>

<span class="comment"># Step 3: Sample QC &mdash; missingness filter</span>
plink2 --bfile variants_qc \\
    --mind 0.1 \\
    --make-bed \\
    --chr-set 27 --allow-extra-chr \\
    --out variants_qc_samples
<span class="result"># 0 samples removed (all 84 passed --mind 0.1)</span>
</pre>

<h3>4. LD Pruning &amp; PCA</h3>
<pre class="code-block"><span class="comment"># LD pruning for PCA</span>
plink2 --bfile variants_qc_samples \\
    --indep-pairwise 50 5 0.2 \\
    --chr-set 27 --allow-extra-chr \\
    --out ld_prune
<span class="result"># 2,933,379 of 3,120,200 variants removed
# 186,821 independent variants retained</span>

<span class="comment"># PCA on LD-pruned variants</span>
plink2 --bfile variants_qc_samples \\
    --extract ld_prune.prune.in \\
    --pca 10 \\
    --chr-set 27 --allow-extra-chr \\
    --out pca
<span class="result"># 10 principal components computed from 186,821 variants, 84 samples</span>
</pre>

<h3>5a. GWAS Association (3 PC covariates)</h3>
<pre class="code-block"><span class="comment"># Logistic regression with Firth fallback</span>
<span class="comment"># Covariates: PC1, PC2, PC3 from pca.eigenvec (columns 3-5)</span>
<span class="comment"># Note: 10 PCs caused massive convergence failures; 3 PCs used instead</span>
plink2 --bfile variants_qc_samples \\
    --pheno phenotypes.txt \\
    --covar pca.eigenvec \\
    --covar-col-nums 3-5 \\
    --glm firth-fallback hide-covar \\
    --chr-set 27 --allow-extra-chr \\
    --out gwas/association_3pc
<span class="result"># Phenotype: 42 cases (DD) + 42 controls (NN)
# 3 covariates loaded (PC1, PC2, PC3)
# 3,116,803 variants with valid p-values
# 3,397 variants skipped (convergence failure)</span>
</pre>

<h3>5b. GWAS v2: PCA without LD pruning</h3>
<pre class="code-block"><span class="comment"># PCA on ALL QC'd variants (no LD pruning)</span>
plink2 --bfile variants_qc_samples \\
    --pca 10 \\
    --chr-set 27 --allow-extra-chr \\
    --out gwas_v2/pca_no_ldprune
<span class="result"># 10 PCs computed from 3,120,200 variants (no LD pruning)</span>

<span class="comment"># GWAS with 3 PCs from non-LD-pruned PCA</span>
plink2 --bfile variants_qc_samples \\
    --pheno phenotypes.txt \\
    --covar gwas_v2/pca_no_ldprune.eigenvec \\
    --covar-col-nums 3-5 \\
    --glm firth-fallback hide-covar \\
    --chr-set 27 --allow-extra-chr \\
    --out gwas_v2/association_v2
<span class="result"># {len(pvals_v2):,} variants with valid p-values</span>
</pre>

<h3>5c. GWAS without covariates</h3>
<pre class="code-block"><span class="comment"># Same model but without any PCA covariates</span>
plink2 --bfile variants_qc_samples \\
    --pheno phenotypes.txt \\
    --glm firth-fallback hide-covar allow-no-covars \\
    --chr-set 27 --allow-extra-chr \\
    --out gwas_no_covar/association_noPC
<span class="result"># {len(pvals_nc):,} variants with valid p-values</span>
</pre>

<h3>6. QTLseq Analysis</h3>
<pre class="code-block"><span class="comment"># Per-pool allele frequency from FreeBayes VCF individual genotypes</span>
<span class="comment"># Pools: NN (42 samples) vs DD (42 samples)</span>
<span class="comment"># Filter: min 10 alleles per pool per site</span>
<span class="comment"># SNP-index = alt_alleles / total_alleles per pool</span>
<span class="comment"># Delta SNP-index = DD_index - NN_index</span>
<span class="comment"># Sliding window: 1 Mb window, 100 kb step, min 5 SNPs per window</span>
python3 qtlseq_analysis.py
<span class="result"># 4,746,690 variants analyzed
# Top window: chr3:11.4 Mb (delta SNP-index = -0.138)</span>
</pre>
</div>

<h2>Manhattan Plot</h2>
<img src="data:image/png;base64,{manhattan_b64}" alt="Manhattan Plot">

<h2>QQ Plot</h2>
<div style="max-width: 500px;">
<img src="data:image/png;base64,{qq_b64}" alt="QQ Plot">
</div>

<h2>PCA (Population Structure)</h2>
<div style="max-width: 700px;">
<img src="data:image/png;base64,{pca_b64}" alt="PCA Plot">
</div>

<div class="highlight">
<strong>PCA outlier cluster ({pca_outlier_count} samples):</strong> Samples with PC1 &lt; -0.1 and PC2 &gt; 0.1.
All are <strong>NN (non-diapause)</strong> individuals. QC investigation (qc1 and qc2 reports) shows
<strong>no quality issues</strong> &mdash; coverage, mapping rates, duplication, and error rates are comparable
to the rest of the dataset. The clustering is driven by <strong>family/genetic structure</strong>:
these samples come from a distinct subset of families ({outlier_family_str}) that share
closer genetic relatedness, while other NN samples come from different families
({non_outlier_family_str}). Cages (C1/C2/C3) are evenly represented (5 each), ruling out
batch effects.
</div>
<table>
<tr><th>Sample</th><th>Group</th><th>Cage</th><th>Family</th><th>Mean Cov</th><th>PC1</th><th>PC2</th></tr>
{pca_outlier_table}
</table>

<h2>Chromosome 3 — GWAS Manhattan</h2>
<div class="highlight">
<strong>Key finding:</strong> Chromosome 3 contains the majority of genome-wide significant hits,
with strong signals spanning ~4-19 Mb. This is consistent with a large structural variant
(likely a deletion or inversion) differentiating NN and DD phenotypes.
</div>
<img src="data:image/png;base64,{chr3_b64}" alt="Chromosome 3 Manhattan">

<h2>QTLseq — Genome-wide Δ SNP-index</h2>
<img src="data:image/png;base64,{qtl_gw_b64}" alt="QTLseq Genome-wide">

<h2>QTLseq — Chromosome 3 Zoom</h2>
<img src="data:image/png;base64,{qtl_chr3_b64}" alt="QTLseq Chromosome 3">

<h2>GWAS Comparison: Original vs v2 vs No Covariates</h2>

<div class="highlight">
<strong>Effect of LD pruning and PCA covariates on GWAS:</strong><br><br>
<strong>Original (LD-pruned PCA, 3 PCs):</strong> The LD pruning step (--indep-pairwise 50 5 0.2) ensures PCA captures
genome-wide population structure rather than local LD. This correctly identifies 22 hits on chr3 with minimal
false positives (1 non-chr3 hit).<br><br>
<strong>v2 (no LD pruning, 3 PCs):</strong> Without LD pruning, the PCA is heavily influenced by the large chr3
LD block that differentiates NN and DD. The PCs therefore <em>absorb part of the chr3 signal itself</em>,
causing the GWAS to overcorrect and lose chr3 power (only {n_sig_v2} significant hits, {(pvals_v2[chroms_v2 == '3'] < 5e-8).sum()} on chr3).
The top hit shifts to chr1, and spurious hits appear across multiple chromosomes.<br><br>
<strong>No covariates:</strong> Without any correction, family structure inflates results massively
(&lambda;={lambda_gc_nc:.1f}), producing {n_sig_nc} spurious genome-wide significant hits across the genome.<br><br>
<strong>Conclusion:</strong> LD pruning before PCA is essential. It prevents the causal locus from being captured
by PCs, preserving signal while still correcting for population structure.
</div>

<div class="summary-box">
<table>
<tr><th>Metric</th><th>Original<br>(LD-pruned PCA)</th><th>v2<br>(no LD prune PCA)</th><th>No Covariates</th></tr>
<tr><td>PCA variants</td><td>186,821</td><td>3,120,200</td><td>&mdash;</td></tr>
<tr><td>Valid p-values</td><td>{len(pvals):,}</td><td>{len(pvals_v2):,}</td><td>{len(pvals_nc):,}</td></tr>
<tr><td>Genomic inflation (&lambda;_GC)</td><td><strong>{lambda_gc:.3f}</strong></td><td>{lambda_gc_v2:.3f}</td><td>{lambda_gc_nc:.3f}</td></tr>
<tr><td>Genome-wide significant (p&lt;5e-8)</td><td><strong>{n_sig}</strong></td><td>{n_sig_v2}</td><td>{n_sig_nc}</td></tr>
<tr><td>Suggestive (p&lt;1e-5)</td><td>{n_sug}</td><td>{n_sug_v2}</td><td>{n_sug_nc}</td></tr>
<tr><td>Chr3 significant</td><td><strong>{(pvals[chroms == '3'] < 5e-8).sum()}</strong></td><td>{(pvals_v2[chroms_v2 == '3'] < 5e-8).sum()}</td><td>{(pvals_nc[chroms_nc == '3'] < 5e-8).sum()}</td></tr>
<tr><td>Non-chr3 significant</td><td>{((pvals < 5e-8) & (chroms != '3')).sum()}</td><td>{((pvals_v2 < 5e-8) & (chroms_v2 != '3')).sum()}</td><td>{((pvals_nc < 5e-8) & (chroms_nc != '3')).sum()}</td></tr>
<tr><td>Top hit</td><td>chr3:6.25 Mb</td><td>chr1:18.20 Mb</td><td>chr3:10.63 Mb</td></tr>
</table>
</div>

<h3>Per-Chromosome Significant Hits (p&lt;5e-8)</h3>
<table>
<tr><th>Chr</th><th>Original (LD-pruned PCA)</th><th>v2 (no LD prune PCA)</th><th>No Covariates</th></tr>
{v2_compare_table}
</table>

<h3>Manhattan Comparison (3 panels)</h3>
<img src="data:image/png;base64,{v2_manhattan_b64}" alt="3-way Manhattan Comparison">

<h3>QQ Plot Comparison (3 panels)</h3>
<img src="data:image/png;base64,{v2_qq_b64}" alt="3-way QQ Comparison">

<h3>Chromosome 3 Comparison (3 panels)</h3>
<img src="data:image/png;base64,{v2_chr3_b64}" alt="3-way Chr3 Comparison">

<h2>Top 30 GWAS Hits (Original: LD-pruned PCA)</h2>
<table>
<tr><th>Chr</th><th>Position</th><th>-log10(p)</th><th>P-value</th></tr>
{top_table}
</table>

<h2>Pipeline Summary</h2>
<table>
<tr><th>Step</th><th>Detail</th></tr>
<tr><td>Variant caller</td><td>FreeBayes v1.3.10 (haplotype-based, per chromosome)</td></tr>
<tr><td>Raw variants</td><td>5,798,154 (4.3M SNPs, 495K MNPs, 833K indels)</td></tr>
<tr><td>After filtering</td><td>4,748,634 (biallelic, QUAL&gt;30, DP&gt;50, normalized)</td></tr>
<tr><td>After QC</td><td>3,120,200 (MAF&gt;0.10, geno&lt;0.2)</td></tr>
<tr><td>LD-pruned for PCA</td><td>186,821 variants</td></tr>
<tr><td>GWAS model</td><td>Logistic regression, Firth fallback, 3 PC covariates</td></tr>
<tr><td>QTLseq variants</td><td>4,746,690 (min 10 alleles per pool)</td></tr>
<tr><td>Top QTLseq window</td><td>chr3:11.4 Mb (Δ SNP-index = -0.138)</td></tr>
</table>

</body>
</html>
"""

with open(f"{PROJ}/summary_report.html", 'w') as f:
    f.write(html)

print(f"Report saved to {PROJ}/summary_report.html")
print("Done!")
