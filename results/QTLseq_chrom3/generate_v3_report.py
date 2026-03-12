#!/usr/bin/env python3
"""
v3 Multi-Evidence Analysis — Chromosome 3 QTL Peak
Uses PHASED VCF from EHH analysis to resolve haplotype-level introgression.
Key advance over v2: haplotype-level tract mapping (168 haplotypes vs 84 diploids).
Includes SparsePainter chromosome painting results.
"""
import gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.stats import fisher_exact
import base64, io, warnings
warnings.filterwarnings('ignore')

# ===== Paths =====
PHASED_VCF = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/EHH_analysis/chr3_phased.vcf.gz"
UNPHASED_VCF = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/variants/filtered_variants.vcf.gz"
BCFTOOLS = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/BayPass.gwas_maf0.1geno0.2/micromamba/envs/freebayes/bin/bcftools"
SAMPLES_TSV = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/samples.tsv"
NN_SAMPLES = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/EHH_analysis/nn_samples.txt"
DD_SAMPLES = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/EHH_analysis/dd_samples.txt"
OUTDIR = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/QTLseq_chrom3"
SPARSEPAINTER_NPZ = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/QTLseq_chrom3/sparsepainter/sparsepainter_parsed.npz"

REGION_START = 9_500_000
REGION_END   = 13_000_000
EXCLUDE_SAMPLES = {'nr59', 'nr60'}  # Same as v2

# ===== Helper =====
def fig_to_b64(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return b64

# ===== 1. Read sample groups =====
print("Reading sample groups...")
with open(NN_SAMPLES) as f:
    nn_set = {line.strip() for line in f if line.strip()} - EXCLUDE_SAMPLES
with open(DD_SAMPLES) as f:
    dd_set = {line.strip() for line in f if line.strip()} - EXCLUDE_SAMPLES

print(f"  NN: {len(nn_set)} samples, DD: {len(dd_set)} samples")
print(f"  Excluded: {EXCLUDE_SAMPLES}")

# ===== 2. Parse phased VCF =====
print("Parsing phased VCF for chr3:9.5–13.0 Mb...")
positions = []
ref_alleles = []
alt_alleles = []
sample_names = []
# haplotypes: list of dicts, one per variant, sample -> (h1, h2)
hap_data = []

with gzip.open(PHASED_VCF, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            parts = line.strip().split('\t')
            sample_names = parts[9:]
            continue
        parts = line.strip().split('\t')
        pos = int(parts[1])
        if pos < REGION_START or pos > REGION_END:
            continue
        ref = parts[3]
        alt = parts[4]
        if ',' in alt:
            continue  # skip multiallelic
        positions.append(pos)
        ref_alleles.append(ref)
        alt_alleles.append(alt)

        gts = {}
        for si, gt_field in enumerate(parts[9:]):
            sname = sample_names[si]
            if sname in EXCLUDE_SAMPLES:
                continue
            gt_str = gt_field.split(':')[0]
            if '|' in gt_str:
                alleles = gt_str.split('|')
            elif '/' in gt_str:
                alleles = gt_str.split('/')
            else:
                continue
            if '.' in alleles:
                gts[sname] = (np.nan, np.nan)
            else:
                gts[sname] = (int(alleles[0]), int(alleles[1]))
        hap_data.append(gts)

positions = np.array(positions)
n_snps = len(positions)
print(f"  Parsed {n_snps:,} biallelic SNPs in region")

# Build ordered sample lists
nn_samples = sorted([s for s in nn_set if s in sample_names])
dd_samples = sorted([s for s in dd_set if s in sample_names])
all_samples = nn_samples + dd_samples
n_nn = len(nn_samples)
n_dd = len(dd_samples)
n_all = n_nn + n_dd
print(f"  NN: {n_nn}, DD: {n_dd}, Total: {n_all}")

# ===== 3. Build haplotype matrix (n_snps × n_all × 2) =====
print("Building haplotype matrix...")
hap_matrix = np.full((n_snps, n_all, 2), np.nan)
for vi in range(n_snps):
    gts = hap_data[vi]
    for si, sname in enumerate(all_samples):
        if sname in gts:
            h1, h2 = gts[sname]
            hap_matrix[vi, si, 0] = h1
            hap_matrix[vi, si, 1] = h2

# Also build diploid genotype matrix (0, 1, 2, nan)
geno_matrix = np.nansum(hap_matrix, axis=2)  # 0+0=0, 0+1=1, 1+0=1, 1+1=2
# Where either haplotype is NaN, set genotype to NaN
any_nan = np.any(np.isnan(hap_matrix), axis=2)
geno_matrix[any_nan] = np.nan

# Per-SNP missingness
miss_per_snp = np.mean(any_nan, axis=1)
good_snps = miss_per_snp < 0.05
n_good = np.sum(good_snps)
print(f"  SNPs with <5% missingness: {n_good:,} / {n_snps:,}")

# ===== 4. Compute per-SNP statistics =====
print("Computing per-SNP statistics (Δ SNP-index, FST, Fisher's p)...")

nn_geno = geno_matrix[:, :n_nn]
dd_geno = geno_matrix[:, n_nn:]

# Allele frequencies
nn_alt_count = np.nansum(hap_matrix[:, :n_nn, :].reshape(n_snps, -1), axis=1)
nn_total = np.sum(~np.isnan(hap_matrix[:, :n_nn, :].reshape(n_snps, -1)), axis=1)
dd_alt_count = np.nansum(hap_matrix[:, n_nn:, :].reshape(n_snps, -1), axis=1)
dd_total = np.sum(~np.isnan(hap_matrix[:, n_nn:, :].reshape(n_snps, -1)), axis=1)

af_nn = np.where(nn_total > 0, nn_alt_count / nn_total, np.nan)
af_dd = np.where(dd_total > 0, dd_alt_count / dd_total, np.nan)
delta_snp = af_dd - af_nn

# FST (Weir & Cockerham approximation)
af_all = np.where((nn_total + dd_total) > 0,
                  (nn_alt_count + dd_alt_count) / (nn_total + dd_total), np.nan)
p_bar = af_all
q_bar = 1 - p_bar
n_mean = (nn_total + dd_total) / 2
# Between-population variance
s2 = ((nn_total * (af_nn - p_bar)**2) + (dd_total * (af_dd - p_bar)**2)) / n_mean
h_bar = (nn_total * 2 * af_nn * (1 - af_nn) / (nn_total - 1) +
         dd_total * 2 * af_dd * (1 - af_dd) / (dd_total - 1)) / 2
numerator = s2 - (1 / (n_mean - 1)) * (p_bar * q_bar - s2 / 2 - h_bar / 4)
denominator = p_bar * q_bar + s2 / 2
fst = np.where(denominator > 0, numerator / denominator, 0)
fst = np.clip(fst, 0, 1)

# Fisher's exact test
log10p = np.zeros(n_snps)
for i in range(n_snps):
    nn_a = int(nn_alt_count[i])
    nn_r = int(nn_total[i] - nn_alt_count[i])
    dd_a = int(dd_alt_count[i])
    dd_r = int(dd_total[i] - dd_alt_count[i])
    if nn_a + nn_r == 0 or dd_a + dd_r == 0:
        log10p[i] = 0
    else:
        try:
            _, p = fisher_exact([[nn_a, nn_r], [dd_a, dd_r]])
            log10p[i] = -np.log10(max(p, 1e-300))
        except:
            log10p[i] = 0

print(f"  Max -log10(p): {np.max(log10p):.1f}")
print(f"  Max FST: {np.max(fst):.4f}")
print(f"  Max |Δ SNP-index|: {np.max(np.abs(delta_snp)):.4f}")

# ===== 5. Smoothed window statistics =====
print("Computing sliding-window smoothed statistics...")
def sliding_mean(pos, vals, window_kb=25):
    half = window_kb * 500  # half-window in bp
    smoothed = np.zeros_like(vals)
    for i in range(len(pos)):
        mask = (pos >= pos[i] - half) & (pos <= pos[i] + half) & ~np.isnan(vals)
        if np.sum(mask) > 0:
            smoothed[i] = np.mean(vals[mask])
    return smoothed

delta_smooth_25 = sliding_mean(positions, delta_snp, 25)
delta_smooth_50 = sliding_mean(positions, delta_snp, 50)
fst_smooth_25 = sliding_mean(positions, fst, 25)
fst_smooth_50 = sliding_mean(positions, fst, 50)
log10p_smooth_25 = sliding_mean(positions, log10p, 25)
log10p_smooth_50 = sliding_mean(positions, log10p, 50)

# ===== 6. Haplotype-level introgression tract mapping =====
print("Mapping haplotype-level introgression tracts...")

# NN consensus: mode allele at each haplotype position across all NN haplotypes
# We use good SNPs only for tract mapping
nn_haps = hap_matrix[good_snps][:, :n_nn, :]  # (n_good, n_nn, 2)
nn_haps_flat = nn_haps.reshape(np.sum(good_snps), -1)  # (n_good, n_nn*2)
nn_consensus = np.zeros(n_good)
for i in range(n_good):
    vals = nn_haps_flat[i]
    vals = vals[~np.isnan(vals)]
    if len(vals) > 0:
        nn_consensus[i] = np.round(np.mean(vals))  # mode ≈ round(mean) for biallelic

good_pos = positions[good_snps]

# Per-haplotype discordance from NN consensus
# For ALL individuals (NN and DD), each has 2 haplotypes
hap_good = hap_matrix[good_snps]  # (n_good, n_all, 2)

# Smooth discordance with 50kb windows, step 10kb
win_size = 50000
win_step = 10000
win_starts = np.arange(REGION_START, REGION_END - win_size + 1, win_step)
n_wins = len(win_starts)

# For each individual, for each haplotype, compute windowed discordance
hap_disc = np.full((n_wins, n_all, 2), np.nan)
for wi in range(n_wins):
    ws = win_starts[wi]
    we = ws + win_size
    mask = (good_pos >= ws) & (good_pos < we)
    if np.sum(mask) < 3:
        continue
    cons = nn_consensus[mask]
    for si in range(n_all):
        for hi in range(2):
            h = hap_good[mask, si, hi]
            valid = ~np.isnan(h)
            if np.sum(valid) > 0:
                hap_disc[wi, si, hi] = np.mean(h[valid] != cons[valid])

# Introgression: contiguous region with >50% discordance
# For each haplotype, find introgression tracts
win_centers = win_starts + win_size // 2

def find_tracts(disc_vec, threshold=0.5):
    """Find contiguous tracts above threshold."""
    tracts = []
    in_tract = False
    for i in range(len(disc_vec)):
        if not np.isnan(disc_vec[i]) and disc_vec[i] > threshold:
            if not in_tract:
                start_i = i
                in_tract = True
        else:
            if in_tract:
                tracts.append((start_i, i - 1))
                in_tract = False
    if in_tract:
        tracts.append((start_i, len(disc_vec) - 1))
    return tracts

all_tracts = {}  # (sample_idx, hap_idx) -> list of (start_mb, end_mb)
for si in range(n_all):
    for hi in range(2):
        tracts = find_tracts(hap_disc[:, si, hi])
        tract_coords = [(win_centers[t[0]] / 1e6, win_centers[t[1]] / 1e6) for t in tracts]
        if tract_coords:
            all_tracts[(si, hi)] = tract_coords

# ===== 7. Haplotype-level coverage and core regions =====
print("Computing haplotype-level introgression coverage...")

# For each window: what fraction of DD haplotypes are introgressed?
dd_hap_introgressed = np.zeros(n_wins)
dd_hap_total = n_dd * 2
for wi in range(n_wins):
    count = 0
    for si in range(n_nn, n_all):  # DD samples
        for hi in range(2):
            d = hap_disc[wi, si, hi]
            if not np.isnan(d) and d > 0.5:
                count += 1
    dd_hap_introgressed[wi] = count / dd_hap_total

# Same for NN haplotypes (should be ~0)
nn_hap_introgressed = np.zeros(n_wins)
nn_hap_total = n_nn * 2
for wi in range(n_wins):
    count = 0
    for si in range(n_nn):
        for hi in range(2):
            d = hap_disc[wi, si, hi]
            if not np.isnan(d) and d > 0.5:
                count += 1
    nn_hap_introgressed[wi] = count / nn_hap_total

# Core regions: >= 90% DD haplotype coverage
def find_core_regions(coverage, threshold=0.90, min_wins=3):
    regions = []
    in_region = False
    for i in range(len(coverage)):
        if coverage[i] >= threshold:
            if not in_region:
                start_i = i
                in_region = True
        else:
            if in_region:
                if i - start_i >= min_wins:
                    regions.append((start_i, i - 1))
                in_region = False
    if in_region and len(coverage) - start_i >= min_wins:
        regions.append((start_i, len(coverage) - 1))
    return regions

core_regions_90 = find_core_regions(dd_hap_introgressed, 0.90)
core_regions_80 = find_core_regions(dd_hap_introgressed, 0.80)

core_table = []
for start_i, end_i in sorted(core_regions_90, key=lambda x: -(x[1] - x[0])):
    s_mb = win_centers[start_i] / 1e6
    e_mb = win_centers[end_i] / 1e6
    size_kb = (e_mb - s_mb) * 1000
    max_cov = np.max(dd_hap_introgressed[start_i:end_i+1])
    mean_cov = np.mean(dd_hap_introgressed[start_i:end_i+1])
    core_table.append((s_mb, e_mb, size_kb, max_cov, mean_cov))

print(f"  Core regions (≥90% DD haplotype coverage): {len(core_table)}")
for i, (s, e, sz, mx, mn) in enumerate(core_table[:5]):
    print(f"    #{i+1}: {s:.2f}–{e:.2f} Mb ({sz:.0f} kb), max={mx:.1%}, mean={mn:.1%}")

# ===== 8. Breakpoint density =====
print("Computing breakpoint density...")
# For each DD haplotype, where does introgression start/end?
bp_positions = []
for si in range(n_nn, n_all):
    for hi in range(2):
        key = (si, hi)
        if key in all_tracts:
            for s, e in all_tracts[key]:
                bp_positions.append(s)
                bp_positions.append(e)

bp_pos_arr = np.array(bp_positions)
bp_bins = np.arange(REGION_START / 1e6, REGION_END / 1e6 + 0.05, 0.05)
bp_hist, _ = np.histogram(bp_pos_arr, bins=bp_bins)
bp_centers = (bp_bins[:-1] + bp_bins[1:]) / 2

# ===== 9. Genotype composition per window =====
print("Computing genotype composition...")
# For each 50kb window: fraction of each genotype class for NN and DD
geno_comp_nn = np.zeros((n_wins, 4))  # 0/0, 0/1, 1/1, miss
geno_comp_dd = np.zeros((n_wins, 4))

geno_good = geno_matrix[good_snps]

for wi in range(n_wins):
    ws = win_starts[wi]
    we = ws + win_size
    mask = (good_pos >= ws) & (good_pos < we)
    if np.sum(mask) == 0:
        continue
    # NN
    nn_g = geno_good[mask][:, :n_nn]
    total = nn_g.size
    if total > 0:
        geno_comp_nn[wi, 0] = np.sum(nn_g == 0) / total
        geno_comp_nn[wi, 1] = np.sum(nn_g == 1) / total
        geno_comp_nn[wi, 2] = np.sum(nn_g == 2) / total
        geno_comp_nn[wi, 3] = np.sum(np.isnan(nn_g)) / total
    # DD
    dd_g = geno_good[mask][:, n_nn:]
    total = dd_g.size
    if total > 0:
        geno_comp_dd[wi, 0] = np.sum(dd_g == 0) / total
        geno_comp_dd[wi, 1] = np.sum(dd_g == 1) / total
        geno_comp_dd[wi, 2] = np.sum(dd_g == 2) / total
        geno_comp_dd[wi, 3] = np.sum(np.isnan(dd_g)) / total

# ===== 10. Count how many haplotypes per DD individual carry introgression =====
print("Counting introgressed haplotypes per individual...")
# At the core region peak, classify each DD individual: 0, 1, or 2 introgressed haplotypes
if core_table:
    peak_s, peak_e = core_table[0][0], core_table[0][1]
else:
    peak_s, peak_e = 10.5, 12.0

peak_mask = (win_centers / 1e6 >= peak_s) & (win_centers / 1e6 <= peak_e)

dd_hap_status = []  # list of (sample, n_introgressed_haps)
for si_rel in range(n_dd):
    si = n_nn + si_rel
    n_intro = 0
    for hi in range(2):
        disc_vals = hap_disc[peak_mask, si, hi]
        mean_disc = np.nanmean(disc_vals) if np.any(~np.isnan(disc_vals)) else 0
        if mean_disc > 0.5:
            n_intro += 1
    dd_hap_status.append((dd_samples[si_rel], n_intro))

nn_hap_status = []
for si in range(n_nn):
    n_intro = 0
    for hi in range(2):
        disc_vals = hap_disc[peak_mask, si, hi]
        mean_disc = np.nanmean(disc_vals) if np.any(~np.isnan(disc_vals)) else 0
        if mean_disc > 0.5:
            n_intro += 1
    nn_hap_status.append((nn_samples[si], n_intro))

dd_0h = sum(1 for _, n in dd_hap_status if n == 0)
dd_1h = sum(1 for _, n in dd_hap_status if n == 1)
dd_2h = sum(1 for _, n in dd_hap_status if n == 2)
nn_0h = sum(1 for _, n in nn_hap_status if n == 0)
nn_1h = sum(1 for _, n in nn_hap_status if n == 1)
nn_2h = sum(1 for _, n in nn_hap_status if n == 2)
print(f"  DD: 0-hap={dd_0h}, 1-hap={dd_1h}, 2-hap={dd_2h}")
print(f"  NN: 0-hap={nn_0h}, 1-hap={nn_1h}, 2-hap={nn_2h}")

# =========================================================================
# PLOTS
# =========================================================================
print("\n=== Generating plots ===")

# ----- Plot 1: QTLseq multi-statistic (separate Fisher exact panels) -----
print("Plot 1: Multi-statistic sliding windows (separate Fisher exact panels)...")
fig, axes = plt.subplots(4, 1, figsize=(14, 13), sharex=True)
pos_mb = positions / 1e6

ax = axes[0]
ax.scatter(pos_mb, delta_snp, s=1, alpha=0.15, c='grey')
ax.plot(pos_mb, delta_smooth_25, '-', color='#e74c3c', linewidth=1.2, label='25 kb', alpha=0.8)
ax.plot(pos_mb, delta_smooth_50, '-', color='#c0392b', linewidth=2.0, label='50 kb')
ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
ax.set_ylabel('Δ(SNP-index)\n(DD − NN)', fontsize=10)
ax.legend(fontsize=8, loc='upper right')
ax.set_title('Approach 1: Multi-Statistic Sliding Windows (phased data)', fontsize=12, fontweight='bold')

ax = axes[1]
ax.scatter(pos_mb, fst, s=1, alpha=0.15, c='grey')
ax.plot(pos_mb, fst_smooth_25, '-', color='#2ecc71', linewidth=1.2, label='25 kb', alpha=0.8)
ax.plot(pos_mb, fst_smooth_50, '-', color='#27ae60', linewidth=2.0, label='50 kb')
ax.set_ylabel('FST', fontsize=10)
ax.legend(fontsize=8, loc='upper right')

# Separate 25 kb Fisher exact panel
ax = axes[2]
ax.scatter(pos_mb, log10p, s=1, alpha=0.15, c='grey')
ax.plot(pos_mb, log10p_smooth_25, '-', color='#3498db', linewidth=2.0, label='25 kb window')
ax.set_ylabel('-log10(p)\nFisher exact', fontsize=10)
ax.set_title('Fisher Exact Test — 25 kb Sliding Window', fontsize=10, fontweight='bold', color='#2563eb')
ax.legend(fontsize=8, loc='upper right')

# Separate 50 kb Fisher exact panel
ax = axes[3]
ax.scatter(pos_mb, log10p, s=1, alpha=0.15, c='grey')
ax.plot(pos_mb, log10p_smooth_50, '-', color='#2980b9', linewidth=2.0, label='50 kb window')
ax.set_ylabel('-log10(p)\nFisher exact', fontsize=10)
ax.set_xlabel('Position on Chr3 (Mb)', fontsize=10)
ax.set_title('Fisher Exact Test — 50 kb Sliding Window', fontsize=10, fontweight='bold', color='#1e40af')
ax.legend(fontsize=8, loc='upper right')

plt.tight_layout()
qtlseq_b64 = fig_to_b64(fig)

# ----- Plot 2: Multi-evidence overlay (normalised) -----
print("Plot 2: Normalised multi-evidence overlay...")
fig, ax = plt.subplots(figsize=(14, 4.5))

def norm01(x):
    mn, mx = np.nanmin(x), np.nanmax(x)
    return (x - mn) / (mx - mn) if mx > mn else x * 0

ax.plot(pos_mb, norm01(delta_smooth_50), '-', color='#e74c3c', linewidth=2, label='Δ SNP-index (50kb)')
ax.plot(pos_mb, norm01(fst_smooth_50), '-', color='#2ecc71', linewidth=2, label='FST (50kb)')
ax.plot(pos_mb, norm01(log10p_smooth_50), '-', color='#3498db', linewidth=2, label='−log₁₀(p) (50kb)')
ax.fill_between(pos_mb, 0, norm01(log10p_smooth_50), alpha=0.08, color='#3498db')
ax.set_xlabel('Position on Chr3 (Mb)', fontsize=10)
ax.set_ylabel('Normalised score (0–1)', fontsize=10)
ax.set_title('Normalised Multi-Evidence Overlay', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
overlay_b64 = fig_to_b64(fig)

# ----- Plot 3: Haplotype-level introgression tracts -----
print("Plot 3: Haplotype-level introgression tracts...")
fig, ax = plt.subplots(figsize=(14, 10))

# Plot each haplotype as a horizontal line
# DD haplotypes first (at top), then NN
y = 0
yticks_dd = []
yticks_nn = []

# DD individuals
for si_rel in range(n_dd):
    si = n_nn + si_rel
    sname = dd_samples[si_rel]
    for hi in range(2):
        key = (si, hi)
        y += 1
        if hi == 0:
            yticks_dd.append((y + 0.5, sname))
        if key in all_tracts:
            for s, e in all_tracts[key]:
                ax.barh(y, e - s, left=s, height=0.8, color='#e74c3c', edgecolor='none', alpha=0.7)
        # Light background
        ax.barh(y, (REGION_END - REGION_START) / 1e6, left=REGION_START / 1e6,
                height=0.8, color='#fee2e2', edgecolor='none', alpha=0.3, zorder=0)

y_sep = y + 1
ax.axhline(y_sep, color='black', linewidth=1.5)

# NN individuals
for si in range(n_nn):
    sname = nn_samples[si]
    for hi in range(2):
        key = (si, hi)
        y += 1
        if hi == 0:
            yticks_nn.append((y + 0.5, sname))
        if key in all_tracts:
            for s, e in all_tracts[key]:
                ax.barh(y, e - s, left=s, height=0.8, color='#3498db', edgecolor='none', alpha=0.7)
        ax.barh(y, (REGION_END - REGION_START) / 1e6, left=REGION_START / 1e6,
                height=0.8, color='#dbeafe', edgecolor='none', alpha=0.3, zorder=0)

ax.set_xlabel('Position on Chr3 (Mb)', fontsize=10)
ax.set_ylabel('Haplotypes (2 per individual)', fontsize=10)
ax.set_title('Haplotype-Level Introgression Tracts (phased)', fontsize=12, fontweight='bold')
ax.text(REGION_START / 1e6, y_sep + 0.5, 'DD haplotypes ↑', fontsize=9, va='bottom', ha='left',
        color='#e74c3c', fontweight='bold')
ax.text(REGION_START / 1e6, y_sep - 0.5, 'NN haplotypes ↓', fontsize=9, va='top', ha='left',
        color='#3498db', fontweight='bold')
ax.set_xlim(REGION_START / 1e6, REGION_END / 1e6)
ax.set_ylim(0.5, y + 1)
ax.invert_yaxis()
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
tracts_b64 = fig_to_b64(fig)

# ----- Plot 4: DD haplotype coverage + breakpoint density -----
print("Plot 4: Coverage & breakpoint density...")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

wc_mb = win_centers / 1e6
ax1.fill_between(wc_mb, 0, dd_hap_introgressed * 100, color='#e74c3c', alpha=0.3)
ax1.plot(wc_mb, dd_hap_introgressed * 100, '-', color='#e74c3c', linewidth=1.5, label='DD haplotypes')
ax1.fill_between(wc_mb, 0, nn_hap_introgressed * 100, color='#3498db', alpha=0.3)
ax1.plot(wc_mb, nn_hap_introgressed * 100, '-', color='#3498db', linewidth=1.5, label='NN haplotypes')
ax1.axhline(90, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
ax1.set_ylabel('% haplotypes\nintrogressed', fontsize=10)
ax1.legend(fontsize=9)
ax1.set_title('Introgression Coverage (haplotype-level)', fontsize=12, fontweight='bold')
ax1.set_ylim(0, 105)

# Highlight core regions
for s_mb, e_mb, _, _, _ in core_table[:5]:
    ax1.axvspan(s_mb, e_mb, color='#fbbf24', alpha=0.15, zorder=0)

ax2.bar(bp_centers, bp_hist, width=0.05, color='#8b5cf6', alpha=0.7, edgecolor='none')
ax2.set_ylabel('Breakpoint\ncount', fontsize=10)
ax2.set_xlabel('Position on Chr3 (Mb)', fontsize=10)
ax2.set_title('Recombination Breakpoint Density (DD haplotypes)', fontsize=11, fontweight='bold')
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)

plt.tight_layout()
coverage_bp_b64 = fig_to_b64(fig)

# ----- Plot 5: Genotype composition -----
print("Plot 5: Genotype composition...")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 5), sharex=True)

ax1.stackplot(wc_mb, geno_comp_nn[:, 2], geno_comp_nn[:, 1], geno_comp_nn[:, 0],
              colors=['#ef4444', '#fbbf24', '#3b82f6'], alpha=0.7,
              labels=['1/1 (alt homo)', '0/1 (het)', '0/0 (ref homo)'])
ax1.set_ylabel('Fraction', fontsize=10)
ax1.set_title(f'NN Genotype Composition (N={n_nn})', fontsize=11, fontweight='bold')
ax1.legend(fontsize=8, loc='upper right')
ax1.set_ylim(0, 1)

ax2.stackplot(wc_mb, geno_comp_dd[:, 2], geno_comp_dd[:, 1], geno_comp_dd[:, 0],
              colors=['#ef4444', '#fbbf24', '#3b82f6'], alpha=0.7)
ax2.set_ylabel('Fraction', fontsize=10)
ax2.set_xlabel('Position on Chr3 (Mb)', fontsize=10)
ax2.set_title(f'DD Genotype Composition (N={n_dd})', fontsize=11, fontweight='bold')
ax2.set_ylim(0, 1)

plt.tight_layout()
genocomp_b64 = fig_to_b64(fig)

# ----- Plot 6: Haplotype heatmap -----
print("Plot 6: Haplotype heatmap (phased)...")
# Use good SNPs, thin if needed for display
thin_step = max(1, n_good // 800)
thin_idx = np.arange(0, n_good, thin_step)
thin_pos = good_pos[thin_idx]

# Build haplotype display matrix: (n_all*2, n_thin_snps)
hap_display = np.full((n_all * 2, len(thin_idx)), np.nan)
hap_good_data = hap_matrix[good_snps]
for si in range(n_all):
    for hi in range(2):
        row = si * 2 + hi
        hap_display[row, :] = hap_good_data[thin_idx, si, hi]

# Reorder: DD at top, NN at bottom
dd_rows = []
for si_rel in range(n_dd):
    si = n_nn + si_rel
    dd_rows.extend([si * 2, si * 2 + 1])
nn_rows = []
for si in range(n_nn):
    nn_rows.extend([si * 2, si * 2 + 1])

ordered_rows = dd_rows + nn_rows
hap_display_ordered = hap_display[ordered_rows, :]

# Custom colormap: 0=blue(ref), 1=red(alt), nan=grey
cmap_hap = ListedColormap(['#3b82f6', '#ef4444'])

fig, ax = plt.subplots(figsize=(16, 10))
display_data = hap_display_ordered.copy()
display_data[np.isnan(display_data)] = -1

im = ax.imshow(display_data, aspect='auto', cmap=cmap_hap, vmin=0, vmax=1,
               extent=[thin_pos[0]/1e6, thin_pos[-1]/1e6, len(ordered_rows), 0],
               interpolation='nearest')

# Mark boundary between DD and NN
ax.axhline(len(dd_rows), color='white', linewidth=2)
ax.text(thin_pos[0]/1e6, len(dd_rows)/2, 'DD haplotypes', fontsize=10, fontweight='bold',
        color='white', ha='left', va='center',
        bbox=dict(boxstyle='round', facecolor='black', alpha=0.5))
ax.text(thin_pos[0]/1e6, len(dd_rows) + len(nn_rows)/2, 'NN haplotypes', fontsize=10, fontweight='bold',
        color='white', ha='left', va='center',
        bbox=dict(boxstyle='round', facecolor='black', alpha=0.5))

ax.set_xlabel('Position on Chr3 (Mb)', fontsize=11)
ax.set_ylabel(f'Haplotypes ({n_dd}×2 DD + {n_nn}×2 NN = {n_all*2})', fontsize=11)
ax.set_title('Phased Haplotype Heatmap — Chr3: 9.5–13.0 Mb', fontsize=13, fontweight='bold')
plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01, ticks=[0, 1],
             label='Allele').ax.set_yticklabels(['REF (0)', 'ALT (1)'])
plt.tight_layout()
hapheat_b64 = fig_to_b64(fig)

# ----- Plot 7: Haplotype dosage summary (per individual: 0, 1, or 2 introgressed haps) -----
print("Plot 7: Per-individual introgression dosage...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

dd_counts = [dd_0h, dd_1h, dd_2h]
nn_counts = [nn_0h, nn_1h, nn_2h]
labels = ['0 haplotypes', '1 haplotype', '2 haplotypes']
colors = ['#94a3b8', '#fbbf24', '#ef4444']

ax1.bar(labels, dd_counts, color=colors, edgecolor='white', linewidth=1.5)
ax1.set_title(f'DD (N={n_dd})', fontsize=11, fontweight='bold')
ax1.set_ylabel('Number of individuals', fontsize=10)
for i, v in enumerate(dd_counts):
    ax1.text(i, v + 0.3, str(v), ha='center', fontsize=11, fontweight='bold')
ax1.set_ylim(0, max(dd_counts) * 1.2)

ax2.bar(labels, nn_counts, color=colors, edgecolor='white', linewidth=1.5)
ax2.set_title(f'NN (N={n_nn})', fontsize=11, fontweight='bold')
for i, v in enumerate(nn_counts):
    ax2.text(i, v + 0.3, str(v), ha='center', fontsize=11, fontweight='bold')
ax2.set_ylim(0, max(max(nn_counts), 1) * 1.2)

plt.suptitle('Introgression Dosage at Core Region', fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
dosage_b64 = fig_to_b64(fig)

# ----- Plot 8: v2 vs v3 core region comparison -----
print("Plot 8: v2 vs v3 comparison...")
# v2 core regions (from the v2 report, diploid-based)
v2_cores = [
    (12.20, 12.38, 180),
    (9.92, 10.07, 150),
    (10.94, 11.09, 150),
    (11.62, 11.73, 110),
    (11.40, 11.46, 60),
]

fig, ax = plt.subplots(figsize=(14, 4))

# v2 regions
for i, (s, e, _) in enumerate(v2_cores):
    ax.barh(1, e - s, left=s, height=0.35, color='#3b82f6', alpha=0.6,
            edgecolor='#1e40af', linewidth=1)
    if i == 0:
        ax.barh(1, e - s, left=s, height=0.35, color='#3b82f6', alpha=0.6,
                edgecolor='#1e40af', linewidth=1, label='v2 (unphased)')

# v3 regions
for i, (s, e, _, _, _) in enumerate(core_table[:10]):
    ax.barh(0, e - s, left=s, height=0.35, color='#ef4444', alpha=0.6,
            edgecolor='#991b1b', linewidth=1)
    if i == 0:
        ax.barh(0, e - s, left=s, height=0.35, color='#ef4444', alpha=0.6,
                edgecolor='#991b1b', linewidth=1, label='v3 (phased)')

ax.set_yticks([0, 1])
ax.set_yticklabels(['v3 (phased)', 'v2 (unphased)'], fontsize=11, fontweight='bold')
ax.set_xlabel('Position on Chr3 (Mb)', fontsize=11)
ax.set_title('Core Introgression Regions — v2 (unphased) vs v3 (phased)', fontsize=13, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.set_xlim(REGION_START / 1e6, REGION_END / 1e6)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
compare_b64 = fig_to_b64(fig)

# ----- Plot 9: SparsePainter chromosome painting -----
print("Plot 9: SparsePainter chromosome painting...")
sp_data = np.load(SPARSEPAINTER_NPZ)
sp_positions = sp_data['positions']
sp_nn_pdd = sp_data['nn_avg_pdd']
sp_dd_pdd = sp_data['dd_avg_pdd']
sp_pos_mb = sp_positions / 1e6

# Smooth the painting probabilities with 50kb window for cleaner plot
def sp_smooth(pos, vals, window_bp=50000):
    half = window_bp // 2
    smoothed = np.zeros_like(vals)
    for i in range(len(pos)):
        mask = (pos >= pos[i] - half) & (pos <= pos[i] + half) & ~np.isnan(vals)
        if np.sum(mask) > 0:
            smoothed[i] = np.mean(vals[mask])
        else:
            smoothed[i] = np.nan
    return smoothed

# Only smooth the focal region for efficiency
focal_sp = (sp_positions >= REGION_START) & (sp_positions <= REGION_END)
sp_nn_focal = sp_nn_pdd[focal_sp]
sp_dd_focal = sp_dd_pdd[focal_sp]
sp_pos_focal = sp_positions[focal_sp]
sp_pos_focal_mb = sp_pos_focal / 1e6

sp_nn_smooth = sp_smooth(sp_pos_focal, sp_nn_focal, 50000)
sp_dd_smooth = sp_smooth(sp_pos_focal, sp_dd_focal, 50000)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

# Panel 1: DD individuals P(DD ancestry)
ax1.scatter(sp_pos_focal_mb, sp_dd_focal, s=0.5, alpha=0.1, c='grey')
ax1.plot(sp_pos_focal_mb, sp_dd_smooth, '-', color='#e74c3c', linewidth=2.0, label='DD individuals')
ax1.axhline(0.5, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
ax1.set_ylabel('P(DD ancestry)', fontsize=10)
ax1.set_title('SparsePainter Chromosome Painting — DD Individuals', fontsize=11, fontweight='bold')
ax1.set_ylim(0, 1.05)
ax1.legend(fontsize=9, loc='lower right')

# Panel 2: NN individuals P(DD ancestry)
ax2.scatter(sp_pos_focal_mb, sp_nn_focal, s=0.5, alpha=0.1, c='grey')
ax2.plot(sp_pos_focal_mb, sp_nn_smooth, '-', color='#3498db', linewidth=2.0, label='NN individuals')
ax2.axhline(0.5, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
ax2.set_ylabel('P(DD ancestry)', fontsize=10)
ax2.set_title('SparsePainter Chromosome Painting — NN Individuals', fontsize=11, fontweight='bold')
ax2.set_ylim(0, 1.05)
ax2.legend(fontsize=9, loc='upper right')

# Panel 3: Difference (DD - NN)
diff_smooth = sp_dd_smooth - sp_nn_smooth
ax3.fill_between(sp_pos_focal_mb, 0, diff_smooth, where=diff_smooth > 0,
                 color='#e74c3c', alpha=0.3)
ax3.fill_between(sp_pos_focal_mb, 0, diff_smooth, where=diff_smooth <= 0,
                 color='#3498db', alpha=0.3)
ax3.plot(sp_pos_focal_mb, diff_smooth, '-', color='#7c3aed', linewidth=2.0,
         label='P(DD|DD) - P(DD|NN)')
ax3.axhline(0, color='black', linestyle='--', linewidth=0.5)
ax3.set_ylabel('Δ P(DD ancestry)', fontsize=10)
ax3.set_xlabel('Position on Chr3 (Mb)', fontsize=10)
ax3.set_title('Differential Painting (DD minus NN)', fontsize=11, fontweight='bold')
ax3.legend(fontsize=9, loc='upper right')

plt.tight_layout()
sparsepainter_b64 = fig_to_b64(fig)

# SparsePainter summary stats for the HTML
sp_peak = (sp_pos_focal >= 10500000) & (sp_pos_focal <= 12000000)
sp_outside = ~((sp_positions >= 9500000) & (sp_positions <= 13000000))
sp_dd_peak_avg = np.nanmean(sp_dd_focal[sp_peak])
sp_nn_peak_avg = np.nanmean(sp_nn_focal[sp_peak])
sp_dd_outside_avg = np.nanmean(sp_dd_pdd[sp_outside])
sp_nn_outside_avg = np.nanmean(sp_nn_pdd[sp_outside])
print(f"  DD at peak: P(DD)={sp_dd_peak_avg:.4f}")
print(f"  NN at peak: P(DD)={sp_nn_peak_avg:.4f}")
print(f"  DD outside: P(DD)={sp_dd_outside_avg:.4f}")
print(f"  NN outside: P(DD)={sp_nn_outside_avg:.4f}")

# ===== Top 30 SNPs by -log10(p) =====
top30_idx = np.argsort(log10p)[::-1][:30]
top30_rows = ""
for rank, idx in enumerate(top30_idx, 1):
    top30_rows += f"<tr><td>{rank}</td><td>{positions[idx]:,}</td><td>{ref_alleles[idx]}</td>"
    top30_rows += f"<td>{alt_alleles[idx]}</td><td>{af_nn[idx]:.3f}</td><td>{af_dd[idx]:.3f}</td>"
    top30_rows += f"<td>{delta_snp[idx]:+.3f}</td><td>{fst[idx]:.4f}</td>"
    top30_rows += f"<td><strong>{log10p[idx]:.1f}</strong></td></tr>\n"

# =========================================================================
# HTML REPORT
# =========================================================================
print("\nWriting HTML report...")

# Core regions table rows
core_rows = ""
for i, (s, e, sz, mx, mn) in enumerate(core_table[:10]):
    style = ' style="background:#eff6ff;font-weight:600;"' if i == 0 else ''
    core_rows += f'<tr{style}><td>{i+1}</td><td>{s:.2f} – {e:.2f} Mb</td>'
    core_rows += f'<td>{sz:.0f} kb</td><td>{mx:.1%}</td><td>{mn:.1%}</td></tr>\n'

# v2 comparison table
v2_rows = ""
for i, (s, e, sz) in enumerate(v2_cores):
    v2_rows += f'<tr><td>{i+1}</td><td>{s:.2f} – {e:.2f} Mb</td><td>{sz:.0f} kb</td></tr>\n'
v3_rows = ""
for i, (s, e, sz, mx, mn) in enumerate(core_table[:5]):
    v3_rows += f'<tr><td>{i+1}</td><td>{s:.2f} – {e:.2f} Mb</td><td>{sz:.0f} kb</td><td>{mn:.1%}</td></tr>\n'

html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Chr3 QTL Peak — Multi-Evidence Analysis (v3, phased)</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:-apple-system,BlinkMacSystemFont,sans-serif; background:#f8fafc;
       color:#1e293b; line-height:1.6; max-width:1300px; margin:0 auto; padding:2rem; }}
h1 {{ font-size:1.6rem; margin-bottom:0.3rem; }}
.subtitle {{ color:#64748b; margin-bottom:2rem; font-size:0.9rem; }}
h2 {{ font-size:1.15rem; margin:2rem 0 0.8rem; border-bottom:2px solid #e2e8f0; padding-bottom:0.3rem; }}
h3 {{ font-size:1rem; margin:1.2rem 0 0.5rem; color:#334155; }}
.card {{ background:#fff; border:1px solid #e2e8f0; border-radius:8px;
         padding:1.2rem; margin-bottom:1rem; box-shadow:0 1px 2px rgba(0,0,0,0.04); }}
table {{ width:100%; border-collapse:collapse; font-size:0.82rem; }}
th {{ background:#f1f5f9; text-align:left; padding:0.5rem 0.7rem; font-weight:600;
      border-bottom:2px solid #e2e8f0; }}
td {{ padding:0.4rem 0.7rem; border-bottom:1px solid #f1f5f9; }}
.highlight {{ color:#2563eb; font-weight:600; }}
.plot-container {{ text-align:center; margin:1rem 0; }}
.plot-container img {{ max-width:100%; border:1px solid #e2e8f0; border-radius:6px; }}
.grid {{ display:grid; grid-template-columns:1fr 1fr; gap:1rem; }}
.grid .card {{ margin-bottom:0; }}
code {{ background:#f1f5f9; padding:0.15rem 0.4rem; border-radius:3px; font-size:0.82rem; }}
pre {{ background:#1e293b; color:#e2e8f0; padding:1rem; border-radius:6px; overflow-x:auto;
       font-size:0.78rem; line-height:1.5; margin:0.5rem 0; }}
.approach {{ border-left:4px solid #3b82f6; padding:0.8rem 1rem; margin:1rem 0;
             background:#f8fafc; }}
.approach-title {{ font-weight:700; color:#1e40af; margin-bottom:0.3rem; }}
.core-box {{ background:linear-gradient(135deg, #eff6ff, #dbeafe); border:2px solid #3b82f6;
             border-radius:10px; padding:1.5rem; margin:1.5rem 0; text-align:center; }}
.core-box h3 {{ color:#1e40af; font-size:1.2rem; margin-bottom:0.5rem; }}
.core-box .coords {{ font-size:1.8rem; font-weight:800; color:#1e3a8a; }}
.core-box .size {{ font-size:1rem; color:#3b82f6; margin-top:0.3rem; }}
.new-badge {{ display:inline-block; background:#10b981; color:white; padding:0.1rem 0.5rem;
              border-radius:4px; font-size:0.75rem; font-weight:700; margin-left:0.5rem; vertical-align:middle; }}
.phased-note {{ background:#fef3c7; border:2px solid #f59e0b; border-radius:10px; padding:1.2rem;
                margin-bottom:1.5rem; }}
</style></head><body>

<h1>Chromosome 3 QTL Peak &mdash; Multi-Evidence Analysis (v3)</h1>
<div class="subtitle"><em>Pararge aegeria</em> &mdash; Haplotype-level introgression mapping
using phased genotypes</div>

<div class="phased-note">
<h3 style="color:#92400e;margin-bottom:0.5rem;">&#x1F9EC; v3: Phased Haplotype Analysis</h3>
<p style="font-size:0.9rem;">This analysis uses the <strong>phased VCF</strong> from the
EHH Haplotype-Based Selection Analysis (SHAPEIT4 phasing). Unlike v2 (unphased diploid genotypes),
v3 resolves each individual into <strong>two separate haplotypes</strong>, enabling:</p>
<ul style="font-size:0.85rem;margin:0.5rem 0 0 1.5rem;">
<li>Haplotype-level introgression tract mapping ({n_all * 2} haplotypes vs {n_all} diploids)</li>
<li>Classification of DD individuals into 0, 1, or 2 introgressed haplotypes</li>
<li>More precise breakpoint localisation at haplotype boundaries</li>
<li>Detection of introgressed haplotypes leaking into NN population</li>
</ul>
<p style="font-size:0.85rem;margin-top:0.5rem;"><strong>Input:</strong> {n_snps:,} phased biallelic
SNPs in chr3: {REGION_START/1e6:.1f}–{REGION_END/1e6:.1f} Mb
(vs 25,679 unphased in v2).</p>
</div>

<div style="background:#fef2f2;border:2px solid #ef4444;border-radius:10px;padding:1.2rem;
            margin-bottom:1.5rem;">
<h3 style="color:#b91c1c;margin-bottom:0.5rem;">Samples Excluded</h3>
<p style="font-size:0.9rem;"><strong>2 sample(s) removed:</strong>
<code>nr59</code>, <code>nr60</code> (same as v2; high missingness in focal region)</p>
</div>

<h2>Summary</h2>
<div class="grid">
<div class="card"><table>
    <tr><td>NN samples / haplotypes</td><td class="highlight">{n_nn} / {n_nn*2}</td></tr>
    <tr><td>DD samples / haplotypes</td><td class="highlight">{n_dd} / {n_dd*2}</td></tr>
    <tr><td>Region analysed</td><td class="highlight">Chr3: 9.5–13.0 Mb</td></tr>
    <tr><td>Phased SNPs in region</td><td class="highlight">{n_snps:,}</td></tr>
    <tr><td>SNPs &lt;5% missing</td><td class="highlight">{n_good:,}</td></tr>
</table></div>
<div class="card"><table>
    <tr><td>Max −log₁₀(p)</td><td class="highlight">{np.max(log10p):.1f}</td></tr>
    <tr><td>SNPs with −log₁₀(p) > 40</td><td class="highlight">{np.sum(log10p > 40):,}</td></tr>
    <tr><td>SNPs with −log₁₀(p) > 50</td><td class="highlight">{np.sum(log10p > 50):,}</td></tr>
    <tr><td>DD: 2 introgressed haps</td><td class="highlight">{dd_2h} / {n_dd}</td></tr>
    <tr><td>DD: 1 introgressed hap</td><td class="highlight">{dd_1h} / {n_dd}</td></tr>
</table></div>
</div>

<div class="core-box">
<h3>Core Introgression Regions (≥90% DD haplotype coverage)</h3>
<div class="coords">#{1 if core_table else '?'}: {core_table[0][0]:.2f} – {core_table[0][1]:.2f} Mb</div>
<div class="size">{len(core_table)} region(s) at ≥90% DD haplotype introgression coverage</div>
</div>
<div class="card">
<table>
<thead><tr><th>Rank</th><th>Coordinates</th><th>Size</th>
<th>Max DD hap coverage</th><th>Mean DD hap coverage</th></tr></thead>
<tbody>{core_rows}</tbody>
</table>
</div>

<h2>Approach 1: Multi-Statistic Sliding Windows</h2>
<div class="approach">
<div class="approach-title">Same statistics, phased allele counts</div>
<p style="font-size:0.85rem;">&Delta;(SNP-index), F<sub>ST</sub>, and Fisher's exact test now computed from
<strong>phased haplotype allele counts</strong> (2N alleles per group). Smoothed with 25 kb and 50 kb windows.
Phasing does not change allele frequencies or population-level statistics, but ensures consistent
allele counting across individuals. The Fisher exact test &minus;log<sub>10</sub>(p) is shown in
<strong>separate panels</strong> for the 25 kb and 50 kb windows to allow detailed comparison of
the peak structure at different smoothing scales.</p>
</div>
<div class="plot-container"><img src="data:image/png;base64,{qtlseq_b64}" alt="QTL-seq"></div>

<h2>Normalised Multi-Evidence Overlay</h2>
<div class="plot-container"><img src="data:image/png;base64,{overlay_b64}" alt="Overlay"></div>

<h2>Approach 2: Haplotype-Level Introgression Tracts <span class="new-badge">NEW in v3</span></h2>
<div class="approach">
<div class="approach-title">Key advance over v2: individual haplotype resolution</div>
<p style="font-size:0.85rem;">Each individual contributes <strong>two haplotypes</strong> (2&times;{n_all} = {n_all*2} total).
The NN consensus allele was determined per SNP from all {n_nn*2} NN haplotypes. Per-haplotype discordance
was computed in 50 kb sliding windows (10 kb step). Introgression = contiguous run of &gt;50% discordance.
This resolves heterozygous individuals into one introgressed + one non-introgressed haplotype.</p>
</div>

<h3>Reading the Tract Plot</h3>
<div class="card">
<p>Each individual is shown as <strong>two horizontal rows</strong> (one per phased haplotype).
The plot is divided by a black line:</p>
<ul style="font-size:0.85rem;margin:0.5rem 0 0 1.2rem;line-height:1.8;">
<li><strong style="color:#e74c3c;">Red blocks (top half)</strong> = DD haplotypes that are
<em>discordant from the NN consensus</em> &mdash; these are the <strong>introgressed tracts</strong>.
A red block means that haplotype carries alleles that differ from the typical NN background,
indicating it derives from a different ancestral source (the introgressed lineage).</li>
<li><strong style="color:#3498db;">Blue blocks (bottom half)</strong> = NN haplotypes that show
discordance from their own consensus. In a clean system these should be mostly absent.
Any blue blocks would indicate NN individuals carrying fragments of the introgressed haplotype.</li>
<li><strong>Light pink/blue backgrounds</strong> = the full extent of the analysed region. Where no
coloured block is drawn, the haplotype matches the NN consensus (i.e., it is <em>not</em> introgressed).</li>
<li><strong>Width of each block</strong> = the physical extent of the introgression tract in Mb.
Wider blocks indicate longer runs of introgressed DNA. The edges of each block are
<strong>recombination breakpoints</strong> where the ancestry switches from introgressed to non-introgressed.</li>
</ul>
</div>

<h3>Step-by-Step Algorithm</h3>
<div class="card">
<p style="font-size:0.85rem;">The algorithm is a simplified form of <strong>local ancestry inference</strong>
that leverages the two-population design (NN = reference, DD = target):</p>

<div style="font-size:0.85rem;margin:0.8rem 0;">
<p><strong>Step 1 &mdash; Build the NN consensus haplotype (reference).</strong><br>
At each SNP position, the <em>modal allele</em> across all {n_nn*2} NN haplotypes is determined.
For biallelic SNPs, this is equivalent to rounding the mean allele value: if &gt;50% of NN haplotypes
carry the ALT allele, the consensus is ALT; otherwise REF. This produces a single &ldquo;typical NN&rdquo;
haplotype vector of length {n_good:,} SNPs.</p>

<p><strong>Step 2 &mdash; Compute per-haplotype discordance in sliding windows.</strong><br>
For every individual (NN and DD alike), each of their two phased haplotypes is compared to the NN
consensus. Within a 50 kb sliding window (stepped every 10 kb), the <em>discordance rate</em> is
calculated as the fraction of non-missing SNPs where the haplotype allele differs from the consensus.
A discordance of 0.0 means perfect match to the NN consensus; a discordance of 1.0 means every
allele differs.</p>

<p><strong>Step 3 &mdash; Call introgression tracts.</strong><br>
A haplotype is classified as <em>introgressed</em> in a given window if discordance &gt; 0.50
(i.e., more alleles differ from the NN consensus than match it). Contiguous runs of introgressed
windows are merged into <strong>introgression tracts</strong>, each defined by a start and end
position in Mb.</p>

<p><strong>Step 4 &mdash; Compute coverage and define core regions.</strong><br>
At each window, the fraction of DD haplotypes that are introgressed is tallied. A <strong>core
introgression region</strong> is a contiguous stretch where &ge;90% of DD haplotypes are introgressed.
These are the regions most likely to contain the causal variant(s).</p>
</div>
</div>

<h3>Code Used</h3>
<div class="card">
<p style="font-size:0.85rem;">The following Python code implements the haplotype-level tract mapping.
The full script is <code>generate_v3_report.py</code>.</p>
<pre>
<span style="color:#6a9955;"># === Step 1: NN consensus haplotype ===</span>
<span style="color:#6a9955;"># For each SNP, determine the modal allele across all NN haplotypes</span>
nn_haps = hap_matrix[good_snps][:, :n_nn, :]      <span style="color:#6a9955;"># shape: (n_snps, n_nn, 2)</span>
nn_haps_flat = nn_haps.reshape(n_good, -1)          <span style="color:#6a9955;"># flatten to (n_snps, n_nn*2)</span>
nn_consensus = np.zeros(n_good)
<span style="color:#569cd6;">for</span> i <span style="color:#569cd6;">in</span> range(n_good):
    vals = nn_haps_flat[i]
    vals = vals[~np.isnan(vals)]
    <span style="color:#569cd6;">if</span> len(vals) > 0:
        nn_consensus[i] = np.round(np.mean(vals))   <span style="color:#6a9955;"># mode ~ round(mean) for biallelic</span>

<span style="color:#6a9955;"># === Step 2: Sliding-window discordance per haplotype ===</span>
win_size = 50000                                     <span style="color:#6a9955;"># 50 kb windows</span>
win_step = 10000                                     <span style="color:#6a9955;"># 10 kb step</span>
win_starts = np.arange(REGION_START, REGION_END - win_size + 1, win_step)

hap_disc = np.full((n_wins, n_all, 2), np.nan)      <span style="color:#6a9955;"># discordance matrix</span>
<span style="color:#569cd6;">for</span> wi <span style="color:#569cd6;">in</span> range(n_wins):
    ws, we = win_starts[wi], win_starts[wi] + win_size
    mask = (good_pos >= ws) & (good_pos < we)
    <span style="color:#569cd6;">if</span> np.sum(mask) < 3:
        <span style="color:#569cd6;">continue</span>
    cons = nn_consensus[mask]                        <span style="color:#6a9955;"># consensus alleles in window</span>
    <span style="color:#569cd6;">for</span> si <span style="color:#569cd6;">in</span> range(n_all):
        <span style="color:#569cd6;">for</span> hi <span style="color:#569cd6;">in</span> range(2):                       <span style="color:#6a9955;"># each haplotype</span>
            h = hap_good[mask, si, hi]
            valid = ~np.isnan(h)
            <span style="color:#569cd6;">if</span> np.sum(valid) > 0:
                hap_disc[wi, si, hi] = np.mean(h[valid] != cons[valid])

<span style="color:#6a9955;"># === Step 3: Call introgression tracts ===</span>
<span style="color:#569cd6;">def</span> <span style="color:#dcdcaa;">find_tracts</span>(disc_vec, threshold=0.5):
    <span style="color:#ce9178;"># Find contiguous runs above threshold.</span>
    tracts = []
    in_tract = False
    <span style="color:#569cd6;">for</span> i <span style="color:#569cd6;">in</span> range(len(disc_vec)):
        <span style="color:#569cd6;">if not</span> np.isnan(disc_vec[i]) <span style="color:#569cd6;">and</span> disc_vec[i] > threshold:
            <span style="color:#569cd6;">if not</span> in_tract:
                start_i = i
                in_tract = True
        <span style="color:#569cd6;">else</span>:
            <span style="color:#569cd6;">if</span> in_tract:
                tracts.append((start_i, i - 1))
                in_tract = False
    <span style="color:#569cd6;">if</span> in_tract:
        tracts.append((start_i, len(disc_vec) - 1))
    <span style="color:#569cd6;">return</span> tracts

<span style="color:#6a9955;"># === Step 4: Coverage = fraction of DD haplotypes introgressed per window ===</span>
dd_hap_introgressed = np.zeros(n_wins)
<span style="color:#569cd6;">for</span> wi <span style="color:#569cd6;">in</span> range(n_wins):
    count = 0
    <span style="color:#569cd6;">for</span> si <span style="color:#569cd6;">in</span> range(n_nn, n_all):             <span style="color:#6a9955;"># DD samples only</span>
        <span style="color:#569cd6;">for</span> hi <span style="color:#569cd6;">in</span> range(2):
            d = hap_disc[wi, si, hi]
            <span style="color:#569cd6;">if not</span> np.isnan(d) <span style="color:#569cd6;">and</span> d > 0.5:
                count += 1
    dd_hap_introgressed[wi] = count / (n_dd * 2)    <span style="color:#6a9955;"># fraction of DD haplotypes</span>
</pre>
</div>

<h3>Relationship to Published Methods</h3>
<div class="card">
<p style="font-size:0.85rem;">This approach is a simplified, two-population form of
<strong>local ancestry inference (LAI)</strong>. Several established methods use conceptually
similar strategies:</p>

<table style="font-size:0.82rem;">
<thead><tr><th>Method</th><th>Approach</th><th>Key similarity</th></tr></thead>
<tbody>
<tr><td><strong>ChromoPainter</strong><br><span style="color:#64748b;">Lawson et al. 2012, <em>PLOS Genetics</em></span></td>
<td>HMM-based chromosome painting; reconstructs each &ldquo;recipient&rdquo; haplotype as a mosaic
of &ldquo;donor&rdquo; haplotype chunks using the Li &amp; Stephens (2003) model.</td>
<td>Both methods compare individual haplotypes against a reference population to identify
segments of different ancestry. ChromoPainter is more sophisticated (HMM, multiple donors),
while our approach uses a single consensus reference and sliding-window discordance.</td></tr>
<tr><td><strong>G<sub>min</sub></strong><br><span style="color:#64748b;">Geneva et al. 2015, <em>PLOS ONE</em></span></td>
<td>Ratio of minimum to mean between-population haplotype distance, computed in sliding
windows. Designed to detect recent introgression tracts.</td>
<td>Both use sliding-window haplotype comparison between two populations. G<sub>min</sub>
focuses on population-level summary statistics rather than per-individual tracts.</td></tr>
<tr><td><strong>IntroUNET</strong><br><span style="color:#64748b;">Durvasula &amp; Sankararaman 2024, <em>PLOS Genetics</em></span></td>
<td>Deep learning (U-Net CNN) that takes a binary genotype tensor as input and predicts
which alleles were introgressed, framing the problem as image segmentation.</td>
<td>Both identify per-individual introgressed segments along the chromosome. IntroUNET uses
neural networks trained on simulated data; our approach uses direct empirical comparison
to an observed reference population.</td></tr>
<tr><td><strong>HAPMIX / LAMP-LD</strong><br><span style="color:#64748b;">Price et al. 2009; Baran et al. 2012</span></td>
<td>HMM-based local ancestry inference that assigns each segment of an admixed individual&rsquo;s
genome to one of two (or more) ancestral populations.</td>
<td>Our approach is functionally equivalent to a two-population LAI with a hard threshold
instead of a probabilistic model. The NN consensus serves as the &ldquo;non-introgressed
reference&rdquo; and high discordance identifies segments from the alternative ancestry.</td></tr>
<tr><td><strong>SparsePainter</strong><br><span style="color:#64748b;">Zhu et al. 2025, <em>Nature Communications</em></span></td>
<td>Fast haplotype-based local ancestry inference using sparse reference panels.
Direct improvement on ChromoPainter with orders-of-magnitude speed gain.</td>
<td>Both use phased haplotypes and a reference panel. SparsePainter is designed for complex
multi-population scenarios; our simplified approach is tailored for the clear
two-population (NN vs DD) design in this experiment.</td></tr>
</tbody></table>

<p style="font-size:0.85rem;margin-top:1rem;">
<strong>Introgression in Lepidoptera:</strong> In the <em>Heliconius</em> butterfly radiation,
Edelman et al. (2019, <em>Science</em>) found that introgression has been a major contributor
of genealogical discordance, with 38% of loci genome-wide showing a history of introgression.
Martin et al. (2019, <em>PLOS Biology</em>) showed that rates of introgression are predicted
by recombination rate variation across butterfly genomes.
For <em>Pararge aegeria</em> specifically, Pruisscher et al. (2018, <em>Molecular Ecology</em>)
identified candidate regions for diapause induction variation including the genes <em>period</em>
and <em>timeless</em>, and recent work has implicated Z-linked inheritance and repeated
colonisation in shaping diapause adaptation across Scandinavian populations.</p>

<p style="font-size:0.85rem;margin-top:0.5rem;">
<strong>Why this simplified approach works here:</strong> In contrast to multi-population admixture
scenarios (humans, <em>Heliconius</em>), our experimental design provides a clean two-population
system where NN (non-diapause) serves as the unambiguous reference and DD (diapause) is the
target population potentially carrying introgressed haplotypes. The crossing design means
the ancestral population structure is known, making a simple consensus-discordance approach
sufficient without requiring the full machinery of HMM-based LAI or simulation-trained neural
networks. The 50% discordance threshold is conservative and biologically motivated: an
introgressed haplotype should differ from the local population consensus at the majority of
linked sites.</p>
</div>
<div class="plot-container"><img src="data:image/png;base64,{tracts_b64}" alt="Haplotype tracts"></div>

<h2>Introgression Dosage <span class="new-badge">NEW in v3</span></h2>
<div class="approach">
<div class="approach-title">How many haplotypes carry the introgression per individual?</div>
<p style="font-size:0.85rem;">At the top core region ({core_table[0][0]:.2f}–{core_table[0][1]:.2f} Mb),
each individual is classified by how many of their two haplotypes are introgressed.
<strong>DD:</strong> {dd_2h} homozygous (2 haps), {dd_1h} heterozygous (1 hap), {dd_0h} non-introgressed (0 haps).
<strong>NN:</strong> {nn_2h} homozygous, {nn_1h} heterozygous, {nn_0h} non-introgressed.</p>
</div>
<div class="plot-container"><img src="data:image/png;base64,{dosage_b64}" alt="Dosage"></div>

<h2>Approach 3: Introgression Coverage &amp; Breakpoint Density</h2>
<div class="approach">
<div class="approach-title">Haplotype-level coverage (vs diploid in v2)</div>
<p style="font-size:0.85rem;">Coverage = fraction of DD haplotypes (out of {n_dd*2}) that are introgressed
at each position. Core regions = contiguous segments with ≥90% coverage. Breakpoint density
shows where recombination has eroded the introgression boundaries.</p>
</div>
<div class="plot-container"><img src="data:image/png;base64,{coverage_bp_b64}" alt="Coverage"></div>

<h2>Approach 4: Genotype Composition</h2>
<div class="approach">
<div class="approach-title">Per-window diploid genotype class frequencies</div>
<p style="font-size:0.85rem;">For each 50 kb sliding window (10 kb step), this plot shows the
fraction of SNP&times;individual observations that fall into each genotype class:
<strong style="color:#3b82f6;">0/0 (REF homozygous)</strong>,
<strong style="color:#fbbf24;">0/1 (heterozygous)</strong>, and
<strong style="color:#ef4444;">1/1 (ALT homozygous)</strong>.
The fractions are computed by pooling all SNPs within the window across all individuals in each
group (NN or DD), so the denominator is
<code>n_snps_in_window &times; n_individuals</code>.</p>
</div>

<h3>Step-by-Step Algorithm</h3>
<div class="card">
<p style="font-size:0.85rem;">For each 50 kb window at step <em>i</em>:</p>
<ol style="font-size:0.85rem;margin:0.5rem 0 0 1.5rem;line-height:1.8;">
<li>Select all SNPs with &lt;5% missingness that fall within the window boundaries.</li>
<li>Extract the diploid genotype sub-matrix for NN samples: shape
<code>(n_snps_in_window, n_nn)</code>. Each cell is 0, 1, or 2 (sum of two haplotype alleles).</li>
<li>Count the fraction of cells that equal 0 (REF homo), 1 (het), or 2 (ALT homo).
The remaining fraction is missing data.</li>
<li>Repeat for DD samples.</li>
</ol>
<pre>
<span style="color:#6a9955;"># Genotype composition per window</span>
<span style="color:#569cd6;">for</span> wi <span style="color:#569cd6;">in</span> range(n_wins):
    ws, we = win_starts[wi], win_starts[wi] + win_size
    mask = (good_pos >= ws) & (good_pos < we)
    <span style="color:#6a9955;"># NN: pool all SNPs x individuals in this window</span>
    nn_g = geno_good[mask][:, :n_nn]        <span style="color:#6a9955;"># shape: (n_snps_in_win, n_nn)</span>
    total = nn_g.size                         <span style="color:#6a9955;"># denominator = n_snps * n_nn</span>
    geno_comp_nn[wi, 0] = np.sum(nn_g == 0) / total   <span style="color:#6a9955;"># frac REF homo</span>
    geno_comp_nn[wi, 1] = np.sum(nn_g == 1) / total   <span style="color:#6a9955;"># frac het</span>
    geno_comp_nn[wi, 2] = np.sum(nn_g == 2) / total   <span style="color:#6a9955;"># frac ALT homo</span>
</pre>
</div>

<h3>Why Genotype Composition Differs from the Haplotype Heatmap</h3>
<div class="card">
<p style="font-size:0.85rem;">At first glance, the genotype composition stacked plot may appear
discordant with the haplotype heatmap. This is because they measure <strong>fundamentally different
things</strong>:</p>
<ul style="font-size:0.85rem;margin:0.5rem 0 0 1.2rem;line-height:1.8;">
<li><strong>Haplotype heatmap</strong> shows the <em>raw allele identity</em> at each SNP (REF=blue, ALT=red)
for each haplotype. In the introgression region, DD haplotypes are red (carry ALT alleles)
because the introgressed lineage is enriched for alternative alleles relative to the reference genome.
The heatmap is <em>allele-agnostic</em> about which allele is &ldquo;introgressed&rdquo; &mdash;
it simply shows REF vs ALT.</li>
<li><strong>Genotype composition</strong> shows the <em>pooled frequency of genotype classes</em> across
many SNPs. Even in a highly differentiated region, a mix of SNPs will be ALT-fixed in DD
(showing 1/1) while others may be REF-fixed in NN (showing 0/0 in NN). The genotype composition
shows the <strong>aggregate pattern</strong>, not per-SNP allele identity.</li>
<li><strong>Key insight:</strong> The genotype composition does not track &ldquo;which allele is the
introgressed allele&rdquo; at each SNP. Instead, it pools REF/ALT identity across all SNPs equally.
Since the reference genome is neither purely NN nor purely DD, some SNPs have the &ldquo;introgressed&rdquo;
allele as REF and others as ALT. This dilutes the genotype composition signal relative to the
haplotype heatmap, which shows the raw pattern clearly.</li>
<li><strong>The tract-based analysis (Approach 2) and the haplotype heatmap are more informative</strong>
for identifying introgression, because they operate at the haplotype level and compare against the
NN consensus rather than against the reference genome.</li>
</ul>
</div>
<div class="plot-container"><img src="data:image/png;base64,{genocomp_b64}" alt="Genotype composition"></div>

<h2>Approach 5: SparsePainter Chromosome Painting <span class="new-badge">NEW</span></h2>
<div class="approach">
<div class="approach-title">Probabilistic local ancestry inference (SparsePainter v1.3.2)</div>
<p style="font-size:0.85rem;"><strong>SparsePainter</strong>
(Yang et al. 2025, <em>Nature Communications</em> 16: 2742) uses the PBWT algorithm to find
sparse haplotype matches and an HMM to assign local ancestry probabilities at each SNP position.
All {n_all} individuals (NN + DD) were painted against a two-population reference panel (NN and DD)
using leave-one-out to avoid self-painting bias. For each individual at each SNP, the tool outputs
<strong>P(NN ancestry)</strong> and <strong>P(DD ancestry)</strong>.</p>
</div>

<h3>Method</h3>
<div class="card">
<p style="font-size:0.85rem;">The SparsePainter analysis was run on the full phased chr3 VCF
(121,472 SNPs, 82 individuals) with the following configuration:</p>
<ul style="font-size:0.85rem;margin:0.5rem 0 0 1.2rem;line-height:1.8;">
<li><strong>Reference file:</strong> Combined phased VCF with all 82 samples (40 NN + 42 DD)</li>
<li><strong>Target file:</strong> Same VCF (self-painting with automatic leave-one-out)</li>
<li><strong>Population labels:</strong> Each sample labeled as &ldquo;NN&rdquo; or &ldquo;DD&rdquo;</li>
<li><strong>Genetic map:</strong> Approximate 1 cM/Mb linear map</li>
<li><strong>Leave-one-out:</strong> Automatic when reffile = targetfile (each individual excluded
from its own population when being painted)</li>
<li><strong>Recombination scaling:</strong> Fixed lambda = 8,822.55 (estimated from data)</li>
</ul>
<pre>
SparsePainter \\
    -reffile all_combined.vcf.gz \\
    -targetfile all_combined.vcf.gz \\
    -popfile popfile_2pop.txt \\      <span style="color:#6a9955;"># 40 NN + 42 DD labels</span>
    -mapfile mapfile.txt \\
    -namefile namefile_all.txt \\
    -out twopop_painting \\
    -prob -aveSNP -aveind -ncores 1
</pre>
</div>

<h3>Results</h3>
<div class="card">
<table style="font-size:0.85rem;">
<thead><tr><th>Group</th><th>Region</th><th>Mean P(DD ancestry)</th></tr></thead>
<tbody>
<tr><td>DD individuals</td><td>Outside focal (rest of chr3)</td><td>{sp_dd_outside_avg:.4f}</td></tr>
<tr><td>DD individuals</td><td>Peak (10.5&ndash;12.0 Mb)</td><td><strong>{sp_dd_peak_avg:.4f}</strong></td></tr>
<tr><td>NN individuals</td><td>Outside focal (rest of chr3)</td><td>{sp_nn_outside_avg:.4f}</td></tr>
<tr><td>NN individuals</td><td>Peak (10.5&ndash;12.0 Mb)</td><td><strong>{sp_nn_peak_avg:.4f}</strong></td></tr>
</tbody>
</table>
<p style="font-size:0.85rem;margin-top:0.8rem;">
<strong>Interpretation:</strong> The SparsePainter results confirm the introgression signal.
At the peak region (10.5&ndash;12.0 Mb), DD individuals show {sp_dd_peak_avg:.1%} DD ancestry
while NN individuals show only {sp_nn_peak_avg:.1%} DD ancestry &mdash; a clean separation.
Outside the focal region, there is still genome-wide differentiation (DD show ~{sp_dd_outside_avg:.0%}
DD ancestry background, reflecting genome-wide population structure),
but the focal region shows an <strong>enhanced signal</strong> above this baseline.</p>
<p style="font-size:0.85rem;margin-top:0.5rem;">
The three panels below show: (top) P(DD ancestry) for DD individuals along the focal region,
(middle) P(DD ancestry) for NN individuals, and (bottom) the difference
&Delta;P = P(DD|DD) &minus; P(DD|NN). The difference plot highlights the <strong>introgression-specific
signal</strong> above the genome-wide background, with the peak region showing maximal differentiation.</p>
</div>
<div class="plot-container"><img src="data:image/png;base64,{sparsepainter_b64}" alt="SparsePainter"></div>

<h2>Phased Haplotype Heatmap <span class="new-badge">NEW in v3</span></h2>
<div class="approach">
<div class="approach-title">Individual haplotypes visible (2 rows per individual)</div>
<p style="font-size:0.85rem;">Each individual is represented by two rows (one per haplotype).
Blue = reference allele, Red = alternative allele. The phased data reveals which specific haplotype
carries the introgressed block, allowing heterozygous individuals to be clearly decomposed into
one NN-like and one DD-like haplotype.</p>
</div>
<div class="plot-container"><img src="data:image/png;base64,{hapheat_b64}" alt="Haplotype heatmap"></div>

<h2>v2 vs v3 &mdash; Core Region Comparison <span class="new-badge">NEW</span></h2>
<div class="approach">
<div class="approach-title">Impact of phasing on core region boundaries</div>
<p style="font-size:0.85rem;">Comparison of core introgression regions identified with unphased
diploid genotypes (v2) vs phased haplotypes (v3). Differences arise because haplotype-level analysis
can distinguish heterozygous carriers (1 introgressed haplotype) from homozygous carriers (2 introgressed
haplotypes), providing more precise boundary resolution.</p>
</div>
<div class="plot-container"><img src="data:image/png;base64,{compare_b64}" alt="v2 vs v3"></div>

<div class="grid">
<div class="card">
<h3>v2 Core Regions (unphased)</h3>
<table>
<thead><tr><th>#</th><th>Coordinates</th><th>Size</th></tr></thead>
<tbody>{v2_rows}</tbody>
</table>
</div>
<div class="card">
<h3>v3 Core Regions (phased)</h3>
<table>
<thead><tr><th>#</th><th>Coordinates</th><th>Size</th><th>Mean cov</th></tr></thead>
<tbody>{v3_rows}</tbody>
</table>
</div>
</div>

<h2>Top 30 Most Significant SNPs</h2>
<div class="card" style="overflow-x:auto;">
<table>
<thead><tr><th>#</th><th>Position</th><th>REF</th><th>ALT</th>
<th>AF(NN)</th><th>AF(DD)</th><th>Δ</th><th>FST</th><th>−log₁₀(p)</th></tr></thead>
<tbody>{top30_rows}</tbody>
</table>
</div>

<h2>Method Details</h2>
<div class="card">
<p><strong>Phased VCF source:</strong> SHAPEIT4 phasing of chr3 biallelic SNPs
(<code>chr3_phased.vcf.gz</code> from EHH analysis). {n_snps:,} SNPs in 9.5–13.0 Mb
(vs 25,679 unphased in v2; fewer because phasing requires biallelic, non-missing sites).</p>
<p><strong>Genotype extraction:</strong> Phased haplotypes (0|0, 0|1, 1|0, 1|1) parsed directly
from VCF GT field. Each diploid individual decomposed into two haplotypes.</p>
<p><strong>Statistics:</strong> Δ(SNP-index) = AF<sub>DD</sub> − AF<sub>NN</sub>;
F<sub>ST</sub> via Weir &amp; Cockerham (1984); Fisher's exact test on 2×2 allele count tables.
All computed from haplotype-level allele counts.</p>
<p><strong>Introgression tract mapping:</strong> NN consensus = modal allele across all {n_nn*2}
NN haplotypes at each SNP. Per-haplotype discordance from consensus smoothed with 50 kb windows
(step 10 kb). Introgression = contiguous region with >50% discordance. Core region = contiguous
segment with ≥90% of DD haplotypes introgressed.</p>
<p><strong>Introgression dosage:</strong> At each core region, individuals classified by count of
introgressed haplotypes (0, 1, or 2).</p>
<p><strong>Excluded samples:</strong> nr59, nr60 (high missingness in focal region).</p>
<p><strong>SparsePainter chromosome painting:</strong> SparsePainter v1.3.2 (Yang et al. 2025,
<em>Nature Communications</em>) was run with both NN (N=40) and DD (N=42) as reference populations.
All 82 individuals were painted as targets using the same VCF as reference, with automatic
leave-one-out for unbiased self-painting. Per-SNP local ancestry probabilities P(NN) and P(DD)
were computed for all 164 haplotypes and averaged per population to produce the chromosome painting plots.
Recombination scaling constant lambda = 8,822.55 (estimated by SparsePainter).</p>
</div>

</body></html>
"""

outpath = f"{OUTDIR}/chr3_QTL_peak_10-12Mb_v3.html"
with open(outpath, 'w') as f:
    f.write(html)

print(f"\nReport saved to {outpath}")
print("Done!")
