#!/usr/bin/env python3
"""
Generate LD heatmap plots for chromosomes 3 and 2.
Shows combined (all 84 samples), NN-only, and DD-only populations.
Chr3 = GWAS signal chromosome; Chr2 = control chromosome for comparison.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import base64
import io

PROJ = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/BayPass.gwas_maf0.1geno0.2"
LDDIR = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/LDplot"

# ===== Helper functions =====
def read_positions(vars_file):
    positions = []
    with open(vars_file) as f:
        for line in f:
            parts = line.strip().split(':')
            positions.append(int(parts[1]))
    return np.array(positions)

def bin_ld(ld_mat, positions, bin_size=10):
    n = len(positions)
    n_bins = n // bin_size
    binned = np.zeros((n_bins, n_bins))
    bin_pos = np.zeros(n_bins)
    for i in range(n_bins):
        i0, i1 = i * bin_size, (i + 1) * bin_size
        bin_pos[i] = positions[i0:i1].mean() / 1e6
        for j in range(i, n_bins):
            j0, j1 = j * bin_size, (j + 1) * bin_size
            binned[i, j] = np.mean(ld_mat[i0:i1, j0:j1])
            binned[j, i] = binned[i, j]
    return binned, bin_pos, n_bins

def compute_decay(ld_mat, positions, max_dist_mb=10.0, n_dist_bins=100, max_pairs=2000):
    dist_edges = np.linspace(0, max_dist_mb, n_dist_bins + 1)
    dist_centers = (dist_edges[:-1] + dist_edges[1:]) / 2
    decay = np.zeros(n_dist_bins)
    counts = np.zeros(n_dist_bins)
    n = len(positions)
    for i in range(n):
        for j in range(i + 1, min(i + max_pairs, n)):
            d = (positions[j] - positions[i]) / 1e6
            if d >= max_dist_mb:
                break
            bi = int(d / max_dist_mb * n_dist_bins)
            if bi >= n_dist_bins:
                bi = n_dist_bins - 1
            decay[bi] += ld_mat[i, j]
            counts[bi] += 1
    mask = counts > 0
    decay[mask] /= counts[mask]
    return dist_centers, decay, mask

def fig_to_b64(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return b64

ld_cmap = LinearSegmentedColormap.from_list('ld',
    ['#FFFFFF', '#FFFFCC', '#FFCC00', '#FF6600', '#CC0000', '#660000'])

# ===== Load chr3 data =====
print("Loading chr3 data...")
pos3 = read_positions(f"{LDDIR}/ld_nn.unphased.vcor2.vars")
ld3_nn = np.nan_to_num(np.loadtxt(f"{LDDIR}/ld_nn.unphased.vcor2"), nan=0.0)
ld3_dd = np.nan_to_num(np.loadtxt(f"{LDDIR}/ld_dd.unphased.vcor2"), nan=0.0)
ld3_all = np.nan_to_num(np.loadtxt(f"{LDDIR}/ld_all_chr3.unphased.vcor2"), nan=0.0)
print(f"  Chr3: {len(pos3)} SNPs, {pos3[0]/1e6:.2f}–{pos3[-1]/1e6:.2f} Mb")

# ===== Load chr2 data =====
print("Loading chr2 data...")
pos2 = read_positions(f"{LDDIR}/ld_nn_chr2.unphased.vcor2.vars")
ld2_nn = np.nan_to_num(np.loadtxt(f"{LDDIR}/ld_nn_chr2.unphased.vcor2"), nan=0.0)
ld2_dd = np.nan_to_num(np.loadtxt(f"{LDDIR}/ld_dd_chr2.unphased.vcor2"), nan=0.0)
ld2_all = np.nan_to_num(np.loadtxt(f"{LDDIR}/ld_all_chr2.unphased.vcor2"), nan=0.0)
print(f"  Chr2: {len(pos2)} SNPs, {pos2[0]/1e6:.2f}–{pos2[-1]/1e6:.2f} Mb")

# ===== Bin matrices =====
print("Binning LD matrices...")
ld3_nn_b, bp3, nb3 = bin_ld(ld3_nn, pos3)
ld3_dd_b, _, _ = bin_ld(ld3_dd, pos3)
ld3_all_b, _, _ = bin_ld(ld3_all, pos3)
ld2_nn_b, bp2, nb2 = bin_ld(ld2_nn, pos2)
ld2_dd_b, _, _ = bin_ld(ld2_dd, pos2)
ld2_all_b, _, _ = bin_ld(ld2_all, pos2)
print(f"  Chr3 binned to {nb3}x{nb3}, Chr2 binned to {nb2}x{nb2}")

# ===== Plot 1: 6-panel LD heatmaps (combined + NN + DD for each chr) =====
print("Generating 6-panel LD heatmaps...")
fig, axes = plt.subplots(2, 3, figsize=(24, 16))

for ax, mat, bp, title in [
    (axes[0, 0], ld3_all_b, bp3, 'Chr3 — All samples (N=84)'),
    (axes[0, 1], ld3_nn_b, bp3, 'Chr3 — NN (N=42)'),
    (axes[0, 2], ld3_dd_b, bp3, 'Chr3 — DD (N=42)'),
    (axes[1, 0], ld2_all_b, bp2, 'Chr2 — All samples (N=84)'),
    (axes[1, 1], ld2_nn_b, bp2, 'Chr2 — NN (N=42)'),
    (axes[1, 2], ld2_dd_b, bp2, 'Chr2 — DD (N=42)'),
]:
    im = ax.imshow(mat, cmap=ld_cmap, vmin=0, vmax=1, aspect='equal',
                   extent=[bp[0], bp[-1], bp[-1], bp[0]], interpolation='nearest')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('Position (Mb)', fontsize=10)
    ax.set_ylabel('Position (Mb)', fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='r²')

plt.suptitle('Pairwise LD (r²) — Chr3 (GWAS signal) vs Chr2 (control)', fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
sixpanel_b64 = fig_to_b64(fig)

# ===== Plot 2: LD difference heatmaps — Combined minus NN, Combined minus DD, NN minus DD =====
print("Generating LD difference plots...")
fig, axes = plt.subplots(2, 3, figsize=(24, 16))

# Chr3 row
diff3_nn_dd = ld3_nn_b - ld3_dd_b
diff3_all_nn = ld3_all_b - ld3_nn_b
diff3_all_dd = ld3_all_b - ld3_dd_b

# Chr2 row
diff2_nn_dd = ld2_nn_b - ld2_dd_b
diff2_all_nn = ld2_all_b - ld2_nn_b
diff2_all_dd = ld2_all_b - ld2_dd_b

for ax, mat, bp, title in [
    (axes[0, 0], diff3_nn_dd, bp3, 'Chr3 — NN − DD'),
    (axes[0, 1], diff3_all_nn, bp3, 'Chr3 — All − NN'),
    (axes[0, 2], diff3_all_dd, bp3, 'Chr3 — All − DD'),
    (axes[1, 0], diff2_nn_dd, bp2, 'Chr2 — NN − DD'),
    (axes[1, 1], diff2_all_nn, bp2, 'Chr2 — All − NN'),
    (axes[1, 2], diff2_all_dd, bp2, 'Chr2 — All − DD'),
]:
    im = ax.imshow(mat, cmap='RdBu_r', vmin=-0.5, vmax=0.5, aspect='equal',
                   extent=[bp[0], bp[-1], bp[-1], bp[0]], interpolation='nearest')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('Position (Mb)', fontsize=10)
    ax.set_ylabel('Position (Mb)', fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Δr²')

plt.suptitle('LD Difference Heatmaps — Chr3 vs Chr2', fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
diff_b64 = fig_to_b64(fig)

# ===== Plot 3: LD decay comparison (3 lines per panel: All, NN, DD) =====
print("Generating LD decay comparison...")
dc3, d3_all, m3_all = compute_decay(ld3_all, pos3)
_, d3_nn, m3_nn = compute_decay(ld3_nn, pos3)
_, d3_dd, m3_dd = compute_decay(ld3_dd, pos3)
dc2, d2_all, m2_all = compute_decay(ld2_all, pos2)
_, d2_nn, m2_nn = compute_decay(ld2_nn, pos2)
_, d2_dd, m2_dd = compute_decay(ld2_dd, pos2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5), sharey=True)

# Chr3
m3 = m3_all & m3_nn & m3_dd
ax1.plot(dc3[m3], d3_all[m3], '-', color='#2ecc71', linewidth=2.0, label='All (N=84)', alpha=0.9)
ax1.plot(dc3[m3], d3_nn[m3], '-', color='#3498db', linewidth=1.5, label='NN (N=42)', alpha=0.8)
ax1.plot(dc3[m3], d3_dd[m3], '-', color='#e74c3c', linewidth=1.5, label='DD (N=42)', alpha=0.8)
ax1.set_title('Chromosome 3 (GWAS signal)', fontsize=12, fontweight='bold')
ax1.set_xlabel('Physical distance (Mb)', fontsize=11)
ax1.set_ylabel('Mean r²', fontsize=11)
ax1.legend(fontsize=9)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_xlim(0, 10)

# Chr2
m2 = m2_all & m2_nn & m2_dd
ax2.plot(dc2[m2], d2_all[m2], '-', color='#2ecc71', linewidth=2.0, label='All (N=84)', alpha=0.9)
ax2.plot(dc2[m2], d2_nn[m2], '-', color='#3498db', linewidth=1.5, label='NN (N=42)', alpha=0.8)
ax2.plot(dc2[m2], d2_dd[m2], '-', color='#e74c3c', linewidth=1.5, label='DD (N=42)', alpha=0.8)
ax2.set_title('Chromosome 2 (control)', fontsize=12, fontweight='bold')
ax2.set_xlabel('Physical distance (Mb)', fontsize=11)
ax2.legend(fontsize=9)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlim(0, 10)

plt.suptitle('LD Decay — Chr3 vs Chr2', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
decay_compare_b64 = fig_to_b64(fig)

# ===== Plot 4: Regional LD decay (signal vs outside vs chr2), 3 lines each =====
print("Generating regional LD decay...")
signal_start, signal_end = 4.0, 19.0
max_dist_mb = 10.0
n_dist_bins = 100
dist_edges = np.linspace(0, max_dist_mb, n_dist_bins + 1)
dist_centers = (dist_edges[:-1] + dist_edges[1:]) / 2

# Chr3 signal region vs outside — for all three groups
all_in = np.zeros(n_dist_bins); nn_in = np.zeros(n_dist_bins); dd_in = np.zeros(n_dist_bins); c_in = np.zeros(n_dist_bins)
all_out = np.zeros(n_dist_bins); nn_out = np.zeros(n_dist_bins); dd_out = np.zeros(n_dist_bins); c_out = np.zeros(n_dist_bins)

for i in range(len(pos3)):
    pi = pos3[i] / 1e6
    for j in range(i + 1, min(i + 2000, len(pos3))):
        pj = pos3[j] / 1e6
        d = pj - pi
        if d >= max_dist_mb:
            break
        bi = min(int(d / max_dist_mb * n_dist_bins), n_dist_bins - 1)
        if signal_start <= pi <= signal_end and signal_start <= pj <= signal_end:
            all_in[bi] += ld3_all[i, j]; nn_in[bi] += ld3_nn[i, j]; dd_in[bi] += ld3_dd[i, j]; c_in[bi] += 1
        elif (pi < signal_start or pi > signal_end) and (pj < signal_start or pj > signal_end):
            all_out[bi] += ld3_all[i, j]; nn_out[bi] += ld3_nn[i, j]; dd_out[bi] += ld3_dd[i, j]; c_out[bi] += 1

mi = c_in > 0; all_in[mi] /= c_in[mi]; nn_in[mi] /= c_in[mi]; dd_in[mi] /= c_in[mi]
mo = c_out > 0; all_out[mo] /= c_out[mo]; nn_out[mo] /= c_out[mo]; dd_out[mo] /= c_out[mo]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 5), sharey=True)

ax1.plot(dist_centers[mi], all_in[mi], '-', color='#2ecc71', linewidth=2.0, label='All (N=84)')
ax1.plot(dist_centers[mi], nn_in[mi], '-', color='#3498db', linewidth=1.5, label='NN (N=42)')
ax1.plot(dist_centers[mi], dd_in[mi], '-', color='#e74c3c', linewidth=1.5, label='DD (N=42)')
ax1.set_title('Chr3: signal region\n(4–19 Mb)', fontsize=11, fontweight='bold')
ax1.set_xlabel('Distance (Mb)', fontsize=10)
ax1.set_ylabel('Mean r²', fontsize=10)
ax1.legend(fontsize=9)
ax1.spines['top'].set_visible(False); ax1.spines['right'].set_visible(False)
ax1.set_xlim(0, 10)

ax2.plot(dist_centers[mo], all_out[mo], '-', color='#2ecc71', linewidth=2.0, label='All (N=84)')
ax2.plot(dist_centers[mo], nn_out[mo], '-', color='#3498db', linewidth=1.5, label='NN (N=42)')
ax2.plot(dist_centers[mo], dd_out[mo], '-', color='#e74c3c', linewidth=1.5, label='DD (N=42)')
ax2.set_title('Chr3: outside signal\n(<4 or >19 Mb)', fontsize=11, fontweight='bold')
ax2.set_xlabel('Distance (Mb)', fontsize=10)
ax2.legend(fontsize=9)
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
ax2.set_xlim(0, 10)

ax3.plot(dc2[m2], d2_all[m2], '-', color='#2ecc71', linewidth=2.0, label='All (N=84)')
ax3.plot(dc2[m2], d2_nn[m2], '-', color='#3498db', linewidth=1.5, label='NN (N=42)')
ax3.plot(dc2[m2], d2_dd[m2], '-', color='#e74c3c', linewidth=1.5, label='DD (N=42)')
ax3.set_title('Chr2: whole chromosome\n(control)', fontsize=11, fontweight='bold')
ax3.set_xlabel('Distance (Mb)', fontsize=10)
ax3.legend(fontsize=9)
ax3.spines['top'].set_visible(False); ax3.spines['right'].set_visible(False)
ax3.set_xlim(0, 10)

plt.suptitle('LD Decay — Chr3 Signal Region vs Background vs Chr2 Control', fontsize=13, fontweight='bold', y=1.03)
plt.tight_layout()
regional_b64 = fig_to_b64(fig)

# ===== Compute summary stats =====
print("Computing LD summary stats...")

def region_mean_r2(mat, bp, start, end):
    idx = np.where((bp >= start) & (bp <= end))[0]
    if len(idx) < 2:
        return np.nan
    return np.mean(mat[np.ix_(idx, idx)][np.triu_indices(len(idx), k=1)])

def whole_mean_r2(mat):
    n = mat.shape[0]
    return np.mean(mat[np.triu_indices(n, k=1)])

# Chr3
all3_sig = region_mean_r2(ld3_all_b, bp3, signal_start, signal_end)
nn3_sig = region_mean_r2(ld3_nn_b, bp3, signal_start, signal_end)
dd3_sig = region_mean_r2(ld3_dd_b, bp3, signal_start, signal_end)
all3_all = whole_mean_r2(ld3_all_b)
nn3_all = whole_mean_r2(ld3_nn_b)
dd3_all = whole_mean_r2(ld3_dd_b)

# Chr2
all2_all = whole_mean_r2(ld2_all_b)
nn2_all = whole_mean_r2(ld2_nn_b)
dd2_all = whole_mean_r2(ld2_dd_b)

print(f"  Chr3 signal (4-19 Mb): All={all3_sig:.4f}, NN={nn3_sig:.4f}, DD={dd3_sig:.4f}")
print(f"  Chr3 whole: All={all3_all:.4f}, NN={nn3_all:.4f}, DD={dd3_all:.4f}")
print(f"  Chr2 whole: All={all2_all:.4f}, NN={nn2_all:.4f}, DD={dd2_all:.4f}")

# ===== Generate HTML =====
print("Writing HTML report...")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>LDplot — Chr3 vs Chr2 LD Analysis</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; background: #fafafa; }}
h1 {{ color: #2c3e50; border-bottom: 3px solid #e67e22; padding-bottom: 10px; }}
h2 {{ color: #34495e; margin-top: 40px; border-bottom: 1px solid #bdc3c7; padding-bottom: 8px; }}
.summary-box {{ background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 20px; margin: 15px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
.stat {{ display: inline-block; margin: 10px 20px 10px 0; padding: 10px 20px; background: #ecf0f1; border-radius: 5px; }}
.stat .value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
.stat .label {{ font-size: 12px; color: #7f8c8d; }}
img {{ max-width: 100%; border: 1px solid #ddd; border-radius: 4px; margin: 10px 0; }}
table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
th {{ background: #e67e22; color: white; }}
tr:nth-child(even) {{ background: #f2f2f2; }}
.highlight {{ background: #fff3cd; border: 1px solid #ffc107; padding: 15px; border-radius: 5px; margin: 15px 0; }}
.method {{ background: #fdf2e9; border-left: 4px solid #e67e22; padding: 15px; margin: 15px 0; }}
pre.code-block {{ background: #1e1e1e; color: #d4d4d4; padding: 16px; border-radius: 6px; overflow-x: auto; font-size: 13px; line-height: 1.5; margin: 10px 0; }}
pre.code-block .comment {{ color: #6a9955; }}
pre.code-block .result {{ color: #9cdcfe; font-style: italic; }}
</style>
</head>
<body>

<h1>LDplot — Chromosome 3 vs Chromosome 2 LD Analysis</h1>
<p><em>Pararge aegeria — Combined (all samples), NN, and DD populations</em></p>

<div class="summary-box">
<div class="stat"><div class="value">{len(pos3):,} + {len(pos2):,}</div><div class="label">SNPs (chr3 + chr2 thinned)</div></div>
<div class="stat"><div class="value">84</div><div class="label">Total samples</div></div>
<div class="stat"><div class="value">42 + 42</div><div class="label">NN + DD samples</div></div>
<div class="stat"><div class="value">Chr3</div><div class="label">GWAS signal chromosome</div></div>
<div class="stat"><div class="value">Chr2</div><div class="label">Control chromosome</div></div>
</div>

<div class="method">
<strong>Methods:</strong> Pairwise r&sup2; (unphased) computed with PLINK2 for three sample groups:
all 84 samples combined, NN only (N=42), and DD only (N=42).
Chr3 ({len(pos3):,} thinned SNPs from 125,313) and chr2 ({len(pos2):,} thinned SNPs from 143,795).
Matrices binned to ~400&times;400 blocks for visualization. Chr3 contains the GWAS signal region (4&ndash;19 Mb);
chr2 serves as a control chromosome with no significant GWAS hits.
</div>

<h2>Analysis Code</h2>
<pre class="code-block"><span class="comment"># Thin SNPs to ~4000 per chromosome</span>
awk '$1 == "3" {{print $2}}' variants_qc_samples.bim | awk 'NR % 31 == 1' > chr3_thin_snps.txt
awk '$1 == "2" {{print $2}}' variants_qc_samples.bim | awk 'NR % 36 == 1' > chr2_thin_snps.txt
<span class="result"># Chr3: 4,043 SNPs; Chr2: 3,995 SNPs</span>

<span class="comment"># Compute pairwise LD — combined (all samples)</span>
for CHR in 3 2; do
  plink2 --bfile variants_qc_samples \\
      --extract chr${{CHR}}_thin_snps.txt \\
      --chr $CHR \\
      --r2-unphased square \\
      --chr-set 27 --allow-extra-chr \\
      --out ld_all_chr${{CHR}}
done

<span class="comment"># Compute pairwise LD — per population (NN, DD)</span>
for CHR in 3 2; do
  for POP in NN DD; do
    plink2 --bfile variants_qc_samples \\
        --keep samples_${{POP}}.txt \\
        --extract chr${{CHR}}_thin_snps.txt \\
        --chr $CHR \\
        --r2-unphased square \\
        --chr-set 27 --allow-extra-chr \\
        --out ld_${{POP,,}}_chr${{CHR}}
  done
done
</pre>

<h2>LD Heatmaps — Combined, NN, and DD</h2>
<div class="highlight">
<strong>Key comparison:</strong> Top row = Chr3 (GWAS-associated chromosome). Bottom row = Chr2 (control).
Left column = all 84 samples combined; middle = NN only; right = DD only.
The combined panel shows the overall LD structure. Large contiguous blocks of high LD that appear in one
population but not the other (or differ between chr3 and chr2) suggest structural variants or inversions.
</div>
<img src="data:image/png;base64,{sixpanel_b64}" alt="6-panel LD Heatmaps">

<h2>LD Difference Heatmaps</h2>
<div class="highlight">
<strong>Interpretation:</strong> Left column: NN &minus; DD difference. Middle: All &minus; NN (shows what DD adds to the combined).
Right: All &minus; DD (shows what NN adds to the combined). Red = higher LD in the first group; Blue = higher in the second.
Chr3 should show structured differential LD in the GWAS signal region (4&ndash;19 Mb) if a structural variant segregates between populations.
Chr2 (bottom row) should show minimal differential LD as a neutral control.
</div>
<img src="data:image/png;base64,{diff_b64}" alt="LD Difference Heatmaps">

<h2>LD Decay — Chr3 vs Chr2</h2>
<div class="highlight">
<strong>Interpretation:</strong> Green = all 84 samples combined; Blue = NN only; Red = DD only.
If LD extends over longer distances in chr3 compared to chr2, this supports extended haplotype blocks
(e.g., from an inversion) in the GWAS region. The combined line shows how population structure
inflates LD relative to within-population estimates.
</div>
<img src="data:image/png;base64,{decay_compare_b64}" alt="LD Decay Comparison">

<h2>LD Decay — Chr3 Signal Region vs Background vs Chr2</h2>
<div class="highlight">
<strong>Three-way comparison:</strong> (Left) LD decay within the chr3 GWAS signal region (4&ndash;19 Mb).
(Middle) LD decay in chr3 flanking regions outside the signal. (Right) LD decay across all of chr2.
All three sample groups shown. In the signal region, the combined (green) line may be elevated above
both within-population lines if the NN/DD haplotype structure creates apparent LD when pooled.
</div>
<img src="data:image/png;base64,{regional_b64}" alt="Regional LD Decay">

<h2>LD Summary Statistics</h2>
<div class="summary-box">
<table>
<tr><th>Region</th><th>All (N=84) mean r&sup2;</th><th>NN (N=42) mean r&sup2;</th><th>DD (N=42) mean r&sup2;</th><th>&Delta; (NN&minus;DD)</th></tr>
<tr><td><strong>Chr3: signal (4&ndash;19 Mb)</strong></td><td>{all3_sig:.4f}</td><td>{nn3_sig:.4f}</td><td>{dd3_sig:.4f}</td><td>{nn3_sig - dd3_sig:+.4f}</td></tr>
<tr><td>Chr3: whole chromosome</td><td>{all3_all:.4f}</td><td>{nn3_all:.4f}</td><td>{dd3_all:.4f}</td><td>{nn3_all - dd3_all:+.4f}</td></tr>
<tr><td>Chr2: whole chromosome</td><td>{all2_all:.4f}</td><td>{nn2_all:.4f}</td><td>{dd2_all:.4f}</td><td>{nn2_all - dd2_all:+.4f}</td></tr>
</table>
</div>

</body>
</html>
"""

with open(f"{LDDIR}/ld_report.html", 'w') as f:
    f.write(html)

print(f"Report saved to {LDDIR}/ld_report.html")
print("Done!")
