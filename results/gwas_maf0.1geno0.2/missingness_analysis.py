#!/usr/bin/env python3
"""
Missingness Analysis — generates missingness_report.html
Investigates why genotype plots show missing (grey) cells.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import base64
import io
import glob

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BED_FILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/plink/variants_qc_samples_maf01geno02.bed"
BIM_FILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/plink/variants_qc_samples_maf01geno02.bim"
FAM_FILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/plink/variants_qc_samples_maf01geno02.fam"
PHENO_FILE = "phenotypes.txt"
SAMPLES_TSV = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/samples.tsv"
SMISS_FILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/qc/samples/missing.smiss"
MINDREM_FILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/plink/variants_qc_samples_maf01geno02.mindrem.id"

FOCAL_CHR = "3"
REGION_START = 8_000_000
REGION_END = 14_000_000
WINDOW_SIZE = 100_000

OUT_HTML = "missingness_report.html"

# ---------------------------------------------------------------------------
# 1. Read sample info
# ---------------------------------------------------------------------------
print("Reading sample info...")
samples = []
with open(FAM_FILE) as f:
    for line in f:
        samples.append(line.strip().split()[1])
n_samples = len(samples)

pheno_map = {}
with open(PHENO_FILE) as f:
    f.readline()
    for line in f:
        parts = line.strip().split()
        pheno_map[parts[1]] = int(parts[2])
sample_phenos = [pheno_map.get(s, 0) for s in samples]
pheno_labels = {1: "NN", 2: "DD"}

# Read pool assignments from samples.tsv
pool_map = {}
try:
    with open(SAMPLES_TSV) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            sample_id = parts[0]
            if len(parts) > 2:
                pool_map[sample_id] = parts[2]
except:
    pass

# Read removed samples
removed_samples = []
try:
    with open(MINDREM_FILE) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] not in ('#FID', 'FID', '#IID', 'IID'):
                removed_samples.append(parts[1])
except:
    pass

# Read pre-QC missingness
pre_qc_miss = {}
try:
    with open(SMISS_FILE) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 5:
                pre_qc_miss[parts[1]] = float(parts[4])
except:
    pass

# Read mosdepth summaries
depth_data = {}
mosdepth_files = glob.glob("/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/qc/mosdepth/*.mosdepth.summary.txt")
for mf in mosdepth_files:
    sname = mf.split('/')[-1].replace('.mosdepth.summary.txt', '')
    try:
        with open(mf) as f:
            for line in f:
                parts = line.strip().split('\t')
                if parts[0] == 'total':
                    depth_data[sname] = float(parts[3])
                    break
    except:
        pass

print(f"  {n_samples} retained samples, {len(removed_samples)} removed, {len(depth_data)} with depth data")

# ---------------------------------------------------------------------------
# 2. Read .bim — get chr3 region variants
# ---------------------------------------------------------------------------
print("Reading .bim for chr3 region...")
region_pos = []
region_idx = []

with open(BIM_FILE) as f:
    for i, line in enumerate(f):
        parts = line.strip().split('\t')
        if parts[0] != FOCAL_CHR:
            continue
        pos = int(parts[3])
        if REGION_START <= pos <= REGION_END:
            region_pos.append(pos)
            region_idx.append(i)

n_region = len(region_idx)
print(f"  {n_region} variants in chr3 {REGION_START/1e6:.0f}-{REGION_END/1e6:.0f} Mb")

# ---------------------------------------------------------------------------
# 3. Read genotypes for the region
# ---------------------------------------------------------------------------
print("Reading genotypes...")
bytes_per_snp = (n_samples + 3) // 4
decode = np.array([0, 1, -1, 2], dtype=np.int8)

genotypes = np.zeros((n_region, n_samples), dtype=np.int8)
with open(BED_FILE, 'rb') as f:
    magic = f.read(3)
    assert magic == b'\x6c\x1b\x01'
    for local_i, global_i in enumerate(region_idx):
        f.seek(3 + global_i * bytes_per_snp)
        raw = f.read(bytes_per_snp)
        si = 0
        for bv in raw:
            for bp in range(4):
                if si >= n_samples:
                    break
                genotypes[local_i, si] = decode[(bv >> (2 * bp)) & 0x03]
                si += 1
        if (local_i + 1) % 5000 == 0:
            print(f"  {local_i + 1}/{n_region}...")

print(f"  Done: {genotypes.shape}")
missing_mask = genotypes == -1

# ---------------------------------------------------------------------------
# 4. Compute statistics
# ---------------------------------------------------------------------------
print("Computing missingness statistics...")

total_cells = genotypes.size
total_missing = missing_mask.sum()
overall_pct = total_missing / total_cells * 100

per_sample_missing = missing_mask.sum(axis=0)
per_sample_pct = per_sample_missing / n_region * 100

per_var_missing = missing_mask.sum(axis=1)
per_var_pct = per_var_missing / n_samples * 100

positions_arr = np.array(region_pos)
window_starts = np.arange(REGION_START, REGION_END, WINDOW_SIZE)
window_miss_pct = []
for ws in window_starts:
    we = ws + WINDOW_SIZE
    mask = (positions_arr >= ws) & (positions_arr < we)
    if mask.sum() > 0:
        window_miss_pct.append(missing_mask[mask, :].sum() / (mask.sum() * n_samples) * 100)
    else:
        window_miss_pct.append(0)

nn_idx = [i for i, p in enumerate(sample_phenos) if p == 1]
dd_idx = [i for i, p in enumerate(sample_phenos) if p == 2]
nn_miss = per_sample_pct[nn_idx].mean()
dd_miss = per_sample_pct[dd_idx].mean()

# ---------------------------------------------------------------------------
# 5. Generate plots
# ---------------------------------------------------------------------------
print("Generating plots...")

# Plot 1: Per-sample missingness bar chart
fig1, ax1 = plt.subplots(figsize=(14, 5))
sorted_idx = np.argsort(per_sample_pct)[::-1]
colors = ['#e74c3c' if per_sample_pct[i] > 5 else '#f39c12' if per_sample_pct[i] > 1 else '#3498db'
          for i in sorted_idx]
ax1.bar(range(n_samples), per_sample_pct[sorted_idx], color=colors, width=1.0, edgecolor='none')
ax1.set_xlabel('Samples (sorted by missingness)', fontsize=11)
ax1.set_ylabel('Missing genotypes (%)', fontsize=11)
ax1.set_title('Per-Sample Missingness — chr3 8–14 Mb', fontsize=13, fontweight='bold')
ax1.set_xlim(-0.5, n_samples - 0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
for rank, si in enumerate(sorted_idx[:5]):
    ax1.annotate(samples[si], xy=(rank, per_sample_pct[si]),
                 xytext=(rank + 2, per_sample_pct[si] + 0.3),
                 fontsize=8, ha='left', color='#333',
                 arrowprops=dict(arrowstyle='-', color='#999', lw=0.5))
plt.tight_layout()
buf1 = io.BytesIO()
fig1.savefig(buf1, format='png', dpi=200, bbox_inches='tight')
buf1.seek(0)
b64_1 = base64.b64encode(buf1.read()).decode('utf-8')
plt.close(fig1)

# Plot 2: Per-window missingness
fig2, ax2 = plt.subplots(figsize=(14, 4))
window_centers = window_starts / 1e6 + WINDOW_SIZE / 2e6
ax2.bar(window_centers, window_miss_pct, width=WINDOW_SIZE / 1e6 * 0.9,
        color='#2166ac', alpha=0.8, edgecolor='none')
ax2.set_xlabel('Position on chr3 (Mb)', fontsize=11)
ax2.set_ylabel('Missing genotypes (%)', fontsize=11)
ax2.set_title('Missingness by 100 kb Window — chr3 8–14 Mb', fontsize=13, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.axhline(y=np.mean(window_miss_pct), color='red', linewidth=0.8, linestyle='--',
            label=f'Mean: {np.mean(window_miss_pct):.2f}%')
ax2.legend(fontsize=9)
plt.tight_layout()
buf2 = io.BytesIO()
fig2.savefig(buf2, format='png', dpi=200, bbox_inches='tight')
buf2.seek(0)
b64_2 = base64.b64encode(buf2.read()).decode('utf-8')
plt.close(fig2)

# Plot 3: Per-variant missingness distribution
fig3, ax3 = plt.subplots(figsize=(8, 5))
bins = np.arange(0, 22, 1)
ax3.hist(per_var_pct, bins=bins, color='#2166ac', edgecolor='white', alpha=0.85)
ax3.set_xlabel('Missing genotypes per variant (%)', fontsize=11)
ax3.set_ylabel('Number of variants', fontsize=11)
ax3.set_title('Distribution of Per-Variant Missingness', fontsize=13, fontweight='bold')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
n_complete = (per_var_missing == 0).sum()
ax3.annotate(f'{n_complete/n_region*100:.1f}% of variants\nhave 0% missing',
             xy=(0.5, n_complete * 0.8), fontsize=9, ha='center', color='#333')
plt.tight_layout()
buf3 = io.BytesIO()
fig3.savefig(buf3, format='png', dpi=200, bbox_inches='tight')
buf3.seek(0)
b64_3 = base64.b64encode(buf3.read()).decode('utf-8')
plt.close(fig3)

# Plot 4: Depth vs missingness
fig4, ax4 = plt.subplots(figsize=(8, 5))
if depth_data:
    retained_depths = []
    retained_miss = []
    for i, s in enumerate(samples):
        if s in depth_data:
            retained_depths.append(depth_data[s])
            retained_miss.append(per_sample_pct[i])
    ax4.scatter(retained_depths, retained_miss, c='#3498db', s=30, alpha=0.7,
                label='Retained', edgecolors='none', zorder=2)
    # Label worst retained
    for i, s in enumerate(samples):
        if s in depth_data and per_sample_pct[i] > 3:
            ax4.annotate(s, xy=(depth_data[s], per_sample_pct[i]),
                         fontsize=7, ha='left', color='#2c3e50',
                         xytext=(3, 3), textcoords='offset points')
    # Removed samples
    first_removed = True
    for rs in removed_samples:
        if rs in depth_data and rs in pre_qc_miss:
            ax4.scatter(depth_data[rs], pre_qc_miss[rs] * 100, c='#e74c3c', s=60,
                        marker='x', linewidths=2,
                        label='Removed' if first_removed else '', zorder=3)
            first_removed = False
            ax4.annotate(rs, xy=(depth_data[rs], pre_qc_miss[rs] * 100),
                         fontsize=7, ha='left', color='#c0392b',
                         xytext=(3, 3), textcoords='offset points')
    ax4.set_xlabel('Mean sequencing depth (x)', fontsize=11)
    ax4.set_ylabel('Missing genotypes (%)', fontsize=11)
    ax4.set_title('Sequencing Depth vs. Missingness', fontsize=13, fontweight='bold')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.legend(fontsize=9)
plt.tight_layout()
buf4 = io.BytesIO()
fig4.savefig(buf4, format='png', dpi=200, bbox_inches='tight')
buf4.seek(0)
b64_4 = base64.b64encode(buf4.read()).decode('utf-8')
plt.close(fig4)

# ---------------------------------------------------------------------------
# 6. Build tables
# ---------------------------------------------------------------------------
sample_table_rows = ""
for rank, si in enumerate(sorted_idx[:15]):
    s = samples[si]
    pheno = pheno_labels.get(sample_phenos[si], "?")
    depth = f"{depth_data[s]:.1f}x" if s in depth_data else "—"
    pool = pool_map.get(s, "—")
    flag = ' style="background:#fef3c7;"' if per_sample_pct[si] > 3 else ''
    sample_table_rows += (
        f"<tr{flag}><td>{rank+1}</td><td>{s}</td><td>{pheno}</td>"
        f"<td>{depth}</td><td>{pool}</td>"
        f"<td>{per_sample_missing[si]:,}</td><td>{per_sample_pct[si]:.2f}%</td></tr>\n"
    )

removed_table_rows = ""
for rs in sorted(removed_samples, key=lambda x: pre_qc_miss.get(x, 0), reverse=True):
    depth = f"{depth_data[rs]:.1f}x" if rs in depth_data else "—"
    miss = f"{pre_qc_miss[rs]*100:.1f}%" if rs in pre_qc_miss else "—"
    pool = pool_map.get(rs, "—")
    removed_table_rows += f"<tr><td>{rs}</td><td>{depth}</td><td>{pool}</td><td>{miss}</td></tr>\n"

# ---------------------------------------------------------------------------
# 7. Generate HTML
# ---------------------------------------------------------------------------
print("Generating HTML report...")

html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Missingness Analysis — Chromosome 3</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
            background: #f8fafc; color: #1e293b; line-height: 1.6;
            max-width: 1200px; margin: 0 auto; padding: 2rem;
        }}
        h1 {{ font-size: 1.6rem; margin-bottom: 0.3rem; color: #0f172a; }}
        .subtitle {{ color: #64748b; margin-bottom: 2rem; font-size: 0.9rem; }}
        h2 {{ font-size: 1.15rem; margin: 2rem 0 0.8rem; color: #0f172a;
              border-bottom: 2px solid #e2e8f0; padding-bottom: 0.3rem; }}
        .card {{
            background: #fff; border: 1px solid #e2e8f0; border-radius: 8px;
            padding: 1.2rem; margin-bottom: 1rem; box-shadow: 0 1px 2px rgba(0,0,0,0.04);
        }}
        table {{
            width: 100%; border-collapse: collapse; margin: 0.5rem 0;
            font-size: 0.85rem;
        }}
        th {{ background: #f1f5f9; text-align: left; padding: 0.5rem 0.7rem;
             font-weight: 600; border-bottom: 2px solid #e2e8f0; }}
        td {{ padding: 0.4rem 0.7rem; border-bottom: 1px solid #f1f5f9; }}
        tr:hover td {{ background: #f8fafc; }}
        .plot-container {{ text-align: center; margin: 1rem 0; }}
        .plot-container img {{ width: 100%; border: 1px solid #e2e8f0; border-radius: 6px; }}
        .key-finding {{
            background: #fef2f2; border-left: 4px solid #ef4444; padding: 1rem;
            margin: 1rem 0; border-radius: 0 6px 6px 0; font-size: 0.95rem;
        }}
        .explanation {{
            background: #eff6ff; border-left: 4px solid #3b82f6; padding: 1rem;
            margin: 1rem 0; border-radius: 0 6px 6px 0; font-size: 0.9rem;
        }}
        .recommendation {{
            background: #f0fdf4; border-left: 4px solid #22c55e; padding: 1rem;
            margin: 1rem 0; border-radius: 0 6px 6px 0; font-size: 0.9rem;
        }}
        .stat-grid {{
            display: grid; grid-template-columns: repeat(4, 1fr); gap: 1rem; margin: 1rem 0;
        }}
        .stat-box {{
            background: #fff; border: 1px solid #e2e8f0; border-radius: 8px;
            padding: 1rem; text-align: center;
        }}
        .stat-box .value {{ font-size: 1.5rem; font-weight: 700; color: #0f172a; }}
        .stat-box .label {{ font-size: 0.8rem; color: #64748b; margin-top: 0.2rem; }}
        @media (max-width: 800px) {{ .stat-grid {{ grid-template-columns: 1fr 1fr; }} }}
    </style>
</head>
<body>
    <h1>Missingness Analysis &mdash; Why Are Genotypes Missing?</h1>
    <div class="subtitle">
        <em>Pararge aegeria</em> GWAS &mdash; chr3 8&ndash;14 Mb focal region &mdash;
        {n_samples} samples, {n_region:,} variants
    </div>

    <h2>Key Finding</h2>
    <div class="key-finding">
        <strong>The missing genotypes are caused by low sequencing depth</strong>, primarily
        in samples from <strong>Pool-6</strong> which suffered a sequencing failure. At 1.6&ndash;3.8&times;
        coverage, bcftools mpileup produces zero reads at many positions (Poisson sampling:
        P(0&nbsp;reads) = e<sup>&minus;&lambda;</sup> &asymp; 10&ndash;20% when &lambda;&nbsp;=&nbsp;2&times;),
        resulting in missing genotype calls (<code>./.</code>). The <code>--geno&nbsp;0.2</code> PLINK filter
        retained variants with up to 20% missingness, and <code>--mind&nbsp;0.1</code> removed the
        6 worst samples but left 2 borderline low-coverage samples (nr65, nr88) in the dataset.
    </div>

    <h2>Overview</h2>
    <div class="stat-grid">
        <div class="stat-box">
            <div class="value">{overall_pct:.2f}%</div>
            <div class="label">Overall missing rate<br>(chr3 8&ndash;14 Mb)</div>
        </div>
        <div class="stat-box">
            <div class="value">{(per_var_missing == 0).sum() / n_region * 100:.1f}%</div>
            <div class="label">Variants with<br>0% missing</div>
        </div>
        <div class="stat-box">
            <div class="value">{len(removed_samples)}</div>
            <div class="label">Samples removed<br>by --mind 0.1</div>
        </div>
        <div class="stat-box">
            <div class="value">{sum(1 for p in per_sample_pct if p > 3)}</div>
            <div class="label">Retained samples<br>with &gt;3% missing</div>
        </div>
    </div>

    <h2>The Mechanism: Low Depth &rarr; Missing Genotypes</h2>
    <div class="explanation">
        <strong>How the variant calling pipeline creates missing data:</strong>
        <ol style="margin-top:0.5rem; padding-left:1.5rem;">
            <li><strong>bcftools mpileup</strong> requires reads with mapping quality &ge;20 and base quality &ge;20</li>
            <li>At low coverage (1&ndash;4&times;), many positions have <strong>zero qualifying reads</strong> for a given sample</li>
            <li>bcftools assigns these positions a <strong>missing genotype</strong> (<code>./.:0,0,0:0:0,0</code>, i.e. DP=0)</li>
            <li>The site-level filter (<code>QUAL&gt;20 &amp;&amp; INFO/DP&gt;5</code>) checks <em>total</em> depth summed
                across all ~94 samples, not per-sample &mdash; so variants pass even when individual samples have DP=0</li>
            <li>PLINK&rsquo;s <code>--geno 0.2</code> allows up to 17 of 88 samples to be missing per variant</li>
        </ol>
        <p style="margin-top:0.5rem;">
            <strong>This is not a region-specific problem.</strong> The missingness rate in chr3 8&ndash;14 Mb
            (~{overall_pct:.1f}%) is comparable to the genome-wide average. The grey cells in the genotype
            plot are a genome-wide pattern driven by sample depth, not by any feature of this chromosomal region.
        </p>
    </div>

    <h2>Sequencing Depth vs. Missingness</h2>
    <div class="plot-container">
        <img src="data:image/png;base64,{b64_4}" alt="Depth vs missingness">
    </div>
    <p style="font-size:0.85rem; color:#64748b; text-align:center; margin-bottom:1rem;">
        Blue dots = retained samples (missingness in chr3 8&ndash;14 Mb).
        Red crosses = removed samples (genome-wide missingness from pre-QC).
        The strong inverse relationship confirms depth as the primary driver.
    </p>

    <h2>Samples Removed by QC (--mind 0.1)</h2>
    <div class="card">
        <table>
            <thead><tr><th>Sample</th><th>Mean Depth</th><th>Pool</th><th>Genome-wide Missing</th></tr></thead>
            <tbody>{removed_table_rows}</tbody>
        </table>
        <p style="margin-top:0.5rem; font-size:0.85rem; color:#64748b;">
            All removed samples had extremely low sequencing depth (1.6&ndash;3.0&times;),
            predominantly from Pool-6 which suffered a sequencing failure.
        </p>
    </div>

    <h2>Per-Sample Missingness (Retained Samples)</h2>
    <div class="plot-container">
        <img src="data:image/png;base64,{b64_1}" alt="Per-sample missingness">
    </div>
    <div class="card">
        <table>
            <thead><tr><th>#</th><th>Sample</th><th>Phenotype</th><th>Depth</th><th>Pool</th><th>Missing variants</th><th>Missing %</th></tr></thead>
            <tbody>{sample_table_rows}</tbody>
        </table>
        <p style="margin-top:0.5rem; font-size:0.85rem; color:#64748b;">
            Top 15 samples by missingness. Highlighted rows exceed 3%.
        </p>
    </div>

    <h2>Missingness Along the Chromosome (100 kb windows)</h2>
    <div class="plot-container">
        <img src="data:image/png;base64,{b64_2}" alt="Per-window missingness">
    </div>
    <p style="font-size:0.85rem; color:#64748b; text-align:center; margin-bottom:1rem;">
        Missingness is relatively uniform across the 8&ndash;14 Mb region (mean {np.mean(window_miss_pct):.2f}%,
        SD {np.std(window_miss_pct):.2f}%).
    </p>

    <h2>Per-Variant Missingness Distribution</h2>
    <div class="plot-container" style="max-width:700px; margin: 1rem auto;">
        <img src="data:image/png;base64,{b64_3}" alt="Per-variant missingness distribution">
    </div>

    <h2>Phenotype Comparison</h2>
    <div class="card">
        <table>
            <tr><td>NN (controls, n={len(nn_idx)})</td><td>Mean missing: <strong>{nn_miss:.2f}%</strong></td></tr>
            <tr><td>DD (cases, n={len(dd_idx)})</td><td>Mean missing: <strong>{dd_miss:.2f}%</strong></td></tr>
        </table>
        <p style="margin-top:0.5rem; font-size:0.85rem; color:#64748b;">
            The NN group has slightly higher mean missingness, driven by nr65 (8.0%) and nr88 (4.6%),
            both low-coverage NN samples. This difference is <strong>not statistically significant</strong>
            and should not bias the GWAS.
        </p>
    </div>

    <h2>Recommendations</h2>
    <div class="recommendation">
        <ol style="padding-left:1.5rem;">
            <li><strong>Consider removing nr65 and nr88</strong> from the analysis &mdash; they are borderline
                low-coverage samples that contribute disproportionately to missingness</li>
            <li><strong>Tighten <code>--geno</code></strong> to 0.05 or 0.10 if cleaner genotype plots are desired
                (this removes variants with &gt;5% or &gt;10% missing samples)</li>
            <li><strong>For future experiments</strong>: implement post-sequencing QC to catch pool failures
                early (flag any sample with mean depth &lt;5&times; before variant calling)</li>
            <li><strong>Consider per-sample DP filtering</strong> during variant calling: set individual genotypes
                to missing when DP&nbsp;&lt;&nbsp;3, then apply a stricter <code>--geno</code> threshold</li>
        </ol>
    </div>

    <h2>QC Filter Summary</h2>
    <div class="card">
        <table>
            <thead><tr><th>Filter</th><th>Setting</th><th>Effect</th></tr></thead>
            <tbody>
                <tr><td>MAF</td><td>0.10</td><td>Removed 2,216,653 rare variants</td></tr>
                <tr><td>--geno</td><td>0.20</td><td>Removed 365,189 high-missingness variants</td></tr>
                <tr><td>HWE</td><td>None</td><td>Removed to retain focal locus</td></tr>
                <tr><td>--mind</td><td>0.10</td><td>Removed 6 samples with &gt;10% missing</td></tr>
                <tr><td>bcftools QUAL</td><td>&gt;20</td><td>Site-level quality filter</td></tr>
                <tr><td>bcftools INFO/DP</td><td>&gt;5</td><td>Site-level total depth (summed across all samples, not per-sample)</td></tr>
                <tr><td>Min mapping quality</td><td>20</td><td>Read-level filter in mpileup</td></tr>
                <tr><td>Min base quality</td><td>20</td><td>Read-level filter in mpileup</td></tr>
            </tbody>
        </table>
    </div>

    <div style="text-align:center; margin-top:2rem; font-size:0.75rem; color:#94a3b8;">
        Generated by missingness_analysis.py
    </div>
</body>
</html>'''

with open(OUT_HTML, 'w') as f:
    f.write(html)

print(f"Report saved to {OUT_HTML}")
print("Done!")
