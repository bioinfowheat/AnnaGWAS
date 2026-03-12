#!/usr/bin/env python3
"""
Genotype Plot — Chromosome 3 — Three regions
Inspired by Paris et al. (2022) Fig. 5a and the genotype_plot R package.

Generates genotype heatmaps for three consecutive 2 Mb windows on chr3,
centered on the most significant GWAS region, with flanking regions
2 Mb upstream and downstream.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import base64
import io
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BED_FILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/plink/variants_qc_samples_maf01geno02.bed"
BIM_FILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/plink/variants_qc_samples_maf01geno02.bim"
FAM_FILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/plink/variants_qc_samples_maf01geno02.fam"
PHENO_FILE = "phenotypes.txt"
GWAS_FILE = "association_results.PHENO1.glm.logistic.hybrid"

FOCAL_CHR = "3"
CENTER = 10_955_000
WINDOW = 1_000_000
SHIFT = 2_000_000

# Three regions: upstream, focal, downstream
REGIONS = [
    {"label": "Upstream",   "center": CENTER - SHIFT},
    {"label": "Focal",      "center": CENTER},
    {"label": "Downstream", "center": CENTER + SHIFT},
]
for r in REGIONS:
    r["start"] = r["center"] - WINDOW
    r["end"]   = r["center"] + WINDOW

MAX_SNPS = 3000  # per region

OUT_HTML = "genotype_plot.html"
OUT_PNG  = "genotype_plot_chr3.png"

# Paris et al. color scheme
COL_HOM_REF = "#FCD225"
COL_HET     = "#C92D59"
COL_HOM_ALT = "#300060"
COL_MISSING = "#CCCCCC"

# ---------------------------------------------------------------------------
# 1. Read .fam
# ---------------------------------------------------------------------------
print("Reading .fam file...")
samples = []
with open(FAM_FILE) as f:
    for line in f:
        samples.append(line.strip().split()[1])
n_samples = len(samples)
print(f"  {n_samples} samples")

# ---------------------------------------------------------------------------
# 2. Read phenotypes
# ---------------------------------------------------------------------------
print("Reading phenotypes...")
pheno_map = {}
with open(PHENO_FILE) as f:
    f.readline()
    for line in f:
        parts = line.strip().split()
        pheno_map[parts[1]] = int(parts[2])

pheno_labels = {1: "NN", 2: "DD"}
sample_phenos = [pheno_map.get(s, 0) for s in samples]

# ---------------------------------------------------------------------------
# 3. Read .bim — collect variant indices for all three regions
# ---------------------------------------------------------------------------
# We also need the full span for GWAS p-values
full_start = min(r["start"] for r in REGIONS)
full_end   = max(r["end"]   for r in REGIONS)

print(f"Reading .bim for chr {FOCAL_CHR}, full span {full_start/1e6:.1f}–{full_end/1e6:.1f} Mb...")

# Store per-region
region_data = [{} for _ in REGIONS]
for rd in region_data:
    rd["pos"] = []
    rd["idx"] = []

with open(BIM_FILE) as f:
    for i, line in enumerate(f):
        parts = line.strip().split('\t')
        if parts[0] != FOCAL_CHR:
            continue
        pos = int(parts[3])
        for j, r in enumerate(REGIONS):
            if r["start"] <= pos <= r["end"]:
                region_data[j]["pos"].append(pos)
                region_data[j]["idx"].append(i)

for j, r in enumerate(REGIONS):
    n = len(region_data[j]["idx"])
    print(f"  {r['label']}: {n} variants in {r['start']/1e6:.1f}–{r['end']/1e6:.1f} Mb")
    # Downsample
    if n > MAX_SNPS:
        step = n / MAX_SNPS
        keep = [int(k * step) for k in range(MAX_SNPS)]
        region_data[j]["pos"] = [region_data[j]["pos"][k] for k in keep]
        region_data[j]["idx"] = [region_data[j]["idx"][k] for k in keep]
        print(f"    Downsampled to {len(region_data[j]['idx'])} SNPs")

# ---------------------------------------------------------------------------
# 4. Read genotypes from .bed for all regions
# ---------------------------------------------------------------------------
print("Reading genotypes from .bed file...")

bytes_per_snp = (n_samples + 3) // 4
decode = np.array([0, 1, -1, 2], dtype=np.int8)

def read_genotypes(bed_path, variant_indices, n_samp):
    n_var = len(variant_indices)
    geno = np.zeros((n_var, n_samp), dtype=np.int8)
    bps = (n_samp + 3) // 4
    with open(bed_path, 'rb') as f:
        magic = f.read(3)
        assert magic == b'\x6c\x1b\x01'
        for local_i, global_i in enumerate(variant_indices):
            f.seek(3 + global_i * bps)
            raw = f.read(bps)
            si = 0
            for bv in raw:
                for bp in range(4):
                    if si >= n_samp:
                        break
                    geno[local_i, si] = decode[(bv >> (2 * bp)) & 0x03]
                    si += 1
            if (local_i + 1) % 1000 == 0:
                print(f"    {local_i + 1}/{n_var}...")
    return geno

for j, r in enumerate(REGIONS):
    print(f"  Reading {r['label']}...")
    region_data[j]["genotypes"] = read_genotypes(
        BED_FILE, region_data[j]["idx"], n_samples
    )
    print(f"    Shape: {region_data[j]['genotypes'].shape}")

# ---------------------------------------------------------------------------
# 5. Cluster samples using FOCAL region genotypes
# ---------------------------------------------------------------------------
print("Clustering samples on focal region genotypes...")

focal_geno = region_data[1]["genotypes"]  # focal = index 1

sample_order = []
for pheno_val in [1, 2]:
    idx = [i for i, p in enumerate(sample_phenos) if p == pheno_val]
    if len(idx) < 3:
        sample_order.extend(idx)
        continue

    sub_geno = focal_geno[:, idx].T.astype(float)
    for col in range(sub_geno.shape[1]):
        mask = sub_geno[:, col] < 0
        if mask.any():
            mean_val = sub_geno[~mask, col].mean() if (~mask).any() else 0
            sub_geno[mask, col] = mean_val

    dist = pdist(sub_geno, metric='hamming')
    Z = linkage(dist, method='ward')
    order = leaves_list(Z)
    sample_order.extend([idx[o] for o in order])

nn_count = sum(1 for p in sample_phenos if p == 1)
dd_count = sum(1 for p in sample_phenos if p == 2)

# ---------------------------------------------------------------------------
# 6. Read GWAS p-values for full span
# ---------------------------------------------------------------------------
print("Reading GWAS p-values...")
gwas_positions = {j: [] for j in range(3)}
gwas_logps     = {j: [] for j in range(3)}

with open(GWAS_FILE) as f:
    header = f.readline().strip().split('\t')
    p_idx = header.index('P')
    for line in f:
        fields = line.strip().split('\t')
        if fields[0] != FOCAL_CHR:
            continue
        pos = int(fields[1])
        for j, r in enumerate(REGIONS):
            if r["start"] <= pos <= r["end"]:
                try:
                    p = float(fields[p_idx])
                    if 0 < p <= 1:
                        gwas_positions[j].append(pos)
                        gwas_logps[j].append(-np.log10(p))
                except ValueError:
                    pass

for j in range(3):
    gwas_positions[j] = np.array(gwas_positions[j])
    gwas_logps[j] = np.array(gwas_logps[j])

# ---------------------------------------------------------------------------
# 7. Create the three-panel genotype plot
# ---------------------------------------------------------------------------
print("Generating three-panel genotype plot...")

cmap = ListedColormap([COL_HOM_REF, COL_HET, COL_HOM_ALT, COL_MISSING])
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
norm = BoundaryNorm(bounds, cmap.N)

fig, axes = plt.subplots(
    6, 1, figsize=(16, 28),
    gridspec_kw={
        'height_ratios': [1, 4, 0.3, 1, 4, 0.3],
        'hspace': 0.05
    }
)

# axes: [gwas0, geno0, spacer0, gwas1, geno1, spacer1] — but we need 3 pairs
# Actually let's do 3×2 = 6 rows + 2 spacers = 8 rows. Simpler: use subplots directly.
plt.close(fig)

fig = plt.figure(figsize=(16, 30))
import matplotlib.gridspec as gridspec

gs = gridspec.GridSpec(
    9, 1,
    height_ratios=[1, 4.5, 0.6, 1, 4.5, 0.6, 1, 4.5, 0.2],
    hspace=0.08
)

region_labels = [
    f"Upstream — {REGIONS[0]['start']/1e6:.1f}–{REGIONS[0]['end']/1e6:.1f} Mb",
    f"Focal — {REGIONS[1]['start']/1e6:.1f}–{REGIONS[1]['end']/1e6:.1f} Mb (most significant)",
    f"Downstream — {REGIONS[2]['start']/1e6:.1f}–{REGIONS[2]['end']/1e6:.1f} Mb",
]

b64_plots = []

for j in range(3):
    ax_gwas = fig.add_subplot(gs[j * 3])
    ax_geno = fig.add_subplot(gs[j * 3 + 1], sharex=ax_gwas)

    geno = region_data[j]["genotypes"][:, sample_order]
    positions = np.array(region_data[j]["pos"])

    # GWAS panel
    if len(gwas_positions[j]) > 0:
        ax_gwas.scatter(gwas_positions[j] / 1e6, gwas_logps[j],
                        s=1.5, alpha=0.5, c='#2166ac', edgecolors='none', rasterized=True)
    ax_gwas.set_ylabel(r'$-\log_{10}(p)$', fontsize=10)
    ax_gwas.set_title(region_labels[j], fontsize=12, fontweight='bold', pad=8)
    ax_gwas.set_ylim(0, 2)
    ax_gwas.spines['top'].set_visible(False)
    ax_gwas.spines['right'].set_visible(False)
    ax_gwas.tick_params(axis='x', labelbottom=False)

    # Genotype heatmap
    plot_data = geno.T.copy().astype(float)
    plot_data[plot_data == -1] = 3

    if len(positions) > 0:
        ax_geno.imshow(
            plot_data, aspect='auto', cmap=cmap, norm=norm,
            interpolation='none',
            extent=[positions[0] / 1e6, positions[-1] / 1e6, n_samples - 0.5, -0.5]
        )

    ax_geno.axhline(y=nn_count - 0.5, color='white', linewidth=2)
    ax_geno.set_yticks([])

    # Phenotype labels
    ax_geno.annotate('NN', xy=(1.01, 1 - (nn_count / 2) / n_samples),
                     xycoords='axes fraction', fontsize=11, fontweight='bold',
                     color='#1a5276', va='center', clip_on=False)
    ax_geno.annotate('DD', xy=(1.01, 1 - (nn_count + dd_count / 2) / n_samples),
                     xycoords='axes fraction', fontsize=11, fontweight='bold',
                     color='#922b21', va='center', clip_on=False)

    if j == 2:
        ax_geno.set_xlabel('Position (Mb)', fontsize=12)
    else:
        ax_geno.tick_params(axis='x', labelbottom=True)

    ax_geno.set_ylabel('Individuals', fontsize=10)

# Legend at bottom of last panel
legend_elements = [
    Patch(facecolor=COL_HOM_REF, edgecolor='grey', label='HOM REF'),
    Patch(facecolor=COL_HET, edgecolor='grey', label='HET'),
    Patch(facecolor=COL_HOM_ALT, edgecolor='grey', label='HOM ALT'),
    Patch(facecolor=COL_MISSING, edgecolor='grey', label='Missing'),
]
fig.legend(handles=legend_elements, loc='lower center', fontsize=10,
           ncol=4, framealpha=0.9, edgecolor='#ccc',
           bbox_to_anchor=(0.5, 0.01))

# Overall title
fig.suptitle(f'Genotype Plot — Chromosome {FOCAL_CHR}', fontsize=16,
             fontweight='bold', y=0.995)

# Save PNG
fig.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
print(f"Saved {OUT_PNG}")

# Save to buffer for HTML
buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=250, bbox_inches='tight')
buf.seek(0)
b64_main = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# ---------------------------------------------------------------------------
# 8. Summary statistics
# ---------------------------------------------------------------------------
all_geno = np.concatenate([region_data[j]["genotypes"] for j in range(3)], axis=0)
n_missing = (all_geno == -1).sum()
missing_pct = n_missing / all_geno.size * 100
total_snps = sum(len(region_data[j]["pos"]) for j in range(3))

# Top GWAS hits per region
def top_hits_table(gpos, glogp, n=5):
    if len(gpos) == 0:
        return "<tr><td colspan='3'>No variants</td></tr>"
    top = np.argsort(glogp)[::-1][:n]
    rows = ""
    for i in top:
        rows += f"<tr><td>{int(gpos[i]):,}</td><td>{10**(-glogp[i]):.4e}</td><td>{glogp[i]:.3f}</td></tr>\n"
    return rows

# ---------------------------------------------------------------------------
# 9. Generate HTML
# ---------------------------------------------------------------------------
print("Generating HTML report...")

region_summaries = ""
for j, r in enumerate(REGIONS):
    n_var = len(region_data[j]["pos"])
    n_gwas = len(gwas_positions[j])
    region_summaries += f"""
    <tr>
        <td><strong>{r['label']}</strong></td>
        <td>{r['start']/1e6:.2f} – {r['end']/1e6:.2f} Mb</td>
        <td>{n_var:,} displayed (of {n_gwas:,})</td>
        <td>{gwas_logps[j].max():.3f}</td>
    </tr>"""

top_hits_sections = ""
for j, r in enumerate(REGIONS):
    top_hits_sections += f"""
    <h3>{r['label']} — {r['start']/1e6:.1f}–{r['end']/1e6:.1f} Mb</h3>
    <div class="card">
        <table>
            <thead><tr><th>Position (bp)</th><th>P-value</th><th>&minus;log<sub>10</sub>(p)</th></tr></thead>
            <tbody>{top_hits_table(gwas_positions[j], gwas_logps[j])}</tbody>
        </table>
    </div>"""

html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Genotype Plot — Chromosome {FOCAL_CHR}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
            background: #f8fafc; color: #1e293b; line-height: 1.6;
            max-width: 1600px; margin: 0 auto; padding: 2rem;
        }}
        h1 {{ font-size: 1.6rem; margin-bottom: 0.3rem; color: #0f172a; }}
        .subtitle {{ color: #64748b; margin-bottom: 2rem; font-size: 0.9rem; }}
        h2 {{ font-size: 1.15rem; margin: 2rem 0 0.8rem; color: #0f172a;
              border-bottom: 2px solid #e2e8f0; padding-bottom: 0.3rem; }}
        h3 {{ font-size: 1rem; margin: 1.2rem 0 0.5rem; color: #334155; }}
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
        .method {{ background: #fffbeb; border-left: 4px solid #f59e0b; padding: 1rem; margin: 1rem 0; border-radius: 0 6px 6px 0; font-size: 0.9rem; }}
        .swatch {{ display: inline-block; width: 14px; height: 14px; border-radius: 2px; vertical-align: middle; margin-right: 4px; border: 1px solid #ccc; }}
    </style>
</head>
<body>
    <h1>Genotype Plot &mdash; Chromosome {FOCAL_CHR}</h1>
    <div class="subtitle">
        <em>Pararge aegeria</em> &mdash; {n_samples} samples ({nn_count} NN &times; {dd_count} DD) &mdash;
        Three 2&nbsp;Mb windows centered &plusmn;2&nbsp;Mb around the most significant GWAS region ({CENTER/1e6:.2f} Mb)
    </div>

    <h2>Method</h2>
    <div class="method">
        <strong>Genotype plot</strong> following
        <a href="https://doi.org/10.1038/s41467-022-28895-4">Paris et al. (2022)</a> and the
        <a href="https://github.com/JimWhiting91/genotype_plot">genotype_plot</a> R package.
        Each column is a SNP position; each row is an individual.
        Genotypes are colour-coded:
        <span class="swatch" style="background:{COL_HOM_REF};"></span>HOM&nbsp;REF &nbsp;
        <span class="swatch" style="background:{COL_HET};"></span>HET &nbsp;
        <span class="swatch" style="background:{COL_HOM_ALT};"></span>HOM&nbsp;ALT &nbsp;
        <span class="swatch" style="background:{COL_MISSING};"></span>Missing.
        <br>
        Samples are grouped by phenotype (NN / DD) and hierarchically clustered within groups
        using Hamming distance and Ward&rsquo;s method on the <strong>focal region</strong> genotypes.
        The same sample ordering is used across all three panels for direct comparison.
        A GWAS &minus;log<sub>10</sub>(<em>p</em>) panel is shown above each heatmap for positional context.
    </div>

    <h2>Region Summary</h2>
    <div class="card">
        <table>
            <thead><tr><th>Region</th><th>Coordinates</th><th>Variants</th><th>Max &minus;log<sub>10</sub>(p)</th></tr></thead>
            <tbody>{region_summaries}</tbody>
        </table>
        <p style="margin-top:0.5rem; font-size:0.85rem; color:#64748b;">
            Total: {total_snps:,} SNPs displayed &bull; {n_samples} samples &bull; {missing_pct:.2f}% missing genotypes
        </p>
    </div>

    <h2>Genotype Plots</h2>
    <div class="plot-container">
        <img src="data:image/png;base64,{b64_main}" alt="Genotype plot chromosome {FOCAL_CHR} three regions">
    </div>

    <h2>Top GWAS Hits per Region</h2>
    {top_hits_sections}

    <div style="text-align:center; margin-top:2rem; font-size:0.75rem; color:#94a3b8;">
        Generated from PLINK binary files &mdash; genotype_plot.py
    </div>
</body>
</html>'''

with open(OUT_HTML, 'w') as f:
    f.write(html)

print(f"HTML report saved to {OUT_HTML}")
print("Done!")
