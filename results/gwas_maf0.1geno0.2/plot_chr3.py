#!/usr/bin/env python3
"""Generate chromosome 3 Manhattan plot and embed as base64 in summary_report.html"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import base64
import io
import re

gwas_file = "association_results.PHENO1.glm.logistic.hybrid"

# Read chr3 data
positions = []
pvalues = []

with open(gwas_file) as f:
    header = f.readline().strip().split('\t')
    chrom_idx = 0   # #CHROM
    pos_idx = 1     # POS
    p_idx = header.index('P')

    for line in f:
        fields = line.strip().split('\t')
        if fields[chrom_idx] != '3':
            continue
        try:
            p = float(fields[p_idx])
            if p > 0 and p <= 1:
                positions.append(int(fields[pos_idx]))
                pvalues.append(p)
        except (ValueError, IndexError):
            continue

positions = np.array(positions)
pvalues = np.array(pvalues)
logp = -np.log10(pvalues)

print(f"Chromosome 3 variants: {len(positions)}")
print(f"Min p-value: {pvalues.min():.2e}")
print(f"Suggestive hits (p < 1e-5): {(pvalues < 1e-5).sum()}")
print(f"Genome-wide significant (p < 5e-8): {(pvalues < 5e-8).sum()}")

# Top 10 hits
top_idx = np.argsort(pvalues)[:10]
print("\nTop 10 associations on chromosome 3:")
for i in top_idx:
    print(f"  Pos: {positions[i]:>12,}  P: {pvalues[i]:.2e}")

# Create Manhattan plot
fig, ax = plt.subplots(figsize=(14, 5), dpi=300)

ax.scatter(positions / 1e6, logp, s=2, alpha=0.6, c='#2166ac',
           edgecolors='none', rasterized=True)

ax.set_xlabel('Position (Mb)', fontsize=12)
ax.set_ylabel(r'$-\log_{10}(p)$', fontsize=12)
ax.set_title('Manhattan Plot — Chromosome 3', fontsize=14, fontweight='bold')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(positions.min() / 1e6 - 0.5, positions.max() / 1e6 + 0.5)
ax.set_ylim(0, 2)

plt.tight_layout()

# Save to PNG buffer
buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
buf.seek(0)
b64 = base64.b64encode(buf.read()).decode('utf-8')
plt.close(fig)

# Also save as file
fig2, ax2 = plt.subplots(figsize=(14, 5), dpi=300)
ax2.scatter(positions / 1e6, logp, s=2, alpha=0.6, c='#2166ac',
            edgecolors='none', rasterized=True)
ax2.set_xlabel('Position (Mb)', fontsize=12)
ax2.set_ylabel(r'$-\log_{10}(p)$', fontsize=12)
ax2.set_title('Manhattan Plot — Chromosome 3', fontsize=14, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlim(positions.min() / 1e6 - 0.5, positions.max() / 1e6 + 0.5)
ax2.set_ylim(0, 2)
plt.tight_layout()
fig2.savefig('manhattan_chr3.png', dpi=300, bbox_inches='tight')
plt.close(fig2)
print("\nSaved manhattan_chr3.png")

# Insert/replace in HTML
html_file = "summary_report.html"
with open(html_file, 'r') as f:
    html = f.read()

chr3_section = f'''    <h2>Manhattan Plot &mdash; Chromosome 3</h2>
    <div class="plot-container">
        <img src="data:image/png;base64,{b64}" alt="Manhattan plot chromosome 3">
    </div>'''

# Replace existing chr3 section if present, otherwise insert before QQ Plot
import re as re_mod
pattern = r'    <h2>Manhattan Plot &mdash; Chromosome 3</h2>\s*<div class="plot-container">\s*<img src="data:image/png;base64,[A-Za-z0-9+/=\s]+" alt="Manhattan plot chromosome 3">\s*</div>'
if re_mod.search(pattern, html):
    html = re_mod.sub(pattern, chr3_section, html)
    print("Replaced existing chr3 section in HTML")
else:
    insertion_point = '    <h2>QQ Plot'
    html = html.replace(insertion_point, chr3_section + '\n\n' + insertion_point)
    print("Inserted new chr3 section in HTML")

with open(html_file, 'w') as f:
    f.write(html)
print(f"Chromosome 3 plot embedded in {html_file}")
