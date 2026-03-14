#!/usr/bin/env python3
"""
Create Annas_phenotypes.tsv and phenotype_summary.html from the xlsx phenotype data.
Merges with existing sample IDs from phenotypes.txt.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import base64
from io import BytesIO
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ── Read data ──
xlsx = pd.read_excel('qc_anna_data_phenotypes.xlsx')
pheno_orig = pd.read_csv('phenotypes.txt', sep='\t')

# Build sample_id from nr column: "nr1", "nr2", etc.
xlsx['sample_id'] = 'nr' + xlsx['nr'].astype(str)

# Verify the sample IDs match
assert set(xlsx['sample_id']) == set(pheno_orig['sample_id']), "Sample ID mismatch!"

# ── Build the phenotypes TSV ──
# Phenotype columns for GWAS (PLINK2 compatible):
#   - Binary traits: 0/1 coding (use --1 flag or --glm cols= to select)
#   - Quantitative traits: numeric, NA for missing
#   - Covariates: numeric coding

out = pd.DataFrame()
out['sample_id'] = xlsx['sample_id']
out['cross'] = xlsx['cross']

# Original diapause phenotype (DD=1, NN=0)
out['diapause'] = (xlsx['cross'] == 'DD').astype(int)

# Binary developmental decisions
out['larval_decision'] = xlsx['Larval decision'].astype(int)
out['pupal_decision'] = xlsx['Pupal decision'].astype(int)

# Sex (M=1, F=2 — PLINK convention)
out['sex'] = xlsx['Sex'].map({'M': 1, 'F': 2})

# Continuous traits
out['pupal_weight_g'] = xlsx['P. Weight']
out['days_to_pupation'] = xlsx['Days until pupation (not exact)']
out['adult_weight_g'] = xlsx['Adult weight']
out['days_to_eclosion'] = xlsx['Days until eclosion']

# Environmental covariate
out['cabinet'] = xlsx['cabinet'].map({'C1': 1, 'C2': 2, 'C3': 3})

# Sort by nr
out = out.sort_values('sample_id', key=lambda x: x.str.extract(r'(\d+)')[0].astype(int))

out.to_csv('Annas_phenotypes.tsv', sep='\t', index=False, na_rep='NA')
print(f"Created Annas_phenotypes.tsv with {len(out)} samples and {len(out.columns)} columns")

# ── Generate HTML report with embedded plots ──

def fig_to_base64(fig, dpi=150):
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return b64

plots = {}

# Color scheme
DD_COLOR = '#e74c3c'
NN_COLOR = '#3498db'
COLORS = [NN_COLOR, DD_COLOR]

dd = out[out['cross'] == 'DD']
nn = out[out['cross'] == 'NN']

# ── Plot 1: Larval & Pupal decisions by cross ──
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Larval decision
ct_larval = pd.crosstab(out['cross'], out['larval_decision'], normalize='index') * 100
ct_larval.plot(kind='bar', ax=axes[0], color=['#2ecc71', '#e67e22'], edgecolor='black', width=0.6)
axes[0].set_title('Larval Decision by Cross Type', fontsize=14, fontweight='bold')
axes[0].set_xlabel('')
axes[0].set_ylabel('Percentage of individuals (%)')
axes[0].legend(['Direct development (0)', 'Diapause (1)'], loc='upper right')
axes[0].set_xticklabels(['DD', 'NN'], rotation=0)
axes[0].set_ylim(0, 110)
for p in axes[0].patches:
    h = p.get_height()
    if h > 3:
        axes[0].annotate(f'{h:.0f}%', (p.get_x() + p.get_width()/2, h), ha='center', va='bottom', fontsize=11)

# Pupal decision
ct_pupal = pd.crosstab(out['cross'], out['pupal_decision'], normalize='index') * 100
ct_pupal.plot(kind='bar', ax=axes[1], color=['#2ecc71', '#e67e22'], edgecolor='black', width=0.6)
axes[1].set_title('Pupal Decision by Cross Type', fontsize=14, fontweight='bold')
axes[1].set_xlabel('')
axes[1].set_ylabel('Percentage of individuals (%)')
axes[1].legend(['Direct development (0)', 'Diapause (1)'], loc='upper right')
axes[1].set_xticklabels(['DD', 'NN'], rotation=0)
axes[1].set_ylim(0, 110)
for p in axes[1].patches:
    h = p.get_height()
    if h > 3:
        axes[1].annotate(f'{h:.0f}%', (p.get_x() + p.get_width()/2, h), ha='center', va='bottom', fontsize=11)

fig.suptitle('Binary Developmental Decisions', fontsize=16, fontweight='bold', y=1.02)
fig.tight_layout()
plots['decisions'] = fig_to_base64(fig)

# ── Plot 2: Larval decision doesn't match cross ──
fig, ax = plt.subplots(figsize=(8, 5))
# Create a combined category
out['decision_match'] = 'Match'
out.loc[(out['cross'] == 'DD') & (out['larval_decision'] == 0), 'decision_match'] = 'DD but direct dev.'
out.loc[(out['cross'] == 'NN') & (out['larval_decision'] == 1), 'decision_match'] = 'NN but diapause'

match_ct = pd.crosstab(out['cross'], out['decision_match'])
match_ct.plot(kind='bar', stacked=True, ax=ax,
              color=['#e74c3c', '#2ecc71', '#3498db'], edgecolor='black', width=0.6)
ax.set_title('Cross Type vs Actual Larval Decision', fontsize=14, fontweight='bold')
ax.set_xlabel('')
ax.set_ylabel('Number of individuals')
ax.set_xticklabels(['DD', 'NN'], rotation=0)
ax.legend(title='Category', bbox_to_anchor=(1.02, 1), loc='upper left')
fig.tight_layout()
plots['mismatch'] = fig_to_base64(fig)

# ── Plot 3: Days to pupation ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Histogram
bins = np.arange(15, 130, 5)
axes[0].hist(nn['days_to_pupation'], bins=bins, alpha=0.7, color=NN_COLOR, label='NN', edgecolor='black')
axes[0].hist(dd['days_to_pupation'], bins=bins, alpha=0.7, color=DD_COLOR, label='DD', edgecolor='black')
axes[0].set_xlabel('Days until pupation')
axes[0].set_ylabel('Count')
axes[0].set_title('Days to Pupation Distribution', fontsize=14, fontweight='bold')
axes[0].legend()
axes[0].axvline(nn['days_to_pupation'].median(), color=NN_COLOR, ls='--', lw=2)
axes[0].axvline(dd['days_to_pupation'].median(), color=DD_COLOR, ls='--', lw=2)

# Box plot
bp = axes[1].boxplot([nn['days_to_pupation'], dd['days_to_pupation']],
                     labels=['NN', 'DD'], patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor(NN_COLOR)
bp['boxes'][1].set_facecolor(DD_COLOR)
for b in bp['boxes']:
    b.set_alpha(0.7)
axes[1].set_ylabel('Days until pupation')
axes[1].set_title('Days to Pupation by Cross', fontsize=14, fontweight='bold')

# Add individual points with jitter
for i, (data, color) in enumerate([(nn['days_to_pupation'], NN_COLOR), (dd['days_to_pupation'], DD_COLOR)]):
    x = np.random.normal(i+1, 0.06, size=len(data))
    axes[1].scatter(x, data, alpha=0.5, color=color, edgecolors='black', linewidth=0.5, s=30, zorder=5)

fig.tight_layout()
plots['pupation'] = fig_to_base64(fig)

# ── Plot 4: Days to eclosion ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

bins = np.arange(8, 35, 1)
axes[0].hist(nn['days_to_eclosion'], bins=bins, alpha=0.7, color=NN_COLOR, label='NN', edgecolor='black')
axes[0].hist(dd['days_to_eclosion'], bins=bins, alpha=0.7, color=DD_COLOR, label='DD', edgecolor='black')
axes[0].set_xlabel('Days until eclosion')
axes[0].set_ylabel('Count')
axes[0].set_title('Days to Eclosion Distribution', fontsize=14, fontweight='bold')
axes[0].legend()
axes[0].axvline(nn['days_to_eclosion'].median(), color=NN_COLOR, ls='--', lw=2)
axes[0].axvline(dd['days_to_eclosion'].median(), color=DD_COLOR, ls='--', lw=2)

bp = axes[1].boxplot([nn['days_to_eclosion'], dd['days_to_eclosion']],
                     labels=['NN', 'DD'], patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor(NN_COLOR)
bp['boxes'][1].set_facecolor(DD_COLOR)
for b in bp['boxes']:
    b.set_alpha(0.7)
axes[1].set_ylabel('Days until eclosion')
axes[1].set_title('Days to Eclosion by Cross', fontsize=14, fontweight='bold')

for i, (data, color) in enumerate([(nn['days_to_eclosion'], NN_COLOR), (dd['days_to_eclosion'], DD_COLOR)]):
    x = np.random.normal(i+1, 0.06, size=len(data))
    axes[1].scatter(x, data, alpha=0.5, color=color, edgecolors='black', linewidth=0.5, s=30, zorder=5)

fig.tight_layout()
plots['eclosion'] = fig_to_base64(fig)

# ── Plot 5: Pupal weight ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

bins = np.arange(0.07, 0.21, 0.005)
axes[0].hist(nn['pupal_weight_g'], bins=bins, alpha=0.7, color=NN_COLOR, label='NN', edgecolor='black')
axes[0].hist(dd['pupal_weight_g'], bins=bins, alpha=0.7, color=DD_COLOR, label='DD', edgecolor='black')
axes[0].set_xlabel('Pupal weight (g)')
axes[0].set_ylabel('Count')
axes[0].set_title('Pupal Weight Distribution', fontsize=14, fontweight='bold')
axes[0].legend()

bp = axes[1].boxplot([nn['pupal_weight_g'], dd['pupal_weight_g']],
                     labels=['NN', 'DD'], patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor(NN_COLOR)
bp['boxes'][1].set_facecolor(DD_COLOR)
for b in bp['boxes']:
    b.set_alpha(0.7)
axes[1].set_ylabel('Pupal weight (g)')
axes[1].set_title('Pupal Weight by Cross', fontsize=14, fontweight='bold')

for i, (data, color) in enumerate([(nn['pupal_weight_g'], NN_COLOR), (dd['pupal_weight_g'], DD_COLOR)]):
    x = np.random.normal(i+1, 0.06, size=len(data))
    axes[1].scatter(x, data, alpha=0.5, color=color, edgecolors='black', linewidth=0.5, s=30, zorder=5)

fig.tight_layout()
plots['pupal_wt'] = fig_to_base64(fig)

# ── Plot 6: Adult weight ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

nn_aw = nn['adult_weight_g'].dropna()
dd_aw = dd['adult_weight_g'].dropna()

bins = np.arange(0.02, 0.12, 0.005)
axes[0].hist(nn_aw, bins=bins, alpha=0.7, color=NN_COLOR, label='NN', edgecolor='black')
axes[0].hist(dd_aw, bins=bins, alpha=0.7, color=DD_COLOR, label='DD', edgecolor='black')
axes[0].set_xlabel('Adult weight (g)')
axes[0].set_ylabel('Count')
axes[0].set_title('Adult Weight Distribution', fontsize=14, fontweight='bold')
axes[0].legend()

bp = axes[1].boxplot([nn_aw, dd_aw], labels=['NN', 'DD'], patch_artist=True, widths=0.5)
bp['boxes'][0].set_facecolor(NN_COLOR)
bp['boxes'][1].set_facecolor(DD_COLOR)
for b in bp['boxes']:
    b.set_alpha(0.7)
axes[1].set_ylabel('Adult weight (g)')
axes[1].set_title('Adult Weight by Cross', fontsize=14, fontweight='bold')

for i, (data, color) in enumerate([(nn_aw, NN_COLOR), (dd_aw, DD_COLOR)]):
    x = np.random.normal(i+1, 0.06, size=len(data))
    axes[1].scatter(x, data, alpha=0.5, color=color, edgecolors='black', linewidth=0.5, s=30, zorder=5)

fig.tight_layout()
plots['adult_wt'] = fig_to_base64(fig)

# ── Plot 7: Sex distribution ──
fig, ax = plt.subplots(figsize=(8, 5))
ct_sex = pd.crosstab(out['cross'], out['sex'].map({1: 'Male', 2: 'Female'}))
ct_sex.plot(kind='bar', ax=ax, color=['#9b59b6', '#1abc9c'], edgecolor='black', width=0.6)
ax.set_title('Sex Distribution by Cross Type', fontsize=14, fontweight='bold')
ax.set_xlabel('')
ax.set_ylabel('Count')
ax.set_xticklabels(['DD', 'NN'], rotation=0)
ax.legend(title='Sex')
for p in ax.patches:
    h = p.get_height()
    if h > 0:
        ax.annotate(f'{int(h)}', (p.get_x() + p.get_width()/2, h), ha='center', va='bottom', fontsize=12)
fig.tight_layout()
plots['sex'] = fig_to_base64(fig)

# ── Plot 8: Correlation heatmap of continuous traits ──
fig, ax = plt.subplots(figsize=(8, 7))
cont_cols = ['days_to_pupation', 'days_to_eclosion', 'pupal_weight_g', 'adult_weight_g']
corr = out[cont_cols].corr()
im = ax.imshow(corr, cmap='RdBu_r', vmin=-1, vmax=1)
ax.set_xticks(range(len(cont_cols)))
ax.set_yticks(range(len(cont_cols)))
labels = ['Days to\npupation', 'Days to\neclosion', 'Pupal\nweight', 'Adult\nweight']
ax.set_xticklabels(labels, fontsize=11)
ax.set_yticklabels(labels, fontsize=11)
for i in range(len(cont_cols)):
    for j in range(len(cont_cols)):
        val = corr.iloc[i, j]
        color = 'white' if abs(val) > 0.5 else 'black'
        ax.text(j, i, f'{val:.2f}', ha='center', va='center', color=color, fontsize=14, fontweight='bold')
plt.colorbar(im, ax=ax, shrink=0.8, label='Pearson r')
ax.set_title('Correlation Between Continuous Traits', fontsize=14, fontweight='bold')
fig.tight_layout()
plots['correlation'] = fig_to_base64(fig)

# ── Plot 9: Scatter — days to pupation vs pupal weight, colored by cross ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].scatter(nn['days_to_pupation'], nn['pupal_weight_g'], color=NN_COLOR, alpha=0.6, edgecolors='black', linewidth=0.5, s=50, label='NN')
axes[0].scatter(dd['days_to_pupation'], dd['pupal_weight_g'], color=DD_COLOR, alpha=0.6, edgecolors='black', linewidth=0.5, s=50, label='DD')
axes[0].set_xlabel('Days to pupation')
axes[0].set_ylabel('Pupal weight (g)')
axes[0].set_title('Days to Pupation vs Pupal Weight', fontsize=13, fontweight='bold')
axes[0].legend()

axes[1].scatter(nn['days_to_pupation'], nn['days_to_eclosion'], color=NN_COLOR, alpha=0.6, edgecolors='black', linewidth=0.5, s=50, label='NN')
axes[1].scatter(dd['days_to_pupation'], dd['days_to_eclosion'], color=DD_COLOR, alpha=0.6, edgecolors='black', linewidth=0.5, s=50, label='DD')
axes[1].set_xlabel('Days to pupation')
axes[1].set_ylabel('Days to eclosion')
axes[1].set_title('Days to Pupation vs Days to Eclosion', fontsize=13, fontweight='bold')
axes[1].legend()

fig.tight_layout()
plots['scatter'] = fig_to_base64(fig)

# ── Plot 10: Cabinet distribution ──
fig, ax = plt.subplots(figsize=(8, 5))
ct_cab = pd.crosstab(out['cross'], out['cabinet'].map({1: 'C1', 2: 'C2', 3: 'C3'}))
ct_cab.plot(kind='bar', ax=ax, color=['#f39c12', '#27ae60', '#8e44ad'], edgecolor='black', width=0.6)
ax.set_title('Cabinet Assignment by Cross Type', fontsize=14, fontweight='bold')
ax.set_xlabel('')
ax.set_ylabel('Count')
ax.set_xticklabels(['DD', 'NN'], rotation=0)
ax.legend(title='Cabinet')
for p in ax.patches:
    h = p.get_height()
    if h > 0:
        ax.annotate(f'{int(h)}', (p.get_x() + p.get_width()/2, h), ha='center', va='bottom', fontsize=11)
fig.tight_layout()
plots['cabinet'] = fig_to_base64(fig)


# ── Build HTML ──

# Summary statistics tables
def describe_by_cross(col_name, label):
    stats = out.groupby('cross')[col_name].describe()
    rows = ''
    for cross in ['NN', 'DD']:
        s = stats.loc[cross]
        rows += f'<tr><td>{cross}</td><td>{int(s["count"])}</td><td>{s["mean"]:.4f}</td><td>{s["std"]:.4f}</td><td>{s["min"]:.4f}</td><td>{s["50%"]:.4f}</td><td>{s["max"]:.4f}</td></tr>\n'
    return f'''<table>
<caption>{label}</caption>
<tr><th>Cross</th><th>N</th><th>Mean</th><th>Std</th><th>Min</th><th>Median</th><th>Max</th></tr>
{rows}</table>'''

# Phenotypes summary table for the TSV file
pheno_table_rows = ''
for _, row in out.head(10).iterrows():
    vals = '</td><td>'.join([str(v) for v in row.values])
    pheno_table_rows += f'<tr><td>{vals}</td></tr>\n'
header_cells = '</th><th>'.join(out.columns)

# Count mismatches
n_dd_direct = len(out[(out['cross'] == 'DD') & (out['larval_decision'] == 0)])
n_nn_diapause = len(out[(out['cross'] == 'NN') & (out['larval_decision'] == 1)])
n_pupal_dia = len(out[out['pupal_decision'] == 1])

html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Phenotype Summary &mdash; P. aegeria GWAS</title>
<style>
    body {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 1100px; margin: 0 auto; padding: 20px; background: #fafafa; color: #333; line-height: 1.6; }}
    h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
    h2 {{ color: #2c3e50; margin-top: 40px; border-left: 4px solid #3498db; padding-left: 12px; }}
    h3 {{ color: #34495e; }}
    .summary-box {{ background: #ecf0f1; border-radius: 8px; padding: 20px; margin: 20px 0; }}
    .warning-box {{ background: #fdf2e9; border-left: 4px solid #e67e22; padding: 15px 20px; margin: 20px 0; border-radius: 4px; }}
    .key-finding {{ background: #eaf2f8; border-left: 4px solid #3498db; padding: 15px 20px; margin: 20px 0; border-radius: 4px; }}
    table {{ border-collapse: collapse; margin: 15px 0; width: 100%; }}
    th, td {{ border: 1px solid #bdc3c7; padding: 8px 12px; text-align: left; }}
    th {{ background: #2c3e50; color: white; }}
    caption {{ font-weight: bold; font-size: 1.1em; margin-bottom: 8px; text-align: left; color: #2c3e50; }}
    tr:nth-child(even) {{ background: #f2f3f4; }}
    img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; margin: 10px 0; }}
    .plot-container {{ background: white; padding: 15px; border-radius: 8px; margin: 20px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
    .two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
    .metric {{ text-align: center; padding: 15px; background: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
    .metric .value {{ font-size: 2em; font-weight: bold; color: #2c3e50; }}
    .metric .label {{ font-size: 0.9em; color: #7f8c8d; }}
    code {{ background: #ecf0f1; padding: 2px 6px; border-radius: 3px; font-size: 0.9em; }}
    .tag-dd {{ background: #e74c3c; color: white; padding: 2px 8px; border-radius: 10px; font-size: 0.85em; }}
    .tag-nn {{ background: #3498db; color: white; padding: 2px 8px; border-radius: 10px; font-size: 0.85em; }}
</style>
</head>
<body>

<h1>Phenotype Summary for <em>P. aegeria</em> GWAS</h1>

<div class="summary-box">
<strong>Source file:</strong> <code>qc_anna_data_phenotypes.xlsx</code><br>
<strong>Output file:</strong> <code>Annas_phenotypes.tsv</code><br>
<strong>Samples:</strong> 94 individuals (47 <span class="tag-dd">DD</span> diapause cross, 47 <span class="tag-nn">NN</span> non-diapause cross)
</div>

<h2>1. Source Data Overview</h2>

<p>The xlsx file contains a single sheet with 94 rows (one per sample) and 16 columns. The columns cover DNA quality metrics, experimental metadata, and biological phenotypes from the diapause experiment.</p>

<h3>Columns in the source xlsx</h3>
<table>
<tr><th>Column</th><th>Type</th><th>Description</th><th>Included in TSV?</th></tr>
<tr><td>Well</td><td>&mdash;</td><td>All empty (96-well plate position, unused)</td><td>No</td></tr>
<tr><td>nr</td><td>Integer</td><td>Sample number (1&ndash;94), maps to sample_id</td><td>Yes (as sample_id)</td></tr>
<tr><td>cabinet</td><td>Categorical</td><td>Growth cabinet assignment (C1, C2, C3)</td><td>Yes (coded 1/2/3)</td></tr>
<tr><td>cross</td><td>Categorical</td><td>Cross type: DD (diapause) or NN (non-diapause)</td><td>Yes</td></tr>
<tr><td>sample</td><td>Integer</td><td>Original sample number within cross</td><td>No (redundant)</td></tr>
<tr><td>qubit</td><td>Numeric</td><td>DNA concentration (Qubit), some with "?" flags</td><td>No (QC metric)</td></tr>
<tr><td>260/280</td><td>Numeric</td><td>DNA purity ratio</td><td>No (QC metric)</td></tr>
<tr><td>gel position</td><td>Mixed</td><td>Gel position, mostly empty, some labelled "pupa"</td><td>No (QC metric)</td></tr>
<tr><td>adult</td><td>Categorical</td><td>Whether individual survived to adult (y/n)</td><td>No (mostly empty)</td></tr>
<tr><td>Larval decision</td><td>Binary</td><td>0 = direct development, 1 = diapause at larval stage</td><td><strong>Yes</strong></td></tr>
<tr><td>Pupal decision</td><td>Binary</td><td>0 = direct development, 1 = diapause at pupal stage</td><td><strong>Yes</strong></td></tr>
<tr><td>Sex</td><td>Categorical</td><td>M or F</td><td><strong>Yes</strong> (M=1, F=2)</td></tr>
<tr><td>P. Weight</td><td>Continuous</td><td>Pupal weight in grams</td><td><strong>Yes</strong></td></tr>
<tr><td>Days until pupation</td><td>Continuous</td><td>Days from experiment start to pupation</td><td><strong>Yes</strong></td></tr>
<tr><td>Adult weight</td><td>Continuous</td><td>Adult weight in grams (3 missing)</td><td><strong>Yes</strong></td></tr>
<tr><td>Days until eclosion</td><td>Continuous</td><td>Days from pupation to adult emergence</td><td><strong>Yes</strong></td></tr>
</table>

<h3>Columns excluded and why</h3>
<ul>
<li><strong>Well</strong> &mdash; entirely empty, no data</li>
<li><strong>sample</strong> &mdash; redundant with <code>nr</code> (the within-cross sample number)</li>
<li><strong>qubit, 260/280, gel position</strong> &mdash; DNA quality metrics, not biological phenotypes. Some have mixed-type data (e.g., "99.2?" with question marks)</li>
<li><strong>adult</strong> &mdash; 84 of 94 values are empty; the remaining entries include non-standard labels ("Data C3", "Data C2")</li>
</ul>

<h2>2. Phenotype File Structure</h2>

<p>The output file <code>Annas_phenotypes.tsv</code> contains {len(out.columns)} columns:</p>

<table>
<tr><th>Column</th><th>Type</th><th>Values</th><th>GWAS usage</th></tr>
<tr><td><code>sample_id</code></td><td>ID</td><td>nr1&ndash;nr94</td><td>Sample identifier</td></tr>
<tr><td><code>cross</code></td><td>Label</td><td>DD, NN</td><td>Reference label</td></tr>
<tr><td><code>diapause</code></td><td>Binary</td><td>0 (NN) / 1 (DD)</td><td>Original GWAS phenotype</td></tr>
<tr><td><code>larval_decision</code></td><td>Binary</td><td>0 / 1</td><td>Binary trait GWAS</td></tr>
<tr><td><code>pupal_decision</code></td><td>Binary</td><td>0 / 1</td><td>Binary trait GWAS</td></tr>
<tr><td><code>sex</code></td><td>Binary</td><td>1 (M) / 2 (F)</td><td>Covariate or sex GWAS</td></tr>
<tr><td><code>pupal_weight_g</code></td><td>Continuous</td><td>0.0802&ndash;0.1963 g</td><td>Quantitative GWAS</td></tr>
<tr><td><code>days_to_pupation</code></td><td>Continuous</td><td>23&ndash;119 days</td><td>Quantitative GWAS</td></tr>
<tr><td><code>adult_weight_g</code></td><td>Continuous</td><td>0.0288&ndash;0.1013 g (3 NA)</td><td>Quantitative GWAS</td></tr>
<tr><td><code>days_to_eclosion</code></td><td>Continuous</td><td>10&ndash;30 days</td><td>Quantitative GWAS</td></tr>
<tr><td><code>cabinet</code></td><td>Categorical</td><td>1 / 2 / 3</td><td>Environmental covariate</td></tr>
</table>

<h3>Preview (first 10 rows)</h3>
<div style="overflow-x:auto;">
<table>
<tr><th>{header_cells}</th></tr>
{pheno_table_rows}
</table>
</div>

<h2>3. Key Observations</h2>

<div class="warning-box">
<strong>Important:</strong> The original <code>diapause</code> phenotype (based on cross type DD/NN) does <strong>not</strong> perfectly match the actual <code>larval_decision</code> phenotype recorded for each individual:
<ul>
<li><strong>{n_dd_direct} DD individuals</strong> actually chose direct development (larval_decision = 0)</li>
<li><strong>{n_nn_diapause} NN individuals</strong> actually entered diapause (larval_decision = 1)</li>
</ul>
This means the cross type is a <em>family-level</em> label reflecting parental genotype, not the individual's own developmental decision. The <code>larval_decision</code> column is the actual individual-level phenotype and may be more informative for GWAS.
</div>

<div class="key-finding">
<strong>Pupal decision</strong> is nearly invariant: only <strong>{n_pupal_dia} individuals</strong> (all from NN crosses) entered pupal diapause. This extremely low case count (7/94) makes it poorly suited as a standalone GWAS phenotype, but it captures a biologically distinct diapause checkpoint.
</div>

<h2>4. Binary Developmental Phenotypes</h2>

<div class="plot-container">
<img src="data:image/png;base64,{plots['decisions']}" alt="Binary decisions by cross type">
<p><strong>Left:</strong> Larval decision is split nearly 50/50 in <em>both</em> DD and NN crosses &mdash; cross type alone does not determine larval fate. <strong>Right:</strong> Pupal diapause occurs only in NN crosses, and even then only in 15% of individuals.</p>
</div>

<div class="plot-container">
<img src="data:image/png;base64,{plots['mismatch']}" alt="Cross vs actual larval decision">
<p>Breakdown showing which individuals match vs. mismatch their cross-type expectation. About half of both DD and NN individuals make the "unexpected" developmental choice, indicating strong individual-level variation beyond family genotype.</p>
</div>

<h2>5. Development Time</h2>

{describe_by_cross('days_to_pupation', 'Days to Pupation by Cross Type')}

<div class="plot-container">
<img src="data:image/png;base64,{plots['pupation']}" alt="Days to pupation">
<p>Both distributions are right-skewed, with a cluster of fast developers (23&ndash;35 days, predominantly direct developers) and a long tail of slow developers (50&ndash;120 days, predominantly diapausing individuals). DD median (41 days) is slightly lower than NN (49 days). The bimodal pattern reflects the two developmental pathways.</p>
</div>

{describe_by_cross('days_to_eclosion', 'Days to Eclosion by Cross Type')}

<div class="plot-container">
<img src="data:image/png;base64,{plots['eclosion']}" alt="Days to eclosion">
<p>Most individuals eclose within 13&ndash;16 days after pupation. NN individuals show a longer tail (up to 30 days), likely driven by the 7 individuals that entered pupal diapause. DD eclosion times are tightly clustered (median 15 days), as none entered pupal diapause.</p>
</div>

<h2>6. Body Mass Traits</h2>

{describe_by_cross('pupal_weight_g', 'Pupal Weight (g) by Cross Type')}

<div class="plot-container">
<img src="data:image/png;base64,{plots['pupal_wt']}" alt="Pupal weight">
<p>Pupal weight distributions are very similar between DD and NN crosses (means ~0.147 vs ~0.149 g). Both show broad, roughly normal distributions. This trait is likely influenced more by larval feeding duration and sex than by cross type.</p>
</div>

{describe_by_cross('adult_weight_g', 'Adult Weight (g) by Cross Type')}

<div class="plot-container">
<img src="data:image/png;base64,{plots['adult_wt']}" alt="Adult weight">
<p>Adult weight (3 missing values, all from NN &mdash; likely individuals that did not survive to adulthood). Distributions overlap almost completely between DD and NN. Adult weight is roughly 40&ndash;45% of pupal weight, reflecting mass lost during metamorphosis.</p>
</div>

<h2>7. Sex Distribution</h2>

<div class="plot-container">
<img src="data:image/png;base64,{plots['sex']}" alt="Sex distribution">
<p>Sex ratios are balanced across both cross types: DD has 26M/21F, NN has 24M/23F. No strong sex bias in either cross. Sex should be included as a covariate in GWAS analyses, and can also serve as a positive control phenotype (should map to the Z chromosome in Lepidoptera).</p>
</div>

<h2>8. Trait Correlations</h2>

<div class="plot-container">
<img src="data:image/png;base64,{plots['correlation']}" alt="Correlation heatmap">
<p>Days to pupation and days to eclosion show moderate positive correlation (r = 0.47), as expected if diapausing individuals are slower at both stages. Body mass traits (pupal and adult weight) are strongly correlated (r = 0.73). Pupation time shows weak positive correlation with pupal weight (r = 0.26), suggesting longer-developing larvae grow somewhat heavier.</p>
</div>

<div class="plot-container">
<img src="data:image/png;base64,{plots['scatter']}" alt="Scatter plots">
<p><strong>Left:</strong> Longer pupation times tend to produce slightly heavier pupae, with DD and NN individuals intermixed across the range. <strong>Right:</strong> The positive correlation between pupation time and eclosion time is driven primarily by the diapausing individuals (long pupation + delayed eclosion).</p>
</div>

<h2>9. Cabinet (Environmental Covariate)</h2>

<div class="plot-container">
<img src="data:image/png;base64,{plots['cabinet']}" alt="Cabinet distribution">
<p>Individuals were reared in three growth cabinets. C1 housed the most individuals from both crosses, while C2 and C3 each had smaller subsets. DD were more evenly spread across cabinets than NN (which were concentrated in C1). Cabinet should be considered as a covariate in GWAS to control for potential environmental batch effects.</p>
</div>

<h2>10. Recommendations for GWAS</h2>

<table>
<tr><th>Phenotype</th><th>Type</th><th>GWAS suitability</th><th>Notes</th></tr>
<tr><td><code>diapause</code></td><td>Binary</td><td>Good (original)</td><td>Already run. Family-level label, not individual decision.</td></tr>
<tr><td><code>larval_decision</code></td><td>Binary</td><td><strong>Excellent</strong></td><td>Individual-level diapause decision. May reveal loci that the cross-level phenotype misses.</td></tr>
<tr><td><code>pupal_decision</code></td><td>Binary</td><td>Poor (underpowered)</td><td>Only 7 cases. Very low statistical power for GWAS.</td></tr>
<tr><td><code>sex</code></td><td>Binary</td><td>Good (positive control)</td><td>Should map to Z chromosome. Useful validation of GWAS pipeline.</td></tr>
<tr><td><code>days_to_pupation</code></td><td>Quantitative</td><td><strong>Excellent</strong></td><td>Continuous trait with broad range. Bimodal distribution reflects diapause biology. Consider log-transform.</td></tr>
<tr><td><code>days_to_eclosion</code></td><td>Quantitative</td><td>Moderate</td><td>Most values clustered at 14&ndash;16 days with a few outliers. May need rank-based normalisation.</td></tr>
<tr><td><code>pupal_weight_g</code></td><td>Quantitative</td><td>Good</td><td>Approximately normal distribution. Body size trait with ecological relevance.</td></tr>
<tr><td><code>adult_weight_g</code></td><td>Quantitative</td><td>Good</td><td>3 missing values. Strongly correlated with pupal weight (r = 0.73).</td></tr>
</table>

<div class="key-finding">
<strong>Suggested covariates for all GWAS runs:</strong> sex, cabinet, and the first N principal components from the PCA. For quantitative trait GWAS, also consider including <code>larval_decision</code> as a covariate to separate diapause-associated body size effects from direct genetic effects on growth.
</div>

<hr>
<p style="color: #95a5a6; font-size: 0.9em;">Generated from <code>qc_anna_data_phenotypes.xlsx</code> &mdash; {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}</p>

</body>
</html>'''

with open('phenotype_summary.html', 'w') as f:
    f.write(html)

print(f"Created phenotype_summary.html")
print("Done!")
