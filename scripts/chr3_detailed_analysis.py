#!/usr/bin/env python3
"""
Detailed chromosome 3 analysis: QTL-seq + genotype heatmap.
Identifies the exact genomic region causing NN vs DD phenotypic differences.

Usage:
    python scripts/chr3_detailed_analysis.py <vcf> <phenotypes> <outdir> [n_cores]
"""

import sys, os, subprocess, base64, gzip
from multiprocessing import Pool
import numpy as np

N_CORES = 20
BCFTOOLS = os.path.expanduser(
    "~/sbl_claudecode/Annas_GWAS/.snakemake/conda/"
    "aedce8cd90aec491ec84dcf2bee0648f_/bin/bcftools"
)
if not os.path.exists(BCFTOOLS):
    BCFTOOLS = "bcftools"

CHROM = "3"
MAX_MISSINGNESS = 0.05  # keep SNPs with < 5% missing genotypes

def parse_phenotypes(pheno_file):
    groups = {}
    with open(pheno_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            groups[parts[0]] = "NN" if int(parts[1]) == 0 else "DD"
    return groups

def extract_chr3_genotypes(vcf, samples_nn, samples_dd):
    """Extract full genotype matrix for chr3."""
    all_samples = samples_nn + samples_dd
    sample_str = ",".join(all_samples)
    n_total = len(all_samples)
    n_nn = len(samples_nn)

    # Single [\\t%GT] iterates over ALL selected samples automatically
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
    cmd = [BCFTOOLS, "query", "-f", fmt, "-s", sample_str, "-r", CHROM, vcf]

    print(f"Extracting chr3 genotypes for {n_total} samples...")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                            bufsize=1048576, stderr=subprocess.DEVNULL)

    positions = []
    # Store genotypes as: 0=hom_ref, 1=het, 2=hom_alt, -1=missing
    genotype_matrix = []
    nn_alt_list, nn_tot_list = [], []
    dd_alt_list, dd_tot_list = [], []
    per_snp_missing = []

    n_parsed = 0
    for line in proc.stdout:
        fields = line.strip().split("\t")
        pos = int(fields[1])
        gts_raw = fields[4:]

        # Parse genotypes into numeric
        geno_row = []
        n_missing = 0
        nn_alt, nn_tot, dd_alt, dd_tot = 0, 0, 0, 0

        for j, gt in enumerate(gts_raw):
            if "." in gt:
                geno_row.append(-1)
                n_missing += 1
                continue
            alleles = gt.replace("|", "/").split("/")
            alt_count = sum(1 for a in alleles if a != "0")
            geno_row.append(alt_count)  # 0, 1, or 2

            if j < n_nn:
                nn_tot += 2
                nn_alt += alt_count
            else:
                dd_tot += 2
                dd_alt += alt_count

        miss_rate = n_missing / n_total
        if nn_tot == 0 or dd_tot == 0:
            continue

        positions.append(pos)
        genotype_matrix.append(geno_row)
        nn_alt_list.append(nn_alt)
        nn_tot_list.append(nn_tot)
        dd_alt_list.append(dd_alt)
        dd_tot_list.append(dd_tot)
        per_snp_missing.append(miss_rate)

        n_parsed += 1
        if n_parsed % 200000 == 0:
            print(f"  {n_parsed:,} SNPs parsed...")

    proc.wait()
    print(f"  Total chr3 SNPs: {n_parsed:,}")

    return (np.array(positions, dtype=np.int64),
            np.array(genotype_matrix, dtype=np.int8),
            np.array(nn_alt_list), np.array(nn_tot_list),
            np.array(dd_alt_list), np.array(dd_tot_list),
            np.array(per_snp_missing))

def compute_stats(nn_alt, nn_tot, dd_alt, dd_tot):
    """Compute Δ(SNP-index) and Fst."""
    af_nn = nn_alt / nn_tot
    af_dd = dd_alt / dd_tot
    delta = af_dd - af_nn

    # Fst (Weir & Cockerham)
    n1, n2 = nn_tot / 2, dd_tot / 2
    n_bar = (n1 + n2) / 2
    p1, p2 = af_nn, af_dd
    p_bar = (n1 * p1 + n2 * p2) / (n1 + n2)
    n_c = (n1 + n2) - (n1**2 + n2**2) / (n1 + n2)
    s2 = (n1 * (p1 - p_bar)**2 + n2 * (p2 - p_bar)**2) / n_bar

    nn_het = np.where(nn_tot > 1, nn_alt * (nn_tot - nn_alt) / (nn_tot * (nn_tot - 1)) * 2, 0)
    dd_het = np.where(dd_tot > 1, dd_alt * (dd_tot - dd_alt) / (dd_tot * (dd_tot - 1)) * 2, 0)
    h_bar = (nn_het + dd_het) / 2

    a = (n_bar / n_c) * (s2 - (1 / (n_bar - 1)) * (p_bar * (1 - p_bar) - s2 / 2 - h_bar / 4))
    b = (n_bar / (n_bar - 1)) * (p_bar * (1 - p_bar) - s2 / 2 - h_bar * (2 * n_bar - 1) / (4 * n_bar))
    c = h_bar / 2
    denom = a + b + c
    fst = np.where(denom > 0, a / denom, 0.0)
    fst = np.clip(fst, -1, 1)

    return delta, af_nn, af_dd, fst

def _fisher_chunk(args):
    from scipy.stats import fisher_exact
    nn_alt, nn_tot, dd_alt, dd_tot = args
    pvals = np.ones(len(nn_alt))
    for i in range(len(nn_alt)):
        nn_ref = nn_tot[i] - nn_alt[i]
        dd_ref = dd_tot[i] - dd_alt[i]
        try:
            _, p = fisher_exact([[nn_ref, nn_alt[i]], [dd_ref, dd_alt[i]]])
            pvals[i] = p
        except:
            pass
    return pvals

def fisher_parallel(nn_alt, nn_tot, dd_alt, dd_tot, n_cores):
    n = len(nn_alt)
    chunk_size = max(1, n // n_cores)
    chunks = []
    for i in range(0, n, chunk_size):
        e = min(i + chunk_size, n)
        chunks.append((nn_alt[i:e], nn_tot[i:e], dd_alt[i:e], dd_tot[i:e]))
    with Pool(n_cores) as pool:
        results = pool.map(_fisher_chunk, chunks)
    return np.concatenate(results)

def sliding_window(positions, values, window_size, step=None):
    if step is None:
        step = window_size // 10
    w_pos, w_vals = [], []
    start, end = int(positions.min()), int(positions.max())
    for w_start in range(start, end, int(step)):
        w_end = w_start + window_size
        mask = (positions >= w_start) & (positions < w_end)
        if mask.sum() < 3:
            continue
        w_pos.append(w_start + window_size // 2)
        w_vals.append(np.nanmean(values[mask]))
    return np.array(w_pos), np.array(w_vals)

def find_differentiated_region(w_pos, w_delta, w_fst, threshold_delta=0.05):
    """Find contiguous region with elevated |Δ(SNP-index)| or Fst."""
    # Use a combination: windows where |delta| is in top 5% or Fst > threshold
    delta_thresh = np.quantile(np.abs(w_delta), 0.95)
    fst_thresh = np.quantile(w_fst, 0.95)

    sig = (np.abs(w_delta) >= delta_thresh) | (w_fst >= fst_thresh)
    if sig.sum() == 0:
        return None, None

    # Find the largest contiguous block
    regions = []
    current_start = None
    for i in range(len(sig)):
        if sig[i]:
            if current_start is None:
                current_start = i
        else:
            if current_start is not None:
                regions.append((current_start, i - 1))
                current_start = None
    if current_start is not None:
        regions.append((current_start, len(sig) - 1))

    # Take largest region
    if not regions:
        return None, None
    best = max(regions, key=lambda r: r[1] - r[0])
    return int(w_pos[best[0]]), int(w_pos[best[1]])

def make_qtlseq_plot(positions, delta, fst, pvals, af_nn, af_dd,
                      w_pos_250k, w_delta_250k, w_fst_250k,
                      w_pos_100k, w_delta_100k, w_fst_100k,
                      region_start, region_end, outdir):
    """4-panel QTL-seq plot for chr3 at fine resolution."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pos_mb = positions / 1e6

    fig, axes = plt.subplots(4, 1, figsize=(16, 14), sharex=True,
                              gridspec_kw={"hspace": 0.08})

    # Panel 1: Δ(SNP-index) with two window sizes
    ax = axes[0]
    ax.scatter(pos_mb, delta, s=0.4, alpha=0.12, c="#67a9cf", linewidths=0, rasterized=True)
    ax.plot(w_pos_250k / 1e6, w_delta_250k, color="#d62728", linewidth=1.0, label="250 kb window")
    ax.plot(w_pos_100k / 1e6, w_delta_100k, color="#ff7f0e", linewidth=0.8, alpha=0.7, label="100 kb window")
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    if region_start and region_end:
        ax.axvspan(region_start / 1e6, region_end / 1e6, alpha=0.1, color="red", label="Candidate region")
    ax.set_ylabel("Δ(SNP-index)\nDD − NN", fontsize=11)
    ax.set_ylim(-1.05, 1.05)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(f"Chromosome 3 — Detailed QTL-seq Analysis", fontsize=13, fontweight="bold", pad=10)

    # Panel 2: Fst
    ax = axes[1]
    ax.scatter(pos_mb, fst, s=0.4, alpha=0.12, c="#67a9cf", linewidths=0, rasterized=True)
    ax.plot(w_pos_250k / 1e6, w_fst_250k, color="#d62728", linewidth=1.0, label="250 kb")
    ax.plot(w_pos_100k / 1e6, w_fst_100k, color="#ff7f0e", linewidth=0.8, alpha=0.7, label="100 kb")
    if region_start and region_end:
        ax.axvspan(region_start / 1e6, region_end / 1e6, alpha=0.1, color="red")
    ax.set_ylabel("Fst", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="upper right", fontsize=8)

    # Panel 3: -log10(p)
    ax = axes[2]
    logp = -np.log10(np.clip(pvals, 1e-300, 1))
    ax.scatter(pos_mb, logp, s=0.4, alpha=0.12, c="#67a9cf", linewidths=0, rasterized=True)
    ax.axhline(-np.log10(5e-8), color="red", linewidth=0.8, linestyle="--", label="5×10⁻⁸")
    ax.axhline(-np.log10(1e-5), color="blue", linewidth=0.6, linestyle=":", label="1×10⁻⁵")
    if region_start and region_end:
        ax.axvspan(region_start / 1e6, region_end / 1e6, alpha=0.1, color="red")
    ax.set_ylabel("-log₁₀(p)\nFisher's exact", fontsize=11)
    ax.legend(loc="upper right", fontsize=8)

    # Panel 4: Allele frequencies
    ax = axes[3]
    ax.scatter(pos_mb, af_nn, s=0.4, alpha=0.08, c="#2166ac", linewidths=0, rasterized=True, label="NN")
    ax.scatter(pos_mb, af_dd, s=0.4, alpha=0.08, c="#d62728", linewidths=0, rasterized=True, label="DD")
    if region_start and region_end:
        ax.axvspan(region_start / 1e6, region_end / 1e6, alpha=0.1, color="red")
    ax.set_ylabel("Alt AF", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Position on chromosome 3 (Mb)", fontsize=12)
    ax.legend(loc="upper right", fontsize=8, markerscale=5)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    path = os.path.join(outdir, "chr3_qtlseq_detail.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, "chr3_qtlseq_detail.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"QTL-seq detail plot saved: {path}")
    return path

def make_genotype_heatmap(positions, genotype_matrix, missingness,
                           samples_nn, samples_dd, outdir,
                           region_start=None, region_end=None,
                           max_snps=5000, tag="full"):
    """
    Genotype heatmap: rows = samples (NN then DD), columns = SNPs.
    Colors: blue=hom_ref(0/0), yellow=het(0/1), red=hom_alt(1/1), grey=missing.
    Only uses SNPs with low missingness.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import matplotlib.patches as mpatches

    n_nn = len(samples_nn)
    n_dd = len(samples_dd)
    n_total = n_nn + n_dd

    # Filter for low missingness
    low_miss = missingness < MAX_MISSINGNESS
    print(f"  [{tag}] SNPs with <{MAX_MISSINGNESS*100:.0f}% missingness: {low_miss.sum():,} / {len(missingness):,}")

    pos_filt = positions[low_miss]
    geno_filt = genotype_matrix[low_miss]

    # Region filter if specified
    if region_start is not None and region_end is not None:
        region_mask = (pos_filt >= region_start) & (pos_filt <= region_end)
        pos_filt = pos_filt[region_mask]
        geno_filt = geno_filt[region_mask]
        print(f"  [{tag}] SNPs in region {region_start:,}-{region_end:,}: {len(pos_filt):,}")

    if len(pos_filt) == 0:
        print(f"  [{tag}] No SNPs passed filters!")
        return None

    # Thin if too many SNPs for visualization
    if len(pos_filt) > max_snps:
        step = len(pos_filt) // max_snps
        idx = np.arange(0, len(pos_filt), step)[:max_snps]
        pos_filt = pos_filt[idx]
        geno_filt = geno_filt[idx]
        print(f"  [{tag}] Thinned to {len(pos_filt):,} SNPs for plotting")

    # Reorder: NN samples first, then DD
    # genotype_matrix columns are already [nn_samples..., dd_samples...]
    # Transpose for heatmap: rows = samples, columns = SNPs
    geno_display = geno_filt.T  # shape: (n_samples, n_snps)

    # Map -1 (missing) to 3 for colormap
    geno_display = np.where(geno_display == -1, 3, geno_display)

    # Custom colormap: 0=hom_ref(blue), 1=het(yellow), 2=hom_alt(red), 3=missing(grey)
    cmap = ListedColormap(["#2166ac", "#ffd700", "#d62728", "#d3d3d3"])

    fig_height = max(6, n_total * 0.12 + 2)
    fig_width = max(12, len(pos_filt) * 0.003 + 3)
    fig_width = min(fig_width, 24)  # cap width

    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))

    im = ax.imshow(geno_display, aspect="auto", cmap=cmap, vmin=0, vmax=3,
                   interpolation="none")

    # Sample labels on y-axis
    sample_labels = samples_nn + samples_dd
    ax.set_yticks(range(n_total))
    ax.set_yticklabels(sample_labels, fontsize=5)

    # Draw line separating NN and DD
    ax.axhline(n_nn - 0.5, color="white", linewidth=2)

    # Group labels
    ax.text(-0.02, (n_nn - 1) / 2 / n_total, "NN", transform=ax.transAxes,
            fontsize=12, fontweight="bold", color="#2166ac", va="center", ha="right")
    ax.text(-0.02, (n_nn + (n_total - 1)) / 2 / n_total, "DD", transform=ax.transAxes,
            fontsize=12, fontweight="bold", color="#d62728", va="center", ha="right")

    # X-axis: positions in Mb
    n_ticks = min(20, len(pos_filt))
    tick_idx = np.linspace(0, len(pos_filt) - 1, n_ticks, dtype=int)
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([f"{pos_filt[i]/1e6:.1f}" for i in tick_idx], fontsize=7, rotation=45)
    ax.set_xlabel("Position on chr3 (Mb)", fontsize=11)

    # Legend
    patches = [mpatches.Patch(color="#2166ac", label="0/0 (hom ref)"),
               mpatches.Patch(color="#ffd700", label="0/1 (het)"),
               mpatches.Patch(color="#d62728", label="1/1 (hom alt)"),
               mpatches.Patch(color="#d3d3d3", label="missing")]
    ax.legend(handles=patches, loc="upper right", fontsize=8, ncol=4,
              bbox_to_anchor=(1, 1.08))

    if tag == "full":
        ax.set_title(f"Genotype Plot — Chromosome 3 (low-missingness SNPs)", fontsize=13, fontweight="bold")
    else:
        r_start_mb = region_start / 1e6 if region_start else 0
        r_end_mb = region_end / 1e6 if region_end else 0
        ax.set_title(f"Genotype Plot — Chr3: {r_start_mb:.1f}–{r_end_mb:.1f} Mb (candidate region)",
                     fontsize=13, fontweight="bold")

    plt.tight_layout()
    path = os.path.join(outdir, f"chr3_genotype_heatmap_{tag}.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, f"chr3_genotype_heatmap_{tag}.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  [{tag}] Genotype heatmap saved: {path}")
    return path

def generate_html(qtlseq_plot, heatmap_full, heatmap_zoom, outdir,
                   n_snps_total, n_snps_lowmiss, n_nn, n_dd,
                   region_start, region_end,
                   top_snps_data, window_stats):
    """Generate comprehensive HTML report."""

    def img_b64(path):
        if path and os.path.exists(path):
            with open(path, "rb") as f:
                return base64.b64encode(f.read()).decode()
        return ""

    qtl_b64 = img_b64(qtlseq_plot)
    hm_full_b64 = img_b64(heatmap_full)
    hm_zoom_b64 = img_b64(heatmap_zoom)

    region_text = ""
    if region_start and region_end:
        size_mb = (region_end - region_start) / 1e6
        region_text = f"""
        <h2>Candidate Region</h2>
        <div class="card">
            <table>
                <tr><td>Chromosome</td><td class="highlight">3</td></tr>
                <tr><td>Start</td><td class="highlight">{region_start:,} ({region_start/1e6:.2f} Mb)</td></tr>
                <tr><td>End</td><td class="highlight">{region_end:,} ({region_end/1e6:.2f} Mb)</td></tr>
                <tr><td>Size</td><td class="highlight">{size_mb:.2f} Mb</td></tr>
            </table>
            <p style="font-size:0.85rem;color:#475569;margin-top:0.8rem;">
                This region was identified as the top 5% of windows by |Δ(SNP-index)| and F<sub>ST</sub>.
                It represents the most likely location of the introgressed segment controlling
                the NN/DD phenotypic difference.
            </p>
        </div>"""

    # Top SNPs table
    top_snps_html = ""
    if top_snps_data:
        rows = ""
        for i, (pos, delta, fst, p, af_nn, af_dd) in enumerate(top_snps_data[:20]):
            rows += f"""<tr><td>{i+1}</td><td>{pos:,}</td><td>{pos/1e6:.3f}</td>
                <td>{delta:.4f}</td><td>{fst:.4f}</td><td>{p:.2e}</td>
                <td>{af_nn:.3f}</td><td>{af_dd:.3f}</td></tr>"""
        top_snps_html = f"""
        <h2>Top 20 Most Differentiated SNPs</h2>
        <div class="card"><table><thead>
            <tr><th>#</th><th>Position</th><th>Mb</th><th>Δ(SNP-index)</th>
                <th>Fst</th><th>p (Fisher)</th><th>AF NN</th><th>AF DD</th></tr>
        </thead><tbody>{rows}</tbody></table></div>"""

    # Window summary
    window_html = ""
    if window_stats:
        rows = ""
        for pos, d, f in window_stats[:15]:
            rows += f"<tr><td>{pos/1e6:.2f}</td><td>{d:.4f}</td><td>{f:.4f}</td></tr>"
        window_html = f"""
        <h2>Top Windows (250 kb, by |Δ SNP-index|)</h2>
        <div class="card"><table><thead>
            <tr><th>Position (Mb)</th><th>Window Δ(SNP-index)</th><th>Window Fst</th></tr>
        </thead><tbody>{rows}</tbody></table></div>"""

    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Chr3 Detailed Analysis</title>
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
table {{ width:100%; border-collapse:collapse; font-size:0.85rem; }}
th {{ background:#f1f5f9; text-align:left; padding:0.5rem 0.7rem; font-weight:600;
      border-bottom:2px solid #e2e8f0; }}
td {{ padding:0.4rem 0.7rem; border-bottom:1px solid #f1f5f9; }}
.highlight {{ color:#2563eb; font-weight:600; }}
.plot-container {{ text-align:center; margin:1rem 0; }}
.plot-container img {{ max-width:100%; border:1px solid #e2e8f0; border-radius:6px; }}
.grid {{ display:grid; grid-template-columns:1fr 1fr; gap:1rem; }}
.grid .card {{ margin-bottom:0; }}
</style></head><body>

<h1>Chromosome 3 — Detailed Analysis</h1>
<div class="subtitle"><em>Pararge aegeria</em> — Identifying the causal region for NN vs DD phenotype</div>

<h2>Rationale</h2>
<div class="card"><p>Genome-wide QTL-seq and GWAS both identified chromosome 3 as harbouring
the strongest differentiation between NN and DD backcross lines. This detailed analysis
focuses exclusively on chromosome 3 at high resolution to:</p>
<ul style="margin:0.5rem 0 0 1.5rem;">
<li>Pinpoint the exact genomic region of the introgressed segment</li>
<li>Visualize individual haplotypes with a genotype plot to confirm the signal</li>
<li>Identify the boundaries of the differentiated block</li>
</ul></div>

<div class="grid">
<div class="card"><table>
    <tr><td>NN samples</td><td class="highlight">{n_nn}</td></tr>
    <tr><td>DD samples</td><td class="highlight">{n_dd}</td></tr>
    <tr><td>Chr3 SNPs (total)</td><td class="highlight">{n_snps_total:,}</td></tr>
    <tr><td>Chr3 SNPs (&lt;5% missing)</td><td class="highlight">{n_snps_lowmiss:,}</td></tr>
</table></div>
<div class="card"><table>
    <tr><td>QTL-seq windows</td><td>250 kb + 100 kb</td></tr>
    <tr><td>Genotype plot filter</td><td>&lt;5% per-SNP missingness</td></tr>
    <tr><td>Genotype encoding</td><td>0/0 (blue), 0/1 (yellow), 1/1 (red)</td></tr>
</table></div>
</div>

{region_text}

<h2>QTL-seq Detail — Chromosome 3</h2>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
Fine-resolution QTL-seq with 250 kb (red) and 100 kb (orange) sliding windows.
The shaded red region marks the candidate introgression. Panels show Δ(SNP-index),
F<sub>ST</sub>, Fisher's exact −log₁₀(p), and per-group allele frequencies.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{qtl_b64}" alt="Chr3 QTL-seq detail">
</div></div>

<h2>Genotype Plot — Full Chromosome 3</h2>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
Each row is one individual (NN above the white line, DD below). Each column is a SNP.
<b>Blue</b> = homozygous reference, <b>red</b> = homozygous alternate, <b>yellow</b> = heterozygote.
Only SNPs with &lt;5% missingness are shown. The introgressed region should appear as a
block where NN and DD have opposite homozygous genotypes.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{hm_full_b64}" alt="Genotype heatmap full chr3">
</div></div>

{"<h2>Genotype Plot — Candidate Region (Zoomed)</h2>" + '''
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
Zoomed view of the candidate region. At this resolution, individual haplotype blocks
should be clearly visible. Look for a contiguous segment where NN samples are predominantly
blue (hom ref) and DD samples are predominantly red (hom alt), or vice versa.</p>
<div class="plot-container">
    <img src="data:image/png;base64,''' + hm_zoom_b64 + '''" alt="Genotype heatmap zoomed">
</div></div>''' if hm_zoom_b64 else ""}

{top_snps_html}

{window_html}

<h2>Method</h2>
<div class="card"><p style="font-size:0.85rem;">
Genotypes for all 94 samples on chromosome 3 were extracted from the filtered VCF (biallelic SNPs,
QUAL&gt;20, DP&gt;5). Per-SNP allele frequencies were computed for each group.
Δ(SNP-index) = AF<sub>DD</sub> − AF<sub>NN</sub> and Weir &amp; Cockerham F<sub>ST</sub>
were smoothed with 250 kb and 100 kb sliding windows. Fisher's exact test was computed per SNP.
The candidate region was defined as contiguous windows in the top 5% of |Δ(SNP-index)| or F<sub>ST</sub>.
For the genotype heatmap, only SNPs with &lt;5% missingness were retained to ensure clean
haplotype visualization. Samples are ordered NN (top) then DD (bottom) to highlight the
introgression pattern.</p></div>

<h2>Output Files</h2>
<div class="card"><table>
<tr><td>QTL-seq plot</td><td><code>results/QTLseq_chrom3/chr3_qtlseq_detail.{{png,pdf}}</code></td></tr>
<tr><td>Genotype heatmap (full)</td><td><code>results/QTLseq_chrom3/chr3_genotype_heatmap_full.{{png,pdf}}</code></td></tr>
<tr><td>Genotype heatmap (zoom)</td><td><code>results/QTLseq_chrom3/chr3_genotype_heatmap_zoom.{{png,pdf}}</code></td></tr>
<tr><td>Per-SNP results</td><td><code>results/QTLseq_chrom3/chr3_per_snp.tsv.gz</code></td></tr>
</table></div>

</body></html>"""

    report_path = os.path.join(outdir, "summary_report.html")
    with open(report_path, "w") as f:
        f.write(html)
    print(f"HTML report saved: {report_path}")

def main():
    if len(sys.argv) < 4:
        print("Usage: python chr3_detailed_analysis.py <vcf> <phenotypes> <outdir> [n_cores]")
        sys.exit(1)

    vcf, pheno_file, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    n_cores = int(sys.argv[4]) if len(sys.argv) > 4 else N_CORES
    os.makedirs(outdir, exist_ok=True)

    groups = parse_phenotypes(pheno_file)
    samples_nn = sorted([s for s, g in groups.items() if g == "NN"])
    samples_dd = sorted([s for s, g in groups.items() if g == "DD"])
    print(f"NN: {len(samples_nn)}, DD: {len(samples_dd)}")

    # Step 1: Extract chr3 genotypes
    print("\n=== Step 1: Extract chr3 genotypes ===")
    positions, geno_matrix, nn_alt, nn_tot, dd_alt, dd_tot, missingness = \
        extract_chr3_genotypes(vcf, samples_nn, samples_dd)
    n_snps_total = len(positions)
    n_snps_lowmiss = (missingness < MAX_MISSINGNESS).sum()

    # Step 2: Compute statistics
    print("\n=== Step 2: Compute statistics ===")
    delta, af_nn, af_dd, fst = compute_stats(nn_alt, nn_tot, dd_alt, dd_tot)

    print("Fisher's exact test...")
    pvals = fisher_parallel(nn_alt, nn_tot, dd_alt, dd_tot, n_cores)

    # Step 3: Sliding windows at multiple resolutions
    print("\n=== Step 3: Sliding windows ===")
    w_pos_250k, w_delta_250k = sliding_window(positions, delta, 250_000, 25_000)
    _, w_fst_250k = sliding_window(positions, fst, 250_000, 25_000)
    w_pos_100k, w_delta_100k = sliding_window(positions, delta, 100_000, 10_000)
    _, w_fst_100k = sliding_window(positions, fst, 100_000, 10_000)
    print(f"  250 kb windows: {len(w_pos_250k):,}")
    print(f"  100 kb windows: {len(w_pos_100k):,}")

    # Step 4: Find candidate region
    print("\n=== Step 4: Identify candidate region ===")
    region_start, region_end = find_differentiated_region(
        w_pos_250k, w_delta_250k, w_fst_250k)
    if region_start and region_end:
        print(f"  Candidate region: {region_start:,} - {region_end:,} "
              f"({(region_end-region_start)/1e6:.2f} Mb)")
    else:
        print("  No clear candidate region identified from windows.")
        # Fall back to region around top SNPs
        top_idx = np.argsort(np.abs(delta))[-100:]
        region_start = int(positions[top_idx].min()) - 500_000
        region_end = int(positions[top_idx].max()) + 500_000
        print(f"  Using top-SNP region: {region_start:,} - {region_end:,}")

    # Step 5: QTL-seq detail plot
    print("\n=== Step 5: QTL-seq detail plot ===")
    qtlseq_plot = make_qtlseq_plot(
        positions, delta, fst, pvals, af_nn, af_dd,
        w_pos_250k, w_delta_250k, w_fst_250k,
        w_pos_100k, w_delta_100k, w_fst_100k,
        region_start, region_end, outdir)

    # Step 6: Genotype heatmaps
    print("\n=== Step 6: Genotype heatmaps ===")
    heatmap_full = make_genotype_heatmap(
        positions, geno_matrix, missingness, samples_nn, samples_dd,
        outdir, max_snps=5000, tag="full")

    # Zoomed heatmap of candidate region — use more SNPs for detail
    heatmap_zoom = make_genotype_heatmap(
        positions, geno_matrix, missingness, samples_nn, samples_dd,
        outdir, region_start=region_start, region_end=region_end,
        max_snps=8000, tag="zoom")

    # Step 7: Top SNPs and window stats
    print("\n=== Step 7: Compile results ===")
    # Top SNPs by |delta|
    top_idx = np.argsort(np.abs(delta))[::-1][:20]
    top_snps_data = [(int(positions[i]), float(delta[i]), float(fst[i]),
                      float(pvals[i]), float(af_nn[i]), float(af_dd[i]))
                     for i in top_idx]

    # Top windows
    w_order = np.argsort(np.abs(w_delta_250k))[::-1][:15]
    window_stats = [(float(w_pos_250k[i]), float(w_delta_250k[i]),
                     float(w_fst_250k[i])) for i in w_order]

    # Save per-SNP results
    print("Saving per-SNP results...")
    with gzip.open(os.path.join(outdir, "chr3_per_snp.tsv.gz"), "wt") as f:
        f.write("pos\taf_nn\taf_dd\tdelta\tfst\tp_fisher\tmissingness\n")
        for i in range(n_snps_total):
            f.write(f"{positions[i]}\t{af_nn[i]:.6f}\t{af_dd[i]:.6f}\t"
                    f"{delta[i]:.6f}\t{fst[i]:.6f}\t{pvals[i]:.6e}\t"
                    f"{missingness[i]:.4f}\n")

    # Step 8: HTML report
    print("\n=== Step 8: Generate HTML report ===")
    generate_html(qtlseq_plot, heatmap_full, heatmap_zoom, outdir,
                   n_snps_total, n_snps_lowmiss, len(samples_nn), len(samples_dd),
                   region_start, region_end, top_snps_data, window_stats)

    print("\n=== Chromosome 3 detailed analysis complete ===")

if __name__ == "__main__":
    main()
