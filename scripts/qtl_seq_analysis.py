#!/usr/bin/env python3
"""
QTL-seq / Bulk Segregant Analysis for introgression mapping.
Parallelized version — uses multiprocessing for genotype extraction and statistics.

Usage:
    python scripts/qtl_seq_analysis.py <vcf> <phenotypes> <outdir> [n_cores]
"""

import sys, os, subprocess, base64, gzip
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import numpy as np

N_CORES = 20  # default, can be overridden via argv

# bcftools path (in Snakemake conda env)
BCFTOOLS = os.path.expanduser(
    "~/sbl_claudecode/Annas_GWAS/.snakemake/conda/"
    "aedce8cd90aec491ec84dcf2bee0648f_/bin/bcftools"
)
if not os.path.exists(BCFTOOLS):
    BCFTOOLS = "bcftools"

def parse_phenotypes(pheno_file):
    groups = {}
    with open(pheno_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            groups[parts[0]] = "NN" if int(parts[1]) == 0 else "DD"
    return groups

def get_chromosomes(vcf):
    """Get list of chromosomes from VCF index."""
    cmd = [BCFTOOLS, "index", "-s", vcf]
    out = subprocess.check_output(cmd, text=True)
    chroms = []
    for line in out.strip().split("\n"):
        if line:
            chroms.append(line.split("\t")[0])
    return chroms

def extract_genotypes_for_region(args):
    """Worker: extract genotypes for one chromosome. Returns numpy arrays."""
    vcf, region, samples_nn, samples_dd = args
    all_samples = samples_nn + samples_dd
    sample_str = ",".join(all_samples)
    n_nn = len(samples_nn)

    fmt = "%CHROM\\t%POS" + "[\\t%GT]" * len(all_samples) + "\\n"
    cmd = [BCFTOOLS, "query", "-f", fmt, "-s", sample_str, "-r", region, vcf]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1048576)

    chroms, positions = [], []
    nn_alt, nn_tot, dd_alt, dd_tot = [], [], [], []

    for line in proc.stdout:
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        gts = fields[2:]

        # Count alt alleles per group
        na, nt = 0, 0
        for i in range(n_nn):
            gt = gts[i]
            if "." in gt:
                continue
            alleles = gt.replace("|", "/").split("/")
            for a in alleles:
                if a != ".":
                    nt += 1
                    if a != "0":
                        na += 1
        da, dt = 0, 0
        for i in range(n_nn, len(gts)):
            gt = gts[i]
            if "." in gt:
                continue
            alleles = gt.replace("|", "/").split("/")
            for a in alleles:
                if a != ".":
                    dt += 1
                    if a != "0":
                        da += 1

        if nt == 0 or dt == 0:
            continue

        chroms.append(chrom)
        positions.append(pos)
        nn_alt.append(na)
        nn_tot.append(nt)
        dd_alt.append(da)
        dd_tot.append(dt)

    proc.wait()
    return (chroms, positions, nn_alt, nn_tot, dd_alt, dd_tot)

def extract_genotypes_parallel(vcf, samples_nn, samples_dd, n_cores):
    """Extract genotypes in parallel by chromosome."""
    chromosomes = get_chromosomes(vcf)
    print(f"Extracting genotypes across {len(chromosomes)} regions using {n_cores} cores...")

    args = [(vcf, chrom, samples_nn, samples_dd) for chrom in chromosomes]

    with Pool(n_cores) as pool:
        results = pool.map(extract_genotypes_for_region, args)

    # Merge results
    all_chroms, all_pos = [], []
    all_nn_alt, all_nn_tot, all_dd_alt, all_dd_tot = [], [], [], []
    for chroms, pos, na, nt, da, dt in results:
        all_chroms.extend(chroms)
        all_pos.extend(pos)
        all_nn_alt.extend(na)
        all_nn_tot.extend(nt)
        all_dd_alt.extend(da)
        all_dd_tot.extend(dt)

    print(f"  Total SNPs: {len(all_chroms):,}")
    return (np.array(all_chroms), np.array(all_pos, dtype=np.int64),
            np.array(all_nn_alt), np.array(all_nn_tot),
            np.array(all_dd_alt), np.array(all_dd_tot))

def compute_delta_snp_index(nn_alt, nn_tot, dd_alt, dd_tot):
    af_nn = nn_alt / nn_tot
    af_dd = dd_alt / dd_tot
    return af_dd - af_nn, af_nn, af_dd

def compute_fst_wc(nn_alt, nn_tot, dd_alt, dd_tot):
    """Weir & Cockerham Fst, vectorized."""
    n1 = nn_tot / 2
    n2 = dd_tot / 2
    n_bar = (n1 + n2) / 2
    r = 2

    p1 = nn_alt / nn_tot
    p2 = dd_alt / dd_tot
    p_bar = (n1 * p1 + n2 * p2) / (n1 + n2)
    n_c = (n1 + n2) - (n1**2 + n2**2) / (n1 + n2)
    s2 = (n1 * (p1 - p_bar)**2 + n2 * (p2 - p_bar)**2) / ((r - 1) * n_bar)

    # Safely compute heterozygosity
    nn_het = np.where(nn_tot > 1, nn_alt * (nn_tot - nn_alt) / (nn_tot * (nn_tot - 1)) * 2, 0)
    dd_het = np.where(dd_tot > 1, dd_alt * (dd_tot - dd_alt) / (dd_tot * (dd_tot - 1)) * 2, 0)
    h_bar = (nn_het + dd_het) / r

    a = (n_bar / n_c) * (s2 - (1 / (n_bar - 1)) * (p_bar * (1 - p_bar) - s2 * (r - 1) / r - h_bar / 4))
    b = (n_bar / (n_bar - 1)) * (p_bar * (1 - p_bar) - s2 * (r - 1) / r - h_bar * (2 * n_bar - 1) / (4 * n_bar))
    c = h_bar / 2

    denom = a + b + c
    fst = np.where(denom > 0, a / denom, 0.0)
    return np.clip(fst, -1, 1)

def _fisher_chunk(args):
    """Worker for parallel Fisher's exact test."""
    from scipy.stats import fisher_exact
    nn_alt, nn_tot, dd_alt, dd_tot = args
    n = len(nn_alt)
    pvals = np.ones(n)
    for i in range(n):
        nn_ref = nn_tot[i] - nn_alt[i]
        dd_ref = dd_tot[i] - dd_alt[i]
        try:
            _, p = fisher_exact([[nn_ref, nn_alt[i]], [dd_ref, dd_alt[i]]])
            pvals[i] = p
        except:
            pvals[i] = 1.0
    return pvals

def fisher_exact_parallel(nn_alt, nn_tot, dd_alt, dd_tot, n_cores):
    """Parallel Fisher's exact test."""
    n = len(nn_alt)
    chunk_size = max(1, n // n_cores)
    chunks = []
    for i in range(0, n, chunk_size):
        end = min(i + chunk_size, n)
        chunks.append((nn_alt[i:end], nn_tot[i:end], dd_alt[i:end], dd_tot[i:end]))

    print(f"  Running Fisher's exact test on {n:,} SNPs across {len(chunks)} chunks...")
    with Pool(n_cores) as pool:
        results = pool.map(_fisher_chunk, chunks)

    return np.concatenate(results)

def sliding_window(chroms, positions, values, window_size=1_000_000, step=100_000):
    result_chroms, result_pos, result_vals = [], [], []
    for chrom in sorted(set(chroms)):
        mask = chroms == chrom
        pos = positions[mask]
        vals = values[mask]
        if len(pos) == 0:
            continue
        start, end = int(pos.min()), int(pos.max())
        for w_start in range(start, end, int(step)):
            w_end = w_start + window_size
            w_mask = (pos >= w_start) & (pos < w_end)
            if w_mask.sum() < 5:
                continue
            result_chroms.append(chrom)
            result_pos.append(w_start + window_size // 2)
            result_vals.append(np.nanmean(vals[w_mask]))
    return np.array(result_chroms), np.array(result_pos), np.array(result_vals)

def _perm_worker(args):
    """Single permutation: hypergeometric reshuffle, compute max windowed |delta|."""
    nn_alt, nn_tot, dd_alt, dd_tot, chroms, positions, window_size, step, seed = args
    rng = np.random.default_rng(seed)
    total_alt = nn_alt + dd_alt
    total_tot = nn_tot + dd_tot
    perm_nn_alt = rng.hypergeometric(total_alt, total_tot - total_alt, nn_tot)
    perm_dd_alt = total_alt - perm_nn_alt
    perm_delta = perm_dd_alt / dd_tot - perm_nn_alt / nn_tot
    _, _, w_vals = sliding_window(chroms, positions, perm_delta, window_size, step)
    return np.max(np.abs(w_vals)) if len(w_vals) > 0 else 0.0

def permutation_thresholds(nn_alt, nn_tot, dd_alt, dd_tot, chroms, positions,
                            n_perm=1000, n_cores=20, window_size=1_000_000,
                            step=100_000, quantiles=[0.95, 0.99]):
    print(f"Running {n_perm} permutations using {n_cores} cores...")
    args = [(nn_alt, nn_tot, dd_alt, dd_tot, chroms, positions,
             window_size, step, i) for i in range(n_perm)]

    with Pool(n_cores) as pool:
        max_deltas = pool.map(_perm_worker, args)

    max_deltas = np.array(max_deltas)
    thresholds = {}
    for q in quantiles:
        thresholds[q] = np.quantile(max_deltas, q)
        print(f"  {q*100:.0f}% threshold: {thresholds[q]:.4f}")
    return thresholds

def make_plots(chroms, positions, delta_snp, fst, pvals, af_nn, af_dd,
               w_chroms, w_pos, w_delta, w_fst, thresholds, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    autosomes = [str(i) for i in range(1, 28)]
    sex_chr = ["Z", "W"]
    all_chr_ordered = autosomes + sex_chr
    scaffolds = sorted(set(chroms) - set(all_chr_ordered))
    chr_order = [c for c in all_chr_ordered if c in set(chroms)] + scaffolds

    chr_cum, cum, gap = {}, 0, 2_000_000
    chr_mids = {}
    for c in chr_order:
        mask = chroms == c
        if mask.sum() == 0:
            continue
        chr_cum[c] = cum
        max_pos = positions[mask].max()
        chr_mids[c] = cum + max_pos // 2
        cum += max_pos + gap

    cum_pos = np.array([positions[i] + chr_cum.get(chroms[i], 0) for i in range(len(chroms))], dtype=np.float64)
    w_cum_pos = np.array([w_pos[i] + chr_cum.get(w_chroms[i], 0) for i in range(len(w_chroms))], dtype=np.float64)

    chr_num = np.zeros(len(chroms), dtype=int)
    for i, c in enumerate(chr_order):
        chr_num[chroms == c] = i
    colors_raw = np.where(chr_num % 2 == 0, "#67a9cf", "#2166ac")

    label_pos = [chr_mids[c] for c in chr_order if c in chr_mids and not c.startswith("NW_")]
    label_names = [c for c in chr_order if c in chr_mids and not c.startswith("NW_")]
    scaff_mids = [chr_mids[c] for c in chr_order if c in chr_mids and c.startswith("NW_")]
    if scaff_mids:
        label_pos.append(np.mean(scaff_mids))
        label_names.append("Scaff")

    fig, axes = plt.subplots(4, 1, figsize=(16, 16), sharex=True,
                              gridspec_kw={"hspace": 0.08})

    # Panel 1: Δ(SNP-index)
    ax = axes[0]
    ax.scatter(cum_pos, delta_snp, s=0.3, alpha=0.15, c=colors_raw, linewidths=0, rasterized=True)
    ax.plot(w_cum_pos, w_delta, color="#d62728", linewidth=1.2, label="1 Mb window")
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    if 0.99 in thresholds:
        ax.axhline(thresholds[0.99], color="#e377c2", linewidth=0.8, linestyle="--",
                   label=f"99% perm ({thresholds[0.99]:.3f})")
        ax.axhline(-thresholds[0.99], color="#e377c2", linewidth=0.8, linestyle="--")
    if 0.95 in thresholds:
        ax.axhline(thresholds[0.95], color="#e377c2", linewidth=0.8, linestyle=":",
                   alpha=0.6, label=f"95% perm ({thresholds[0.95]:.3f})")
        ax.axhline(-thresholds[0.95], color="#e377c2", linewidth=0.8, linestyle=":", alpha=0.6)
    ax.set_ylabel("Δ(SNP-index)\nDD − NN", fontsize=11)
    ax.set_ylim(-1.05, 1.05)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.set_title("QTL-seq Analysis — Pararge aegeria (NN vs DD backcross)", fontsize=13, fontweight="bold", pad=10)

    # Panel 2: Fst
    ax = axes[1]
    ax.scatter(cum_pos, fst, s=0.3, alpha=0.15, c=colors_raw, linewidths=0, rasterized=True)
    ax.plot(w_cum_pos, w_fst, color="#d62728", linewidth=1.2, label="1 Mb window")
    ax.set_ylabel("Fst\n(Weir & Cockerham)", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # Panel 3: -log10(p)
    ax = axes[2]
    logp = -np.log10(np.clip(pvals, 1e-300, 1))
    ax.scatter(cum_pos, logp, s=0.3, alpha=0.15, c=colors_raw, linewidths=0, rasterized=True)
    ax.axhline(-np.log10(5e-8), color="red", linewidth=0.8, linestyle="--", label="5×10⁻⁸")
    ax.axhline(-np.log10(1e-5), color="blue", linewidth=0.8, linestyle=":", alpha=0.6, label="1×10⁻⁵")
    ax.set_ylabel("-log₁₀(p)\nFisher's exact", fontsize=11)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # Panel 4: Allele frequencies
    ax = axes[3]
    ax.scatter(cum_pos, af_nn, s=0.3, alpha=0.1, c="#2166ac", linewidths=0, rasterized=True, label="NN")
    ax.scatter(cum_pos, af_dd, s=0.3, alpha=0.1, c="#d62728", linewidths=0, rasterized=True, label="DD")
    ax.set_ylabel("Alt allele\nfrequency", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9, markerscale=5)

    axes[3].set_xticks(label_pos)
    axes[3].set_xticklabels(label_names, fontsize=8)
    axes[3].set_xlabel("Chromosome", fontsize=12)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    png_path = os.path.join(outdir, "qtl_seq_genome_plot.png")
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, "qtl_seq_genome_plot.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plots saved: {png_path}")
    return png_path

def find_qtl_peaks(w_chroms, w_pos, w_delta, w_fst, thresholds, window_size=1_000_000):
    thresh = thresholds.get(0.95, 0.1)
    sig_mask = np.abs(w_delta) > thresh
    if sig_mask.sum() == 0:
        return []

    peaks, current = [], None
    for i in range(len(w_chroms)):
        if sig_mask[i]:
            if current is None or w_chroms[i] != current["chrom"] or \
               w_pos[i] - (current["end"] - window_size // 2) > window_size * 2:
                if current is not None:
                    peaks.append(current)
                current = {"chrom": w_chroms[i], "start": w_pos[i] - window_size // 2,
                           "end": w_pos[i] + window_size // 2,
                           "max_delta": w_delta[i], "max_fst": w_fst[i], "n_windows": 1}
            else:
                current["end"] = w_pos[i] + window_size // 2
                current["n_windows"] += 1
                if abs(w_delta[i]) > abs(current["max_delta"]):
                    current["max_delta"] = w_delta[i]
                if w_fst[i] > current["max_fst"]:
                    current["max_fst"] = w_fst[i]
        else:
            if current is not None:
                peaks.append(current)
                current = None
    if current is not None:
        peaks.append(current)

    peaks.sort(key=lambda p: abs(p["max_delta"]), reverse=True)
    return peaks

def generate_html_report(peaks, thresholds, plot_path, n_snps, n_nn, n_dd, outdir):
    with open(plot_path, "rb") as f:
        plot_b64 = base64.b64encode(f.read()).decode()

    if peaks:
        rows = ""
        for i, p in enumerate(peaks[:20]):
            size_mb = (p["end"] - p["start"]) / 1_000_000
            rows += f"""<tr><td>{i+1}</td><td class="highlight">{p["chrom"]}</td>
                <td>{p["start"]:,} – {p["end"]:,}</td><td>{size_mb:.1f} Mb</td>
                <td>{p["max_delta"]:.4f}</td><td>{p["max_fst"]:.4f}</td>
                <td>{p["n_windows"]}</td></tr>"""
        peaks_html = f"""<h2>Significant QTL Regions (95% permutation threshold)</h2>
        <div class="card"><table><thead>
            <tr><th>#</th><th>Chr</th><th>Region</th><th>Size</th>
                <th>Max Δ(SNP-index)</th><th>Max Fst</th><th>Windows</th></tr>
        </thead><tbody>{rows}</tbody></table></div>"""
    else:
        peaks_html = """<h2>Significant QTL Regions</h2>
        <div class="card"><p>No regions exceeded the 95% permutation threshold.</p></div>"""

    t95 = f"{thresholds[0.95]:.4f}" if 0.95 in thresholds else "N/A"
    t99 = f"{thresholds[0.99]:.4f}" if 0.99 in thresholds else "N/A"

    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>QTL-seq Report</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:-apple-system,BlinkMacSystemFont,sans-serif; background:#f8fafc;
       color:#1e293b; line-height:1.6; max-width:1200px; margin:0 auto; padding:2rem; }}
h1 {{ font-size:1.6rem; margin-bottom:0.3rem; }}
.subtitle {{ color:#64748b; margin-bottom:2rem; font-size:0.9rem; }}
h2 {{ font-size:1.15rem; margin:2rem 0 0.8rem; border-bottom:2px solid #e2e8f0; padding-bottom:0.3rem; }}
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
<h1>QTL-seq Analysis Report</h1>
<div class="subtitle"><em>Pararge aegeria</em> — Backcross introgression mapping (NN vs DD)</div>

<h2>Experimental Design</h2>
<div class="card"><p>Backcross design with selection for alternative homozygous states at a
focal diapause locus. After several generations of backcrossing, only the introgressed region
should show allele frequency divergence.</p></div>

<div class="grid">
<div class="card"><table>
    <tr><td>NN samples</td><td class="highlight">{n_nn}</td></tr>
    <tr><td>DD samples</td><td class="highlight">{n_dd}</td></tr>
    <tr><td>Total SNPs</td><td class="highlight">{n_snps:,}</td></tr>
</table></div>
<div class="card"><table>
    <tr><td>Window size</td><td>1 Mb (100 kb step)</td></tr>
    <tr><td>Permutations</td><td>1,000</td></tr>
    <tr><td>95% threshold</td><td>{t95}</td></tr>
    <tr><td>99% threshold</td><td>{t99}</td></tr>
</table></div>
</div>

<h2>Genome-wide Plot</h2>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
<b>Panel 1:</b> Δ(SNP-index) = AF<sub>DD</sub> − AF<sub>NN</sub>. Red = 1 Mb window. Dashed = permutation thresholds.<br>
<b>Panel 2:</b> Weir &amp; Cockerham F<sub>ST</sub>.<br>
<b>Panel 3:</b> Fisher's exact test (−log₁₀ p).<br>
<b>Panel 4:</b> Allele frequencies (blue = NN, red = DD).</p>
<div class="plot-container">
    <img src="data:image/png;base64,{plot_b64}" alt="QTL-seq genome plot">
</div></div>

{peaks_html}

<h2>Method</h2>
<div class="card"><p style="font-size:0.85rem;">
Per-SNP allele frequencies were calculated in each group from the filtered VCF (biallelic SNPs,
QUAL&gt;20, DP&gt;5). Δ(SNP-index) was smoothed with 1 Mb sliding windows (100 kb step, min 5
SNPs/window). F<sub>ST</sub> estimated per Weir &amp; Cockerham (1984). Significance thresholds
from 1,000 permutations using hypergeometric resampling. Fisher's exact test per SNP.</p></div>

<h2>Output Files</h2>
<div class="card"><table>
<tr><td>Genome plot (PNG)</td><td><code>results/qtl_seq/qtl_seq_genome_plot.png</code></td></tr>
<tr><td>Genome plot (PDF)</td><td><code>results/qtl_seq/qtl_seq_genome_plot.pdf</code></td></tr>
<tr><td>Per-SNP results</td><td><code>results/qtl_seq/qtl_seq_per_snp.tsv.gz</code></td></tr>
<tr><td>Window results</td><td><code>results/qtl_seq/qtl_seq_windows.tsv</code></td></tr>
<tr><td>QTL peaks</td><td><code>results/qtl_seq/qtl_peaks.tsv</code></td></tr>
</table></div>
</body></html>"""

    with open(os.path.join(outdir, "summary_report.html"), "w") as f:
        f.write(html)
    print("HTML report saved.")

def main():
    if len(sys.argv) < 4:
        print("Usage: python qtl_seq_analysis.py <vcf> <phenotypes> <outdir> [n_cores]")
        sys.exit(1)

    vcf, pheno_file, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    n_cores = int(sys.argv[4]) if len(sys.argv) > 4 else N_CORES
    os.makedirs(outdir, exist_ok=True)

    groups = parse_phenotypes(pheno_file)
    samples_nn = sorted([s for s, g in groups.items() if g == "NN"])
    samples_dd = sorted([s for s, g in groups.items() if g == "DD"])
    print(f"NN: {len(samples_nn)}, DD: {len(samples_dd)}, Cores: {n_cores}")

    # Step 1: Extract genotypes (parallel by chromosome)
    chroms, positions, nn_alt, nn_tot, dd_alt, dd_tot = \
        extract_genotypes_parallel(vcf, samples_nn, samples_dd, n_cores)
    n_snps = len(chroms)

    # Step 2: Δ(SNP-index)
    print("Computing Δ(SNP-index)...")
    delta_snp, af_nn, af_dd = compute_delta_snp_index(nn_alt, nn_tot, dd_alt, dd_tot)

    # Step 3: Fst
    print("Computing Fst...")
    valid = (nn_tot >= 4) & (dd_tot >= 4)
    fst = np.zeros(n_snps)
    fst[valid] = compute_fst_wc(nn_alt[valid], nn_tot[valid], dd_alt[valid], dd_tot[valid])

    # Step 4: Fisher's exact test (parallel)
    print("Computing Fisher's exact test...")
    pvals = fisher_exact_parallel(nn_alt, nn_tot, dd_alt, dd_tot, n_cores)

    # Step 5: Sliding windows
    print("Sliding windows...")
    w_chroms, w_pos, w_delta = sliding_window(chroms, positions, delta_snp)
    _, _, w_fst = sliding_window(chroms, positions, fst)

    # Step 6: Permutation thresholds (parallel)
    thresholds = permutation_thresholds(nn_alt, nn_tot, dd_alt, dd_tot,
                                         chroms, positions, n_perm=1000, n_cores=n_cores)

    # Step 7: Find peaks
    print("Finding QTL peaks...")
    peaks = find_qtl_peaks(w_chroms, w_pos, w_delta, w_fst, thresholds)
    print(f"Found {len(peaks)} significant regions")
    for p in peaks[:5]:
        print(f"  Chr {p['chrom']}: {p['start']:,}-{p['end']:,} "
              f"(Δ={p['max_delta']:.4f}, Fst={p['max_fst']:.4f})")

    # Step 8: Save results
    print("Saving results...")
    with gzip.open(os.path.join(outdir, "qtl_seq_per_snp.tsv.gz"), "wt") as f:
        f.write("chrom\tpos\taf_nn\taf_dd\tdelta_snp_index\tfst\tp_fisher\n")
        for i in range(n_snps):
            f.write(f"{chroms[i]}\t{positions[i]}\t{af_nn[i]:.6f}\t{af_dd[i]:.6f}\t"
                    f"{delta_snp[i]:.6f}\t{fst[i]:.6f}\t{pvals[i]:.6e}\n")

    with open(os.path.join(outdir, "qtl_seq_windows.tsv"), "w") as f:
        f.write("chrom\tpos\twindow_delta\twindow_fst\n")
        for i in range(len(w_chroms)):
            f.write(f"{w_chroms[i]}\t{w_pos[i]}\t{w_delta[i]:.6f}\t{w_fst[i]:.6f}\n")

    with open(os.path.join(outdir, "qtl_peaks.tsv"), "w") as f:
        f.write("chrom\tstart\tend\tsize_bp\tmax_delta\tmax_fst\tn_windows\n")
        for p in peaks:
            f.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\t{p['end']-p['start']}\t"
                    f"{p['max_delta']:.6f}\t{p['max_fst']:.6f}\t{p['n_windows']}\n")

    # Step 9: Plots + HTML
    print("Generating plots...")
    plot_path = make_plots(chroms, positions, delta_snp, fst, pvals, af_nn, af_dd,
                            w_chroms, w_pos, w_delta, w_fst, thresholds, outdir)
    generate_html_report(peaks, thresholds, plot_path, n_snps,
                          len(samples_nn), len(samples_dd), outdir)

    print("\n=== QTL-seq analysis complete ===")

if __name__ == "__main__":
    main()
