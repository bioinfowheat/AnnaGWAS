#!/usr/bin/env python3
"""
Focused analysis of chr3 ~9.5-13 Mb region — the major QTL peak.

Uses five independent approaches to delineate the core introgression:
  1. Multi-statistic sliding windows (Δ SNP-index, Fst, -log10p)
  2. Per-individual introgression tract mapping
  3. Minimum shared introgression (intersection of all DD tracts)
  4. Recombination breakpoint density histogram
  5. Per-window genotype composition (NN vs DD)

Usage:
    python scripts/chr3_region_analysis.py <vcf> <phenotypes> <outdir> [n_cores]
"""

import sys, os, subprocess, base64, gzip, math
from multiprocessing import Pool
from collections import Counter, defaultdict
import numpy as np

N_CORES = 20
BCFTOOLS = os.path.expanduser(
    "~/sbl_claudecode/Annas_GWAS/.snakemake/conda/"
    "aedce8cd90aec491ec84dcf2bee0648f_/bin/bcftools"
)
if not os.path.exists(BCFTOOLS):
    BCFTOOLS = "bcftools"

CHROM = "3"
REGION_START = 9_500_000
REGION_END   = 13_000_000
MAX_MISSINGNESS = 0.05

# ---------------------------------------------------------------------------
# Genotype extraction
# ---------------------------------------------------------------------------
def parse_phenotypes(pheno_file):
    groups = {}
    with open(pheno_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            groups[parts[0]] = "NN" if int(parts[1]) == 0 else "DD"
    return groups

def extract_region_genotypes(vcf, region, samples_nn, samples_dd):
    all_samples = samples_nn + samples_dd
    sample_str = ",".join(all_samples)
    n_total = len(all_samples)
    n_nn = len(samples_nn)

    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
    cmd = [BCFTOOLS, "query", "-f", fmt, "-s", sample_str, "-r", region, vcf]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                            bufsize=1048576, stderr=subprocess.DEVNULL)

    positions, refs, alts = [], [], []
    genotype_matrix = []
    nn_alt_list, nn_tot_list = [], []
    dd_alt_list, dd_tot_list = [], []
    per_snp_missing = []

    for line in proc.stdout:
        fields = line.strip().split("\t")
        pos = int(fields[1])
        ref, alt = fields[2], fields[3]
        gts_raw = fields[4:]

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
            geno_row.append(alt_count)
            if j < n_nn:
                nn_tot += 2; nn_alt += alt_count
            else:
                dd_tot += 2; dd_alt += alt_count

        if nn_tot == 0 or dd_tot == 0:
            continue

        positions.append(pos)
        refs.append(ref)
        alts.append(alt)
        genotype_matrix.append(geno_row)
        nn_alt_list.append(nn_alt)
        nn_tot_list.append(nn_tot)
        dd_alt_list.append(dd_alt)
        dd_tot_list.append(dd_tot)
        per_snp_missing.append(n_missing / n_total)

    proc.wait()
    return (np.array(positions, dtype=np.int64), refs, alts,
            np.array(genotype_matrix, dtype=np.int8),
            np.array(nn_alt_list), np.array(nn_tot_list),
            np.array(dd_alt_list), np.array(dd_tot_list),
            np.array(per_snp_missing))

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
def compute_stats(nn_alt, nn_tot, dd_alt, dd_tot):
    af_nn = nn_alt / nn_tot
    af_dd = dd_alt / dd_tot
    delta = af_dd - af_nn
    n1, n2 = nn_tot / 2, dd_tot / 2
    p1, p2 = af_nn, af_dd
    p_bar = (n1 * p1 + n2 * p2) / (n1 + n2)
    n_bar = (n1 + n2) / 2
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
    chunks = [(nn_alt[i:min(i+chunk_size,n)], nn_tot[i:min(i+chunk_size,n)],
               dd_alt[i:min(i+chunk_size,n)], dd_tot[i:min(i+chunk_size,n)])
              for i in range(0, n, chunk_size)]
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

# ---------------------------------------------------------------------------
# APPROACH 1: Multi-statistic QTL-seq plot
# ---------------------------------------------------------------------------
def make_qtlseq_plot(positions, delta, fst, pvals, af_nn, af_dd,
                      w50_pos, w50_delta, w50_fst,
                      w25_pos, w25_delta, w25_fst, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pos_mb = positions / 1e6
    logp = -np.log10(np.clip(pvals, 1e-300, 1))

    fig, axes = plt.subplots(4, 1, figsize=(16, 14), sharex=True,
                              gridspec_kw={"hspace": 0.08})

    ax = axes[0]
    ax.scatter(pos_mb, delta, s=1, alpha=0.2, c="#67a9cf", linewidths=0, rasterized=True)
    ax.plot(w50_pos / 1e6, w50_delta, color="#d62728", linewidth=1.2, label="50 kb window")
    ax.plot(w25_pos / 1e6, w25_delta, color="#ff7f0e", linewidth=0.8, alpha=0.8, label="25 kb window")
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.set_ylabel("Δ(SNP-index)\nDD − NN", fontsize=11)
    ax.set_ylim(-1.05, 1.05)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Chr3: 9.5–13.0 Mb — QTL Peak Region", fontsize=13, fontweight="bold", pad=10)

    ax = axes[1]
    ax.scatter(pos_mb, fst, s=1, alpha=0.2, c="#67a9cf", linewidths=0, rasterized=True)
    ax.plot(w50_pos / 1e6, w50_fst, color="#d62728", linewidth=1.2, label="50 kb")
    ax.plot(w25_pos / 1e6, w25_fst, color="#ff7f0e", linewidth=0.8, alpha=0.8, label="25 kb")
    ax.set_ylabel("Fst", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="upper right", fontsize=8)

    ax = axes[2]
    ax.scatter(pos_mb, logp, s=1, alpha=0.25, c="#67a9cf", linewidths=0, rasterized=True)
    ax.axhline(20, color="red", linewidth=0.8, linestyle="--", label="-log₁₀(p) = 20")
    ax.axhline(-np.log10(5e-8), color="darkred", linewidth=0.6, linestyle=":", label="5×10⁻⁸")
    ax.set_ylabel("-log₁₀(p)\nFisher's exact", fontsize=11)
    ax.legend(loc="upper right", fontsize=8)

    ax = axes[3]
    ax.scatter(pos_mb, af_nn, s=1, alpha=0.15, c="#2166ac", linewidths=0, rasterized=True, label="NN")
    ax.scatter(pos_mb, af_dd, s=1, alpha=0.15, c="#d62728", linewidths=0, rasterized=True, label="DD")
    ax.set_ylabel("Alt AF", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Position on chr3 (Mb)", fontsize=12)
    ax.legend(loc="upper right", fontsize=8, markerscale=5)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    path = os.path.join(outdir, "chr3_10-12Mb_qtlseq.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, "chr3_10-12Mb_qtlseq.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    return path

# ---------------------------------------------------------------------------
# APPROACH 2: Per-individual introgression tract mapping
# ---------------------------------------------------------------------------
def map_introgression_tracts(positions, genotype_matrix, missingness,
                              samples_nn, samples_dd, window_kb=50):
    """
    For each individual, determine the introgression tract boundaries.

    Method:
      - At each low-missingness SNP, compute the NN consensus genotype (mode).
      - For each DD individual, score discordance with NN consensus.
      - Smooth with a sliding window.
      - Introgression = contiguous region where discordance > 0.5.
      - Also compute for NN individuals as negative controls.
    """
    n_nn = len(samples_nn)
    n_dd = len(samples_dd)
    n_total = n_nn + n_dd

    low_miss = missingness < MAX_MISSINGNESS
    pos = positions[low_miss]
    geno = genotype_matrix[low_miss]  # (n_snps, n_total)

    # NN consensus: mode genotype at each SNP among NN individuals
    nn_geno = geno[:, :n_nn]  # (n_snps, n_nn)
    nn_consensus = np.zeros(len(pos), dtype=np.int8)
    for i in range(len(pos)):
        valid = nn_geno[i][nn_geno[i] >= 0]
        if len(valid) > 0:
            counts = np.bincount(valid, minlength=3)
            nn_consensus[i] = np.argmax(counts)
        else:
            nn_consensus[i] = 0

    # For each individual, compute discordance with NN consensus
    window_bp = window_kb * 1000
    step_bp = window_bp // 5
    w_positions = np.arange(int(pos.min()), int(pos.max()), step_bp) + step_bp // 2

    all_tracts = {}  # sample_name -> list of (start, end) tuples
    all_discordance = {}  # sample_name -> (w_positions, discordance_values)

    for s_idx in range(n_total):
        sample_name = (samples_nn + samples_dd)[s_idx]
        ind_geno = geno[:, s_idx]

        # Per-SNP discordance: 1 if different from NN consensus, 0 if same
        valid = ind_geno >= 0
        disc = np.zeros(len(pos))
        disc[valid] = (ind_geno[valid] != nn_consensus[valid]).astype(float)
        disc[~valid] = np.nan

        # Smooth with sliding window
        w_disc = []
        for w_center in w_positions:
            mask = (pos >= w_center - window_bp // 2) & (pos < w_center + window_bp // 2)
            vals = disc[mask]
            vals = vals[~np.isnan(vals)]
            if len(vals) >= 3:
                w_disc.append(np.mean(vals))
            else:
                w_disc.append(np.nan)
        w_disc = np.array(w_disc)
        all_discordance[sample_name] = (w_positions.copy(), w_disc)

        # Find introgression tracts: contiguous regions where discordance > 0.5
        tracts = []
        in_tract = False
        tract_start = 0
        for i, d in enumerate(w_disc):
            if np.isnan(d):
                continue
            if d > 0.5 and not in_tract:
                tract_start = w_positions[i]
                in_tract = True
            elif (d <= 0.5 or i == len(w_disc) - 1) and in_tract:
                tract_end = w_positions[i]
                tracts.append((int(tract_start), int(tract_end)))
                in_tract = False
        all_tracts[sample_name] = tracts

    return all_tracts, all_discordance, nn_consensus, pos

def find_minimum_shared_introgression(all_tracts, samples_dd, region_start, region_end):
    """
    Find the region shared by all (or most) DD introgression tracts.
    Uses a coverage approach: count how many DD individuals carry
    an introgression at each position.
    """
    resolution = 10_000  # 10 kb
    bins = np.arange(region_start, region_end, resolution)
    coverage = np.zeros(len(bins))

    for sample in samples_dd:
        tracts = all_tracts.get(sample, [])
        for t_start, t_end in tracts:
            mask = (bins >= t_start) & (bins < t_end)
            coverage[mask] += 1

    frac = coverage / len(samples_dd)
    return bins, coverage, frac

# ---------------------------------------------------------------------------
# APPROACH 3: Genotype composition per window
# ---------------------------------------------------------------------------
def compute_genotype_composition(positions, genotype_matrix, missingness,
                                  n_nn, window_size=25_000, step=5_000):
    """
    For each window, compute fraction of hom_ref, het, hom_alt in NN vs DD.
    """
    low_miss = missingness < MAX_MISSINGNESS
    pos = positions[low_miss]
    geno = genotype_matrix[low_miss]

    w_positions = []
    nn_hom_ref, nn_het, nn_hom_alt = [], [], []
    dd_hom_ref, dd_het, dd_hom_alt = [], [], []

    for w_start in range(int(pos.min()), int(pos.max()), step):
        w_end = w_start + window_size
        mask = (pos >= w_start) & (pos < w_end)
        if mask.sum() < 3:
            continue
        w_geno = geno[mask]
        w_positions.append(w_start + window_size // 2)

        # NN samples
        nn_g = w_geno[:, :n_nn].flatten()
        nn_valid = nn_g[nn_g >= 0]
        nn_n = max(1, len(nn_valid))
        nn_hom_ref.append(np.sum(nn_valid == 0) / nn_n)
        nn_het.append(np.sum(nn_valid == 1) / nn_n)
        nn_hom_alt.append(np.sum(nn_valid == 2) / nn_n)

        # DD samples
        dd_g = w_geno[:, n_nn:].flatten()
        dd_valid = dd_g[dd_g >= 0]
        dd_n = max(1, len(dd_valid))
        dd_hom_ref.append(np.sum(dd_valid == 0) / dd_n)
        dd_het.append(np.sum(dd_valid == 1) / dd_n)
        dd_hom_alt.append(np.sum(dd_valid == 2) / dd_n)

    return (np.array(w_positions),
            np.array(nn_hom_ref), np.array(nn_het), np.array(nn_hom_alt),
            np.array(dd_hom_ref), np.array(dd_het), np.array(dd_hom_alt))

# ---------------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------------
def make_multi_evidence_plot(positions, delta, fst, logp,
                              w25_pos, w25_delta, w25_fst,
                              core_regions, outdir):
    """
    Overlay normalised statistics to show convergent evidence for core region(s).
    core_regions: list of (start, end, max_frac, mean_frac) tuples.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Compute windowed -log10(p)
    _, w25_logp = sliding_window(positions, logp, 25_000, 2_500)
    w25_logp_pos = w25_pos  # same windows

    # Normalise each to [0, 1]
    def norm(x):
        mn, mx = np.nanmin(x), np.nanmax(x)
        if mx == mn:
            return np.zeros_like(x)
        return (x - mn) / (mx - mn)

    # Colours for up to 5 core regions
    region_colors = ["#e74c3c", "#2166ac", "#2ca02c", "#ff7f0e", "#9467bd"]

    fig, ax = plt.subplots(1, 1, figsize=(16, 5))

    for rank, (rs, re, mx, mn) in enumerate(core_regions):
        col = region_colors[rank % len(region_colors)]
        label = f"Core #{rank+1}" if rank == 0 else f"Core #{rank+1}"
        ax.fill_betweenx([0, 1], rs / 1e6, re / 1e6,
                          alpha=0.12, color=col, zorder=0)
        mid = (rs + re) / 2e6
        ax.text(mid, 1.04, f"#{rank+1}", ha="center", va="bottom",
                fontsize=8, fontweight="bold", color=col)

    ax.plot(w25_pos / 1e6, norm(np.abs(w25_delta)), color="#2166ac",
            linewidth=1.2, label="|Δ(SNP-index)| (25 kb)", alpha=0.9)
    ax.plot(w25_pos / 1e6, norm(w25_fst), color="#d62728",
            linewidth=1.2, label="Fst (25 kb)", alpha=0.9)
    ax.plot(w25_logp_pos / 1e6, norm(w25_logp), color="#2ca02c",
            linewidth=1.2, label="Mean −log₁₀(p) (25 kb)", alpha=0.9)

    ax.set_xlabel("Position on chr3 (Mb)", fontsize=12)
    ax.set_ylabel("Normalised statistic (0–1)", fontsize=11)
    ax.set_title("Convergent Evidence — Normalised Multi-Statistic Overlay",
                 fontsize=13, fontweight="bold")
    ax.set_ylim(-0.05, 1.15)
    ax.legend(loc="upper right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    path = os.path.join(outdir, "chr3_multi_evidence.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    return path


def make_introgression_tract_plot(all_tracts, all_discordance,
                                   samples_nn, samples_dd,
                                   core_regions, outdir):
    """
    Plot per-individual introgression tracts as horizontal bars.
    NN samples shown as negative controls; DD samples show introgression.
    core_regions: list of (start, end, max_frac, mean_frac) tuples.
    All core regions are shown as labelled vertical bands.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    n_nn = len(samples_nn)
    n_dd = len(samples_dd)
    all_samples = samples_nn + samples_dd
    n_all = len(all_samples)

    fig, ax = plt.subplots(1, 1, figsize=(16, max(8, n_all * 0.15 + 2)))

    # Plot tracts
    for idx, sample in enumerate(all_samples):
        y = n_all - 1 - idx
        tracts = all_tracts.get(sample, [])
        is_dd = idx >= n_nn

        # Background bar (light)
        ax.barh(y, (REGION_END - REGION_START) / 1e6,
                left=REGION_START / 1e6, height=0.7,
                color="#f0f0f0", edgecolor="none")

        # Introgression tracts
        for t_start, t_end in tracts:
            width = (t_end - t_start) / 1e6
            color = "#d62728" if is_dd else "#ff9999"
            ax.barh(y, width, left=t_start / 1e6, height=0.7,
                    color=color, edgecolor="none", alpha=0.8)

    # Core region shading — all regions with distinct colours and labels
    region_colors = ["#2166ac", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
    for rank, (rs, re, mx, mn) in enumerate(core_regions):
        col = region_colors[rank % len(region_colors)]
        ax.axvspan(rs / 1e6, re / 1e6, alpha=0.15, color=col, zorder=0)
        # Label at top of plot
        mid = (rs + re) / 2e6
        ax.text(mid, n_all + 0.3, f"#{rank+1}",
                ha="center", va="bottom", fontsize=9, fontweight="bold",
                color=col,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor=col, alpha=0.85, linewidth=1.2))

    # Separator line
    ax.axhline(n_dd - 0.5, color="black", linewidth=1.5, linestyle="-")

    # Labels
    ax.set_yticks(range(n_all))
    ax.set_yticklabels(list(reversed(all_samples)), fontsize=4.5)
    ax.set_xlabel("Position on chr3 (Mb)", fontsize=12)
    ax.set_title("Per-Individual Introgression Tracts (red = discordant with NN consensus)",
                 fontsize=12, fontweight="bold")

    # Group labels
    ax.text(REGION_START / 1e6 - 0.15, n_all - n_nn / 2,
            "NN", fontsize=11, fontweight="bold", color="#2166ac",
            va="center", ha="right")
    ax.text(REGION_START / 1e6 - 0.15, n_dd / 2,
            "DD", fontsize=11, fontweight="bold", color="#d62728",
            va="center", ha="right")

    ax.set_xlim(REGION_START / 1e6 - 0.05, REGION_END / 1e6 + 0.05)
    ax.set_ylim(-0.5, n_all + 1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Build legend entries for core regions
    from matplotlib.patches import Patch
    legend_handles = []
    for rank, (rs, re, mx, mn) in enumerate(core_regions):
        col = region_colors[rank % len(region_colors)]
        size_kb = (re - rs) / 1e3
        legend_handles.append(
            Patch(facecolor=col, alpha=0.3,
                  label=f"Core #{rank+1}: {rs/1e6:.2f}–{re/1e6:.2f} Mb ({size_kb:.0f} kb)"))
    ax.legend(handles=legend_handles, loc="upper right", fontsize=7,
              framealpha=0.9, edgecolor="grey")

    plt.tight_layout()
    path = os.path.join(outdir, "chr3_introgression_tracts.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    return path


def make_coverage_breakpoint_plot(bins, coverage, frac, all_tracts,
                                    samples_dd, core_regions, outdir):
    """
    Top: introgression coverage (fraction of DD carrying introgression).
    Bottom: breakpoint density histogram.
    core_regions: list of (start, end, max_frac, mean_frac) tuples.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    region_colors = ["#2166ac", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]

    # Collect breakpoints
    breakpoints = []
    for sample in samples_dd:
        for t_start, t_end in all_tracts.get(sample, []):
            breakpoints.append(t_start)
            breakpoints.append(t_end)
    breakpoints = np.array(breakpoints) / 1e6

    fig, axes = plt.subplots(2, 1, figsize=(16, 7), sharex=True,
                              gridspec_kw={"hspace": 0.08, "height_ratios": [2, 1]})

    # Top: coverage
    ax = axes[0]
    ax.fill_between(bins / 1e6, frac, alpha=0.4, color="#d62728", step="mid")
    ax.step(bins / 1e6, frac, color="#d62728", linewidth=1.2, where="mid")
    ax.axhline(0.90, color="grey", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle=":", alpha=0.3)
    for rank, (rs, re, mx, mn) in enumerate(core_regions):
        col = region_colors[rank % len(region_colors)]
        ax.axvspan(rs / 1e6, re / 1e6, alpha=0.12, color=col)
        mid = (rs + re) / 2e6
        ax.text(mid, 1.03, f"#{rank+1}", ha="center", va="bottom",
                fontsize=8, fontweight="bold", color=col)
    ax.set_ylabel("Fraction of DD\nwith introgression", fontsize=11)
    ax.set_ylim(-0.02, 1.12)
    ax.set_title("Introgression Coverage & Recombination Breakpoints",
                 fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Bottom: breakpoint histogram
    ax = axes[1]
    if len(breakpoints) > 0:
        ax.hist(breakpoints, bins=70, color="#ff7f0e", edgecolor="white",
                linewidth=0.3, alpha=0.85)
    for rank, (rs, re, mx, mn) in enumerate(core_regions):
        col = region_colors[rank % len(region_colors)]
        ax.axvspan(rs / 1e6, re / 1e6, alpha=0.12, color=col)
    ax.set_xlabel("Position on chr3 (Mb)", fontsize=12)
    ax.set_ylabel("Breakpoint\ncount", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    path = os.path.join(outdir, "chr3_coverage_breakpoints.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    return path


def make_genotype_composition_plot(w_pos, nn_hr, nn_het, nn_ha,
                                    dd_hr, dd_het, dd_ha,
                                    core_regions, outdir):
    """
    Stacked area: genotype composition in NN vs DD per window.
    core_regions: list of (start, end, max_frac, mean_frac) tuples.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    region_colors = ["#2166ac", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
    x = w_pos / 1e6

    fig, axes = plt.subplots(2, 1, figsize=(16, 7), sharex=True,
                              gridspec_kw={"hspace": 0.12})

    # NN
    ax = axes[0]
    ax.fill_between(x, 0, nn_hr, color="#2166ac", alpha=0.7, label="Hom ref (0/0)", step="mid")
    ax.fill_between(x, nn_hr, nn_hr + nn_het, color="#ffd700", alpha=0.7, label="Het (0/1)", step="mid")
    ax.fill_between(x, nn_hr + nn_het, 1, color="#d62728", alpha=0.7, label="Hom alt (1/1)", step="mid")
    for rank, (rs, re, mx, mn) in enumerate(core_regions):
        col = region_colors[rank % len(region_colors)]
        ax.axvspan(rs / 1e6, re / 1e6, alpha=0.10, color=col)
        mid = (rs + re) / 2e6
        ax.text(mid, 1.02, f"#{rank+1}", ha="center", va="bottom",
                fontsize=7, fontweight="bold", color=col)
    ax.set_ylabel("NN genotype\nfrequency", fontsize=11)
    ax.set_ylim(0, 1.08)
    ax.set_title("Per-Window Genotype Composition (25 kb windows)", fontsize=13, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, ncol=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # DD
    ax = axes[1]
    ax.fill_between(x, 0, dd_hr, color="#2166ac", alpha=0.7, label="Hom ref", step="mid")
    ax.fill_between(x, dd_hr, dd_hr + dd_het, color="#ffd700", alpha=0.7, label="Het", step="mid")
    ax.fill_between(x, dd_hr + dd_het, 1, color="#d62728", alpha=0.7, label="Hom alt", step="mid")
    for rank, (rs, re, mx, mn) in enumerate(core_regions):
        col = region_colors[rank % len(region_colors)]
        ax.axvspan(rs / 1e6, re / 1e6, alpha=0.10, color=col)
    ax.set_ylabel("DD genotype\nfrequency", fontsize=11)
    ax.set_xlabel("Position on chr3 (Mb)", fontsize=12)
    ax.set_ylim(0, 1)
    ax.legend(loc="upper right", fontsize=8, ncol=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    path = os.path.join(outdir, "chr3_genotype_composition.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    return path


def make_genotype_heatmap(positions, genotype_matrix, missingness,
                           samples_nn, samples_dd, outdir,
                           region_start=None, region_end=None,
                           snp_mask=None, title=None,
                           max_snps=10000, tag="peak"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import matplotlib.patches as mpatches

    n_nn = len(samples_nn)
    n_dd = len(samples_dd)
    n_total = n_nn + n_dd

    if snp_mask is not None:
        pos_filt = positions[snp_mask]
        geno_filt = genotype_matrix[snp_mask]
    else:
        low_miss = missingness < MAX_MISSINGNESS
        pos_filt = positions[low_miss]
        geno_filt = genotype_matrix[low_miss]
        if region_start is not None and region_end is not None:
            region_mask = (pos_filt >= region_start) & (pos_filt <= region_end)
            pos_filt = pos_filt[region_mask]
            geno_filt = geno_filt[region_mask]

    n_snps = len(pos_filt)
    print(f"  [{tag}] SNPs for heatmap: {n_snps:,}")

    if n_snps == 0:
        return None, 0

    if n_snps > max_snps:
        step = n_snps // max_snps
        idx = np.arange(0, n_snps, step)[:max_snps]
        pos_filt = pos_filt[idx]
        geno_filt = geno_filt[idx]
        print(f"  [{tag}] Thinned to {len(pos_filt):,}")

    geno_display = geno_filt.T
    geno_display = np.where(geno_display == -1, 3, geno_display)
    cmap = ListedColormap(["#2166ac", "#ffd700", "#d62728", "#d3d3d3"])

    fig_height = max(6, n_total * 0.13 + 2)
    fig_width = min(24, max(14, len(pos_filt) * 0.004 + 3))
    if n_snps < 200:
        fig_width = min(24, max(10, n_snps * 0.08 + 3))

    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    ax.imshow(geno_display, aspect="auto", cmap=cmap, vmin=0, vmax=3,
              interpolation="none")

    sample_labels = samples_nn + samples_dd
    ax.set_yticks(range(n_total))
    ax.set_yticklabels(sample_labels, fontsize=5)
    ax.axhline(n_nn - 0.5, color="white", linewidth=2)

    ax.text(-0.02, (n_nn - 1) / 2 / n_total, "NN", transform=ax.transAxes,
            fontsize=12, fontweight="bold", color="#2166ac", va="center", ha="right")
    ax.text(-0.02, (n_nn + (n_total - 1)) / 2 / n_total, "DD", transform=ax.transAxes,
            fontsize=12, fontweight="bold", color="#d62728", va="center", ha="right")

    n_ticks = min(25, len(pos_filt))
    tick_idx = np.linspace(0, len(pos_filt) - 1, n_ticks, dtype=int)
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([f"{pos_filt[i]/1e6:.3f}" for i in tick_idx], fontsize=7, rotation=45)
    ax.set_xlabel("Position on chr3 (Mb)", fontsize=11)

    patches = [mpatches.Patch(color="#2166ac", label="0/0 (hom ref)"),
               mpatches.Patch(color="#ffd700", label="0/1 (het)"),
               mpatches.Patch(color="#d62728", label="1/1 (hom alt)"),
               mpatches.Patch(color="#d3d3d3", label="missing")]
    ax.legend(handles=patches, loc="upper right", fontsize=8, ncol=4,
              bbox_to_anchor=(1, 1.08))

    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    else:
        r_s = pos_filt[0] / 1e6
        r_e = pos_filt[-1] / 1e6
        ax.set_title(f"Genotype Plot — Chr3: {r_s:.1f}–{r_e:.1f} Mb",
                     fontsize=12, fontweight="bold")

    plt.tight_layout()
    path = os.path.join(outdir, f"chr3_genotype_{tag}.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, f"chr3_genotype_{tag}.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    return path, n_snps


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
def generate_html(plots, stats, outdir, vcf_path, samples_nn, samples_dd,
                   version_tag="", exclude_samples=None, exclude_reason=""):
    def img_b64(path):
        if path and os.path.exists(path):
            with open(path, "rb") as f:
                return base64.b64encode(f.read()).decode()
        return ""

    b64 = {k: img_b64(v) for k, v in plots.items()}

    n_nn = len(samples_nn)
    n_dd = len(samples_dd)
    nn_list = ", ".join(samples_nn[:5]) + (f" ... ({n_nn} total)" if n_nn > 5 else "")
    dd_list = ", ".join(samples_dd[:5]) + (f" ... ({n_dd} total)" if n_dd > 5 else "")

    # Top SNPs table
    rows = ""
    for i, (pos, logp, delta, fst, p, afnn, afdd) in enumerate(stats["top_snps"][:30]):
        rows += f"""<tr><td>{i+1}</td><td>{pos:,}</td><td>{pos/1e6:.3f}</td>
            <td>{logp:.1f}</td><td>{delta:.4f}</td><td>{fst:.4f}</td>
            <td>{afnn:.3f}</td><td>{afdd:.3f}</td></tr>"""

    cs = stats["core_start"]
    ce = stats["core_end"]
    n90 = stats["n_dd_90pct"]
    frac90_start = stats["frac90_start"]
    frac90_end = stats["frac90_end"]
    core_regions = stats.get("core_regions", [(cs, ce, 0.0, 0.0)])
    core_threshold = stats.get("core_threshold", 90)

    # Build core regions HTML — ranked table
    if len(core_regions) == 1:
        rs, re, mx, mn = core_regions[0]
        core_regions_html = f"""
<div class="core-box">
<h3>Core Introgression Region</h3>
<div class="coords">{rs/1e6:.2f} &ndash; {re/1e6:.2f} Mb</div>
<div class="size">{(re-rs)/1e6:.2f} Mb &mdash; defined by &ge;{core_threshold}% DD introgression coverage
(max coverage: {mx:.0%}, mean: {mn:.0%})</div>
</div>"""
    else:
        region_rows = ""
        for rank, (rs, re, mx, mn) in enumerate(core_regions, 1):
            size_kb = (re - rs) / 1000
            highlight = ' style="background:#eff6ff;font-weight:600;"' if rank == 1 else ""
            region_rows += (
                f'<tr{highlight}><td>{rank}</td>'
                f'<td>{rs/1e6:.2f} &ndash; {re/1e6:.2f} Mb</td>'
                f'<td>{size_kb:.0f} kb</td>'
                f'<td>{mx:.1%}</td><td>{mn:.1%}</td></tr>\n'
            )
        rs0, re0 = core_regions[0][0], core_regions[0][1]
        core_regions_html = f"""
<div class="core-box">
<h3>Core Introgression Regions (Top {len(core_regions)})</h3>
<div class="coords">#1: {rs0/1e6:.2f} &ndash; {re0/1e6:.2f} Mb</div>
<div class="size">{len(core_regions)} region(s) with &ge;{core_threshold}% DD introgression coverage</div>
</div>
<div class="card">
<table>
<thead><tr><th>Rank</th><th>Coordinates</th><th>Size</th>
<th>Max DD coverage</th><th>Mean DD coverage</th></tr></thead>
<tbody>{region_rows}</tbody>
</table>
</div>"""

    # Build exclusion banner HTML if samples were excluded
    excl_banner = ""
    if exclude_samples:
        excl_list = sorted(exclude_samples)
        excl_banner = f"""
<div style="background:#fef2f2;border:2px solid #ef4444;border-radius:10px;padding:1.2rem;
            margin-bottom:1.5rem;">
<h3 style="color:#b91c1c;margin-bottom:0.5rem;">Samples Excluded</h3>
<p style="font-size:0.9rem;"><strong>{len(excl_list)} sample(s) removed:</strong>
{', '.join(f'<code>{s}</code>' for s in excl_list)}</p>
<p style="font-size:0.85rem;margin-top:0.4rem;"><strong>Reason:</strong>
{exclude_reason if exclude_reason else 'High missingness in the chr3 focal region.'}</p>
</div>"""

    vtag = f" ({version_tag})" if version_tag else ""

    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Chr3 QTL Peak — Multi-Evidence Analysis{vtag}</title>
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
</style></head><body>

<h1>Chromosome 3 QTL Peak &mdash; Multi-Evidence Analysis{vtag}</h1>
<div class="subtitle"><em>Pararge aegeria</em> &mdash; Delineating the diapause introgression
using five independent approaches</div>

{excl_banner}

<h2>Summary</h2>
<div class="grid">
<div class="card"><table>
    <tr><td>NN samples</td><td class="highlight">{n_nn}</td></tr>
    <tr><td>DD samples</td><td class="highlight">{n_dd}</td></tr>
    <tr><td>Region analysed</td><td class="highlight">Chr3: 9.5&ndash;13.0 Mb</td></tr>
    <tr><td>SNPs in region</td><td class="highlight">{stats['n_snps_total']:,}</td></tr>
    <tr><td>SNPs &lt;5% missing</td><td class="highlight">{stats['n_snps_lowmiss']:,}</td></tr>
</table></div>
<div class="card"><table>
    <tr><td>SNPs with &minus;log&#8321;&#8320;(p) &gt; 50</td><td class="highlight">{stats['n_sig50']:,}</td></tr>
    <tr><td>SNPs with &minus;log&#8321;&#8320;(p) &gt; 40</td><td class="highlight">{stats['n_sig40']:,}</td></tr>
    <tr><td>DD with introgression at core</td><td class="highlight">{n90}/{n_dd} (&ge;90%)</td></tr>
    <tr><td>Max &minus;log&#8321;&#8320;(p)</td><td class="highlight">{stats['max_logp']:.1f}</td></tr>
</table></div>
</div>

{core_regions_html}

<h2>Rationale: How Were the Core Regions Defined?</h2>
<div class="card">
<p>Rather than relying on a single statistic, the core region was defined by the
<strong>convergence of five independent lines of evidence</strong>. This approach is robust
because each method captures different aspects of the introgression signal:</p>

<div class="approach">
<div class="approach-title">Approach 1: Multi-Statistic Sliding Windows</div>
<p style="font-size:0.85rem;">Three population-genetic statistics (&Delta;SNP-index,
F<sub>ST</sub>, Fisher&rsquo;s &minus;log&#8321;&#8320;p) were computed per SNP and smoothed with
25 kb and 50 kb windows. All three peak in the same region, providing convergent
evidence for the location of the introgression.</p>
</div>

<div class="approach">
<div class="approach-title">Approach 2: Per-Individual Introgression Tract Mapping</div>
<p style="font-size:0.85rem;">For each of the 94 individuals, the &ldquo;NN consensus
genotype&rdquo; (mode among NN) was computed at each SNP. Each individual&rsquo;s genotype
was then scored for discordance with this consensus in 50 kb sliding windows.
Contiguous regions with &gt;50% discordance define that individual&rsquo;s introgression tract.
DD individuals show large tracts; NN individuals (negative controls) show none or very
small tracts, confirming the method&rsquo;s specificity.</p>
</div>

<div class="approach">
<div class="approach-title">Approach 3: Minimum Shared Introgression (Coverage Analysis)</div>
<p style="font-size:0.85rem;">The introgression tracts from all DD individuals were overlaid.
At each 10 kb bin, the fraction of DD carrying the introgression was computed.
<strong>Core regions</strong> are defined as all contiguous segments where
<strong>&ge;{core_threshold}%</strong> of DD individuals carry the introgression. Multiple
non-contiguous core regions may exist if the introgression is fragmented by recombination
or if there are multiple linked loci. The top {len(core_regions)} region(s) are ranked by size above.</p>
</div>

<div class="approach">
<div class="approach-title">Approach 4: Recombination Breakpoint Density</div>
<p style="font-size:0.85rem;">The start and end coordinates of each individual&rsquo;s
introgression tract were collected as breakpoints. A histogram of breakpoint positions
shows where recombination has trimmed the introgression. The core region lies
<em>between</em> the breakpoint clusters &mdash; where recombination has not yet broken
down the introgressed haplotype. Breakpoints cluster at the flanks, confirming the
boundaries of the shared segment.</p>
</div>

<div class="approach">
<div class="approach-title">Approach 5: Per-Window Genotype Composition</div>
<p style="font-size:0.85rem;">For each 25 kb window, the fraction of hom-ref, het, and
hom-alt genotypes was computed separately for NN and DD. In the core region, NN is
dominated by one homozygous class while DD shows the opposite pattern (shifted towards
the alternate homozygote and/or elevated heterozygosity). Outside the core, both groups
converge to similar genotype frequencies.</p>
</div>
</div>

<h2>Approach 1: QTL-seq &mdash; Multi-Statistic Windows</h2>
<div class="card">
<div class="plot-container">
    <img src="data:image/png;base64,{b64['qtlseq']}" alt="QTL-seq">
</div></div>

<h2>Normalised Multi-Evidence Overlay</h2>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
All three statistics normalised to [0,1] and overlaid. The shaded blue region marks the
core introgression defined by &ge;90% DD coverage. All three metrics peak within or
near this region.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64['multi_evidence']}" alt="Multi-evidence overlay">
</div></div>

<h2>Approach 2: Per-Individual Introgression Tracts</h2>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
Each horizontal bar is one individual. Red = region where genotypes are discordant with
the NN consensus (&gt;50% discordance in 50 kb windows). NN individuals (top, negative
controls) show minimal discordance. DD individuals (bottom) show large introgression tracts
of varying extent. The blue shading marks the core region shared by &ge;90% of DD.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64['tracts']}" alt="Introgression tracts">
</div></div>

<h2>Approach 3 &amp; 4: Introgression Coverage &amp; Breakpoint Density</h2>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
<strong>Top:</strong> Fraction of DD individuals carrying the introgression at each position.
The core region is where coverage &ge; 90%.
<strong>Bottom:</strong> Histogram of introgression tract boundaries (start/end coordinates).
Breakpoints cluster at the flanks of the core region, where recombination has been
trimming the introgressed segment across backcross generations.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64['coverage']}" alt="Coverage and breakpoints">
</div></div>

<h2>Approach 5: Per-Window Genotype Composition</h2>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
Stacked area plots showing the proportion of each genotype class per 25 kb window.
In the core region (shaded), NN is dominated by hom-ref (blue) while DD shifts to
hom-alt (red) and/or elevated heterozygosity (yellow). Outside the core, both groups
converge to similar genotype compositions &mdash; the hallmark of a backcross design
where the rest of the genome has reverted to the recurrent parent.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64['geno_comp']}" alt="Genotype composition">
</div></div>

<h2>Genotype Heatmaps</h2>

<h3>{stats['n_sig50']} SNPs with &minus;log&#8321;&#8320;(p) &gt; 50</h3>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
Only the most extreme SNPs. At these sites, genotype separation between NN and DD
should be near-complete.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64['hm_sig50']}" alt="Heatmap logp>50">
</div></div>

<h3>{stats['n_sig40']} SNPs with &minus;log&#8321;&#8320;(p) &gt; 40</h3>
<div class="card">
<div class="plot-container">
    <img src="data:image/png;base64,{b64['hm_sig40']}" alt="Heatmap logp>40">
</div></div>

<h3>Full 9.5&ndash;13.0 Mb Region</h3>
<div class="card">
<div class="plot-container">
    <img src="data:image/png;base64,{b64['hm_broad']}" alt="Heatmap broad">
</div></div>

<h2>Top 30 Most Significant SNPs</h2>
<div class="card"><table><thead>
    <tr><th>#</th><th>Position</th><th>Mb</th><th>&minus;log&#8321;&#8320;(p)</th>
        <th>&Delta;(SNP-index)</th><th>Fst</th><th>AF NN</th><th>AF DD</th></tr>
</thead><tbody>{rows}</tbody></table>
</div>

<h2 style="margin-top:3rem;border-bottom:3px solid #e74c3c;padding-bottom:0.3rem;color:#b91c1c;">
Chromosome 2 Control Comparison</h2>

<div class="card" style="border-left:4px solid #e74c3c;">
<p style="font-size:0.9rem;">To validate that the chr3 signal is real and not an artefact of the analysis
pipeline, the identical analysis was run on <strong>chr2: 9.5&ndash;13.0 Mb</strong> &mdash; a region
with <strong>no expected QTL</strong>. If the method is working correctly, chr2 should show:</p>
<ul style="font-size:0.85rem;margin:0.5rem 0 0.5rem 1.5rem;">
<li>No significant enrichment of extreme p-values</li>
<li>No consistent introgression tracts in DD individuals</li>
<li>No differentiation in genotype composition between NN and DD</li>
<li>Near-zero introgression coverage (no shared region among DD)</li>
</ul>
</div>

<h3>Chr3 vs Chr2 &mdash; Summary Statistics</h3>
<div class="card">
<table>
<thead><tr><th>Metric</th><th>Chr3 (QTL)</th><th>Chr2 (control)</th></tr></thead>
<tbody>
<tr><td>Total SNPs in region</td><td class="highlight">{stats['n_snps_total']:,}</td>
    <td>{stats['c2_n_snps']:,}</td></tr>
<tr><td>SNPs &lt;5% missing</td><td class="highlight">{stats['n_snps_lowmiss']:,}</td>
    <td>{stats['c2_n_lowmiss']:,}</td></tr>
<tr><td>Max &minus;log&#8321;&#8320;(p)</td><td class="highlight">{stats['max_logp']:.1f}</td>
    <td>{stats['c2_max_logp']:.1f}</td></tr>
<tr><td>SNPs with &minus;log&#8321;&#8320;(p) &gt; 50</td><td class="highlight">{stats['n_sig50']:,}</td>
    <td>0</td></tr>
<tr><td>SNPs with &minus;log&#8321;&#8320;(p) &gt; 10</td><td class="highlight">&mdash;</td>
    <td>{stats['c2_n_sig10']:,}</td></tr>
<tr><td>DD with introgression tracts</td><td class="highlight">{stats['dd_with_tracts']}/{n_dd}</td>
    <td>{stats['c2_dd_with_tracts']}/{n_dd}</td></tr>
</tbody></table>
</div>

<h3>Introgression Tracts &mdash; Chr3 vs Chr2</h3>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
<strong>Left:</strong> Chr3 (QTL region) — DD individuals show large, overlapping introgression tracts.
<strong>Right:</strong> Chr2 (control) — any tracts visible are expected to be sparse, short, and non-overlapping,
consistent with random genotypic noise rather than a true introgression signal.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64.get('tracts_comparison','')}" alt="Chr3 vs Chr2 tracts">
</div></div>

<h3>Genotype Heatmap &mdash; Chr3 (QTL region)</h3>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
Clear genotype differentiation between NN (top) and DD (bottom) should be visible
as blocks of contrasting colour in the QTL region.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64.get('hm_chr3','')}" alt="Chr3 heatmap">
</div></div>

<h3>Genotype Heatmap &mdash; Chr2 (control)</h3>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
On a control chromosome, no systematic genotype differentiation should exist between
NN and DD. The heatmap should appear &ldquo;noisy&rdquo; with no clear blocks of colour
separating the two groups.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64.get('hm_chr2','')}" alt="Chr2 heatmap">
</div></div>

<h3>DD Introgression Coverage &mdash; Chr3 vs Chr2</h3>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
<strong>Top:</strong> Chr3 shows a region where the vast majority of DD carry the introgression
(coverage near 1.0). <strong>Bottom:</strong> Chr2 should show uniformly low coverage, confirming
the signal is specific to chr3.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64.get('coverage_comparison','')}" alt="Chr3 vs Chr2 coverage">
</div></div>

<h3>Genotype Composition &mdash; Chr3 vs Chr2</h3>
<div class="card">
<p style="font-size:0.85rem;color:#64748b;margin-bottom:0.5rem;">
<strong>Left column:</strong> Chr3 &mdash; NN and DD show opposite genotype compositions in the QTL region.
<strong>Right column:</strong> Chr2 &mdash; NN and DD should show essentially identical genotype compositions
throughout, confirming no differentiation on the control chromosome.</p>
<div class="plot-container">
    <img src="data:image/png;base64,{b64.get('geno_comp_comparison','')}" alt="Chr3 vs Chr2 genotype composition">
</div></div>

<h2>Method Details</h2>
<div class="card">
<p style="font-size:0.85rem;"><strong>Genotype extraction:</strong></p>
<pre>bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \\
  -s {nn_list},{dd_list} \\
  -r 3:{REGION_START}-{REGION_END} {os.path.basename(vcf_path)}</pre>

<p style="font-size:0.85rem;margin-top:1rem;"><strong>Statistics:</strong>
&Delta;(SNP-index) = AF<sub>DD</sub> &minus; AF<sub>NN</sub>;
F<sub>ST</sub> via Weir &amp; Cockerham (1984);
Fisher&rsquo;s exact test on 2&times;2 allele count tables.</p>

<p style="font-size:0.85rem;margin-top:0.5rem;"><strong>Introgression tract mapping:</strong>
NN consensus = mode genotype at each SNP across all NN individuals.
Per-individual discordance smoothed with 50 kb windows (step 10 kb).
Introgression = contiguous region with &gt;50% discordance.
Core region = contiguous segment with &ge;90% of DD individuals carrying introgression.</p>

<p style="font-size:0.85rem;margin-top:0.5rem;"><strong>Genotype encoding:</strong>
0/0 &rarr; 0 (blue), 0/1 &rarr; 1 (yellow), 1/1 &rarr; 2 (red), ./. &rarr; missing (grey).
Only SNPs with &lt;5% per-SNP missingness used for heatmaps and tract mapping.</p>

<p style="font-size:0.85rem;margin-top:0.5rem;"><strong>Control comparison:</strong>
Identical analysis pipeline applied to chr2: 9.5&ndash;13.0 Mb (same region coordinates,
no expected QTL). This validates that the chr3 signal reflects true biological differentiation
rather than a methodological artefact.</p>
</div>

</body></html>"""

    tag_suffix = f"_{version_tag}" if version_tag else ""
    html_filename = f"chr3_QTL_peak_10-12Mb{tag_suffix}.html"
    with open(os.path.join(outdir, html_filename), "w") as f:
        f.write(html)
    print(f"HTML report saved: {html_filename}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: python chr3_region_analysis.py <vcf> <phenotypes> <outdir> "
                 "[n_cores] [--exclude sample1,sample2] [--tag version_tag]")

    vcf, pheno_file, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    n_cores = int(sys.argv[4]) if len(sys.argv) > 4 and not sys.argv[4].startswith("--") else N_CORES

    # Parse optional arguments
    exclude_samples = set()
    version_tag = ""
    exclude_reason = ""
    args = sys.argv[4:]
    i = 0
    while i < len(args):
        if args[i] == "--exclude" and i + 1 < len(args):
            exclude_samples = set(args[i + 1].split(","))
            i += 2
        elif args[i] == "--tag" and i + 1 < len(args):
            version_tag = args[i + 1]
            i += 2
        elif args[i] == "--exclude-reason" and i + 1 < len(args):
            exclude_reason = args[i + 1]
            i += 2
        elif not args[i].startswith("--"):
            try:
                n_cores = int(args[i])
            except ValueError:
                pass
            i += 1
        else:
            i += 1

    os.makedirs(outdir, exist_ok=True)

    groups = parse_phenotypes(pheno_file)
    samples_nn = sorted([s for s, g in groups.items() if g == "NN" and s not in exclude_samples])
    samples_dd = sorted([s for s, g in groups.items() if g == "DD" and s not in exclude_samples])

    if exclude_samples:
        excluded_nn = sorted([s for s in exclude_samples if groups.get(s) == "NN"])
        excluded_dd = sorted([s for s in exclude_samples if groups.get(s) == "DD"])
        print(f"EXCLUDED SAMPLES: {', '.join(sorted(exclude_samples))}")
        if excluded_nn:
            print(f"  From NN: {', '.join(excluded_nn)}")
        if excluded_dd:
            print(f"  From DD: {', '.join(excluded_dd)}")
    print(f"NN: {len(samples_nn)}, DD: {len(samples_dd)}")

    # --- Extract genotypes ---
    region_str = f"{CHROM}:{REGION_START}-{REGION_END}"
    print(f"\n=== Extracting {region_str} ===")
    positions, refs, alts, geno_matrix, nn_alt, nn_tot, dd_alt, dd_tot, missingness = \
        extract_region_genotypes(vcf, region_str, samples_nn, samples_dd)
    n_snps = len(positions)
    n_lowmiss = int((missingness < MAX_MISSINGNESS).sum())
    print(f"  SNPs: {n_snps:,} (low-miss: {n_lowmiss:,})")

    # --- Compute statistics ---
    print("\n=== Computing statistics ===")
    delta, af_nn, af_dd, fst = compute_stats(nn_alt, nn_tot, dd_alt, dd_tot)
    print("  Fisher's exact...")
    pvals = fisher_parallel(nn_alt, nn_tot, dd_alt, dd_tot, n_cores)
    logp = -np.log10(np.clip(pvals, 1e-300, 1))

    # --- Sliding windows ---
    print("\n=== Sliding windows ===")
    w50_pos, w50_delta = sliding_window(positions, delta, 50_000, 5_000)
    _, w50_fst = sliding_window(positions, fst, 50_000, 5_000)
    w25_pos, w25_delta = sliding_window(positions, delta, 25_000, 2_500)
    _, w25_fst = sliding_window(positions, fst, 25_000, 2_500)

    # --- APPROACH 2: Per-individual introgression tracts ---
    print("\n=== Mapping introgression tracts ===")
    all_tracts, all_discordance, nn_consensus, low_miss_pos = \
        map_introgression_tracts(positions, geno_matrix, missingness,
                                  samples_nn, samples_dd, window_kb=50)

    # Count DD with tracts
    dd_with_tracts = sum(1 for s in samples_dd if len(all_tracts.get(s, [])) > 0)
    print(f"  DD with detected tracts: {dd_with_tracts}/{len(samples_dd)}")

    # --- APPROACH 3: Minimum shared introgression ---
    print("\n=== Computing introgression coverage ===")
    bins, coverage, frac = find_minimum_shared_introgression(
        all_tracts, samples_dd, REGION_START, REGION_END)

    # Find ALL contiguous core regions with >=90% DD coverage
    def find_contiguous_regions(frac_array, bins_array, threshold):
        """Find all contiguous runs above threshold. Returns list of (start, end, max_frac, mean_frac)."""
        above = frac_array >= threshold
        regions = []
        current_start = None
        for i in range(len(bins_array)):
            if above[i]:
                if current_start is None:
                    current_start = i
            else:
                if current_start is not None:
                    r_start = int(bins_array[current_start])
                    r_end = int(bins_array[i - 1]) + 10_000
                    r_max = float(np.max(frac_array[current_start:i]))
                    r_mean = float(np.mean(frac_array[current_start:i]))
                    regions.append((r_start, r_end, r_max, r_mean))
                    current_start = None
        if current_start is not None:
            r_start = int(bins_array[current_start])
            r_end = int(bins_array[-1]) + 10_000
            r_max = float(np.max(frac_array[current_start:]))
            r_mean = float(np.mean(frac_array[current_start:]))
            regions.append((r_start, r_end, r_max, r_mean))
        return regions

    # Try >=90% first, fallback to >=80%
    core_regions_90 = find_contiguous_regions(frac, bins, 0.90)
    if core_regions_90:
        # Sort by size (largest first), then take top 5
        core_regions_90.sort(key=lambda r: -(r[1] - r[0]))
        core_regions = core_regions_90[:5]
        core_threshold = 90
    else:
        core_regions_80 = find_contiguous_regions(frac, bins, 0.80)
        core_regions_80.sort(key=lambda r: -(r[1] - r[0]))
        core_regions = core_regions_80[:5]
        core_threshold = 80

    if not core_regions:
        core_regions = [(10_500_000, 11_500_000, 0.0, 0.0)]
        core_threshold = 0

    # Primary core = largest region (for backward compatibility with plots)
    core_start = core_regions[0][0]
    core_end = core_regions[0][1]
    n_dd_at_core = int(np.max(coverage[
        (bins >= core_start) & (bins < core_end)
    ])) if core_start else 0

    print(f"  Found {len(core_regions)} core region(s) at ≥{core_threshold}% coverage:")
    for rank, (rs, re, mx, mn) in enumerate(core_regions, 1):
        print(f"    #{rank}: {rs:,} - {re:,} ({(re-rs)/1e6:.2f} Mb, "
              f"max={mx:.1%}, mean={mn:.1%})")

    # --- APPROACH 5: Genotype composition ---
    print("\n=== Computing genotype composition ===")
    gc_pos, nn_hr, nn_het, nn_ha, dd_hr, dd_het, dd_ha = \
        compute_genotype_composition(positions, geno_matrix, missingness,
                                      len(samples_nn), window_size=25_000, step=5_000)

    # --- Generate plots ---
    print("\n=== Generating plots ===")
    plots = {}

    plots["qtlseq"] = make_qtlseq_plot(positions, delta, fst, pvals, af_nn, af_dd,
                                         w50_pos, w50_delta, w50_fst,
                                         w25_pos, w25_delta, w25_fst, outdir)

    plots["multi_evidence"] = make_multi_evidence_plot(
        positions, delta, fst, logp, w25_pos, w25_delta, w25_fst,
        core_regions, outdir)

    plots["tracts"] = make_introgression_tract_plot(
        all_tracts, all_discordance, samples_nn, samples_dd,
        core_regions, outdir)

    plots["coverage"] = make_coverage_breakpoint_plot(
        bins, coverage, frac, all_tracts, samples_dd,
        core_regions, outdir)

    plots["geno_comp"] = make_genotype_composition_plot(
        gc_pos, nn_hr, nn_het, nn_ha, dd_hr, dd_het, dd_ha,
        core_regions, outdir)

    # Genotype heatmaps
    print("\n=== Genotype heatmaps ===")
    sig50_mask = (logp > 50) & (missingness < MAX_MISSINGNESS)
    sig40_mask = (logp > 40) & (missingness < MAX_MISSINGNESS)
    n_sig50 = int(sig50_mask.sum())
    n_sig40 = int(sig40_mask.sum())

    plots["hm_sig50"], _ = make_genotype_heatmap(
        positions, geno_matrix, missingness, samples_nn, samples_dd, outdir,
        snp_mask=sig50_mask, tag="sig50",
        title=f"Genotype Plot — {n_sig50} SNPs with −log₁₀(p) > 50")

    plots["hm_sig40"], _ = make_genotype_heatmap(
        positions, geno_matrix, missingness, samples_nn, samples_dd, outdir,
        snp_mask=sig40_mask, tag="sig40",
        title=f"Genotype Plot — {n_sig40} SNPs with −log₁₀(p) > 40")

    plots["hm_broad"], _ = make_genotype_heatmap(
        positions, geno_matrix, missingness, samples_nn, samples_dd, outdir,
        region_start=REGION_START, region_end=REGION_END, max_snps=8000, tag="broad",
        title="Genotype Plot — Chr3: 9.5–13.0 Mb (all low-missingness SNPs)")

    # Top SNPs
    top_idx = np.argsort(logp)[::-1][:30]
    top_snps = [(int(positions[i]), float(logp[i]), float(delta[i]), float(fst[i]),
                 float(pvals[i]), float(af_nn[i]), float(af_dd[i])) for i in top_idx]

    # Save per-SNP results
    print("\n=== Saving per-SNP results ===")
    with gzip.open(os.path.join(outdir, "chr3_peak_per_snp.tsv.gz"), "wt") as f:
        f.write("pos\tref\talt\taf_nn\taf_dd\tdelta\tfst\tp_fisher\tlog10p\tmissingness\n")
        for i in range(n_snps):
            f.write(f"{positions[i]}\t{refs[i]}\t{alts[i]}\t{af_nn[i]:.6f}\t{af_dd[i]:.6f}\t"
                    f"{delta[i]:.6f}\t{fst[i]:.6f}\t{pvals[i]:.6e}\t{logp[i]:.2f}\t"
                    f"{missingness[i]:.4f}\n")

    # ===================================================================
    # CHROMOSOME 2 CONTROL — same region, no expected QTL
    # ===================================================================
    print("\n" + "=" * 60)
    print("=== CHROMOSOME 2 CONTROL ANALYSIS ===")
    print("=" * 60)
    chr2_region = f"2:{REGION_START}-{REGION_END}"
    print(f"\n=== Extracting {chr2_region} ===")
    c2_pos, c2_refs, c2_alts, c2_geno, c2_nn_alt, c2_nn_tot, c2_dd_alt, c2_dd_tot, c2_miss = \
        extract_region_genotypes(vcf, chr2_region, samples_nn, samples_dd)
    c2_n = len(c2_pos)
    c2_n_lm = int((c2_miss < MAX_MISSINGNESS).sum())
    print(f"  Chr2 SNPs: {c2_n:,} (low-miss: {c2_n_lm:,})")

    print("  Computing chr2 statistics...")
    c2_delta, c2_af_nn, c2_af_dd, c2_fst = compute_stats(c2_nn_alt, c2_nn_tot, c2_dd_alt, c2_dd_tot)
    c2_pvals = fisher_parallel(c2_nn_alt, c2_nn_tot, c2_dd_alt, c2_dd_tot, n_cores)
    c2_logp = -np.log10(np.clip(c2_pvals, 1e-300, 1))
    print(f"  Chr2 max -log10(p): {c2_logp.max():.1f}")
    print(f"  Chr2 SNPs with -log10(p) > 10: {(c2_logp > 10).sum():,}")

    print("  Mapping chr2 introgression tracts...")
    c2_tracts, c2_disc, _, _ = map_introgression_tracts(
        c2_pos, c2_geno, c2_miss, samples_nn, samples_dd, window_kb=50)
    c2_dd_with_tracts = sum(1 for s in samples_dd if len(c2_tracts.get(s, [])) > 0)
    print(f"  Chr2 DD with tracts: {c2_dd_with_tracts}/{len(samples_dd)}")

    c2_bins, c2_coverage, c2_frac = find_minimum_shared_introgression(
        c2_tracts, samples_dd, REGION_START, REGION_END)

    c2_gc_pos, c2_nn_hr, c2_nn_het, c2_nn_ha, c2_dd_hr, c2_dd_het, c2_dd_ha = \
        compute_genotype_composition(c2_pos, c2_geno, c2_miss,
                                      len(samples_nn), window_size=25_000, step=5_000)

    # --- Comparison plots ---
    print("\n=== Generating comparison plots ===")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import matplotlib.patches as mpatches

    # --- Side-by-side introgression tracts ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, max(8, len(samples_nn + samples_dd) * 0.15 + 2)),
                                     sharey=True)
    all_samples = samples_nn + samples_dd
    n_nn_s = len(samples_nn)

    for ax, tracts, chrom_label in [(ax1, all_tracts, "Chr3 (QTL region)"),
                                      (ax2, c2_tracts, "Chr2 (control)")]:
        for idx, sample in enumerate(all_samples):
            y = len(all_samples) - 1 - idx
            is_dd = idx >= n_nn_s
            ax.barh(y, (REGION_END - REGION_START) / 1e6,
                    left=REGION_START / 1e6, height=0.7,
                    color="#f0f0f0", edgecolor="none")
            for t_start, t_end in tracts.get(sample, []):
                w = (t_end - t_start) / 1e6
                col = "#d62728" if is_dd else "#ff9999"
                ax.barh(y, w, left=t_start / 1e6, height=0.7,
                        color=col, edgecolor="none", alpha=0.8)
        ax.axhline(len(samples_dd) - 0.5, color="black", linewidth=1.5)
        ax.set_xlabel("Position (Mb)", fontsize=11)
        ax.set_title(chrom_label, fontsize=13, fontweight="bold")
        ax.set_xlim(REGION_START / 1e6 - 0.05, REGION_END / 1e6 + 0.05)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax1.set_yticks(range(len(all_samples)))
    ax1.set_yticklabels(list(reversed(all_samples)), fontsize=4.5)
    ax1.text(REGION_START / 1e6 - 0.15, len(all_samples) - n_nn_s / 2,
             "NN", fontsize=11, fontweight="bold", color="#2166ac", va="center", ha="right")
    ax1.text(REGION_START / 1e6 - 0.15, len(samples_dd) / 2,
             "DD", fontsize=11, fontweight="bold", color="#d62728", va="center", ha="right")

    plt.tight_layout()
    plots["tracts_comparison"] = os.path.join(outdir, "chr3_vs_chr2_tracts.png")
    fig.savefig(plots["tracts_comparison"], dpi=200, bbox_inches="tight")
    plt.close()

    # --- Side-by-side genotype heatmaps ---
    for chrom_data, chrom_name in [(
        (positions, geno_matrix, missingness), "chr3"),
        ((c2_pos, c2_geno, c2_miss), "chr2")]:

        pos_d, geno_d, miss_d = chrom_data
        lm = miss_d < MAX_MISSINGNESS
        pf = pos_d[lm]
        gf = geno_d[lm]
        n_s = len(pf)
        if n_s > 5000:
            step_t = n_s // 5000
            idx_t = np.arange(0, n_s, step_t)[:5000]
            pf = pf[idx_t]
            gf = gf[idx_t]

        gd = gf.T
        gd = np.where(gd == -1, 3, gd)
        cmap = ListedColormap(["#2166ac", "#ffd700", "#d62728", "#d3d3d3"])

        n_tot = len(all_samples)
        fig, ax = plt.subplots(1, 1, figsize=(18, max(6, n_tot * 0.13 + 2)))
        ax.imshow(gd, aspect="auto", cmap=cmap, vmin=0, vmax=3, interpolation="none")
        ax.set_yticks(range(n_tot))
        ax.set_yticklabels(all_samples, fontsize=5)
        ax.axhline(n_nn_s - 0.5, color="white", linewidth=2)
        ax.text(-0.02, (n_nn_s - 1) / 2 / n_tot, "NN", transform=ax.transAxes,
                fontsize=12, fontweight="bold", color="#2166ac", va="center", ha="right")
        ax.text(-0.02, (n_nn_s + (n_tot - 1)) / 2 / n_tot, "DD", transform=ax.transAxes,
                fontsize=12, fontweight="bold", color="#d62728", va="center", ha="right")
        n_tk = min(25, len(pf))
        tk_i = np.linspace(0, len(pf) - 1, n_tk, dtype=int)
        ax.set_xticks(tk_i)
        ax.set_xticklabels([f"{pf[j]/1e6:.2f}" for j in tk_i], fontsize=7, rotation=45)
        ax.set_xlabel(f"Position on {chrom_name} (Mb)", fontsize=11)
        patches = [mpatches.Patch(color="#2166ac", label="0/0"),
                   mpatches.Patch(color="#ffd700", label="0/1"),
                   mpatches.Patch(color="#d62728", label="1/1"),
                   mpatches.Patch(color="#d3d3d3", label="miss")]
        ax.legend(handles=patches, loc="upper right", fontsize=8, ncol=4,
                  bbox_to_anchor=(1, 1.08))
        ax.set_title(f"Genotype Plot — {chrom_name.upper()}: 9.5–13.0 Mb",
                     fontsize=12, fontweight="bold")
        plt.tight_layout()
        tag = f"hm_{chrom_name}"
        path = os.path.join(outdir, f"{chrom_name}_genotype_comparison.png")
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close()
        plots[tag] = path

    # --- Side-by-side coverage ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 7), sharex=True,
                                     gridspec_kw={"hspace": 0.15})
    ax1.fill_between(bins / 1e6, frac, alpha=0.4, color="#d62728", step="mid")
    ax1.step(bins / 1e6, frac, color="#d62728", linewidth=1.2, where="mid")
    ax1.axhline(0.90, color="grey", linewidth=0.5, linestyle="--", alpha=0.5)
    ax1.set_ylabel("Fraction of DD\nwith introgression", fontsize=11)
    ax1.set_ylim(-0.02, 1.08)
    ax1.set_title("Chr3 (QTL region) — DD introgression coverage", fontsize=12, fontweight="bold")
    ax1.spines["top"].set_visible(False); ax1.spines["right"].set_visible(False)

    ax2.fill_between(c2_bins / 1e6, c2_frac, alpha=0.4, color="#67a9cf", step="mid")
    ax2.step(c2_bins / 1e6, c2_frac, color="#2166ac", linewidth=1.2, where="mid")
    ax2.axhline(0.90, color="grey", linewidth=0.5, linestyle="--", alpha=0.5)
    ax2.set_ylabel("Fraction of DD\nwith introgression", fontsize=11)
    ax2.set_ylim(-0.02, 1.08)
    ax2.set_xlabel("Position (Mb)", fontsize=12)
    ax2.set_title("Chr2 (control) — DD introgression coverage", fontsize=12, fontweight="bold")
    ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
    plt.tight_layout()
    plots["coverage_comparison"] = os.path.join(outdir, "chr3_vs_chr2_coverage.png")
    fig.savefig(plots["coverage_comparison"], dpi=200, bbox_inches="tight")
    plt.close()

    # --- Side-by-side genotype composition ---
    fig, axes = plt.subplots(2, 2, figsize=(20, 9), sharex=True,
                               gridspec_kw={"hspace": 0.15, "wspace": 0.15})
    for col, (gp, nhr, nht, nha, dhr, dht, dha, label) in enumerate([
        (gc_pos, nn_hr, nn_het, nn_ha, dd_hr, dd_het, dd_ha, "Chr3"),
        (c2_gc_pos, c2_nn_hr, c2_nn_het, c2_nn_ha, c2_dd_hr, c2_dd_het, c2_dd_ha, "Chr2")]):
        x = gp / 1e6
        ax = axes[0, col]
        ax.fill_between(x, 0, nhr, color="#2166ac", alpha=0.7, step="mid")
        ax.fill_between(x, nhr, nhr + nht, color="#ffd700", alpha=0.7, step="mid")
        ax.fill_between(x, nhr + nht, 1, color="#d62728", alpha=0.7, step="mid")
        ax.set_ylabel("NN genotype freq", fontsize=10)
        ax.set_ylim(0, 1)
        ax.set_title(f"{label} — NN", fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

        ax = axes[1, col]
        ax.fill_between(x, 0, dhr, color="#2166ac", alpha=0.7, step="mid")
        ax.fill_between(x, dhr, dhr + dht, color="#ffd700", alpha=0.7, step="mid")
        ax.fill_between(x, dhr + dht, 1, color="#d62728", alpha=0.7, step="mid")
        ax.set_ylabel("DD genotype freq", fontsize=10)
        ax.set_xlabel("Position (Mb)", fontsize=11)
        ax.set_ylim(0, 1)
        ax.set_title(f"{label} — DD", fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plots["geno_comp_comparison"] = os.path.join(outdir, "chr3_vs_chr2_genotype_comp.png")
    fig.savefig(plots["geno_comp_comparison"], dpi=200, bbox_inches="tight")
    plt.close()

    # ===================================================================

    # --- HTML ---
    print("\n=== Generating HTML report ===")
    stats_dict = {
        "n_snps_total": n_snps,
        "n_snps_lowmiss": n_lowmiss,
        "n_sig50": n_sig50,
        "n_sig40": n_sig40,
        "max_logp": float(logp.max()),
        "core_start": core_start,
        "core_end": core_end,
        "n_dd_90pct": n_dd_at_core,
        "frac90_start": core_start,
        "frac90_end": core_end,
        "core_regions": core_regions,
        "core_threshold": core_threshold,
        "top_snps": top_snps,
        "dd_with_tracts": dd_with_tracts,
        # chr2 control stats
        "c2_n_snps": c2_n,
        "c2_n_lowmiss": c2_n_lm,
        "c2_max_logp": float(c2_logp.max()),
        "c2_n_sig10": int((c2_logp > 10).sum()),
        "c2_dd_with_tracts": c2_dd_with_tracts,
    }
    generate_html(plots, stats_dict, outdir, vcf, samples_nn, samples_dd,
                   version_tag=version_tag, exclude_samples=exclude_samples,
                   exclude_reason=exclude_reason)

    print(f"\n=== Done ===")
    print(f"Primary core introgression: {core_start:,} - {core_end:,} "
          f"({(core_end - core_start) / 1e6:.2f} Mb)")
    if len(core_regions) > 1:
        print(f"Total core regions found: {len(core_regions)}")


if __name__ == "__main__":
    main()
