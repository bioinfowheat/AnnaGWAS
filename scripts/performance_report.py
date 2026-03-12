#!/usr/bin/env python3
"""
Pipeline Performance Report Generator
Parses Snakemake log files and file sizes to produce a standalone HTML report
on BAM and VCF generation run times.

Usage:
    python performance_report.py <snakemake_log> <results_dir> <output_html>
"""

import sys
import re
import os
import glob
from datetime import datetime, timedelta
from pathlib import Path

# ── CLI args ──────────────────────────────────────────────────────────────────
if len(sys.argv) != 4:
    sys.exit("Usage: python performance_report.py <snakemake_log> <results_dir> <output_html>")

log_path = sys.argv[1]
results_dir = sys.argv[2]
output_html = sys.argv[3]

# ── Parse Snakemake log ──────────────────────────────────────────────────────
with open(log_path) as fh:
    log_lines = fh.readlines()

# Regex for timestamp lines: [Thu Mar  5 10:02:04 2026]
ts_re = re.compile(r'^\[(.+?)\]')
# Rule start line:  localrule <rulename>:
rule_re = re.compile(r'^\s*localrule\s+(\w+):')
# Sample wildcard:  wildcards: sample=nr1
sample_re = re.compile(r'^\s*wildcards:\s+sample=(\S+)')
# Finished line:  Finished job XX.
finish_re = re.compile(r'^\s*Finished job (\d+)\.')
# Jobid line:  jobid: 25
jobid_re = re.compile(r'^\s*jobid:\s+(\d+)')

def parse_ts(s):
    """Parse 'Thu Mar  5 10:02:04 2026' → datetime."""
    return datetime.strptime(s.strip(), "%a %b %d %H:%M:%S %Y")

# Build events: for each job, record rule, sample, start, end
# Pass 1: collect job blocks — a block starts at a timestamp+rule and contains
# jobid + optional wildcards. We collect all fields first, then store by jobid.
job_starts = {}   # jobid → {rule, sample, start_ts}
current_ts = None
current_rule = None
current_sample = None
current_jobid = None

for line in log_lines:
    m = ts_re.match(line)
    if m:
        # Before starting a new block, save any pending block
        if current_jobid is not None and current_rule:
            job_starts[current_jobid] = {
                "rule": current_rule,
                "sample": current_sample,
                "start": current_ts,
            }
        current_ts = parse_ts(m.group(1))
        current_rule = None
        current_sample = None
        current_jobid = None
        continue

    m = rule_re.match(line)
    if m:
        current_rule = m.group(1)
        current_sample = None
        current_jobid = None
        continue

    m = jobid_re.match(line)
    if m and current_rule:
        current_jobid = int(m.group(1))
        continue

    m = sample_re.match(line)
    if m and current_rule:
        current_sample = m.group(1)
        # Now we have both jobid and sample — save immediately
        if current_jobid is not None:
            job_starts[current_jobid] = {
                "rule": current_rule,
                "sample": current_sample,
                "start": current_ts,
            }
        continue

# Save final pending block if any
if current_jobid is not None and current_rule:
    job_starts[current_jobid] = {
        "rule": current_rule,
        "sample": current_sample,
        "start": current_ts,
    }

# Pass 2: find job finishes
job_events = []  # list of {rule, sample, start, end, duration}
current_ts = None
for line in log_lines:
    m = ts_re.match(line)
    if m:
        current_ts = parse_ts(m.group(1))
        continue

    m = finish_re.match(line)
    if m and current_ts:
        jid = int(m.group(1))
        if jid in job_starts:
            info = job_starts[jid]
            dur = current_ts - info["start"]
            job_events.append({
                "rule": info["rule"],
                "sample": info["sample"],
                "start": info["start"],
                "end": current_ts,
                "duration": dur,
            })

# ── Collect file sizes ────────────────────────────────────────────────────────
def file_size_gb(path):
    """Return file size in GB, or None if missing."""
    try:
        return os.path.getsize(path) / (1024**3)
    except OSError:
        return None

def file_size_human(gb):
    if gb is None:
        return "—"
    if gb >= 1.0:
        return f"{gb:.1f} GB"
    return f"{gb*1024:.0f} MB"

# Per-sample FASTQ sizes (R1+R2 combined)
samples_seen = sorted(set(e["sample"] for e in job_events if e["sample"] and e["rule"] == "bwa_mem2_align"))
fastq_sizes = {}
for s in samples_seen:
    r1 = glob.glob(os.path.join(results_dir, f"trimmed/{s}_R1.trimmed.fastq.gz"))
    r2 = glob.glob(os.path.join(results_dir, f"trimmed/{s}_R2.trimmed.fastq.gz"))
    total = 0.0
    for f in r1 + r2:
        sz = file_size_gb(f)
        if sz:
            total += sz
    fastq_sizes[s] = total

# Per-sample dedup BAM sizes
bam_sizes = {}
for s in samples_seen:
    bam = os.path.join(results_dir, f"aligned/{s}.dedup.bam")
    bam_sizes[s] = file_size_gb(bam)

# VCF sizes
raw_vcf_size = file_size_gb(os.path.join(results_dir, "variants/raw_variants.vcf.gz"))
filt_vcf_size = file_size_gb(os.path.join(results_dir, "variants/filtered_variants.vcf.gz"))

# ── Organise per-sample timing ────────────────────────────────────────────────
rules_of_interest = ["bwa_mem2_align", "mark_duplicates", "index_bam"]

per_sample = {}  # sample → {rule → duration}
for e in job_events:
    if e["rule"] in rules_of_interest and e["sample"]:
        per_sample.setdefault(e["sample"], {})[e["rule"]] = e["duration"]

# Variant calling (bcftools_call) — many regions, aggregate
vc_events = [e for e in job_events if e["rule"] == "bcftools_call"]
if vc_events:
    vc_start = min(e["start"] for e in vc_events)
    vc_end = max(e["end"] for e in vc_events)
    vc_wall = vc_end - vc_start
    vc_total_cpu = sum((e["duration"] for e in vc_events), timedelta())
    vc_n_regions = len(vc_events)
    vc_avg_per_region = vc_total_cpu / vc_n_regions
else:
    vc_wall = vc_total_cpu = timedelta()
    vc_n_regions = 0
    vc_avg_per_region = timedelta()

# Concat + filter
concat_events = [e for e in job_events if e["rule"] == "concat_vcfs"]
filter_events = [e for e in job_events if e["rule"] == "filter_variants"]
concat_dur = concat_events[0]["duration"] if concat_events else timedelta()
filter_dur = filter_events[0]["duration"] if filter_events else timedelta()

# Pipeline wall time
all_starts = [e["start"] for e in job_events]
all_ends = [e["end"] for e in job_events]
pipeline_start = min(all_starts) if all_starts else datetime.now()
pipeline_end = max(all_ends) if all_ends else datetime.now()
pipeline_wall = pipeline_end - pipeline_start

# ── Helper: format timedelta ──────────────────────────────────────────────────
def fmt_td(td):
    total_sec = int(td.total_seconds())
    if total_sec < 0:
        return "—"
    h, rem = divmod(total_sec, 3600)
    m, s = divmod(rem, 60)
    if h > 0:
        return f"{h}h {m:02d}m {s:02d}s"
    return f"{m}m {s:02d}s"

# ── Build HTML ────────────────────────────────────────────────────────────────
def make_bar(value_sec, max_sec, color="#3498db"):
    if max_sec == 0:
        return ""
    pct = min(value_sec / max_sec * 100, 100)
    return f'<div style="background:{color};height:18px;width:{pct:.1f}%;border-radius:3px;"></div>'

# Per-sample data rows for BAM table
bam_rows = []
align_durations = []
for s in samples_seen:
    d = per_sample.get(s, {})
    align = d.get("bwa_mem2_align", timedelta())
    dedup = d.get("mark_duplicates", timedelta())
    index = d.get("index_bam", timedelta())
    total_bam = align + dedup + index
    align_durations.append(align.total_seconds())
    bam_rows.append({
        "sample": s,
        "fastq_size": fastq_sizes.get(s, 0),
        "align": align,
        "dedup": dedup,
        "index": index,
        "total": total_bam,
        "bam_size": bam_sizes.get(s),
    })

# Averages
n = len(bam_rows)
if n > 0:
    avg_align = timedelta(seconds=sum(r["align"].total_seconds() for r in bam_rows) / n)
    avg_dedup = timedelta(seconds=sum(r["dedup"].total_seconds() for r in bam_rows) / n)
    avg_index = timedelta(seconds=sum(r["index"].total_seconds() for r in bam_rows) / n)
    avg_total = timedelta(seconds=sum(r["total"].total_seconds() for r in bam_rows) / n)
    avg_fastq = sum(r["fastq_size"] for r in bam_rows) / n
    avg_bam_vals = [r["bam_size"] for r in bam_rows if r["bam_size"] is not None]
    avg_bam = sum(avg_bam_vals) / len(avg_bam_vals) if avg_bam_vals else None
else:
    avg_align = avg_dedup = avg_index = avg_total = timedelta()
    avg_fastq = 0
    avg_bam = None

max_align_sec = max(align_durations) if align_durations else 1

# BAM table rows HTML
bam_table_rows = ""
for r in bam_rows:
    bar = make_bar(r["align"].total_seconds(), max_align_sec)
    bam_table_rows += f"""
    <tr>
        <td>{r['sample']}</td>
        <td>{file_size_human(r['fastq_size'])}</td>
        <td>{fmt_td(r['align'])}</td>
        <td>{fmt_td(r['dedup'])}</td>
        <td>{fmt_td(r['index'])}</td>
        <td><strong>{fmt_td(r['total'])}</strong></td>
        <td>{file_size_human(r['bam_size'])}</td>
        <td style="width:200px">{bar}</td>
    </tr>"""

bam_table_rows += f"""
    <tr class="avg-row">
        <td><strong>Average</strong></td>
        <td>{file_size_human(avg_fastq)}</td>
        <td><strong>{fmt_td(avg_align)}</strong></td>
        <td><strong>{fmt_td(avg_dedup)}</strong></td>
        <td><strong>{fmt_td(avg_index)}</strong></td>
        <td><strong>{fmt_td(avg_total)}</strong></td>
        <td>{file_size_human(avg_bam)}</td>
        <td></td>
    </tr>"""

# VCF total wall time
vcf_total_wall = vc_wall + concat_dur + filter_dur

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Pipeline Performance Report</title>
<style>
    body {{
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        max-width: 1100px;
        margin: 30px auto;
        padding: 0 20px;
        color: #2c3e50;
        background: #f8f9fa;
    }}
    h1 {{
        color: #2c3e50;
        border-bottom: 3px solid #3498db;
        padding-bottom: 10px;
    }}
    h2 {{
        color: #34495e;
        margin-top: 35px;
        border-bottom: 1px solid #bdc3c7;
        padding-bottom: 5px;
    }}
    .summary-box {{
        background: #fff;
        border: 1px solid #ddd;
        border-left: 4px solid #3498db;
        border-radius: 4px;
        padding: 15px 20px;
        margin: 20px 0;
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
        gap: 10px;
    }}
    .summary-item {{
        text-align: center;
    }}
    .summary-item .label {{
        font-size: 0.85em;
        color: #7f8c8d;
        text-transform: uppercase;
    }}
    .summary-item .value {{
        font-size: 1.4em;
        font-weight: 700;
        color: #2c3e50;
    }}
    table {{
        width: 100%;
        border-collapse: collapse;
        background: #fff;
        border-radius: 4px;
        overflow: hidden;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        margin: 15px 0;
    }}
    th {{
        background: #34495e;
        color: #fff;
        padding: 10px 12px;
        text-align: left;
        font-weight: 600;
        font-size: 0.9em;
    }}
    td {{
        padding: 8px 12px;
        border-bottom: 1px solid #ecf0f1;
        font-size: 0.9em;
    }}
    tr:hover {{
        background: #f1f8ff;
    }}
    .avg-row {{
        background: #eaf2f8 !important;
        border-top: 2px solid #3498db;
    }}
    .avg-row:hover {{
        background: #d4e6f1 !important;
    }}
    .section-note {{
        font-size: 0.85em;
        color: #7f8c8d;
        margin: 5px 0 15px 0;
        font-style: italic;
    }}
    .footer {{
        margin-top: 40px;
        padding-top: 15px;
        border-top: 1px solid #ddd;
        font-size: 0.8em;
        color: #95a5a6;
        text-align: center;
    }}
</style>
</head>
<body>

<h1>Pipeline Performance Report</h1>
<p>Run-time assessment for mapping (BAM generation) and variant calling (VCF generation) steps.</p>

<div class="summary-box">
    <div class="summary-item">
        <div class="label">Total Wall Time</div>
        <div class="value">{fmt_td(pipeline_wall)}</div>
    </div>
    <div class="summary-item">
        <div class="label">Samples</div>
        <div class="value">{len(samples_seen)}</div>
    </div>
    <div class="summary-item">
        <div class="label">Avg BAM Time</div>
        <div class="value">{fmt_td(avg_total)}</div>
    </div>
    <div class="summary-item">
        <div class="label">VCF Generation</div>
        <div class="value">{fmt_td(vcf_total_wall)}</div>
    </div>
</div>

<h2>BAM Generation (per sample)</h2>
<p class="section-note">
    Pipeline: bwa-mem2 mem | samtools sort &rarr; Picard MarkDuplicates &rarr; samtools index.
    Settings: 7 threads per alignment, 3 concurrent, samtools sort -m 10G.
</p>

<table>
<thead>
<tr>
    <th>Sample</th>
    <th>FASTQ Size</th>
    <th>Alignment</th>
    <th>Mark Dups</th>
    <th>Indexing</th>
    <th>Total</th>
    <th>BAM Size</th>
    <th>Alignment Time</th>
</tr>
</thead>
<tbody>
{bam_table_rows}
</tbody>
</table>

<h2>VCF Generation</h2>
<p class="section-note">
    Pipeline: bcftools mpileup | bcftools call (per region) &rarr; bcftools concat &rarr; bcftools filter + norm.
</p>

<table>
<thead>
<tr>
    <th>Step</th>
    <th>Details</th>
    <th>Wall Time</th>
</tr>
</thead>
<tbody>
<tr>
    <td>Variant calling</td>
    <td>{vc_n_regions} regions, avg {fmt_td(vc_avg_per_region)} per region</td>
    <td>{fmt_td(vc_wall)}</td>
</tr>
<tr>
    <td>Concatenate VCFs</td>
    <td>Merge {vc_n_regions} region VCFs &rarr; raw_variants.vcf.gz ({file_size_human(raw_vcf_size)})</td>
    <td>{fmt_td(concat_dur)}</td>
</tr>
<tr>
    <td>Filter variants</td>
    <td>QUAL&gt;20, DP&gt;5, biallelic SNPs &rarr; filtered_variants.vcf.gz ({file_size_human(filt_vcf_size)})</td>
    <td>{fmt_td(filter_dur)}</td>
</tr>
<tr class="avg-row">
    <td><strong>Total VCF generation</strong></td>
    <td></td>
    <td><strong>{fmt_td(vcf_total_wall)}</strong></td>
</tr>
</tbody>
</table>

<div class="footer">
    Generated from Snakemake log: {os.path.basename(log_path)}<br>
    Pipeline run: {pipeline_start.strftime('%Y-%m-%d %H:%M')} &mdash; {pipeline_end.strftime('%Y-%m-%d %H:%M')}
</div>

</body>
</html>
"""

os.makedirs(os.path.dirname(output_html), exist_ok=True)
with open(output_html, "w") as fh:
    fh.write(html)

print(f"Performance report written to {output_html}")
print(f"Pipeline wall time: {fmt_td(pipeline_wall)}")
print(f"Average BAM generation: {fmt_td(avg_total)} per sample")
print(f"VCF generation: {fmt_td(vcf_total_wall)}")
