#!/usr/bin/env python3
"""Pipeline monitor — generates an auto-refreshing HTML dashboard every hour."""
import subprocess, glob, time, os, re
from datetime import datetime

OUTFILE = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/reports/pipeline_dashboard.html"
PIPELINE_LOG = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/pipeline_run2.log"
ALIGNED_DIR = "/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/aligned"
TOTAL_SAMPLES = 94
DISK_TOTAL = 926

def get_status():
    dedup = len(glob.glob(os.path.join(ALIGNED_DIR, "*.dedup.bam")))

    jobs_done, jobs_total = 0, 368
    try:
        with open(PIPELINE_LOG) as f:
            for line in f:
                m = re.search(r'(\d+) of (\d+) steps', line)
                if m:
                    jobs_done, jobs_total = int(m.group(1)), int(m.group(2))
    except:
        pass

    try:
        out = subprocess.check_output(["df", "-g", "/Users/stockholmbutterflylab"], text=True)
        disk_free = int(out.strip().split("\n")[-1].split()[3])
    except:
        disk_free = 0
    disk_used = DISK_TOTAL - disk_free

    ps = subprocess.run(["pgrep", "-f", "snakemake"], capture_output=True)
    running = ps.returncode == 0

    return dedup, jobs_done, jobs_total, disk_free, disk_used, running

def write_html(dedup, jobs_done, jobs_total, disk_free, disk_used, running):
    bam_pct = dedup * 100 // TOTAL_SAMPLES
    jobs_pct = jobs_done * 100 // jobs_total if jobs_total > 0 else 0
    disk_pct = disk_used * 100 // DISK_TOTAL

    if running:
        status_text, status_color = "Running", "#22c55e"
    elif jobs_done >= jobs_total:
        status_text, status_color = "Complete", "#3b82f6"
    else:
        status_text, status_color = "Stopped", "#f59e0b"

    disk_bar_color = "#ef4444, #f87171" if disk_pct > 80 else "#22c55e, #4ade80"

    now = datetime.now()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="refresh" content="3600">
    <title>Pipeline Monitor</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', sans-serif;
            background: #0f172a;
            color: #e2e8f0;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
            padding: 2rem;
        }}
        .dashboard {{ width: 100%; max-width: 640px; }}
        .header {{ text-align: center; margin-bottom: 2.5rem; }}
        .header h1 {{ font-size: 1.5rem; font-weight: 600; color: #f8fafc; letter-spacing: -0.025em; }}
        .header .subtitle {{ font-size: 0.8rem; color: #f1f5f9; margin-top: 0.4rem; }}
        .status-badge {{
            display: inline-flex; align-items: center; gap: 0.4rem;
            padding: 0.3rem 0.8rem; border-radius: 999px;
            font-size: 0.75rem; font-weight: 600; margin-top: 0.8rem;
            background: {status_color}20; color: {status_color};
            border: 1px solid {status_color}40;
        }}
        .status-dot {{ width: 6px; height: 6px; border-radius: 50%; background: {status_color}; }}
        .card {{
            background: #1e293b; border-radius: 12px; padding: 1.5rem;
            margin-bottom: 1rem; border: 1px solid #334155;
        }}
        .metric-row {{ display: flex; justify-content: space-between; align-items: baseline; margin-bottom: 0.6rem; }}
        .metric-label {{ font-size: 0.85rem; color: #f1f5f9; font-weight: 500; }}
        .metric-value {{ font-size: 0.85rem; color: #ffffff; font-weight: 600; font-variant-numeric: tabular-nums; }}
        .bar-container {{
            width: 100%; height: 28px; background: #0f172a;
            border-radius: 8px; overflow: hidden; position: relative;
        }}
        .bar-fill {{
            height: 100%; border-radius: 8px;
            transition: width 1s ease; position: relative;
            display: flex; align-items: center; justify-content: flex-end; padding-right: 10px;
        }}
        .bar-label {{
            font-size: 0.7rem; font-weight: 700; color: rgba(255,255,255,0.9);
        }}
        .bar-bam {{ background: linear-gradient(90deg, #3b82f6, #60a5fa); }}
        .bar-jobs {{ background: linear-gradient(90deg, #8b5cf6, #a78bfa); }}
        .bar-disk {{ background: linear-gradient(90deg, {disk_bar_color}); }}
        .footer {{ text-align: center; margin-top: 1.5rem; font-size: 0.7rem; color: #e2e8f0; }}
    </style>
</head>
<body>
    <div class="dashboard">
        <div class="header">
            <h1>GWAS Pipeline</h1>
            <div class="subtitle">Pararge aegeria &mdash; 94 samples</div>
            <div class="status-badge">
                <span class="status-dot"></span>
                {status_text}
            </div>
        </div>

        <div class="card">
            <div class="metric-row">
                <span class="metric-label">Samples mapped</span>
                <span class="metric-value">{dedup} / {TOTAL_SAMPLES}</span>
            </div>
            <div class="bar-container">
                <div class="bar-fill bar-bam" style="width: {max(bam_pct, 3)}%">
                    <span class="bar-label">{bam_pct}%</span>
                </div>
            </div>
        </div>

        <div class="card">
            <div class="metric-row">
                <span class="metric-label">Pipeline jobs</span>
                <span class="metric-value">{jobs_done} / {jobs_total}</span>
            </div>
            <div class="bar-container">
                <div class="bar-fill bar-jobs" style="width: {max(jobs_pct, 3)}%">
                    <span class="bar-label">{jobs_pct}%</span>
                </div>
            </div>
        </div>

        <div class="card">
            <div class="metric-row">
                <span class="metric-label">Disk used</span>
                <span class="metric-value">{disk_used}G / {DISK_TOTAL}G ({disk_free}G free)</span>
            </div>
            <div class="bar-container">
                <div class="bar-fill bar-disk" style="width: {max(disk_pct, 3)}%">
                    <span class="bar-label">{disk_pct}%</span>
                </div>
            </div>
        </div>

        <div class="footer">
            Last updated: {now.strftime('%Y-%m-%d %H:%M:%S')} &mdash; refreshes every hour
        </div>
    </div>
</body>
</html>"""

    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
    with open(OUTFILE, "w") as f:
        f.write(html)

def main():
    while True:
        dedup, jobs_done, jobs_total, disk_free, disk_used, running = get_status()
        write_html(dedup, jobs_done, jobs_total, disk_free, disk_used, running)
        print(f"[{datetime.now().strftime('%H:%M:%S')}] BAMs: {dedup}/94 | Jobs: {jobs_done}/{jobs_total} | Disk free: {disk_free}G | Running: {running}")
        if not running:
            print("Pipeline no longer running. Final dashboard written. Exiting.")
            break
        time.sleep(3600)

if __name__ == "__main__":
    main()
