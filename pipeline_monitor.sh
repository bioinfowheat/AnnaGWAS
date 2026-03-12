#!/bin/bash
# Pipeline monitor - generates an auto-refreshing HTML dashboard
OUTFILE="/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/reports/pipeline_dashboard.html"
PIPELINE_LOG="/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/pipeline_run2.log"
ALIGNED_DIR="/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/aligned"
TOTAL_SAMPLES=94
DISK_TOTAL=926

while true; do
    DEDUP=$(ls "$ALIGNED_DIR"/*.dedup.bam 2>/dev/null | wc -l | tr -d ' ')
    PROGRESS_LINE=$(grep "steps (" "$PIPELINE_LOG" 2>/dev/null | tail -1)
    JOBS_DONE=$(echo "$PROGRESS_LINE" | grep -oE '^[0-9]+' || echo "0")
    JOBS_TOTAL=$(echo "$PROGRESS_LINE" | grep -oE 'of [0-9]+' | grep -oE '[0-9]+' || echo "368")
    DISK_FREE=$(df -g /Users/stockholmbutterflylab | tail -1 | awk '{print $4}')
    DISK_USED=$((DISK_TOTAL - DISK_FREE))
    TIMESTAMP=$(date '+%H:%M:%S')
    DATE=$(date '+%Y-%m-%d')
    SNAKEMAKE_RUNNING=$(ps aux | grep "[s]nakemake" | wc -l | tr -d ' ')

    if [ "$SNAKEMAKE_RUNNING" -gt 0 ]; then
        STATUS_TEXT="Running"
        STATUS_COLOR="#22c55e"
        STATUS_EMOJI="&#9654;"
    else
        # Check if it completed or failed
        if grep -q "Error\|error\|failed" "$PIPELINE_LOG" 2>/dev/null | tail -5 | grep -q -i "error"; then
            STATUS_TEXT="Failed"
            STATUS_COLOR="#ef4444"
            STATUS_EMOJI="&#10060;"
        else
            STATUS_TEXT="Stopped"
            STATUS_COLOR="#f59e0b"
            STATUS_EMOJI="&#9209;"
        fi
    fi

    BAM_PCT=$((DEDUP * 100 / TOTAL_SAMPLES))
    JOBS_PCT=$((JOBS_DONE * 100 / JOBS_TOTAL))
    DISK_PCT=$((DISK_USED * 100 / DISK_TOTAL))

    cat > "$OUTFILE" << HTMLEOF
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="refresh" content="3600">
    <title>Pipeline Monitor</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', sans-serif;
            background: #0f172a;
            color: #e2e8f0;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
            padding: 2rem;
        }
        .dashboard {
            width: 100%;
            max-width: 640px;
        }
        .header {
            text-align: center;
            margin-bottom: 2.5rem;
        }
        .header h1 {
            font-size: 1.5rem;
            font-weight: 600;
            color: #f8fafc;
            letter-spacing: -0.025em;
        }
        .header .subtitle {
            font-size: 0.8rem;
            color: #64748b;
            margin-top: 0.4rem;
        }
        .status-badge {
            display: inline-flex;
            align-items: center;
            gap: 0.4rem;
            padding: 0.3rem 0.8rem;
            border-radius: 999px;
            font-size: 0.75rem;
            font-weight: 600;
            margin-top: 0.8rem;
            background: ${STATUS_COLOR}20;
            color: ${STATUS_COLOR};
            border: 1px solid ${STATUS_COLOR}40;
        }
        .status-dot {
            width: 6px;
            height: 6px;
            border-radius: 50%;
            background: ${STATUS_COLOR};
        }
        .card {
            background: #1e293b;
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 1rem;
            border: 1px solid #334155;
        }
        .metric-row {
            display: flex;
            justify-content: space-between;
            align-items: baseline;
            margin-bottom: 0.6rem;
        }
        .metric-label {
            font-size: 0.85rem;
            color: #94a3b8;
            font-weight: 500;
        }
        .metric-value {
            font-size: 0.85rem;
            color: #cbd5e1;
            font-weight: 600;
            font-variant-numeric: tabular-nums;
        }
        .bar-container {
            width: 100%;
            height: 28px;
            background: #0f172a;
            border-radius: 8px;
            overflow: hidden;
            position: relative;
        }
        .bar-fill {
            height: 100%;
            border-radius: 8px;
            transition: width 1s ease;
            position: relative;
        }
        .bar-fill::after {
            content: attr(data-label);
            position: absolute;
            right: 10px;
            top: 50%;
            transform: translateY(-50%);
            font-size: 0.7rem;
            font-weight: 700;
            color: rgba(255,255,255,0.9);
        }
        .bar-bam { background: linear-gradient(90deg, #3b82f6, #60a5fa); }
        .bar-jobs { background: linear-gradient(90deg, #8b5cf6, #a78bfa); }
        .bar-disk {
            background: linear-gradient(90deg,
                ${DISK_PCT > 80 ? "#ef4444, #f87171" : "#22c55e, #4ade80"});
        }
        .footer {
            text-align: center;
            margin-top: 1.5rem;
            font-size: 0.7rem;
            color: #475569;
        }
    </style>
</head>
<body>
    <div class="dashboard">
        <div class="header">
            <h1>GWAS Pipeline</h1>
            <div class="subtitle">Pararge aegeria &mdash; 94 samples</div>
            <div class="status-badge">
                <span class="status-dot"></span>
                ${STATUS_TEXT}
            </div>
        </div>

        <div class="card">
            <div class="metric-row">
                <span class="metric-label">Samples mapped</span>
                <span class="metric-value">${DEDUP} / ${TOTAL_SAMPLES}</span>
            </div>
            <div class="bar-container">
                <div class="bar-fill bar-bam" style="width: ${BAM_PCT}%" data-label="${BAM_PCT}%"></div>
            </div>
        </div>

        <div class="card">
            <div class="metric-row">
                <span class="metric-label">Pipeline jobs</span>
                <span class="metric-value">${JOBS_DONE} / ${JOBS_TOTAL}</span>
            </div>
            <div class="bar-container">
                <div class="bar-fill bar-jobs" style="width: ${JOBS_PCT}%" data-label="${JOBS_PCT}%"></div>
            </div>
        </div>

        <div class="card">
            <div class="metric-row">
                <span class="metric-label">Disk used</span>
                <span class="metric-value">${DISK_USED}G / ${DISK_TOTAL}G (${DISK_FREE}G free)</span>
            </div>
            <div class="bar-container">
                <div class="bar-fill bar-disk" style="width: ${DISK_PCT}%" data-label="${DISK_PCT}%"></div>
            </div>
        </div>

        <div class="footer">
            Last updated: ${DATE} ${TIMESTAMP} &mdash; refreshes every hour
        </div>
    </div>
</body>
</html>
HTMLEOF

    # If snakemake is no longer running, write one final update and exit
    if [ "$SNAKEMAKE_RUNNING" -eq 0 ]; then
        echo "$(date): Pipeline no longer running. Monitor exiting." >> /Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/pipeline_status.log
        exit 0
    fi

    sleep 3600
done
