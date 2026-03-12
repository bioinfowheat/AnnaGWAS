#!/bin/bash
LOGFILE="/Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/pipeline_status.log"
while true; do
    echo "========== $(date '+%Y-%m-%d %H:%M:%S') ==========" >> "$LOGFILE"
    DEDUP=$(ls /Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/results/aligned/*.dedup.bam 2>/dev/null | wc -l | tr -d ' ')
    PROGRESS=$(grep "steps (" /Users/stockholmbutterflylab/sbl_claudecode/Annas_GWAS/pipeline_run2.log 2>/dev/null | tail -1)
    DISK=$(df -h /Users/stockholmbutterflylab | tail -1 | awk '{print $4}')
    RUNNING=$(ps aux | grep snakemake | grep -v grep | wc -l | tr -d ' ')
    echo "Dedup BAMs: ${DEDUP}/94" >> "$LOGFILE"
    echo "Progress: ${PROGRESS}" >> "$LOGFILE"
    echo "Disk free: ${DISK}" >> "$LOGFILE"
    echo "Snakemake running: $( [ "$RUNNING" -gt 0 ] && echo 'yes' || echo 'NO - may have finished or failed')" >> "$LOGFILE"
    echo "" >> "$LOGFILE"
    # Exit if snakemake is no longer running
    if [ "$RUNNING" -eq 0 ]; then
        echo "Pipeline no longer running. Monitor exiting." >> "$LOGFILE"
        exit 0
    fi
    sleep 3600
done
