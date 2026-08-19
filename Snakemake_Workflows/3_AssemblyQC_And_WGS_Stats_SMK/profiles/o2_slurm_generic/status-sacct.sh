#!/usr/bin/env bash
# Job-status script for the Snakemake cluster-generic executor plugin.
# Takes one SLURM jobid as $1, prints one of: running / success / failed.
# Simplified re-implementation of DeepVarMtb's workflow/profiles/o2_slurm/slurm-status.py
# for the new (Snakemake >=8) cluster-generic executor interface, which expects a
# single positional jobid argument and one of these three literal strings on stdout.
set -euo pipefail

jobid="$1"

state=$(sacct -j "$jobid" -P -b -n 2>/dev/null | head -n1 | cut -d'|' -f2 || true)

if [ -z "$state" ]; then
    echo "running"
    exit 0
fi

case "$state" in
    COMPLETED)
        echo "success"
        ;;
    BOOT_FAIL|CANCELLED*|DEADLINE|FAILED|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|TIMEOUT)
        echo "failed"
        ;;
    *)
        echo "running"
        ;;
esac
