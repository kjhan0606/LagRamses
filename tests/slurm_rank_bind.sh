#!/bin/bash
# [RESIZABLE] Bind one Slurm local rank to its deterministic allocation CPU.
# PHASE4_CPU_LIST must be a canonical CSV list produced by expand_cpuset.sh.

set -euo pipefail

: "${PHASE4_CPU_LIST:?PHASE4_CPU_LIST is required}"
: "${SLURM_LOCALID:?SLURM_LOCALID is required}"
[ "$#" -gt 0 ] || { echo "rank command is required" >&2; exit 64; }
[[ $SLURM_LOCALID =~ ^[0-9]+$ ]] || exit 96

IFS=, read -r -a cpu_ids <<< "$PHASE4_CPU_LIST"
if [ "$SLURM_LOCALID" -ge "${#cpu_ids[@]}" ]; then
  exit 96
fi
cpu=${cpu_ids[$SLURM_LOCALID]}
[[ $cpu =~ ^[0-9]+$ ]] || exit 97

exec /usr/bin/taskset --cpu-list "$cpu" "$@"
