#!/bin/bash
# [RESIZABLE] Give each Slurm task a complete, non-interleaved runtime log.
set -euo pipefail

prefix=${1:?log prefix is required}
shift
[ "$#" -gt 0 ] || { echo 'rank command is required' >&2; exit 64; }
: "${SLURM_PROCID:?SLURM_PROCID is required}"
[[ $SLURM_PROCID =~ ^[0-9]+$ ]] || exit 65
exec "$@" >"${prefix}_${SLURM_PROCID}.log" 2>&1
