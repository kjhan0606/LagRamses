#!/bin/bash
set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
lock_dir=${root}/.manual_prepare_lageunha.lock
log=${root}/prepare_zoom_lageunha.log
pid_file=${root}/prepare_zoom_lageunha.pid

if [ "$(hostname)" != "LagEunha" ]; then
    printf 'This manual fallback must run on LagEunha, not %s\n' \
           "$(hostname)" >&2
    exit 2
fi

if ! mkdir "${lock_dir}" 2>/dev/null; then
    printf 'Manual preparation lock already exists: %s\n' "${lock_dir}" >&2
    exit 3
fi

cleanup() {
    rmdir "${lock_dir}" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

printf '%s\n' "$$" > "${pid_file}"
export SLURM_CPUS_PER_TASK=64

printf '%s Manual zoom preparation starting on %s pid=%s\n' \
       "$(date --iso-8601=seconds)" "$(hostname)" "$$" >> "${log}"

set +e
"${root}/parent_run/prepare_zoom_l11_grammar.sh" >> "${log}" 2>&1
rc=$?
set -e

printf '%s Manual zoom preparation finished rc=%s\n' \
       "$(date --iso-8601=seconds)" "${rc}" >> "${log}"
exit "${rc}"
