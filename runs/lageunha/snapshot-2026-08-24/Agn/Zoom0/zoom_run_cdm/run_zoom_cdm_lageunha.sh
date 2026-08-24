#!/bin/bash
set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
run_dir=${root}/zoom_run_cdm
ic_dir=${root}/zoom_l11/ic_zoom_l11
lock_dir=${run_dir}/.manual_lageunha.lock
pid_file=${run_dir}/zoom_cdm_lageunha.pid
log=${run_dir}/zoom_cdm_lageunha.log

if [ "$(hostname)" != "LagEunha" ]; then
    printf 'This manual fallback must run on LagEunha, not %s\n' \
           "$(hostname)" >&2
    exit 2
fi

if ! mkdir "${lock_dir}" 2>/dev/null; then
    printf 'Manual CDM zoom lock already exists: %s\n' "${lock_dir}" >&2
    exit 3
fi

cleanup() {
    rmdir "${lock_dir}" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

cd "${run_dir}"
printf '%s\n' "$$" > "${pid_file}"

for level in 008 009 010 011; do
    if [ ! -d "${ic_dir}/level_${level}" ]; then
        printf 'Missing zoom initial condition level_%s\n' "${level}" >&2
        exit 4
    fi
done

if [ ! -x ./ramses_zoom3d ]; then
    printf 'Missing zoom RAMSES executable\n' >&2
    exit 5
fi

if find . -maxdepth 1 -type d -name 'output_*' -print -quit |
   grep -q .; then
    printf 'Existing CDM zoom output found; refusing a duplicate start\n' >&2
    exit 6
fi

export OMP_NUM_THREADS=2
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact
export I_MPI_DEBUG=0
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

printf '' > jobcontrol.txt
printf '%s Manual CDM zoom pid=%s starting on %s\n' \
       "$(date --iso-8601=seconds)" "$$" "$(hostname)" >> "${log}"
printf '%s MPI ranks=32 OpenMP threads=2\n' \
       "$(date --iso-8601=seconds)" >> "${log}"

set +e
mpirun -np 32 ./ramses_zoom3d zoom_cdm.nml >> "${log}" 2>&1
rc=$?
set -e

printf '%s Manual CDM zoom finished rc=%s\n' \
       "$(date --iso-8601=seconds)" "${rc}" >> "${log}"
exit "${rc}"
