#!/usr/bin/env bash
set -euo pipefail

run_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_sidm03
ic_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_l11/ic_zoom_l11
exe=${run_dir}/ramses_zoom3d
nml=zoom_sidm03.nml

cd "${run_dir}"

for level in 008 009 010 011; do
    if [ ! -d "${ic_dir}/level_${level}" ]; then
        printf 'Missing zoom initial condition level_%s\n' "${level}" >&2
        exit 2
    fi
done

if [ ! -x "${exe}" ]; then
    printf 'Missing executable: %s\n' "${exe}" >&2
    exit 3
fi

if pgrep -u "${USER}" -af "ramses_zoom3d ${nml}" >/dev/null 2>&1; then
    printf 'SIDM03 zoom appears to be already running; refusing duplicate start.\n' >&2
    pgrep -u "${USER}" -af "ramses_zoom3d ${nml}" >&2 || true
    exit 4
fi

mpi_ranks=24
omp_threads=2

# LagEunha socket0 / NUMA0:
#   physical cores 0-31, SMT siblings 64-95.
# Use only cores 0-23 plus siblings 64-87, leaving 24-31 and 88-95 free.
cpu_list=
pin_list=
for core in $(seq 0 23); do
    sibling=$((core + 64))
    if [ -z "${cpu_list}" ]; then
        cpu_list="${core},${sibling}"
        pin_list="${core},${sibling}"
    else
        cpu_list="${cpu_list},${core},${sibling}"
        pin_list="${pin_list},${core},${sibling}"
    fi
done

export OMP_NUM_THREADS=${omp_threads}
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close
export OMP_PLACES=threads
unset KMP_AFFINITY

export I_MPI_FABRICS=shm
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact
export I_MPI_PIN_PROCESSOR_LIST="${pin_list}"
export I_MPI_PIN_RESPECT_CPUSET=1
export I_MPI_DEBUG=${I_MPI_DEBUG:-4}

export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

ulimit -s unlimited || true
printf '' > jobcontrol.txt

printf 'SIDM03 manual zoom run started on %s at %s\n' "$(hostname)" "$(date --iso-8601=seconds)"
printf 'run_dir=%s\n' "${run_dir}"
printf 'model=constant isotropic SIDM, sigma/m=0.3 cm^2/g\n'
printf 'mpi_ranks=%s omp_threads=%s logical_cpus=%s\n' "${mpi_ranks}" "${omp_threads}" "${cpu_list}"
printf 'memory_policy=interleave:all\n'
printf 'executable=%s\n' "${exe}"

set +e
numactl --physcpubind="${cpu_list}" --interleave=all \
    mpirun -np "${mpi_ranks}" "${exe}" "${nml}"
rc=$?
set -e

printf 'SIDM03 manual zoom run finished on %s at %s rc=%s\n' \
    "$(hostname)" "$(date --iso-8601=seconds)" "${rc}"
printf '%s\n' "${rc}" > zoom_sidm03_lageunha_manual.rc
exit "${rc}"
