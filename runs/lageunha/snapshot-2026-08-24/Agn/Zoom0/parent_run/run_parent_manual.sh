#!/bin/bash
set -euo pipefail

run_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/parent_run
cd "${run_dir}"

exec 9>parent_manual.lock
if ! flock -n 9; then
    printf 'Another parent manual run holds the lock\n' >&2
    exit 3
fi

export OMP_NUM_THREADS=2
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export I_MPI_FABRICS=shm
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=omp
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

if [ ! -d ../parent/ic_parent_lv8/level_008 ]; then
    printf 'Missing parent initial conditions\n' >&2
    exit 2
fi

printf '' > jobcontrol.txt
printf 'Manual parent started on %s at %s\n' "$(hostname)" "$(date)"
printf 'MPI ranks=16 OpenMP threads=%s\n' "${OMP_NUM_THREADS}"

/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun \
    -np 16 ./ramses_final3d parent.nml
rc=$?

printf 'Manual parent finished on %s at %s with rc=%s\n' \
       "$(hostname)" "$(date)" "${rc}"
exit "${rc}"
