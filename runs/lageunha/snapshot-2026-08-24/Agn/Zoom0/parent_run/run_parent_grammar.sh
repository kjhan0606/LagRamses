#!/bin/bash
#SBATCH --job-name=sidm_parent
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --mem=256G
#SBATCH --time=14-00:00:00
#SBATCH --output=grammar_%j.log
#SBATCH --error=grammar_%j.err

set -euo pipefail

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

export UCX_RC_TIMEOUT=30s
export UCX_RC_RETRY_COUNT=7
export UCX_LOG_LEVEL=warn
export UCX_IB_REG_METHODS=direct
export UCX_RNDV_THRESH=65536

cd /gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/parent_run

if [ ! -d ../parent/ic_parent_lv8/level_008 ]; then
    printf 'Missing parent initial conditions\n' >&2
    exit 2
fi

printf '' > jobcontrol.txt
printf 'CPU parent job %s started on %s at %s\n' \
       "${SLURM_JOB_ID}" "$(hostname)" "$(date)"
printf 'MPI ranks=%s OpenMP threads=%s\n' \
       "${SLURM_NTASKS}" "${OMP_NUM_THREADS}"

srun --mpi=pmi2 ./ramses_parent3d parent.nml
rc=$?

printf 'CPU parent job %s finished at %s with rc=%s\n' \
       "${SLURM_JOB_ID}" "$(date)" "${rc}"
exit "${rc}"
