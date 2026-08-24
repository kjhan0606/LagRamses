#!/bin/bash
#SBATCH --job-name=sidm_zoom_cdm
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --mem=400G
#SBATCH --time=14-00:00:00
#SBATCH --output=zoom_cdm_%j.log
#SBATCH --error=zoom_cdm_%j.err

set -euo pipefail

run_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_cdm
ic_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_l11/ic_zoom_l11
cd "${run_dir}"

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

for level in 008 009 010 011; do
    if [ ! -d "${ic_dir}/level_${level}" ]; then
        printf 'Missing zoom initial condition level_%s\n' "${level}" >&2
        exit 2
    fi
done

if [ ! -x ./ramses_zoom3d ]; then
    printf 'Missing zoom RAMSES executable\n' >&2
    exit 3
fi

printf '' > jobcontrol.txt
printf 'CDM zoom job %s started on %s at %s\n' \
       "${SLURM_JOB_ID}" "$(hostname)" "$(date)"
printf 'MPI ranks=%s OpenMP threads=%s\n' \
       "${SLURM_NTASKS}" "${OMP_NUM_THREADS}"

srun --mpi=pmi2 ./ramses_zoom3d zoom_cdm.nml
rc=$?

printf 'CDM zoom job %s finished at %s with rc=%s\n' \
       "${SLURM_JOB_ID}" "$(date)" "${rc}"
exit "${rc}"
