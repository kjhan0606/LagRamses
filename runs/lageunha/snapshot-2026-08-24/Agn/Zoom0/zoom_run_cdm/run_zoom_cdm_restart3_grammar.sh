#!/bin/bash
#SBATCH --job-name=sidm_zoom_cdm_r3
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --mem=450G
#SBATCH --time=14-00:00:00
#SBATCH --output=zoom_cdm_restart3_%j.log
#SBATCH --error=zoom_cdm_restart3_%j.err

set -euo pipefail

run_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_cdm
restart_dir=${run_dir}/output_00003
cd "${run_dir}"

for prefix in amr grav part; do
    count=$(find "${restart_dir}" -maxdepth 1 -type f \
        -name "${prefix}_00003.out*" | wc -l)
    if [ "${count}" -ne 32 ]; then
        printf 'Incomplete restart: expected 32 %s files, found %s\n' \
            "${prefix}" "${count}" >&2
        exit 2
    fi
done

for metadata in info_00003.txt header_00003.txt; do
    if [ ! -s "${restart_dir}/${metadata}" ]; then
        printf 'Incomplete restart: missing %s\n' "${metadata}" >&2
        exit 3
    fi
done

if [ ! -x ./ramses_zoom3d ]; then
    printf 'Missing zoom RAMSES executable\n' >&2
    exit 4
fi

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

printf '' > jobcontrol.txt
printf 'CDM zoom restart job %s started on %s at %s\n' \
       "${SLURM_JOB_ID}" "$(hostname)" "$(date)"
printf 'Restart=3, ngridtot=240000000, nparttot=200000000\n'
printf 'MPI ranks=%s OpenMP threads=%s\n' \
       "${SLURM_NTASKS}" "${OMP_NUM_THREADS}"

srun --mpi=pmi2 ./ramses_zoom3d zoom_cdm_restart3.nml
rc=$?

printf 'CDM zoom restart job %s finished at %s with rc=%s\n' \
       "${SLURM_JOB_ID}" "$(date)" "${rc}"
exit "${rc}"
