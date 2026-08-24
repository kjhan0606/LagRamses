#!/bin/bash
#SBATCH --job-name=sidm_cdm_z63_test
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --mem=450G
#SBATCH --time=1-00:00:00
#SBATCH --output=zoom_cdm_z63_test_%j.log
#SBATCH --error=zoom_cdm_z63_test_%j.err

set -euo pipefail

run_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_cdm_z63_test
ic_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_l11_z63/ic_zoom_l11_z63
cd "${run_dir}"

for level in 008 009 010 011; do
    if [ ! -d "${ic_dir}/level_${level}" ]; then
        printf 'Missing zstart=63 IC level_%s\n' "${level}" >&2
        exit 2
    fi
done
if [ ! -x ./ramses_zoom3d ]; then
    printf 'Missing RAMSES executable\n' >&2
    exit 3
fi
if find . -maxdepth 1 -type d -name 'output_*' -print -quit | grep -q .; then
    printf 'Refusing to overwrite an existing test output\n' >&2
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
printf 'zstart=63 short CDM test %s started on %s at %s\n' \
    "${SLURM_JOB_ID}" "$(hostname)" "$(date --iso-8601=seconds)"
printf 'Targets: z=31 and z=15; production zstart=100 run is untouched\n'

srun --mpi=pmi2 ./ramses_zoom3d zoom_cdm_z63_test.nml
