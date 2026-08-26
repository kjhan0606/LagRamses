#!/bin/bash
#SBATCH --job-name=sidm_z4_l17
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-1
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/z4_force_l17/slurm_%A_%a.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/z4_force_l17/slurm_%A_%a.err

set -euo pipefail

source_root=/home/kjhan/BACKUP/lagRamses-SIDM
validation_root=/gpfs/kjhan/Hydro/Sidm/validation/z4_force_l17

case "${SLURM_ARRAY_TASK_ID:?}" in
  0)
    label=cdm
    input=z4_force_l17_cdm.nml
    executable=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_cdm/ramses_zoom3d
    ;;
  1)
    label=sidm3
    input=z4_force_l17_sidm3.nml
    executable=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_sidm3/ramses_zoom3d
    ;;
  *)
    printf 'unsupported array task: %s\n' "${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
    ;;
esac

case_root=${validation_root}/${label}_${SLURM_ARRAY_JOB_ID}
if [[ -e "${case_root}" ]]; then
  printf 'refusing to reuse existing case directory: %s\n' "${case_root}" >&2
  exit 3
fi
mkdir -p "${case_root}/snapshot"

if [[ ! -x "${executable}" ]]; then
  printf 'missing executable: %s\n' "${executable}" >&2
  exit 4
fi

git -C "${source_root}" status --short >"${case_root}/snapshot/git_status.txt"
git -C "${source_root}" diff --binary >"${case_root}/snapshot/source.diff"
cp -p "${source_root}/tests/sidm/${input}" "${case_root}/run.nml"
cp -p "${source_root}/tests/sidm/run_z4_force_l17.sbatch" \
  "${case_root}/snapshot/"
cp -p "${case_root}/run.nml" "${case_root}/snapshot/"
cp -p "${executable}" "${case_root}/ramses_z4_l17"
sha256sum "${case_root}/ramses_z4_l17" >"${case_root}/snapshot/binary.sha256"
sha256sum "${executable}" >"${case_root}/snapshot/source_binary.sha256"

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
export I_MPI_FABRICS=shm

: >"${case_root}/jobcontrol.txt"
(
  cd "${case_root}"
  mpirun -np "${SLURM_NTASKS}" ./ramses_z4_l17 run.nml >run.log 2>&1
)
