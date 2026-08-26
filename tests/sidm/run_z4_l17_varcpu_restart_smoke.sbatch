#!/bin/bash
# Validate the production binary's Hilbert restart from the 12-MPI L17
# SIDM3 force-convergence snapshot to the 24-MPI continuation layout.
# Example:
#   sbatch --export=ALL,SEED_JOB_ID=323771 run_z4_l17_varcpu_restart_smoke.sbatch
#SBATCH --job-name=sidm_l17_varcpu_smoke
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/slurm_l17_varcpu_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/slurm_l17_varcpu_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
seed_job_id=${SEED_JOB_ID:?SEED_JOB_ID must name the 12-MPI L17 SIDM3 run}
validation_root=/gpfs/kjhan/Hydro/Sidm/validation
source_root=${validation_root}/z4_force_l17/sidm3_${seed_job_id}
source_output=${source_root}/output_00001
source_binary=${source_root}/ramses_z4_l17
case_root=${validation_root}/z4_varcpu_restart/sidm3_${SLURM_JOB_ID}
input=${code}/tests/sidm/z4_l17_varcpu_restart_sidm3.nml

test -d "${source_output}"
test -f "${source_output}/info_00001.txt"
test -x "${source_binary}"
if [[ -e "${case_root}" ]]; then
  printf 'refusing to reuse case directory: %s\n' "${case_root}" >&2
  exit 3
fi
mkdir -p "${case_root}/snapshot"

git -C "${code}" status --short >"${case_root}/snapshot/git_status.txt"
git -C "${code}" diff --binary >"${case_root}/snapshot/source.diff"
cp -p "${input}" "${case_root}/run.nml"
cp -p "${code}/tests/sidm/run_z4_l17_varcpu_restart_smoke.sbatch" \
  "${case_root}/snapshot/"
cp -p "${case_root}/run.nml" "${case_root}/snapshot/"
cp -p "${source_binary}" "${case_root}/ramses_varcpu_restart"
sha256sum "${source_binary}" >"${case_root}/snapshot/source_binary.sha256"
sha256sum "${case_root}/ramses_varcpu_restart" >"${case_root}/snapshot/binary.sha256"
ln -s "${source_output}" "${case_root}/output_00001"
printf '0 -1\n' >"${case_root}/jobcontrol.txt"

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

(
  cd "${case_root}"
  mpirun -np "${SLURM_NTASKS}" ./ramses_varcpu_restart run.nml >run.log 2>&1
)

test -d "${case_root}/output_00002"
grep -Fq 'Variable-ncpu restart: ncpu_old=   12 ncpu_new=   24' "${case_root}/run.log"
grep -Fq 'Binary variable-ncpu AMR restore done.' "${case_root}/run.log"
grep -Fq 'Starting time integration' "${case_root}/run.log"
grep -Fq 'Job control: stop requested' "${case_root}/run.log"
grep -Fq 'Run completed' "${case_root}/run.log"

source_aexp=$(awk '/^aexp[[:space:]]*=/ {print $3; exit}' \
  "${source_output}/info_00001.txt")
next_aexp=$(awk '/^aexp[[:space:]]*=/ {print $3; exit}' \
  "${case_root}/output_00002/info_00002.txt")
awk -v source_aexp="${source_aexp}" -v next_aexp="${next_aexp}" \
  'BEGIN {exit !(next_aexp > source_aexp)}'

{
  printf 'source_output=output_00001\n'
  printf 'source_ncpu=12\n'
  printf 'restart_ncpu=24\n'
  printf 'source_aexp=%s\n' "${source_aexp}"
  printf 'next_aexp=%s\n' "${next_aexp}"
  printf '%s\n' 'L17_VARCPU_RESTART_SMOKE_RESULT=passed'
} | tee "${case_root}/restart_smoke_summary.txt"
