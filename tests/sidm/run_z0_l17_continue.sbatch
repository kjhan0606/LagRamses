#!/bin/bash
# Submit only after the paired L17 z≈4 force-convergence run and its analysis
# have completed successfully.  Example:
#   sbatch --export=ALL,CASE=sidm3,SEED_RUN_ID=323683 run_z0_l17_continue.sbatch
#
# Each invocation stops cleanly before the walltime limit.  Resubmit the same
# command to continue from the newest completed output.
#SBATCH --job-name=sidm_z0_l17
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/slurm_z0_l17_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/slurm_z0_l17_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
case_name=${CASE:?CASE must be cdm or sidm3}
seed_run_id=${SEED_RUN_ID:?SEED_RUN_ID must name the completed L17 array run}
validation_root=/gpfs/kjhan/Hydro/Sidm/validation
force_root=${validation_root}/z4_force_l17
run_root=${validation_root}/z0_production_l17/${case_name}
stop_after_seconds=${STOP_AFTER_SECONDS:-79200}
lock_dir=${run_root}/.continuation_lock
lock_acquired=0

case "${case_name}" in
  cdm)
    seed_output=${force_root}/cdm_${seed_run_id}/output_00002
    source_binary=${force_root}/cdm_${seed_run_id}/ramses_z4_l17
    template=${code}/tests/sidm/z0_l17_continue_cdm.nml.in
    ;;
  sidm3)
    seed_output=${force_root}/sidm3_${seed_run_id}/output_00002
    source_binary=${force_root}/sidm3_${seed_run_id}/ramses_z4_l17
    template=${code}/tests/sidm/z0_l17_continue_sidm3.nml.in
    ;;
  *)
    printf 'unsupported CASE: %s\n' "${case_name}" >&2
    exit 2
    ;;
esac

cleanup() {
  if [[ -n ${watchdog_pid:-} ]]; then
    kill "${watchdog_pid}" 2>/dev/null || true
  fi
  if [[ ${lock_acquired} -eq 1 ]]; then
    rmdir "${lock_dir}" 2>/dev/null || true
  fi
}
trap cleanup EXIT

test -d "${seed_output}"
test -f "${seed_output}/info_00002.txt"
test -x "${source_binary}"
mkdir -p "${run_root}/snapshot"
if ! mkdir "${lock_dir}"; then
  printf 'continuation lock already exists: %s\n' "${lock_dir}" >&2
  exit 4
fi
lock_acquired=1

if ! find -L "${run_root}" -maxdepth 1 -type d -name 'output_*' -print -quit | grep -q .; then
  ln -s "${seed_output}" "${run_root}/output_00002"
fi
if [[ ! -x "${run_root}/ramses_z0_l17" ]]; then
  cp -p "${source_binary}" "${run_root}/ramses_z0_l17"
fi

latest_output=$(find -L "${run_root}" -maxdepth 1 -type d -name 'output_*' -printf '%f\n' | sort -V | tail -1)
restart_number=${latest_output#output_}
restart_number=$((10#${restart_number}))
next_output=$((restart_number + 1))

git -C "${code}" status --short >"${run_root}/snapshot/git_status_${SLURM_JOB_ID}.txt"
git -C "${code}" diff --binary >"${run_root}/snapshot/source_${SLURM_JOB_ID}.diff"
cp -p "${template}" "${run_root}/snapshot/"
cp -p "${code}/tests/sidm/run_z0_l17_continue.sbatch" "${run_root}/snapshot/"
sha256sum "${source_binary}" >"${run_root}/snapshot/source_binary_${SLURM_JOB_ID}.sha256"
sha256sum "${run_root}/ramses_z0_l17" >"${run_root}/snapshot/binary_${SLURM_JOB_ID}.sha256"
sed "s/__NRESTART__/${restart_number}/" "${template}" >"${run_root}/run.nml"
cp -p "${run_root}/run.nml" "${run_root}/snapshot/run_${SLURM_JOB_ID}.nml"
: >"${run_root}/jobcontrol.txt"

(
  sleep "${stop_after_seconds}"
  printf '0 -1\n' >"${run_root}/jobcontrol.txt"
) &
watchdog_pid=$!

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
  cd "${run_root}"
  mpirun -np "${SLURM_NTASKS}" ./ramses_z0_l17 run.nml >"run_${SLURM_JOB_ID}.log" 2>&1
)

test -d "${run_root}/output_$(printf '%05d' "${next_output}")"
grep -Fq 'Job control: stop requested' "${run_root}/run_${SLURM_JOB_ID}.log"
grep -Fq 'Run completed' "${run_root}/run_${SLURM_JOB_ID}.log"
