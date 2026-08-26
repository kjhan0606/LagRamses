#!/bin/bash
#SBATCH --job-name=sidm3_z0_l18
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/z0_production_l18/sidm3/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/z0_production_l18/sidm3/slurm_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
run_root=/gpfs/kjhan/Hydro/Sidm/validation/z0_production_l18/sidm3
seed_output=/gpfs/kjhan/Hydro/Sidm/validation/z0_restart_smoke/sidm3_323693/output_00010
source_binary=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_sidm3/ramses_zoom3d
template=$code/tests/sidm/z0_l18_continue_sidm3.nml.in
stop_after_seconds=${STOP_AFTER_SECONDS:-79200}
lock_dir=${run_root}/.continuation_lock
lock_acquired=0

cleanup() {
  if [[ -n ${watchdog_pid:-} ]]; then
    kill "$watchdog_pid" 2>/dev/null || true
  fi
  if [[ $lock_acquired -eq 1 ]]; then
    rmdir "$lock_dir" 2>/dev/null || true
  fi
}
trap cleanup EXIT

test -d "$seed_output"
test -x "$source_binary"
mkdir -p "$run_root/snapshot"
if ! mkdir "$lock_dir"; then
  echo "continuation lock already exists: $lock_dir" >&2
  exit 4
fi
lock_acquired=1

if ! find -L "$run_root" -maxdepth 1 -type d -name 'output_*' -print -quit | grep -q .; then
  ln -s "$seed_output" "$run_root/output_00010"
fi
if [[ ! -x "$run_root/ramses_z0_l18" ]]; then
  cp -p "$source_binary" "$run_root/ramses_z0_l18"
fi

latest_output=$(find -L "$run_root" -maxdepth 1 -type d -name 'output_*' -printf '%f\n' | sort -V | tail -1)
restart_number=${latest_output#output_}
restart_number=$((10#$restart_number))
next_output=$((restart_number + 1))

git -C "$code" status --short >"$run_root/snapshot/git_status_${SLURM_JOB_ID}.txt"
git -C "$code" diff --binary >"$run_root/snapshot/source_${SLURM_JOB_ID}.diff"
cp -p "$template" "$run_root/snapshot/"
cp -p "$code/tests/sidm/run_z0_l18_continue_sidm3.sbatch" "$run_root/snapshot/"
sha256sum "$source_binary" >"$run_root/snapshot/source_binary_${SLURM_JOB_ID}.sha256"
sha256sum "$run_root/ramses_z0_l18" >"$run_root/snapshot/binary_${SLURM_JOB_ID}.sha256"
sed "s/__NRESTART__/${restart_number}/" "$template" >"$run_root/run.nml"
cp -p "$run_root/run.nml" "$run_root/snapshot/run_${SLURM_JOB_ID}.nml"
: >"$run_root/jobcontrol.txt"

(
  sleep "$stop_after_seconds"
  printf '0 -1\n' >"$run_root/jobcontrol.txt"
) &
watchdog_pid=$!

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export I_MPI_FABRICS=shm

(
  cd "$run_root"
  mpirun -np "${SLURM_NTASKS}" ./ramses_z0_l18 run.nml >"run_${SLURM_JOB_ID}.log" 2>&1
)

test -d "$run_root/output_$(printf '%05d' "$next_output")"
grep -Fq 'Job control: stop requested' "$run_root/run_${SLURM_JOB_ID}.log"
grep -Fq 'Run completed' "$run_root/run_${SLURM_JOB_ID}.log"
