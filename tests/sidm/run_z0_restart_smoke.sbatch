#!/bin/bash
#SBATCH --job-name=sidm_z0_restart_smoke
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --array=0-1
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/z0_restart_smoke/slurm_%A_%a.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/z0_restart_smoke/slurm_%A_%a.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
validation_root=/gpfs/kjhan/Hydro/Sidm/validation/z0_restart_smoke

case "${SLURM_ARRAY_TASK_ID:?}" in
  0)
    label=sidm1
    source_run=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_sidm1
    input=z0_restart_smoke_sidm1.nml
    restart_number=10
    next_output=11
    mpi_ranks=32
    ;;
  1)
    label=sidm3
    source_run=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_sidm3
    input=z0_restart_smoke_sidm3.nml
    restart_number=9
    next_output=10
    mpi_ranks=24
    ;;
  *)
    echo "unsupported array task: ${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
    ;;
esac

source_output=${source_run}/output_$(printf '%05d' "$restart_number")
source_binary=${source_run}/ramses_zoom3d
case_root=${validation_root}/${label}_${SLURM_ARRAY_JOB_ID}
snapshot=${case_root}/snapshot

test -d "$source_output"
test -x "$source_binary"
if test -e "$case_root"; then
  echo "refusing to reuse existing case directory: $case_root" >&2
  exit 3
fi
mkdir -p "$snapshot"

git -C "$code" status --short >"$snapshot/git_status.txt"
git -C "$code" diff --binary >"$snapshot/source.diff"
cp -p "$code/tests/sidm/${input}" "$case_root/run.nml"
cp -p "$code/tests/sidm/run_z0_restart_smoke.sbatch" "$snapshot/"
cp -p "$case_root/run.nml" "$snapshot/"
cp -p "$source_binary" "$case_root/ramses_z0_restart"
sha256sum "$source_binary" >"$snapshot/source_binary.sha256"
sha256sum "$case_root/ramses_z0_restart" >"$snapshot/binary.sha256"
ln -s "$source_output" "$case_root/output_$(printf '%05d' "$restart_number")"
printf '0 -1\n' >"$case_root/jobcontrol.txt"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export I_MPI_FABRICS=shm

(
  cd "$case_root"
  mpirun -np "$mpi_ranks" ./ramses_z0_restart run.nml >run.log 2>&1
)

test -d "$case_root/output_$(printf '%05d' "$next_output")"
grep -Fq 'Starting time integration' "$case_root/run.log"
grep -Fq 'Job control: stop requested' "$case_root/run.log"
grep -Fq 'Run completed' "$case_root/run.log"

source_aexp=$(awk '/^aexp[[:space:]]*=/ {print $3; exit}' \
  "$source_output/info_$(printf '%05d' "$restart_number").txt")
next_aexp=$(awk '/^aexp[[:space:]]*=/ {print $3; exit}' \
  "$case_root/output_$(printf '%05d' "$next_output")/info_$(printf '%05d' "$next_output").txt")
awk -v source_aexp="$source_aexp" -v next_aexp="$next_aexp" \
  'BEGIN {exit !(next_aexp > source_aexp)}'

{
  echo "source_output=$(basename "$source_output")"
  echo "next_output=output_$(printf '%05d' "$next_output")"
  echo "source_aexp=${source_aexp}"
  echo "next_aexp=${next_aexp}"
  echo 'Z0_RESTART_SMOKE_RESULT=same-level-restart-and-clean-checkpoint-passed'
} | tee "$case_root/restart_smoke_summary.txt"
