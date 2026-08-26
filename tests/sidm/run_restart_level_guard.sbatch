#!/bin/bash
#SBATCH --job-name=sidm_restart_guard
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/restart_level_guard/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/restart_level_guard/slurm_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
source_output=/gpfs/kjhan/Hydro/Sidm/validation/amr_estimator/323678/levelmax_4/output_00001
hdf5_source=/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/test_ksection/run_hdf5_test/output_00001
root=/gpfs/kjhan/Hydro/Sidm/validation/restart_level_guard/${SLURM_JOB_ID:?}
case_directory=${root}/binary_case
hdf5_case_directory=${root}/hdf5_case
snapshot=${root}/snapshot

test -d "$source_output"
test -f "$hdf5_source/data_00001.h5"
mkdir -p "$case_directory" "$hdf5_case_directory" "$snapshot"

git -C "$code" status --short >"$snapshot/git_status.txt"
git -C "$code" diff --binary >"$snapshot/source.diff"
cp -p "$code/tests/sidm/restart_level_guard.nml" "$snapshot/"
cp -p "$code/tests/sidm/restart_hdf5_level_guard.nml" "$snapshot/"

"$code/bin/build_lageunha.sh"
cp -p "$code/bin/build_lageunha.log" "$snapshot/build_lageunha.log"
grep -Fq '=== make exit=0' "$snapshot/build_lageunha.log"
cp -p "$code/bin/ramses_final3d" "$case_directory/ramses_final3d"
sha256sum "$case_directory/ramses_final3d" >"$snapshot/binary.sha256"

ln -s "$source_output" "$case_directory/output_00001"
cp -p "$snapshot/restart_level_guard.nml" "$case_directory/restart_level_guard.nml"
: >"$case_directory/jobcontrol.txt"

set +e
(
  cd "$case_directory"
  I_MPI_FABRICS=shm mpirun -np 1 \
    ./ramses_final3d restart_level_guard.nml >run.log 2>&1
)
ramses_status=$?
set -e

grep -Fq 'Restart level guard: checkpoint levelmax=   4 exceeds requested levelmax=   3' \
  "$case_directory/run.log"
if grep -Fq 'Starting time integration' "$case_directory/run.log"; then
  echo 'Restart level guard failed: time integration was reached.' >&2
  exit 1
fi

{
  echo "ramses_exit_status=${ramses_status}"
  echo 'RESTART_LEVEL_GUARD_RESULT=downlevel-restart-rejected-before-array-read'
} | tee "$root/restart_level_guard_summary.txt"

ln -s "$hdf5_source" "$hdf5_case_directory/output_00001"
cp -p "$snapshot/restart_hdf5_level_guard.nml" \
  "$hdf5_case_directory/restart_hdf5_level_guard.nml"
: >"$hdf5_case_directory/jobcontrol.txt"

set +e
(
  cd "$hdf5_case_directory"
  I_MPI_FABRICS=shm mpirun -np 1 \
    "$case_directory/ramses_final3d" restart_hdf5_level_guard.nml >run.log 2>&1
)
hdf5_ramses_status=$?
set -e

grep -Fq 'HDF5 restart level guard: checkpoint levelmax=  10 exceeds requested levelmax=   9' \
  "$hdf5_case_directory/run.log"
if grep -Fq 'Starting time integration' "$hdf5_case_directory/run.log"; then
  echo 'HDF5 restart level guard failed: time integration was reached.' >&2
  exit 1
fi

{
  echo "hdf5_ramses_exit_status=${hdf5_ramses_status}"
  echo 'HDF5_RESTART_LEVEL_GUARD_RESULT=downlevel-restart-rejected-before-array-read'
} | tee "$root/restart_hdf5_level_guard_summary.txt"
