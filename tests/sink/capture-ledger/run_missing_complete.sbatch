#!/bin/bash
# Negative restart gate: a copied checkpoint without COMPLETE must be refused.
#SBATCH --job-name=smbh_capture_no_complete
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_missing_complete_%j.log

set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export I_MPI_FABRICS=shm

source_run=${CAPTURE_SOURCE_RUN:?set CAPTURE_SOURCE_RUN to the completed smoke directory}
test -f "$source_run/output_00001/COMPLETE"
test -f "$source_run/restart.nml"
run=/gpfs/kjhan/lagramses_capture_missing_complete_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
mkdir "$run"
cp -a "$source_run/output_00001" "$run/output_00001"
cp "$source_run/restart.nml" "$run/effective.nml"
cp "$source_run/legacy_yields.dat" "$run/legacy_yields.dat"
cp "$source_run/ramses_final3d" "$run/ramses_final3d"
# Delete only the marker in this new test copy, never in the source checkpoint.
rm -- "$run/output_00001/COMPLETE"
echo "purpose=negative SMBH capture restart completeness gate"
echo "source_checkpoint=$source_run/output_00001"
echo "effective_namelist=$run/effective.nml"
sha256sum "$run/effective.nml" "$run/ramses_final3d"
cd "$run"
if srun --kill-on-bad-exit=1 -n 1 ./ramses_final3d effective.nml >run.log 2>&1; then
    echo 'FATAL: restart without COMPLETE was accepted'
    exit 1
fi
grep -Fq 'FATAL: SMBH capture restart lacks COMPLETE output 1' run.log
test ! -e capture.jsonl
echo 'MISSING_COMPLETE_RESTART_REJECTED_PASS'
