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
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export I_MPI_PIN=0

source_run=${CAPTURE_SOURCE_RUN:?set CAPTURE_SOURCE_RUN to the completed smoke directory}
test -f "$source_run/output_00001/COMPLETE"
test -f "$source_run/restart.nml"
source_ncpu=$(awk '$1 == "ncpu" {print $3; exit}' "$source_run/output_00001/info_00001.txt")
if [ "$source_ncpu" != "$SLURM_NTASKS" ]; then
    echo "FATAL: source checkpoint ncpu=$source_ncpu differs from Slurm tasks=$SLURM_NTASKS"
    exit 2
fi
run=/gpfs/kjhan/lagramses_capture_missing_complete_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
mkdir "$run"
cp -a "$source_run/output_00001" "$run/output_00001"
cp "$source_run/restart.nml" "$run/effective.nml"
cp "$source_run/capture.jsonl" "$run/capture.jsonl"
cp "$source_run/legacy_yields.dat" "$run/legacy_yields.dat"
cp "$source_run/ramses_final3d" "$run/ramses_final3d"
# Delete only the marker in this new test copy, never in the source checkpoint.
rm -- "$run/output_00001/COMPLETE"
echo "purpose=negative SMBH capture restart completeness gate"
echo "source_checkpoint=$source_run/output_00001"
echo "effective_namelist=$run/effective.nml"
sha256sum "$run/effective.nml" "$run/ramses_final3d"
cd "$run"
ledger_sha_before=$(sha256sum capture.jsonl | awk '{print $1}')
if srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none \
    -n "${SLURM_NTASKS}" ./ramses_final3d effective.nml >run.log 2>&1; then
    echo 'FATAL: restart without COMPLETE was accepted'
    exit 1
fi
grep -Fq 'FATAL: SMBH capture restart lacks COMPLETE output 1' run.log
ledger_sha_after=$(sha256sum capture.jsonl | awk '{print $1}')
test "$ledger_sha_before" = "$ledger_sha_after"
echo "unchanged_capture_ledger_sha256=$ledger_sha_after"
echo 'MISSING_COMPLETE_RESTART_REJECTED_PASS'
