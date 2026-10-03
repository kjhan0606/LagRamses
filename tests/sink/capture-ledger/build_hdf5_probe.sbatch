#!/bin/bash
# Isolated-name HDF5 executable with a bounded capture restart census.
#SBATCH --job-name=smbh_capture_probe_build
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_probe_build_%j.log

set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

code=/home/kjhan/BACKUP/lagRamses-capture-ledger
run=/gpfs/kjhan/lagramses_capture_probe_build_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing build directory $run"
    exit 2
fi
mkdir "$run"
echo 'purpose=diagnostic-only HDF5 build for evolved capture restart cloud census'
echo "source_commit=$(git -C "$code" rev-parse HEAD)"
echo "source_status=$(git -C "$code" status --porcelain --untracked-files=no)"
sha256sum "$code/patch/lagRamses/amr_step.jaehyun.f90" "$code/bin/Makefile" "$code/tests/sink/capture-ledger/build_hdf5_probe.sbatch"
make -C "$code/bin" -B -j1 HDF5=1 CAPTURE_RESTART_DIAG=1 EXEC=ramses_capture_probe ramses
test -x "$code/bin/ramses_capture_probe3d"
cp "$code/bin/ramses_capture_probe3d" "$run/ramses_capture_probe3d"
sha256sum "$run/ramses_capture_probe3d"
ldd "$run/ramses_capture_probe3d" | grep -F libhdf5_fortran
echo "probe_binary=$run/ramses_capture_probe3d"
echo 'CAPTURE_RESTART_PROBE_BUILD_PASS'
