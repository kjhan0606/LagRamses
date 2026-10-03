#!/bin/bash
# Bounded CPU build for the capture-ledger HDF5 restart smoke.
#SBATCH --job-name=smbh_capture_hdf5_build
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_hdf5_build_%j.log

set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

code=/home/kjhan/BACKUP/lagRamses-capture-ledger
test -f /home/kjhan/local/hdf5/include/hdf5.mod
test -f /home/kjhan/local/hdf5/lib/libhdf5_fortran.so
echo "purpose=HDF5 CPU build for bounded capture-ledger restart structural smoke"
echo "source_commit=$(git -C "$code" rev-parse HEAD)"
echo "source_status=$(git -C "$code" status --porcelain --untracked-files=no)"
sha256sum "$code/patch/lagRamses/backup_hdf5.f90" \
    "$code/patch/lagRamses/sink_particle.kjhan.f90" \
    "$code/tests/sink/capture-ledger/build_hdf5.sbatch"
# The historical Makefile has incomplete Fortran-module dependency edges;
# forcing a parallel full rebuild can read modules while they are replaced.
make -C "$code/bin" -B -j1 HDF5=1 ramses
test -x "$code/bin/ramses_final3d"
sha256sum "$code/bin/ramses_final3d"
ldd "$code/bin/ramses_final3d" | grep -F libhdf5_fortran
echo 'CAPTURE_HDF5_BUILD_PASS'
