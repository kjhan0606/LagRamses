#!/usr/bin/env bash
# CPU chemistry smoke; one GPU requested only to enter a shared GPU partition.
#SBATCH --job-name=grackle-neq-smoke
#SBATCH --partition=h200
#SBATCH --exclude=syn102
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-%j.log
set -euo pipefail
module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
root=/gpfs/kjhan/chimes-cvode-gpu-20261001
build=$root/.grackle-neq-build-$SLURM_JOB_ID
mkdir "$build"
for folder in patch hydro mhd pm amr poisson aton utils; do ln -s "$root/$folder" "$build/$folder"; done
mkdir "$build/bin" "$build/run"
ln -s "$root/bin/Makefile" "$build/bin/Makefile"
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
cd "$build/bin"
# Existing base Makefile lacks complete Fortran module dependencies; a fresh
# build must be serial rather than racing module producers and consumers.
make -j1 GRACKLE_NEQ=1 HDF5=1 USE_CUDA=0 SNRT=0 CHIMES=0 DUST_LIVE=0 \
  EXEC=ramses_grackle USE_FFTW=0 NVECTOR=32
sha256sum ramses_grackle3d
cd "$build/run"
cp "$root/simulation/snrt/config/grackle_neq_smoke.nml" effective.nml
export OMP_NUM_THREADS=2 OMP_STACKSIZE=64M OMP_PROC_BIND=close
timeout 120s "$build/bin/ramses_grackle3d" effective.nml > ramses.log 2>&1
cat ramses.log
if rg -n 'ERROR:|Aborting|Too many errors|failed' ramses.log; then exit 1; fi
# NaN_CHK is also printed for successful runs with zero counters.
if rg -n 'NaN_CHK.*=[[:space:]]*[1-9][0-9]*' ramses.log; then exit 1; fi
rg 'Grackle NEQ: nine species' ramses.log
rg 'Run completed' ramses.log
echo 'GRACKLE_NEQ_RAMSES_SMOKE_COMPLETED'
