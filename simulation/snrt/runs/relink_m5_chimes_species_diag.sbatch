#!/usr/bin/env bash
# Relink an isolated M5+CHIMES executable with one-cell negative-species diagnostics.
#SBATCH --job-name=m5-chimes-species-build
#SBATCH --partition=a100,h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-chimes-species-build-%j.sbatch.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly BUILD=$ROOT/.m5-cosmo-dust-chimes-build-419104
readonly BIN=$BUILD/bin/ramses_m5_chimes_species_diag_v2_3d

test -d "$BUILD/bin"
test -f "$BUILD/bin/snrt_chimes_runtime.o.pre_reconcile_diag"
test ! -e "$BIN"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2

cd "$BUILD/bin"
make -j1 \
  GITBRANCH="$(git -C "$ROOT" branch --show-current)" \
  GITHASH="$(git -C "$ROOT" rev-parse HEAD)" \
  GITREPO="$(git -C "$ROOT" remote get-url origin)" \
  EXEC="$BUILD/bin/ramses_m5_chimes_species_diag_v2" \
  SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 USE_CUDA=1 USE_FFTW=1 \
  NENER=0 NVECTOR=32 CHIMES_DIR="$BUILD/chimes-install" \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release \
  ramses > "$ROOT/.m5-chimes-species-build-${SLURM_JOB_ID}.log" 2>&1

test -x "$BIN"
if ldd "$BIN" | rg 'not found'; then
  echo 'Unresolved shared library dependency in diagnostic executable' >&2
  exit 3
fi
rg -F '../patch/lagRamses/snrt_chimes_runtime.f90' "$ROOT/.m5-chimes-species-build-${SLURM_JOB_ID}.log"
sha256sum "$BIN" "$BUILD/bin/snrt_chimes_runtime.o" \
  > "$ROOT/.m5-chimes-species-build-${SLURM_JOB_ID}.sha256"
