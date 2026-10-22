#!/usr/bin/env bash
# Relink an isolated M5+CHIMES executable with first-cell CHIMES preparation diagnostics.
#SBATCH --job-name=m5-chimes-reconcile-build
#SBATCH --partition=a100,h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-chimes-reconcile-build-%j.sbatch.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly BUILD=$ROOT/.m5-cosmo-dust-chimes-build-419104
readonly BIN=$BUILD/bin/ramses_m5_chimes_reconcile_diag3d
readonly OBJ=$BUILD/bin/snrt_chimes_runtime.o
readonly OBJ_BACKUP=$BUILD/bin/snrt_chimes_runtime.o.pre_reconcile_diag
readonly CHIMES_OBJ=$BUILD/bin/snrt_chimes.o
readonly CHIMES_OBJ_BACKUP=$BUILD/bin/snrt_chimes.o.pre_stage_diag
readonly PHOTO_OBJ=$BUILD/bin/snrt_chimes_photo.o
readonly PHOTO_OBJ_BACKUP=$BUILD/bin/snrt_chimes_photo.o.pre_subnormal_diag

test -f "$OBJ"
test -f "$CHIMES_OBJ"
test -f "$BUILD/bin/snrt_chimes_runtime.o.pre_stage_diag"
test -f "$CHIMES_OBJ_BACKUP"
test -f "$PHOTO_OBJ"
test -f "$BUILD/bin/snrt_chimes_photo.o.pre_input_detail_diag"
test -f "$PHOTO_OBJ_BACKUP"
test ! -e "$OBJ_BACKUP"
test ! -e "$BIN"
cp -p -- "$OBJ" "$OBJ_BACKUP"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2

cd "$BUILD/bin"
make -j1 \
  GITBRANCH="$(git -C "$ROOT" branch --show-current)" \
  GITHASH="$(git -C "$ROOT" rev-parse HEAD)" \
  GITREPO="$(git -C "$ROOT" remote get-url origin)" \
  EXEC="$BUILD/bin/ramses_m5_chimes_reconcile_diag" \
  SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 USE_CUDA=1 USE_FFTW=1 \
  NENER=0 NVECTOR=32 CHIMES_DIR="$BUILD/chimes-install" \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release \
  ramses > "$ROOT/.m5-chimes-reconcile-build-${SLURM_JOB_ID}.log" 2>&1

test -x "$BIN"
if ldd "$BIN" | rg 'not found'; then
  echo 'Unresolved shared-library dependency in diagnostic executable' >&2
  exit 3
fi
rg -F '../patch/lagRamses/snrt_chimes_runtime.f90' "$ROOT/.m5-chimes-reconcile-build-${SLURM_JOB_ID}.log"
sha256sum "$BIN" "$OBJ" "$OBJ_BACKUP" "$CHIMES_OBJ" "$CHIMES_OBJ_BACKUP" "$PHOTO_OBJ" "$PHOTO_OBJ_BACKUP" \
  > "$ROOT/.m5-chimes-reconcile-build-${SLURM_JOB_ID}.sha256"
