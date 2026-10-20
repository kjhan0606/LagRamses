#!/usr/bin/env bash
# Isolated incremental M5+CHIMES build for RAMSES-state momentum closure checks.
# The baseline build tree and executable are kept untouched.
#SBATCH --job-name=m5-momentum-build
#SBATCH --partition=a100,h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-momentum-build-%j.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly BASE=${M5_BUILD_BASE:-$ROOT/.m5-cosmo-dust-chimes-build-419104}
readonly BUILD=$ROOT/.m5-ramses-momentum-build-$SLURM_JOB_ID
readonly BIN=$BUILD/bin/ramses_m5_ramses_momentum_closure_v13d
readonly EXEC=$BUILD/bin/ramses_m5_ramses_momentum_closure_v1

test -d "$BASE/bin"
test -f "$BASE/chimes-install/libchimes.so"
test ! -e "$BUILD"
cp -a -- "$BASE" "$BUILD"
git -C "$ROOT" show HEAD:bin/Makefile > "$BUILD/bin/Makefile"
git -C "$ROOT" diff --quiet HEAD -- \
  patch/lagRamses/load_balance.kjhan.f90 \
  patch/lagRamses/snrt_moment_amr.f90 \
  patch/lagRamses/snrt_moment_band.cpp \
  patch/lagRamses/snrt_moment_ramses.f90 \
  patch/lagRamses/snrt_regrid.f90

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2

cd "$BUILD/bin"
make -j1 \
  GITBRANCH="$(git -C "$ROOT" branch --show-current)" \
  GITHASH="$(git -C "$ROOT" rev-parse HEAD)" \
  GITREPO="$(git -C "$ROOT" remote get-url origin)" \
  EXEC="$EXEC" \
  SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 USE_CUDA=1 USE_FFTW=1 \
  NENER=0 NVECTOR=32 CHIMES_DIR="$BUILD/chimes-install" \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release \
  ramses > "$ROOT/.m5-momentum-build-$SLURM_JOB_ID.build.log" 2>&1

test -x "$BIN"
if ldd "$BIN" | rg 'not found'; then
  echo 'Unresolved shared library dependency in momentum-closure executable' >&2
  exit 3
fi
strings "$BIN" | rg -F 'M_N radiation momentum finite check failed'
strings "$BIN" | rg -F 'M_N RAMSES momentum closure failed'
sha256sum "$BIN" "$BUILD/bin/snrt_moment_ramses.o" > "$BUILD/build.sha256"
printf 'M5_MOMENTUM_CLOSURE_BUILD_PASS\nBINARY=%s\n' "$BIN" > "$BUILD/build-result.txt"
