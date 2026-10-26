#!/usr/bin/env bash
# Isolated incremental M5+CHIMES build for RAMSES-state momentum closure checks.
# The baseline build tree and executable are kept untouched.
#SBATCH --job-name=m5-momentum-build
# CUDA toolkit modules are available on syntax's GPU nodes, not grammar normal.
# The compiler does not execute kernels, but reserve one typed GPU for the build.
#SBATCH --partition=h100
#SBATCH --gres=gpu:H100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-momentum-build-%j.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly SOURCE_COMMIT=$(git -C "$ROOT" rev-parse HEAD)
readonly BASE=${M5_BUILD_BASE:-$ROOT/.m5-cosmo-dust-chimes-build-419104}
readonly SUNDIALS_CUDA_SOURCE=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials
readonly SUNDIALS_CUDA_ROOT=$ROOT/.cvode-gpu/sundials-cuda-vector-413666
readonly BUILD=$ROOT/.m5-ramses-momentum-build-$SLURM_JOB_ID
readonly SOURCE_TREE=$BUILD/committed-source
readonly BIN=$BUILD/bin/ramses_m5_ramses_momentum_closure_v13d
readonly EXEC=$BUILD/bin/ramses_m5_ramses_momentum_closure_v1

test -d "$BASE/bin"
test -f "$BASE/chimes-install/libchimes.so"
test -f "$SUNDIALS_CUDA_SOURCE/include/nvector/nvector_cuda.h"
test -f "$SUNDIALS_CUDA_ROOT/install/include/sundials/sundials_config.h"
test -f "$SUNDIALS_CUDA_ROOT/build/src/nvector/cuda/libsundials_nveccuda.so"
test ! -e "$BUILD"
cp -a -- "$BASE" "$BUILD"
git -C "$ROOT" worktree add --detach "$SOURCE_TREE" "$SOURCE_COMMIT"
test -z "$(git -C "$SOURCE_TREE" status --porcelain --untracked-files=normal)"
for source_dir in patch amr aton hydro pm poisson; do
  test -L "$BUILD/$source_dir"
  unlink "$BUILD/$source_dir"
  ln -s "$SOURCE_TREE/$source_dir" "$BUILD/$source_dir"
done
git -C "$SOURCE_TREE" show "$SOURCE_COMMIT:bin/Makefile" > "$BUILD/bin/Makefile"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2

cd "$BUILD/bin"
make -j1 \
  GITBRANCH="$(git -C "$ROOT" branch --show-current)" \
  GITHASH="$SOURCE_COMMIT" \
  GITREPO="$(git -C "$ROOT" remote get-url origin)" \
  EXEC="$EXEC" \
  SNRT=1 DUST_LIVE=1 CHIMES=1 CHIMES_CUDA=1 HDF5=1 USE_CUDA=1 USE_FFTW=1 \
  NENER=0 NVECTOR=32 CHIMES_DIR="$BUILD/chimes-install" \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release \
  CHIMES_CUDA_SUNDIALS_SOURCE_DIR="$SUNDIALS_CUDA_SOURCE" \
  CHIMES_CUDA_SUNDIALS_BUILD_DIR="$SUNDIALS_CUDA_ROOT/build" \
  CHIMES_CUDA_SUNDIALS_INSTALL_DIR="$SUNDIALS_CUDA_ROOT/install" \
  ramses > "$ROOT/.m5-momentum-build-$SLURM_JOB_ID.build.log" 2>&1

test -x "$BIN"
test -f "$BUILD/bin/snrt_chimes_cvode_cuda.o"
test -f "$BUILD/bin/snrt_chimes_rhs_cuda.o"
if ldd "$BIN" | rg 'not found'; then
  echo 'Unresolved shared library dependency in momentum-closure executable' >&2
  exit 3
fi
strings "$BIN" | rg -F 'hybrid (host CVODE/LU, asynchronous device RHS/cooling)'
strings "$BIN" | rg -F 'M_N radiation momentum finite check failed'
strings "$BIN" | rg -F 'M_N RAMSES momentum closure failed'
sha256sum "$BIN" "$BUILD/bin/snrt_moment_ramses.o" \
  "$BUILD/bin/snrt_chimes_bridge.o" "$BUILD/bin/snrt_chimes_cvode_cuda.o" \
  "$BUILD/bin/snrt_chimes_rhs_cuda.o" > "$BUILD/build.sha256"
sha256sum "$SOURCE_TREE/bin/Makefile" \
  "$SOURCE_TREE/patch/lagRamses/snrt_chimes_bridge.c" \
  "$SOURCE_TREE/patch/lagRamses/snrt_chimes_cvode_cuda.cu" \
  "$SOURCE_TREE/patch/lagRamses/snrt_chimes_rhs_cuda.cu" \
  "$SOURCE_TREE/patch/lagRamses/snrt_chimes_rhs_cuda.h" \
  "$SOURCE_TREE/patch/lagRamses/load_balance.kjhan.f90" \
  "$SOURCE_TREE/patch/lagRamses/snrt_moment_amr.f90" \
  "$SOURCE_TREE/patch/lagRamses/snrt_moment_band.cpp" \
  "$SOURCE_TREE/patch/lagRamses/snrt_moment_ramses.f90" \
  "$SOURCE_TREE/patch/lagRamses/snrt_regrid.f90" > "$BUILD/source.sha256"
printf 'M5_MOMENTUM_CLOSURE_BUILD_PASS\nCHIMES_CUDA=1\nSOURCE_COMMIT=%s\nSOURCE_TREE=%s\nBINARY=%s\n' \
  "$SOURCE_COMMIT" "$SOURCE_TREE" "$BIN" > "$BUILD/build-result.txt"
