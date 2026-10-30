#!/usr/bin/env bash
# Isolated incremental M5+CHIMES build for RAMSES-state momentum closure checks.
# The baseline build tree and executable are kept untouched.
# M5_CHIMES_LIBRARY_MODE=trace opts into the separately qualified dark-trace
# CVODE recovery; the default reference library remains unchanged.
#SBATCH --job-name=m5-momentum-build
# CUDA toolkit modules are available on syntax's GPU nodes, not grammar normal.
# The compiler does not execute kernels, but reserve one typed GPU for the build.
# CPU-only build: submit on grammar with --partition=normal --gres=none and M5_BUILD_CUDA=0.
#SBATCH --partition=h100
#SBATCH --gres=gpu:H100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-momentum-build-%j.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly SOURCE_COMMIT=$(git -C "$ROOT" rev-parse HEAD)
readonly BASE=${M5_BUILD_BASE:-$ROOT/.m5-cosmo-dust-chimes-build-419104}
readonly SUNDIALS_CUDA_SOURCE=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials
readonly SUNDIALS_CUDA_ROOT=$ROOT/.cvode-gpu/sundials-cuda-vector-413666
readonly SUNDIALS_CPU_ROOT=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
readonly BUILD_CUDA=${M5_BUILD_CUDA:-1}
readonly BUILD=$ROOT/.m5-ramses-momentum-build-$SLURM_JOB_ID
readonly SOURCE_TREE=$BUILD/committed-source
readonly BIN=$BUILD/bin/ramses_m5_ramses_momentum_closure_v13d
readonly EXEC=$BUILD/bin/ramses_m5_ramses_momentum_closure_v1
readonly CHIMES_LIBRARY_MODE=${M5_CHIMES_LIBRARY_MODE:-reference}
readonly CHIMES_SOURCE_ROOT=$ROOT/.cvode-gpu/chimes
readonly CHIMES_LIBRARY_BUILDER=$ROOT/simulation/snrt/tools/build_chimes_cpu_library.sh
readonly CHIMES_TRACE_PATCH=$ROOT/simulation/snrt/data/chimes_trace_solver.patch
readonly HDF5_CPU_ROOT=/home/kjhan/local/hdf5

case "$BUILD_CUDA" in
  0|1) ;;
  *) echo 'M5_BUILD_CUDA must be 0 or 1' >&2; exit 2 ;;
esac
case "$CHIMES_LIBRARY_MODE" in
  reference|trace) ;;
  *) echo 'M5_CHIMES_LIBRARY_MODE must be reference or trace' >&2; exit 2 ;;
esac
test -d "$BASE/bin"
test -f "$BASE/chimes-install/libchimes.so"
test -f "$SUNDIALS_CPU_ROOT/include/cvode/cvode.h"
test -f "$SUNDIALS_CPU_ROOT/lib64/libsundials_cvode.so"
if [[ "$BUILD_CUDA" == 1 ]]; then
  test -f "$SUNDIALS_CUDA_SOURCE/include/nvector/nvector_cuda.h"
  test -f "$SUNDIALS_CUDA_ROOT/install/include/sundials/sundials_config.h"
  test -f "$SUNDIALS_CUDA_ROOT/build/src/nvector/cuda/libsundials_nveccuda.so"
fi
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
module load intel/compiler/2025.3.0 intel/mpi/2021.17
if [[ "$BUILD_CUDA" == 1 ]]; then module load cuda/13.0.2; fi

if [[ "$CHIMES_LIBRARY_MODE" == trace ]]; then
  test -f "$CHIMES_SOURCE_ROOT/src/chimes.c"
  test -f "$CHIMES_TRACE_PATCH"
  test -f "$HDF5_CPU_ROOT/include/hdf5.h"
  expected_chimes_files=$(printf '%s\n' \
    src/chimes.c src/chimes_cooling.c src/chimes_proto.h src/init_chimes.c \
    src/rate_equations.c src/update_rates.c | sort)
  actual_chimes_files=$(git -C "$CHIMES_SOURCE_ROOT" diff --name-only | sort)
  if [[ "$actual_chimes_files" != "$expected_chimes_files" ]]; then
    echo 'Trace mode requires the six provenance-recorded CHIMES source modifications' >&2
    printf 'Expected:\n%s\nActual:\n%s\n' "$expected_chimes_files" "$actual_chimes_files" >&2
    exit 4
  fi
  readonly CHIMES_SOURCE_COMMIT=$(git -C "$CHIMES_SOURCE_ROOT" rev-parse HEAD)
  readonly CHIMES_SOURCE_DIFF_SHA256=$(git -C "$CHIMES_SOURCE_ROOT" diff --binary -- \
    src/chimes.c src/chimes_cooling.c src/chimes_proto.h src/init_chimes.c \
    src/rate_equations.c src/update_rates.c | sha256sum | awk '{print $1}')
  readonly CHIMES_TRACE_PATCH_SHA256=$(sha256sum "$CHIMES_TRACE_PATCH" | awk '{print $1}')
  readonly CHIMES_TRACE_STAGE=$BUILD/chimes-trace-library-stage
  CC=mpicc bash "$CHIMES_LIBRARY_BUILDER" "$CHIMES_SOURCE_ROOT" "$CHIMES_TRACE_STAGE" \
    "$SUNDIALS_CPU_ROOT" "$HDF5_CPU_ROOT" trace \
    > "$ROOT/.m5-momentum-build-$SLURM_JOB_ID.chimes-trace-build.log" 2>&1
  mv -- "$BUILD/chimes-install/libchimes.so" "$BUILD/chimes-install/libchimes-reference.so"
  cp -- "$CHIMES_TRACE_STAGE/libchimes.so" "$BUILD/chimes-install/libchimes.so"
  nm -D --defined-only "$BUILD/chimes-install/libchimes-reference.so" \
    | awk '{print $NF}' | sort -u > "$BUILD/chimes-reference-symbols.txt"
  nm -D --defined-only "$BUILD/chimes-install/libchimes.so" \
    | awk '{print $NF}' | sort -u > "$BUILD/chimes-trace-symbols.txt"
  if ! diff -u "$BUILD/chimes-reference-symbols.txt" "$BUILD/chimes-trace-symbols.txt" \
    > "$BUILD/chimes-symbol-diff.txt"; then
    echo 'Trace CHIMES library changed the exported symbol ABI' >&2
    cat "$BUILD/chimes-symbol-diff.txt" >&2
    exit 4
  fi
  export LD_LIBRARY_PATH="$SUNDIALS_CPU_ROOT/lib64:$SUNDIALS_CPU_ROOT/lib:$HDF5_CPU_ROOT/lib:${LD_LIBRARY_PATH:-}"
  ldd "$BUILD/chimes-install/libchimes.so" > "$BUILD/chimes-trace-library.ldd"
  if grep -F 'not found' "$BUILD/chimes-trace-library.ldd"; then
    echo 'Trace CHIMES library has an unresolved runtime dependency' >&2
    exit 4
  fi
  {
    printf 'CHIMES_LIBRARY_MODE=trace\nCHIMES_SOURCE_COMMIT=%s\n' "$CHIMES_SOURCE_COMMIT"
    printf 'CHIMES_SOURCE_DIFF_SHA256=%s\nCHIMES_TRACE_PATCH_SHA256=%s\n' \
      "$CHIMES_SOURCE_DIFF_SHA256" "$CHIMES_TRACE_PATCH_SHA256"
    sha256sum "$BUILD/chimes-install/libchimes-reference.so" "$BUILD/chimes-install/libchimes.so"
  } > "$BUILD/chimes-library-provenance.txt"
fi

cd "$BUILD/bin"
if [[ "$BUILD_CUDA" == 0 ]]; then make clean; fi
cuda_flags=("CHIMES_CUDA=$BUILD_CUDA" "USE_CUDA=$BUILD_CUDA")
if [[ "$BUILD_CUDA" == 1 ]]; then
  cuda_flags+=(
    "CHIMES_CUDA_SUNDIALS_SOURCE_DIR=$SUNDIALS_CUDA_SOURCE"
    "CHIMES_CUDA_SUNDIALS_BUILD_DIR=$SUNDIALS_CUDA_ROOT/build"
    "CHIMES_CUDA_SUNDIALS_INSTALL_DIR=$SUNDIALS_CUDA_ROOT/install"
  )
fi
make -j1 \
  GITBRANCH="$(git -C "$ROOT" branch --show-current)" \
  GITHASH="$SOURCE_COMMIT" \
  GITREPO="$(git -C "$ROOT" remote get-url origin)" \
  EXEC="$EXEC" \
  SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 USE_FFTW=1 \
  NENER=0 NVECTOR=32 CHIMES_DIR="$BUILD/chimes-install" \
  SUNDIALS_DIR="$SUNDIALS_CPU_ROOT" \
  "${cuda_flags[@]}" \
  ramses > "$ROOT/.m5-momentum-build-$SLURM_JOB_ID.build.log" 2>&1

test -x "$BIN"
if [[ "$BUILD_CUDA" == 1 ]]; then
  test -f "$BUILD/bin/snrt_chimes_cvode_cuda.o"
  test -f "$BUILD/bin/snrt_chimes_rhs_cuda.o"
fi
if ldd "$BIN" | grep -E 'not found'; then
  echo 'Unresolved shared library dependency in momentum-closure executable' >&2
  exit 3
fi
if [[ "$BUILD_CUDA" == 1 ]]; then
  strings "$BIN" | grep -F 'hybrid (host CVODE/LU, asynchronous device RHS/cooling)'
else
  if ldd "$BIN" | grep -E 'lib(cudart|cuda)\.so'; then
    echo 'CPU-only momentum-closure executable unexpectedly links CUDA runtime libraries' >&2
    exit 3
  fi
fi
strings "$BIN" | grep -F 'M_N radiation momentum finite check failed'
strings "$BIN" | grep -F 'M_N RAMSES momentum closure failed'
build_artifacts=("$BIN" "$BUILD/chimes-install/libchimes.so" \
  "$BUILD/bin/snrt_moment_ramses.o" "$BUILD/bin/snrt_chimes_bridge.o")
if [[ "$CHIMES_LIBRARY_MODE" == trace ]]; then
  build_artifacts+=("$BUILD/chimes-install/libchimes-reference.so")
fi
if [[ "$BUILD_CUDA" == 1 ]]; then
  build_artifacts+=("$BUILD/bin/snrt_chimes_cvode_cuda.o" "$BUILD/bin/snrt_chimes_rhs_cuda.o")
else
  build_artifacts+=("$BUILD/bin/snrt_chimes.o" "$BUILD/bin/snrt_chimes_photo.o")
fi
sha256sum "${build_artifacts[@]}" > "$BUILD/build.sha256"
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
printf 'M5_MOMENTUM_CLOSURE_BUILD_PASS\nCHIMES_CUDA=%s\nSOURCE_COMMIT=%s\nSOURCE_TREE=%s\nBINARY=%s\n' \
  "$BUILD_CUDA" "$SOURCE_COMMIT" "$SOURCE_TREE" "$BIN" > "$BUILD/build-result.txt"
printf 'CHIMES_LIBRARY_MODE=%s\n' "$CHIMES_LIBRARY_MODE" >> "$BUILD/build-result.txt"
if [[ "$CHIMES_LIBRARY_MODE" == trace ]]; then
  cat "$BUILD/chimes-library-provenance.txt" >> "$BUILD/build-result.txt"
fi
