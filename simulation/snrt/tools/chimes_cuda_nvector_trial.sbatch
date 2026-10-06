#!/bin/bash
#SBATCH --job-name=chimes-gpu-cvode
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/chimes-gpu-cvode-%j.log

set -euo pipefail
module purge
module load gnu13 intel/mpi/2021.17 cuda/13.0.2

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly SUN_BUILD=$ROOT/.cvode-gpu/sundials-cuda-vector-413666/build
readonly SUN_SOURCE=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials
readonly SUN_ROOT=$SUN_BUILD
readonly CHIMES_ROOT=$ROOT/.cvode-gpu/chimes
readonly HDF_ROOT=/home/kjhan/local
readonly WORK=$ROOT/.cvode-gpu/chimes-cuda-nvector-$SLURM_JOB_ID
readonly TEST_NCELL=${SNRT_CHIMES_TEST_NCELL:-24}
mkdir -p "$WORK/chimes"
cd "$WORK"

cmake -S "$SUN_SOURCE" -B "$SUN_BUILD" "-DCMAKE_CUDA_ARCHITECTURES=80;86;90"
cmake --build "$SUN_BUILD" --target sundials_cvode_shared sundials_nvecserial_shared \
  sundials_nveccuda_shared sundials_sunlinsoldense_shared sundials_sunlinsolspgmr_shared \
  sundials_sunmatrixdense_shared -j4
test -f "$SUN_ROOT/include/sundials/sundials_config.h"
test -f "$SUN_SOURCE/include/nvector/nvector_cuda.h"

# Patch a private build copy only; the pinned CHIMES source and the sibling
# LRD_JWST SUNDIALS source/install are not modified.
cp -a "$CHIMES_ROOT/src" "$WORK/chimes/src"
patch -p1 -d "$WORK/chimes" < "$ROOT/simulation/snrt/data/chimes_cvode_cuda_nvector_trial.patch"

export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
sunlibdirs=("$SUN_BUILD/src/cvode" "$SUN_BUILD/src/nvector/serial" \
  "$SUN_BUILD/src/nvector/cuda" "$SUN_BUILD/src/sunlinsol/dense" \
  "$SUN_BUILD/src/sunlinsol/spgmr" \
  "$SUN_BUILD/src/sunmatrix/dense")
sunpaths=()
for path in "${sunlibdirs[@]}"; do sunpaths+=("-L$path" "-Wl,-rpath,$path"); done
export LD_LIBRARY_PATH="$(IFS=:; printf '%s' "${sunlibdirs[*]}"):$WORK:$HDF_ROOT/lib:/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}"
export SNRT_CHIMES_CVODE_BACKEND=cpu SNRT_CHIMES_CVODE_NVECTOR=cpu SNRT_CHIMES_RHS_BACKEND=cpu
unset SNRT_CHIMES_COMPUTE_BACKEND SNRT_CHIMES_RHS_VERIFY SNRT_CHIMES_GPU_RESULT

cuda_root=$(dirname "$(dirname "$(command -v nvcc)")")
common=(-O3 -ffp-contract=off -fPIC -DCHIMES_USE_DOUBLE_PRECISION \
  -I"$SUN_BUILD/include" -I"$SUN_SOURCE/include" -I"$cuda_root/include" \
  -I"$HDF_ROOT/include" -I"$WORK/chimes/src")
sunlibs=("${sunpaths[@]}" -lsundials_cvode -lsundials_nvecserial -lsundials_nveccuda \
  -lsundials_sunlinsoldense -lsundials_sunlinsolspgmr -lsundials_sunmatrixdense)
for unit in chimes chimes_cooling init_chimes rate_equations update_rates; do
  extra=()
  if [[ "$unit" == chimes || "$unit" == rate_equations ]]; then
    extra=(-DSNRT_CHIMES_CVODE_CUDA_NVECTOR)
  fi
  mpicc "${common[@]}" "${extra[@]}" -c "$WORK/chimes/src/$unit.c" -o "$unit.o"
done
mpicc -shared chimes.o chimes_cooling.o init_chimes.o rate_equations.o update_rates.o \
  "${sunlibs[@]}" -L"$HDF_ROOT/lib" -Wl,-rpath,"$HDF_ROOT/lib" -lhdf5 -lm -o libchimes.so

nvflags=(-O3 -std=c++17 --fmad=false -Xcompiler -fopenmp \
  -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 \
  -gencode arch=compute_90,code=sm_90)
mpicc "${common[@]}" -DSNRT_CHIMES_CUDA -I"$ROOT/patch/lagRamses" \
  -c "$ROOT/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
nvcc "${nvflags[@]}" -I"$SUN_BUILD/include" -I"$SUN_SOURCE/include" \
  -c "$ROOT/patch/lagRamses/snrt_chimes_cvode_cuda.cu" -o backend.o
nvcc "${nvflags[@]}" -DCHIMES_USE_DOUBLE_PRECISION -I"$SUN_BUILD/include" -I"$SUN_SOURCE/include" \
  -I"$HDF_ROOT/include" -I"$WORK/chimes/src" \
  -c "$ROOT/patch/lagRamses/snrt_chimes_rhs_cuda.cu" -o rhs.o
nvcc "${nvflags[@]}" -c "$ROOT/patch/cuRamses/cuda_stream_pool.cu" -o pool.o
mpicxx "${common[@]}" -std=c++17 -fopenmp -I"$ROOT/patch/lagRamses" \
  -DSNRT_CHIMES_TEST_NCELL="$TEST_NCELL" \
  -c "$ROOT/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o test.o
mpicxx -fopenmp -Wl,--export-dynamic test.o bridge.o backend.o rhs.o pool.o \
  -L"$WORK" -lchimes "${sunlibs[@]}" -L"$cuda_root/lib64" -lcudart -lcublas \
  -L"$HDF_ROOT/lib" -lhdf5 -lcrypto -ldl -lm -o test

ldd ./test | rg 'sundials_nvec(cuda|serial)|not found'
if ldd ./test | rg -q 'not found'; then exit 3; fi
sha256sum ./test ./libchimes.so "$SUN_BUILD/src/nvector/cuda/libsundials_nveccuda.so.5.8.0"
timeout 240 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 cuda-direct-solve
