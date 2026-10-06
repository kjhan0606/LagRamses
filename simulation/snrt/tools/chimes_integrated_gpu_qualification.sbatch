#!/bin/bash
#SBATCH --job-name=chimes-integrated-a100
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/integrated-a100-%j.log
# Bounded native chemistry fixture only; no RAMSES evolution or snapshots.
set -euo pipefail
module purge
module load gnu13 intel/mpi/2021.17 cuda/12.8.1
root=/gpfs/kjhan/chimes-cvode-gpu-20261001
sun=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
chimes=$root/.cvode-gpu/chimes
table=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
build=$root/.cvode-gpu/integrated-$SLURM_JOB_ID
mkdir "$build"
cd "$build"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
unset LD_PRELOAD SNRT_CHIMES_COMPUTE_BACKEND SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND \
 SNRT_CHIMES_RHS_VERIFY SNRT_CHIMES_RHS_DIAGNOSTIC SNRT_CHIMES_CVODE_DIAGNOSTICS \
 SNRT_CHIMES_TEST_REFERENCE SNRT_CHIMES_TEST_STATE_FILE SNRT_CHIMES_TEST_COMPOSITION \
 SNRT_CHIMES_BUNDLE_PROFILE
gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader)
echo "$gpu"
test_arch=sm_80
case "$gpu" in
  *A100*) ;;
  *H100*|*H200*)
    [[ "${SNRT_CHIMES_TILE_REGRESSION_ONLY:-0}" == 1 ]] || exit 2
    test_arch=sm_90 ;;
  *) exit 2 ;;
esac
bash "$root/simulation/snrt/tools/build_chimes_cpu_library.sh" "$chimes" "$build/library" "$sun" /home/kjhan/local trace
export LD_LIBRARY_PATH="$build/library:$sun/lib64:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}"
common=(-O3 -ffp-contract=off -fPIC -DCHIMES_USE_DOUBLE_PRECISION -I"$sun/include" -I/home/kjhan/local/include -I"$chimes/src" -I"$root/patch/lagRamses")
sunlibs=(-L"$sun/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense)
mpicc "${common[@]}" -DSNRT_CHIMES_CUDA -c "$root/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
nvcc -O3 -std=c++17 -arch="$test_arch" --fmad=false -Xcompiler -fopenmp -I"$sun/include" -c "$root/patch/lagRamses/snrt_chimes_cvode_cuda.cu" -o backend.o
nvcc -O3 -std=c++17 -arch="$test_arch" --fmad=false -Xcompiler -fopenmp -DSNRT_CHIMES_RHS_TESTING -DCHIMES_USE_DOUBLE_PRECISION -I"$sun/include" -I/home/kjhan/local/include -I"$chimes/src" -c "$root/patch/lagRamses/snrt_chimes_rhs_cuda.cu" -o rhs.o
nvcc -O3 -std=c++17 -arch="$test_arch" -c "$root/patch/cuRamses/cuda_stream_pool.cu" -o pool.o
mpicxx "${common[@]}" -std=c++17 -fopenmp -DSNRT_CHIMES_RHS_TESTING -c "$root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o test.o
cuda_root=$(dirname "$(dirname "$(command -v nvcc)")")
mpicxx -fopenmp -Wl,--export-dynamic test.o bridge.o backend.o rhs.o pool.o -L"$build/library" -lchimes "${sunlibs[@]}" -L"$cuda_root/lib64" -lcudart -lcublas -L/home/kjhan/local/lib -lhdf5 -lcrypto -ldl -lm -o test
ldd ./test
if ldd ./test | rg 'not found'; then exit 3; fi
sha256sum ./test ./library/libchimes.so

if [ "${SNRT_CHIMES_TILE_REGRESSION_ONLY:-0}" = 1 ]; then
  unset SNRT_CHIMES_COMPUTE_BACKEND
  timeout 180 ./test "$table" tile > tile-regression.log 2>&1
  rg 'PHYSICS|STATUSES|TILE_STAGED|PARITY_PASS|FAIL' tile-regression.log
  exit 0
fi

# Reuse the existing forced-midcell-drop fixture for bounded adapter fixes.
# This is correctness regression only, not a substitute for a RAMSES run or
# a performance promotion. Default qualification coverage is unchanged.
if [ "${SNRT_CHIMES_RETRY_REGRESSION_ONLY:-0}" = 1 ]; then
  unset SNRT_CHIMES_COMPUTE_BACKEND
  timeout 180 ./test "$table" rhs > retry-regression.log 2>&1
  rg 'PHYSICS|RHS_MIDCELL_RETRIES|RHS_CPU_S|PARITY_PASS|FAIL' retry-regression.log
  # Exercise contention between LU and RHS as well as isolated RHS leases.
  SNRT_CHIMES_COMPUTE_BACKEND=cuda_integrated timeout 180 ./test "$table" check-integrated > retry-integrated.log 2>&1
  rg 'PHYSICS|INTEGRATED_COLD|BUNDLE_GPU_WORK|BUNDLE_PASS|FAIL' retry-integrated.log
  exit 0
fi

failures=0
for composition in default mixed; do
  if [ "$composition" = mixed ]; then export SNRT_CHIMES_TEST_COMPOSITION=mixed; else unset SNRT_CHIMES_TEST_COMPOSITION; fi
  for factor in 0.1 1; do
    rc=0
    SNRT_CHIMES_COMPUTE_BACKEND=cuda_integrated SNRT_CHIMES_TEST_DT_FACTOR=$factor \
      timeout 180 ./test "$table" check-integrated > "check-$composition-$factor.log" 2>&1 || rc=$?
    echo "INTEGRATED_QUALIFICATION composition=$composition dt_factor=$factor rc=$rc"
    rg 'integrated compute backend|PHYSICS|STATUSES|INTEGRATED_COLD|BUNDLE_GPU_WORK|FAIL' "check-$composition-$factor.log" || true
    if [ "$rc" -ne 0 ]; then failures=1; fi
  done
done

if [ "$failures" -eq 0 ]; then
  rc=0
  timeout 180 ./test "$table" check-lu > check-lu.log 2>&1 || rc=$?
  echo "LU_MULTISTREAM_QUALIFICATION rc=$rc"
  rg 'PHYSICS|STATUSES|BUNDLE_GPU_WORK|FAIL' check-lu.log || true
  if [ "$rc" -ne 0 ]; then failures=1; fi
fi

if [ "$failures" -eq 0 ]; then
  for composition in default mixed; do
    if [ "$composition" = mixed ]; then export SNRT_CHIMES_TEST_COMPOSITION=mixed; else unset SNRT_CHIMES_TEST_COMPOSITION; fi
    for factor in 0.1 1; do
      SNRT_CHIMES_COMPUTE_BACKEND=cuda_integrated SNRT_CHIMES_TEST_DT_FACTOR=$factor \
        timeout 240 ./test "$table" bench-integrated > "bench-$composition-$factor.log" 2>&1
      echo "INTEGRATED_TIMING composition=$composition dt_factor=$factor"
      rg 'BUNDLE_TIME|BUNDLE_GPU_WORK' "bench-$composition-$factor.log"
    done
  done
fi
echo "INTEGRATED_QUALIFICATION_COMPLETE failures=$failures"
exit "$failures"
