#!/bin/bash
#SBATCH --job-name=chimes-requalify
#SBATCH --partition=h100,h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/requalify-%j.log
# No RAMSES evolution or snapshots: bounded native chemistry fixtures only.
set -euo pipefail
module purge
module load gnu13 intel/mpi/2021.17 cuda/12.8.1
root=/gpfs/kjhan/chimes-cvode-gpu-20261001
sun=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
chimes=$root/.cvode-gpu/chimes
table=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
build=$root/.cvode-gpu/requalify-$SLURM_JOB_ID
mkdir "$build"
cd "$build"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
unset LD_PRELOAD SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND SNRT_CHIMES_RHS_VERIFY \
 SNRT_CHIMES_RHS_DIAGNOSTIC SNRT_CHIMES_CVODE_DIAGNOSTICS SNRT_CHIMES_TEST_REFERENCE \
 SNRT_CHIMES_TEST_STATE_FILE SNRT_CHIMES_TEST_COMPOSITION SNRT_CHIMES_BUNDLE_PROFILE
gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader)
echo "$gpu"
case "$gpu" in *H100*|*H200*) ;; *) exit 2;; esac
bash "$root/simulation/snrt/tools/build_chimes_cpu_library.sh" "$chimes" "$build/library" "$sun" /home/kjhan/local trace
export LD_LIBRARY_PATH="$build/library:$sun/lib64:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}"
common=(-O3 -ffp-contract=off -fPIC -DCHIMES_USE_DOUBLE_PRECISION -I"$sun/include" -I/home/kjhan/local/include -I"$chimes/src" -I"$root/patch/lagRamses")
sunlibs=(-L"$sun/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense)
mpicc "${common[@]}" -DSNRT_CHIMES_CUDA -c "$root/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
nvcc -O3 -std=c++17 -arch=sm_90 --fmad=false -Xcompiler -fopenmp -I"$sun/include" -c "$root/patch/lagRamses/snrt_chimes_cvode_cuda.cu" -o backend.o
nvcc -O3 -std=c++17 -arch=sm_90 --fmad=false -Xcompiler -fopenmp -DSNRT_CHIMES_RHS_TESTING -DCHIMES_USE_DOUBLE_PRECISION -I"$sun/include" -I/home/kjhan/local/include -I"$chimes/src" -c "$root/patch/lagRamses/snrt_chimes_rhs_cuda.cu" -o rhs.o
nvcc -O3 -std=c++17 -arch=sm_90 -c "$root/patch/cuRamses/cuda_stream_pool.cu" -o pool.o
mpicxx "${common[@]}" -std=c++17 -fopenmp -DSNRT_CHIMES_RHS_TESTING -c "$root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o test.o
cuda_root=$(dirname "$(dirname "$(command -v nvcc)")")
mpicxx -fopenmp -Wl,--export-dynamic test.o bridge.o backend.o rhs.o pool.o -L"$build/library" -lchimes "${sunlibs[@]}" -L"$cuda_root/lib64" -lcudart -lcublas -L/home/kjhan/local/lib -lhdf5 -lcrypto -ldl -lm -o test
ldd ./test
sha256sum ./test ./library/libchimes.so
failures=0
for backend in lu rhs; do
 passed=1
 for composition in default mixed; do
  if [ "$composition" = mixed ]; then export SNRT_CHIMES_TEST_COMPOSITION=mixed; else unset SNRT_CHIMES_TEST_COMPOSITION; fi
  for factor in 0.1 1; do
   rc=0
   SNRT_CHIMES_TEST_DT_FACTOR=$factor timeout 120 ./test "$table" "check-$backend" > "$backend-$composition-$factor.log" 2>&1 || rc=$?
   echo "QUALIFICATION backend=$backend composition=$composition dt_factor=$factor rc=$rc"
   grep -E 'PHYSICS|STATUSES|BUNDLE_GPU_WORK|FAIL' "$backend-$composition-$factor.log" || true
   if [ "$rc" -ne 0 ]; then passed=0; failures=1; fi
  done
 done
 if [ "$passed" -ne 1 ]; then echo "NO_PERFORMANCE_PROMOTION backend=$backend"; continue; fi
 unset SNRT_CHIMES_TEST_COMPOSITION
 export SNRT_CHIMES_TEST_DT_FACTOR=0.1
 timeout 120 ./test "$table" "bench-$backend" > "$backend-bench.log" 2>&1
 grep -E 'BUNDLE_TIME|BUNDLE_GPU_WORK' "$backend-bench.log"
 # Profile only a warmed GPU solve, separately from unprofiled timing.
 SNRT_CHIMES_BUNDLE_PROFILE=1 timeout 180 nsys profile --trace=cuda,osrt --sample=none --cpuctxsw=none \
  --capture-range=cudaProfilerApi --capture-range-end=stop --output="$build/profile-$backend" \
  ./test "$table" "bench-$backend" > "$backend-profile.log" 2>&1
 nsys stats --report cuda_gpu_kern_sum,cuda_gpu_mem_time_sum,cuda_api_sum,osrt_sum --format csv \
  "$build/profile-$backend.nsys-rep" > "$backend-profile-stats.csv" 2>&1
done
echo "REQUALIFICATION_COMPLETE failures=$failures"
exit "$failures"
