#!/bin/bash
#SBATCH --job-name=chimes-rhs-gpu
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/rhs-%j.log
set -euo pipefail
module purge
module load gnu13 intel/mpi/2021.17 cuda/12.8.1
task_root=/gpfs/kjhan/chimes-cvode-gpu-20261001
project_root=/gpfs/kjhan/LRD_JWST
sun_root=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
chimes_root=$task_root/.cvode-gpu/chimes
spgmr_source_root=$project_root/.chimes-transition.YdVRvD/chimes/src
spgmr_chimes_source=$project_root/.chimes-spgmr-pilot.jQ3hgB/patch-final.WRVIc1/src/chimes.c
hdf_root=/home/kjhan/local
task_build=$task_root/.cvode-gpu/rhs-build-$SLURM_JOB_ID
mkdir "$task_build"
cd "$task_build"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores
export LD_LIBRARY_PATH="$sun_root/lib64:$task_build:$hdf_root/lib:${LD_LIBRARY_PATH:-}"
gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader)
echo "$gpu_name"
case "$gpu_name" in *H100*) ;; *) exit 2;; esac
rg '^#define SUNDIALS_(VERSION |DOUBLE_PRECISION|INT64_T)' "$sun_root/include/sundials/sundials_config.h"
common=(-O3 -ffp-contract=off -fPIC -DCHIMES_USE_DOUBLE_PRECISION -I"$sun_root/include" -I"$hdf_root/include" -I"$chimes_root/src")
sunlibs=(-L"$sun_root/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense)
for unit in chimes chimes_cooling init_chimes rate_equations update_rates; do
    mpicc "${common[@]}" -c "$chimes_root/src/$unit.c" -o "$unit.o"
done
mpicc -shared chimes.o chimes_cooling.o init_chimes.o rate_equations.o update_rates.o "${sunlibs[@]}" -L"$hdf_root/lib" -lhdf5 -lm -o libchimes.so
mpicc "${common[@]}" -DSNRT_CHIMES_CUDA -DSNRT_CHIMES_CVODE_DIAGNOSTICS -c "$task_root/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
nvcc -O3 -std=c++17 -arch=sm_90 --fmad=false -Xcompiler -fopenmp -I"$sun_root/include" -c "$task_root/patch/lagRamses/snrt_chimes_cvode_cuda.cu" -o backend.o
nvcc -O3 -std=c++17 -arch=sm_90 --fmad=false -Xcompiler -fopenmp -DSNRT_CHIMES_RHS_TESTING -DCHIMES_USE_DOUBLE_PRECISION -I"$sun_root/include" -I"$hdf_root/include" -I"$chimes_root/src" -c "$task_root/patch/lagRamses/snrt_chimes_rhs_cuda.cu" -o rhs.o
nvcc -O3 -std=c++17 -arch=sm_90 -c "$task_root/patch/cuRamses/cuda_stream_pool.cu" -o pool.o
mpicxx "${common[@]}" -std=c++17 -fopenmp -DSNRT_CHIMES_RHS_TESTING -I"$task_root/patch/lagRamses" -c "$task_root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o test.o
cuda_root=$(dirname "$(dirname "$(command -v nvcc)")")
mpicxx -fopenmp -Wl,--export-dynamic test.o bridge.o backend.o rhs.o pool.o -L. -lchimes "${sunlibs[@]}" -L"$cuda_root/lib64" -lcudart -lcublas -L"$hdf_root/lib" -lhdf5 -lcrypto -ldl -lm -o test
ldd ./test
sha256sum ./test ./libchimes.so
unset SNRT_COST_EXACT_CACHE SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_CVODE_DIAGNOSTICS \
    SNRT_CHIMES_RHS_BACKEND SNRT_CHIMES_RHS_VERIFY SNRT_CHIMES_TEST_DT_FACTOR
diagnostic_failures=0
if [ "${SNRT_CHIMES_RUN_RHS_DIAGNOSTICS:-1}" = "1" ]; then
for dt_factor in 1 0.1; do
    for trajectory in cpu gpu; do
        echo "BEGIN_PAIRED_RHS_DIAGNOSTIC dt_factor=$dt_factor trajectory=$trajectory"
        test_rc=0
        SNRT_CHIMES_RHS_DIAGNOSTIC=1 SNRT_CHIMES_RHS_VERIFY="$trajectory" \
            SNRT_CHIMES_TEST_DT_FACTOR="$dt_factor" SNRT_CHIMES_CVODE_BACKEND=cpu \
            timeout 240 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 rhs || test_rc=$?
        echo "END_PAIRED_RHS_DIAGNOSTIC dt_factor=$dt_factor trajectory=$trajectory rc=$test_rc"
        if [ "$test_rc" -ne 0 ]; then diagnostic_failures=1; fi
    done
done
echo "BEGIN_CPU_CVODE_STATUS_DETAIL"
SNRT_CHIMES_CVODE_DIAGNOSTICS=1 SNRT_CHIMES_TEST_DT_FACTOR=1 \
    timeout 120 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 profile
echo "END_CPU_CVODE_STATUS_DETAIL"
fi
for dt_factor in 1 0.1; do
    echo "DENSE_PROFILE_DT_FACTOR=$dt_factor"
    SNRT_CHIMES_TEST_DT_FACTOR="$dt_factor" timeout 120 ./test \
        /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 profile
done
echo "BEGIN_ADAPTIVE_SPGMR_FAILURE_FIXTURE"
test -r "$spgmr_chimes_source"
mpicc "${common[@]}" -DSNRT_CHIMES_ADAPTIVE_SPGMR -c "$spgmr_chimes_source" -o spgmr_chimes.o
for unit in chimes_cooling init_chimes rate_equations update_rates; do
    mpicc "${common[@]}" -c "$spgmr_source_root/$unit.c" -o "spgmr_$unit.o"
done
mpicc -shared spgmr_chimes.o spgmr_chimes_cooling.o spgmr_init_chimes.o \
    spgmr_rate_equations.o spgmr_update_rates.o "${sunlibs[@]}" \
    -L"$sun_root/lib64" -lsundials_sunlinsolspgmr -L"$hdf_root/lib" -lhdf5 -lm -o libchimes_spgmr.so
mpicxx "${common[@]}" -std=c++17 -fopenmp -DSNRT_CHIMES_RHS_TESTING \
    -DSNRT_CHIMES_CPU_ONLY_TEST_HOOKS -I"$task_root/patch/lagRamses" \
    -c "$task_root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o spgmr_test.o
mpicxx -fopenmp -Wl,--export-dynamic spgmr_test.o bridge.o backend.o rhs.o pool.o \
    -L. -l:libchimes_spgmr.so "${sunlibs[@]}" -L"$cuda_root/lib64" -lcudart -lcublas \
    -L"$hdf_root/lib" -lhdf5 -lcrypto -ldl -lm -o test-spgmr
ldd ./test-spgmr
if ldd ./test-spgmr | rg 'not found'; then exit 3; fi
sha256sum ./test-spgmr ./libchimes_spgmr.so "$spgmr_chimes_source"
unset SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND SNRT_CHIMES_RHS_VERIFY \
    SNRT_CHIMES_CVODE_DIAGNOSTICS SNRT_CHIMES_TEST_DT_FACTOR
for dt_factor in 1 0.1; do
    echo "SPGMR_PROFILE_DT_FACTOR=$dt_factor"
    SNRT_CHIMES_TEST_DT_FACTOR="$dt_factor" timeout 120 ./test-spgmr \
        /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 profile
done
echo "END_ADAPTIVE_SPGMR_FAILURE_FIXTURE"
test "$diagnostic_failures" -eq 0
