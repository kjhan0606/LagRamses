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
sun_root=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
chimes_root=$task_root/.cvode-gpu/chimes
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
mpicc "${common[@]}" -DSNRT_CHIMES_CUDA -c "$task_root/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
nvcc -O3 -std=c++17 -arch=sm_90 --fmad=false -Xcompiler -fopenmp -I"$sun_root/include" -c "$task_root/patch/lagRamses/snrt_chimes_cvode_cuda.cu" -o backend.o
nvcc -O3 -std=c++17 -arch=sm_90 --fmad=false -Xcompiler -fopenmp -DSNRT_CHIMES_RHS_TESTING -DCHIMES_USE_DOUBLE_PRECISION -I"$sun_root/include" -I"$hdf_root/include" -I"$chimes_root/src" -c "$task_root/patch/lagRamses/snrt_chimes_rhs_cuda.cu" -o rhs.o
nvcc -O3 -std=c++17 -arch=sm_90 -c "$task_root/patch/cuRamses/cuda_stream_pool.cu" -o pool.o
mpicxx "${common[@]}" -std=c++17 -fopenmp -DSNRT_CHIMES_RHS_TESTING -I"$task_root/patch/lagRamses" -c "$task_root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o test.o
cuda_root=$(dirname "$(dirname "$(command -v nvcc)")")
mpicxx -fopenmp -Wl,--export-dynamic test.o bridge.o backend.o rhs.o pool.o -L. -lchimes "${sunlibs[@]}" -L"$cuda_root/lib64" -lcudart -lcublas -L"$hdf_root/lib" -lhdf5 -lcrypto -ldl -lm -o test
ldd ./test
sha256sum ./test ./libchimes.so
unset SNRT_COST_EXACT_CACHE SNRT_CHIMES_RHS_BACKEND
SNRT_CHIMES_CVODE_BACKEND=cpu timeout 120 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 rhs
