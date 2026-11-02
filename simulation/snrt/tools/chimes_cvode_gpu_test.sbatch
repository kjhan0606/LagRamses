#!/bin/bash
#SBATCH --job-name=chimes-cvode-gpu
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H100:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/job-%j.log
# Other supported typed requests: --partition=h200 --gres=gpu:H200:1,
# a100/A100, a40/A40, or a10/A10. Submit with sbatch overrides.
set -euo pipefail
module purge
module load gnu13 intel/mpi/2021.17 cuda/12.8.1
task_root=/gpfs/kjhan/chimes-cvode-gpu-20261001
sun_root=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
# SUNDIALS 5.8's installed CPU package does not install the optional CUDA
# N_Vector header. The adapter queries the vector ID, so compile against the
# matching 5.8 source header while retaining the installed 5.8 libraries.
sun_source=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials
chimes_root=$task_root/.cvode-gpu/chimes
hdf_root=/home/kjhan/local
task_build=$task_root/.cvode-gpu/build-$SLURM_JOB_ID
mkdir "$task_build"
cd "$task_build"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores
export LD_LIBRARY_PATH="$sun_root/lib64:$task_build:$hdf_root/lib:${LD_LIBRARY_PATH:-}"
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader)
case "$gpu_name" in
    *A100*) cuda_arch=sm_80 ;;
    *A10*|*A40*) cuda_arch=sm_86 ;;
    *H100*|*H200*) cuda_arch=sm_90 ;;
    *) echo "Unsupported GPU for CHIMES parity test: $gpu_name" >&2; exit 2 ;;
esac
printf 'GPU_MODEL=%s CUDA_ARCH=%s\n' "$gpu_name" "$cuda_arch"
git -C "$task_root" rev-parse HEAD | sed 's/^/LAGRAMSES_SOURCE_HEAD=/'
git -C "$chimes_root" rev-parse HEAD | sed 's/^/CHIMES_SOURCE_HEAD=/'
sha256sum "$chimes_root/src/chimes.c" "$chimes_root/src/chimes_cooling.c" \
    "$chimes_root/src/init_chimes.c" "$chimes_root/src/rate_equations.c" \
    "$chimes_root/src/update_rates.c" "$task_root/patch/lagRamses/snrt_chimes_bridge.c" \
    "$task_root/patch/lagRamses/snrt_chimes_cvode_cuda.cu" \
    "$task_root/patch/lagRamses/snrt_chimes_rhs_cuda.cu" \
    "$task_root/patch/cuRamses/cuda_stream_pool.cu" \
    "$task_root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp"
grep -E 'SUNDIALS_VERSION |SUNDIALS_(DOUBLE_PRECISION|INT64_T)' "$sun_root/include/sundials/sundials_config.h"
test -f "$sun_source/include/nvector/nvector_cuda.h"
grep -Eq 'set\(PACKAGE_STRING "SUNDIALS 5\.8\.0"\)' "$sun_source/CMakeLists.txt"
common=(-O3 -ffp-contract=off -fPIC -DCHIMES_USE_DOUBLE_PRECISION -I"$sun_root/include" -I"$hdf_root/include" -I"$chimes_root/src")
sunlibs=(-L"$sun_root/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense)
for unit in chimes chimes_cooling init_chimes rate_equations update_rates; do
    mpicc "${common[@]}" -c "$chimes_root/src/$unit.c" -o "$unit.o"
done
mpicc -shared chimes.o chimes_cooling.o init_chimes.o rate_equations.o update_rates.o "${sunlibs[@]}" -L"$hdf_root/lib" -lhdf5 -lm -o libchimes.so
mpicc "${common[@]}" -DSNRT_CHIMES_CUDA -c "$task_root/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
nvcc -O3 -std=c++17 -arch="$cuda_arch" --fmad=false -Xcompiler -fopenmp -I"$sun_root/include" -I"$sun_source/include" -c "$task_root/patch/lagRamses/snrt_chimes_cvode_cuda.cu" -o backend.o
nvcc -O3 -std=c++17 -arch="$cuda_arch" --fmad=false -Xcompiler -fopenmp -DCHIMES_USE_DOUBLE_PRECISION -I"$sun_root/include" -I"$hdf_root/include" -I"$chimes_root/src" -c "$task_root/patch/lagRamses/snrt_chimes_rhs_cuda.cu" -o rhs.o
nvcc -O3 -std=c++17 -arch="$cuda_arch" -c "$task_root/patch/cuRamses/cuda_stream_pool.cu" -o pool.o
# The level-queue smoke includes four brokers at eight OpenMP workers. Its
# dispatch floor is 16 cells per broker plus 16 CPU cells, hence 80 total.
mpicxx "${common[@]}" -std=c++17 -fopenmp -I"$task_root/patch/lagRamses" -DSNRT_CHIMES_TEST_NCELL=80 -c "$task_root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o test.o
cuda_root=$(dirname "$(dirname "$(command -v nvcc)")")
mpicxx -fopenmp -Wl,--export-dynamic test.o bridge.o backend.o rhs.o pool.o -L. -lchimes "${sunlibs[@]}" -L"$cuda_root/lib64" -lcudart -lcublas -L"$hdf_root/lib" -lhdf5 -lcrypto -ldl -lm -o test
ldd ./test
if ldd ./test | grep -q 'not found'; then exit 3; fi
sha256sum ./test ./libchimes.so "$sun_root/lib64/libsundials_cvode.so.5.8.0"
SNRT_CHIMES_CVODE_BACKEND=cuda_batched_lu timeout 30 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 config
config_rc=0
SNRT_CHIMES_CVODE_BACKEND=invalid timeout 30 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 config || config_rc=$?
test "$config_rc" -eq 5
unset SNRT_CHIMES_CVODE_BACKEND
timeout 240 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
# Isolate the cell RHS backend math and measure the direct CPU/GPU arms.
timeout 240 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 rhs
# Exercise the production level_queue CPU/GPU brokers, busy-pool fallback, and
# compact/full result parity. The rhs mode above does not enter this scheduler.
timeout 240 ./test /gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 tile
