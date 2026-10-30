#!/bin/bash
# Bounded reproduction of the CHIMES tile/GPU-RHS status discrepancy.
# The paired diagnostic defaults to six copies of logical cell 5. The
# level-queue qualification uses a larger independent-cell batch so it crosses
# the runtime's minimum-work threshold and actually leases GPU brokers.
# Paired RHS verification computes CPU rates at each exact GPU-CVODE state,
# then restores the GPU rates so the tested trajectory remains the GPU path.
#SBATCH --job-name=chimes-rhs-tile
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H100:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/tile-rhs-%j.log
set -euo pipefail
module purge
module load gnu13 intel/mpi/2021.17 cuda/13.0.2

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
sun_build=$root/.cvode-gpu/sundials-cuda-vector-413666/build
sun_source=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials
chimes=$root/.cvode-gpu/chimes
table=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
build=$root/.cvode-gpu/tile-rhs-build-$SLURM_JOB_ID
diagnostic_mode=${SNRT_CHIMES_DIAGNOSTIC_MODE:-paired}
if [[ -n ${SNRT_CHIMES_TEST_NCELL:-} ]]; then
    test_ncell=$SNRT_CHIMES_TEST_NCELL
elif [[ "$diagnostic_mode" == level-qualify ]]; then
    # The level scheduler intentionally keeps short levels on CPU; qualify
    # the broker path with enough independent cells to cross its 16*OMP floor.
    test_ncell=512
else
    test_ncell=6
fi
if [[ "$diagnostic_mode" == level-qualify ]] && (( test_ncell < 128 )); then
    echo "level-qualify needs at least 128 cells for 8 OMP workers (16*workers floor)" >&2
    exit 2
fi
solver_mode=${SNRT_CHIMES_DIAGNOSTIC_LIBRARY:-trace}
case "$solver_mode" in reference|trace) ;; *) echo "Unknown solver mode: $solver_mode" >&2; exit 2;; esac
test -r "$sun_build/include/sundials/sundials_config.h"
test -r "$sun_source/include/nvector/nvector_cuda.h"
test -r "$chimes/src/chimes.c"
test -r "$table"
mkdir "$build"
cd "$build"

export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
cuda_root=$(dirname "$(dirname "$(command -v nvcc)")")
sunlibdirs=("$sun_build/src/cvode" "$sun_build/src/nvector/serial" \
    "$sun_build/src/nvector/cuda" "$sun_build/src/sunlinsol/dense" \
    "$sun_build/src/sunlinsol/spgmr" "$sun_build/src/sunmatrix/dense")
sunpaths=()
for path in "${sunlibdirs[@]}"; do sunpaths+=("-L$path" "-Wl,-rpath,$path"); done
export LD_LIBRARY_PATH="$(IFS=:; printf '%s' "${sunlibdirs[*]}"):$cuda_root/lib64:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}"
gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader)
echo "GPU=$gpu"
case "$gpu" in *H100*) ;; *) echo "This sm_90 diagnostic requires H100, got: $gpu" >&2; exit 2;; esac
rg 'SUNDIALS_VERSION |SUNDIALS_(DOUBLE_PRECISION|INT64_T)' "$sun_build/include/sundials/sundials_config.h"

common=(-O3 -ffp-contract=off -fPIC -DCHIMES_USE_DOUBLE_PRECISION \
    -I"$sun_build/include" -I"$sun_source/include" -I"$cuda_root/include" \
    -I/home/kjhan/local/include -I"$chimes/src")
sunlibs=("${sunpaths[@]}" -lsundials_cvode -lsundials_nvecserial -lsundials_nveccuda \
    -lsundials_sunlinsoldense -lsundials_sunlinsolspgmr -lsundials_sunmatrixdense)
for unit in chimes chimes_cooling init_chimes rate_equations update_rates; do
    input="$chimes/src/$unit.c"
    if [[ "$unit" == chimes && "$solver_mode" == trace ]]; then
        patch --batch --fuzz=0 --output=chimes.c "$input" "$root/simulation/snrt/data/chimes_trace_solver.patch"
        input=chimes.c
    fi
    mpicc "${common[@]}" -c "$input" -o "$unit.o"
done
echo "CHIMES_DIAGNOSTIC_LIBRARY=$solver_mode"
mpicc -shared chimes.o chimes_cooling.o init_chimes.o rate_equations.o update_rates.o \
    "${sunlibs[@]}" -L/home/kjhan/local/lib -lhdf5 -lm -o libchimes.so
mpicc "${common[@]}" -DSNRT_CHIMES_CUDA -DSNRT_CHIMES_CVODE_DIAGNOSTICS \
    -c "$root/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
nvcc -O3 -std=c++17 -arch=sm_90 --fmad=false -Xcompiler -fopenmp \
    -I"$sun_build/include" -I"$sun_source/include" -I"$cuda_root/include" \
    -c "$root/patch/lagRamses/snrt_chimes_cvode_cuda.cu" -o backend.o
nvcc -O3 -std=c++17 -arch=sm_90 --fmad=false -Xcompiler -fopenmp -DSNRT_CHIMES_RHS_TESTING \
    -DCHIMES_USE_DOUBLE_PRECISION -I"$sun_build/include" -I"$sun_source/include" \
    -I"$cuda_root/include" -I/home/kjhan/local/include -I"$chimes/src" \
    -c "$root/patch/lagRamses/snrt_chimes_rhs_cuda.cu" -o rhs.o
nvcc -O3 -std=c++17 -arch=sm_90 -c "$root/patch/cuRamses/cuda_stream_pool.cu" -o pool.o
mpicxx "${common[@]}" -std=c++17 -fopenmp -DSNRT_CHIMES_RHS_TESTING -DSNRT_CHIMES_FIXED_DIAGNOSTIC \
    -DSNRT_CHIMES_TEST_NCELL="$test_ncell" -I"$root/patch/lagRamses" \
    -c "$root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o test.o
mpicxx -fopenmp -Wl,--export-dynamic test.o bridge.o backend.o rhs.o pool.o -L. -lchimes \
    "${sunlibs[@]}" -L"$cuda_root/lib64" -lcudart -lcublas -L/home/kjhan/local/lib \
    -lhdf5 -lcrypto -ldl -lm -o test
if ldd ./test | rg 'not found'; then exit 3; fi
sha256sum ./test ./libchimes.so

unset SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND
export SNRT_CHIMES_RHS_DIAGNOSTIC=1
if [[ "$diagnostic_mode" == level-qualify ]]; then
    unset SNRT_CHIMES_RHS_DIAGNOSTIC SNRT_CHIMES_RHS_VERIFY SNRT_CHIMES_CVODE_DIAGNOSTICS \
        SNRT_CHIMES_TEST_CELL_OFFSET SNRT_CHIMES_TEST_CELL_STRIDE \
        SNRT_CHIMES_TEST_TOLERANCE_FACTOR SNRT_CHIMES_TEST_DT_FACTOR
    export SNRT_CHIMES_GPU_BATCH_CELLS=${SNRT_CHIMES_LEVEL_GPU_BATCH_CELLS:-64}
    export SNRT_CHIMES_GPU_RESULT=compact SNRT_CHIMES_GPU_COOLING=1
    export SNRT_CHIMES_GPU_BROKERS=${SNRT_CHIMES_LEVEL_GPU_BROKERS:-4}
    for composition in default mixed; do
        if [[ "$composition" == mixed ]]; then export SNRT_CHIMES_TEST_COMPOSITION=mixed; else unset SNRT_CHIMES_TEST_COMPOSITION; fi
        echo "LEVEL_QUALIFICATION composition=$composition cells=$test_ncell"
        timeout 240 ./test "$table" bench-level
    done
elif [[ "$diagnostic_mode" == fixed-route ]]; then
    unset SNRT_CHIMES_RHS_VERIFY SNRT_CHIMES_CVODE_DIAGNOSTICS \
        SNRT_CHIMES_TEST_COMPOSITION SNRT_CHIMES_TEST_DT_FACTOR
    export SNRT_CHIMES_TEST_CELL_OFFSET=5 SNRT_CHIMES_TEST_CELL_STRIDE=12
    export SNRT_CHIMES_GPU_WORKER=any
    timeout 240 ./test "$table" fixed-route-diagnostic
elif [[ "$diagnostic_mode" == tile-route ]]; then
    unset SNRT_CHIMES_RHS_VERIFY SNRT_CHIMES_CVODE_DIAGNOSTICS \
        SNRT_CHIMES_TEST_CELL_OFFSET SNRT_CHIMES_TEST_CELL_STRIDE \
        SNRT_CHIMES_TEST_COMPOSITION SNRT_CHIMES_TEST_DT_FACTOR
    export SNRT_CHIMES_GPU_BATCH_CELLS=64 SNRT_CHIMES_GPU_RESULT=compact
    timeout 240 ./test "$table" tile-route-diagnostic
elif [[ "$diagnostic_mode" == paired ]]; then
    export SNRT_CHIMES_TEST_CELL_OFFSET=5 SNRT_CHIMES_TEST_CELL_STRIDE=12
    export SNRT_CHIMES_CVODE_DIAGNOSTICS=1
    for trajectory in cpu gpu; do
        export SNRT_CHIMES_RHS_VERIFY="$trajectory"
        echo "BEGIN_TILE_PAIRED_RHS_TRAJECTORY=$trajectory"
        timeout 240 ./test "$table" tile-rhs-diagnostic
        echo "END_TILE_PAIRED_RHS_TRAJECTORY=$trajectory"
    done
else
    echo "Unknown SNRT_CHIMES_DIAGNOSTIC_MODE=$diagnostic_mode; expected paired|tile-route|fixed-route|level-qualify" >&2
    exit 4
fi
