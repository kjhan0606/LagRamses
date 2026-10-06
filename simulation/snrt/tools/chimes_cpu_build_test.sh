#!/bin/bash
# Run on an allocated compute node or the authorized LagEunha server.
set -euo pipefail
test "$#" -eq 1
task_root=$(cd "$(dirname "$0")/../../.." && pwd)
sun=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
chimes=$task_root/.cvode-gpu/chimes
hdf=/home/kjhan/local
table=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
build=$(realpath -m "$1")
candidate=${SNRT_CHIMES_CPU_CANDIDATE:-lto}
case "$candidate" in lto|inactive|zero) ;; *) exit 2;; esac
mkdir "$build"
cd "$build"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
export LD_LIBRARY_PATH="$sun/lib64:$hdf/lib:${LD_LIBRARY_PATH:-}"
unset SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND SNRT_CHIMES_RHS_VERIFY \
    SNRT_CHIMES_RHS_DIAGNOSTIC SNRT_CHIMES_CVODE_DIAGNOSTICS SNRT_COST_EXACT_CACHE \
    SNRT_CHIMES_TEST_COMPOSITION SNRT_CHIMES_TEST_REFERENCE SNRT_CHIMES_TEST_STATE_FILE \
    SNRT_COST_JACOBIAN_AGE SNRT_COST_JACOBIAN_COLD_ONLY
for arm in reference "$candidate"; do
    bash "$task_root/simulation/snrt/tools/build_chimes_cpu_library.sh" \
        "$chimes" "$build/$arm" "$sun" "$hdf" "$arm"
done
common=(-O3 -ffp-contract=off -fPIC -DCHIMES_USE_DOUBLE_PRECISION -I"$sun/include" -I"$hdf/include" -I"$chimes/src")
sunlibs=(-L"$sun/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense)
mpicc "${common[@]}" -c "$task_root/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
mpicxx "${common[@]}" -std=c++17 -fopenmp -DSNRT_CHIMES_CPU_BUILD_ONLY \
    -I"$task_root/patch/lagRamses" -c "$task_root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" -o test.o
mpicxx -fopenmp -Wl,--export-dynamic test.o bridge.o -L"$build/reference" -lchimes \
    "${sunlibs[@]}" -L"$hdf/lib" -lhdf5 -lcrypto -ldl -lm -o test
sha256sum ./test "$table"
for arm in reference "$candidate"; do
    LD_LIBRARY_PATH="$build/$arm:$LD_LIBRARY_PATH" ldd ./test
    LD_LIBRARY_PATH="$build/$arm:$LD_LIBRARY_PATH" SNRT_CHIMES_TEST_DT_FACTOR=1 \
        timeout 120 ./test "$table" cpu-state > "$arm-state.log" 2>&1
    grep -E '^CPU_STATE' "$arm-state.log" > "$arm-state.txt"
done
cmp reference-state.txt "$candidate-state.txt"
for sample in 1 2 3 4 5; do
    arms=(reference "$candidate")
    if (( sample % 2 == 0 )); then arms=("$candidate" reference); fi
    for arm in "${arms[@]}"; do
        LD_LIBRARY_PATH="$build/$arm:$LD_LIBRARY_PATH" SNRT_CHIMES_TEST_DT_FACTOR=0.1 \
            timeout 120 ./test "$table" cpu-bench > "$arm-$sample.log" 2>&1
        grep -E '^CPU_STATE' "$arm-$sample.log" > "$arm-$sample-state.txt"
        echo "CPU_BUILD_SAMPLE arm=$arm sample=$sample"
        grep -E '^CPU_(STATE|BUILD_TIME)' "$arm-$sample.log"
    done
    cmp "reference-$sample-state.txt" "$candidate-$sample-state.txt"
done
for arm in reference "$candidate"; do
    LD_LIBRARY_PATH="$build/$arm:$LD_LIBRARY_PATH" SNRT_CHIMES_TEST_DT_FACTOR=0.1 \
        SNRT_CHIMES_TEST_COMPOSITION=mixed timeout 120 ./test "$table" cpu-state \
        > "$arm-mixed.log" 2>&1
    grep -E '^CPU_STATE' "$arm-mixed.log" > "$arm-mixed-state.txt"
    echo "CPU_BUILD_MIXED arm=$arm"
    cat "$arm-mixed-state.txt"
done
cmp reference-mixed-state.txt "$candidate-mixed-state.txt"
echo CHIMES_CPU_BUILD_EXACT_STATE_PASS
