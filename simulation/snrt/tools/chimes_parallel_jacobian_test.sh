#!/bin/bash
# Bounded CPU experiment; run on authorized LagEunha, not the login node.
set -euo pipefail
test "$#" -ge 1 && test "$#" -le 2
experiment=${2:-parallel}
case "$experiment" in
 parallel) fixture=chimes_parallel_jacobian; workers=(1 2 4 8); timed=(0 2 4 8); reverse=(8 4 2 0); marker=PARALLEL_JACOBIAN;;
 cooling_reuse|reaction_reuse) fixture=chimes_jacobian_cooling_reuse; workers=(1); timed=(0 1); reverse=(1 0); marker=JACOBIAN_COOLING_REUSE;;
 sparse_solve|sparse_lu) fixture=chimes_jacobian_cooling_reuse; workers=(1); timed=(0 1); reverse=(1 0); marker='JACOBIAN_COOLING_REUSE|SPARSE_TRIANGULAR|SPARSE_LU';;
 lapack) fixture=chimes_jacobian_cooling_reuse; workers=(1); timed=(0 1); reverse=(1 0); marker='JACOBIAN_COOLING_REUSE|LAPACK_DENSE';;
 semianalytic) fixture=chimes_semianalytic_jacobian; workers=(1); timed=(0 1); reverse=(1 0); marker=SEMIANALYTIC_JAC;;
 semianalytic_sparse) fixture=chimes_semianalytic_jacobian; workers=(1); timed=(0 1); reverse=(1 0); marker='SEMIANALYTIC_JAC|SPARSE_TRIANGULAR|SPARSE_LU';;
 semianalytic_thermal|semianalytic_thermal_sparse|semianalytic_thermal_lapack) fixture=chimes_semianalytic_jacobian; workers=(1); timed=(0 1); reverse=(1 0); marker='SEMIANALYTIC_JAC|SEMIANALYTIC_THERMAL|SPARSE_TRIANGULAR|SPARSE_LU|LAPACK_DENSE';;
 *) exit 2;;
esac
compare_replay() {
 if [[ "$experiment" != semianalytic* ]]; then cmp "$1" "$2"; return; fi
 # Analytic derivatives are not bitwise DQ derivatives. Use the existing
 # C++ fixture endpoint bounds, without changing any solver tolerance.
 awk '
 function number(s) {
   gsub(/D/,"E",s)
   if(s !~ /[Ee]/) sub(/[-+][0-9][0-9][0-9]$/,"E&",s)
   if(s !~ /^ *[-+]?[0-9]+[.][0-9]+[Ee][-+][0-9]+ *$/) bad=1
   return s+0
 }
 function abs(x){return x<0?-x:x}
 {if(substr($0,1,13)!="REPLAY_STATE " || length($0)!=13+158*25)bad=1
  for(i=1;i<=158;i++) {
   b=number(substr($0,14+25*(i-1),25))
   if(NR==1)a[i]=b
   else if(i==1){trel=abs(b-a[i])/abs(a[i]);if(trel>1e-6)bad=1}
   else {e=abs(b-a[i])/(1e-17+1e-6*abs(a[i]));if(e>max)max=e;if(e>1)bad=1}
  }
 }
 END{printf "REPLAY_PARITY temperature_rel=%.17g species_scaled_error=%.17g\n",trel,max
     exit (bad || NR!=2)}' "$1" "$2"
}
root=$(cd "$(dirname "$0")/../../.." && pwd)
build=$(realpath -m "$1")
sun=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
sun_src=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials
chimes=$root/.cvode-gpu/chimes
reference=$root/.cvode-gpu/cpu-build-lageunha-20261001-zero/reference
table=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
mkdir "$build"
cd "$build"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
export LD_LIBRARY_PATH="$reference:$sun/lib64:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}"
unset LD_PRELOAD SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND \
 SNRT_CHIMES_TEST_REFERENCE SNRT_CHIMES_TEST_STATE_FILE SNRT_CHIMES_TEST_COMPOSITION \
 SNRT_CHIMES_JAC_THREADS SNRT_CHIMES_JAC_VERIFY SNRT_CHIMES_RHS_DIAGNOSTIC SNRT_CHIMES_TEST_CAPTURE \
 SNRT_CHIMES_REPLAY_REPEATS SNRT_CHIMES_JAC_REACTION_REUSE SNRT_CHIMES_SPARSE_FACTOR SNRT_CHIMES_JAC_THERMAL_CHAIN
if [[ "$experiment" = semianalytic_thermal* ]]; then export SNRT_CHIMES_JAC_THERMAL_CHAIN=1; fi
if [ "$experiment" = reaction_reuse ] || [ "$experiment" = lapack ] || [[ "$experiment" = sparse_* ]]; then export SNRT_CHIMES_JAC_REACTION_REUSE=1; fi
if [ "$experiment" = sparse_lu ] || [[ "$experiment" = semianalytic*_sparse ]]; then export SNRT_CHIMES_SPARSE_FACTOR=1; fi
flags=(-O3 -ffp-contract=off -fno-fast-math -fPIC -DCHIMES_USE_DOUBLE_PRECISION
 -I"$chimes/src" -I"$sun/include" -I/home/kjhan/local/include -I"$root/patch/lagRamses")
libs=(-L"$reference" -lchimes -L"$sun/lib64" -lsundials_cvode
 -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense
 -L/home/kjhan/local/lib -lhdf5 -lcrypto -ldl -lm)
mpicc "${flags[@]}" -DSNRT_CHIMES_JACOBIAN_TESTING \
 -c "$root/patch/lagRamses/snrt_chimes_bridge.c" -o bridge.o
mpicxx "${flags[@]}" -std=c++17 -fopenmp -DSNRT_CHIMES_CPU_BUILD_ONLY \
 "$root/simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp" bridge.o \
 -Wl,--export-dynamic "${libs[@]}" -o test
extra=()
if [[ "$experiment" = sparse_* ]] || [[ "$experiment" = semianalytic*_sparse ]]; then extra+=("$root/simulation/snrt/tests/fixtures/phase0/chimes_sparse_triangular.c"); fi
if [ "$experiment" = lapack ] || [ "$experiment" = semianalytic_thermal_lapack ]; then
 mkl=/opt/ohpc/pub/intel/oneapi/mkl/2025.3
 export MKL_INTERFACE_LAYER=LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1
 export LD_LIBRARY_PATH="$mkl/lib:$LD_LIBRARY_PATH"
 flags+=(-I"$mkl/include")
 libs+=(-L"$mkl/lib" -lmkl_rt)
 extra+=("$root/simulation/snrt/tests/fixtures/phase0/chimes_lapack_dense.c")
fi
mpicc "${flags[@]}" -fopenmp -shared -I"$sun_src/src/cvode" "${extra[@]}" \
 "$root/simulation/snrt/tests/fixtures/phase0/$fixture.c" \
 "${libs[@]}" -o parallel-jacobian.so
sha256sum test parallel-jacobian.so "$reference/libchimes.so"
for composition in default mixed; do
 for factor in 1 0.1; do
  tag=$composition-$factor
  export SNRT_CHIMES_TEST_DT_FACTOR=$factor
  if [ "$composition" = mixed ]; then export SNRT_CHIMES_TEST_COMPOSITION=mixed;
  else unset SNRT_CHIMES_TEST_COMPOSITION; fi
  SNRT_CHIMES_TEST_STATE_FILE="$build/$tag.state" timeout 120 ./test "$table" cpu-state > "$tag-reference.log" 2>&1
  for threads in "${workers[@]}"; do
   SNRT_CHIMES_TEST_REFERENCE="$build/$tag.state" SNRT_CHIMES_JAC_THREADS=$threads \
    SNRT_CHIMES_JAC_VERIFY=1 timeout 120 env LD_PRELOAD="$build/parallel-jacobian.so" \
    ./test "$table" cpu-state > "$tag-verify-$threads.log" 2>&1
   grep -E "cpu-build-library|$marker" "$tag-verify-$threads.log"
   grep -Eq 'matrices_checked=[1-9][0-9]*' "$tag-verify-$threads.log"
  done
 done
done
unset SNRT_CHIMES_TEST_COMPOSITION
export SNRT_CHIMES_TEST_DT_FACTOR=0.1
for sample in 1 2 3; do
 arms=("${timed[@]}")
 if [ "$sample" = 2 ]; then arms=("${reverse[@]}"); fi
 for threads in "${arms[@]}"; do
  if [ "$threads" = 0 ]; then
   timeout 120 ./test "$table" cpu-bench > "bench-$threads-$sample.log" 2>&1
  else
   SNRT_CHIMES_TEST_REFERENCE="$build/default-0.1.state" SNRT_CHIMES_JAC_THREADS=$threads \
    timeout 120 env LD_PRELOAD="$build/parallel-jacobian.so" ./test "$table" cpu-bench > "bench-$threads-$sample.log" 2>&1
  fi
  echo "JACOBIAN_SAMPLE threads=$threads sample=$sample"
  grep -E "CPU_BUILD_TIME|$marker" "bench-$threads-$sample.log"
 done
done
objects=/gpfs/kjhan/LRD_JWST/.m5-throughput-20260927.riy1Oo
mpiifx -qopenmp -O3 -no-ftz -I"$objects" \
 "$root/simulation/snrt/tests/fixtures/phase0/chimes_dark_mixed_test.f90" \
 "$objects/snrt_thermochemistry.o" "$objects/snrt_agn_source.o" bridge.o \
 "$objects/snrt_chimes_spectrum.o" "$objects/snrt_chimes_photo.o" "$objects/snrt_chimes.o" \
 -Wl,--export-dynamic -Wl,--no-as-needed "${libs[@]}" -lstdc++ -o replay
source /gpfs/kjhan/LRD_JWST/.chimes-sources.r1TaD8/environment.sh
export LD_LIBRARY_PATH="$reference:$sun/lib64:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
export SNRT_CHIMES_REPLAY_FILE=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/live-cell-20260925/live-403698.log
unset SNRT_CHIMES_TARGET_SUBDT_S SNRT_CHIMES_TARGET_REPEATS
timeout 120 ./replay > replay-reference.log 2>&1
grep '^REPLAY_STATE ' replay-reference.log > replay-reference.state
for threads in "${workers[@]}"; do
 SNRT_CHIMES_JAC_THREADS=$threads SNRT_CHIMES_JAC_VERIFY=1 \
  timeout 120 env LD_PRELOAD="$build/parallel-jacobian.so" ./replay > "replay-verify-$threads.log" 2>&1
 grep '^REPLAY_STATE ' "replay-verify-$threads.log" > "replay-verify-$threads.state"
 compare_replay replay-reference.state "replay-verify-$threads.state"
 grep -E "CHIMES_REPLAY_PASS|$marker" "replay-verify-$threads.log"
 grep -Eq 'matrices_checked=[1-9][0-9]*' "replay-verify-$threads.log"
done
export SNRT_CHIMES_REPLAY_REPEATS=100
for sample in 1 2 3 4 5; do
 arms=("${timed[@]}")
 if ((sample%2==0)); then arms=("${reverse[@]}"); fi
 for threads in "${arms[@]}"; do
  log="replay-$threads-$sample.log"
  if [ "$threads" = 0 ]; then timeout 120 ./replay > "$log" 2>&1;
  else SNRT_CHIMES_JAC_THREADS=$threads timeout 120 env LD_PRELOAD="$build/parallel-jacobian.so" ./replay > "$log" 2>&1; fi
  grep '^REPLAY_STATE ' "$log" > "$log.state"
  compare_replay replay-reference.state "$log.state"
  echo "JACOBIAN_REPLAY threads=$threads sample=$sample"
  grep -E "CHIMES_REPLAY_PASS|CHIMES_REPLAY_WARM|$marker" "$log"
 done
done
echo "${marker}_TEST_PASS"
