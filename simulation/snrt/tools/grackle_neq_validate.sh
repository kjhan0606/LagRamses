#!/usr/bin/env bash
# Execute on lageunha or inside a compute allocation, never on syntax.
set -euo pipefail
case $(hostname -s) in syntax*) echo 'Compute host required' >&2; exit 2;; esac
root=/gpfs/kjhan/chimes-cvode-gpu-20261001
binary=$root/.grackle-neq-build-lageunha-20261006T1943/bin/ramses_grackle3d
destination=${1:?new validation directory required}
test ! -e "$destination"
mkdir "$destination"
module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_STACKSIZE=64M OMP_PROC_BIND=close
sha256sum "$binary" > "$destination/binary.sha256"
for name in omp1 omp2 mpi2 checkpoint restart; do
  mkdir "$destination/$name"
  sed -e 's/nstepmax=1/nstepmax=4/' -e 's/foutput=1000000/foutput=1/' \
    "$root/simulation/snrt/config/grackle_neq_smoke.nml" > "$destination/$name/effective.nml"
done
sed -i 's/nstepmax=4/nstepmax=2/' "$destination/checkpoint/effective.nml"
sed -i 's/nstepmax=4/nstepmax=4, nrestart=2/' "$destination/restart/effective.nml"
for name in omp1 omp2 mpi2 checkpoint; do
  export OMP_NUM_THREADS=2
  ranks=1
  if [[ $name == omp1 ]]; then export OMP_NUM_THREADS=1; fi
  if [[ $name == mpi2 ]]; then ranks=2; fi
  (
    cd "$destination/$name"
    timeout 120s mpiexec -n "$ranks" "$binary" effective.nml > run.log 2>&1
    grep -E 'Run completed' run.log
    if grep -E 'ERROR:|Aborting|Too many errors|failed|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' run.log; then exit 1; fi
  )
done
cp -a "$destination/checkpoint/output_00002" "$destination/restart/output_00002"
(
  cd "$destination/restart"
  export OMP_NUM_THREADS=2
  timeout 120s mpiexec -n 1 "$binary" effective.nml > run.log 2>&1
  grep -E 'Run completed' run.log
  if grep -E 'ERROR:|Aborting|Too many errors|failed|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' run.log; then exit 1; fi
)
echo 'GRACKLE_NEQ_VALIDATION_RUNS_COMPLETED'
