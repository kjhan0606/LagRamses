#!/bin/bash
# [RESIZABLE] Small same-rank CPU HDF5 continuous-versus-restart smoke.
#SBATCH --job-name=rz_h5_restart
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --chdir=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
#SBATCH --output=/home/kjhan/BACKUP/lagramses_qa/hdf5_restart_jobs/h5_%j.slurm

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
qa=/home/kjhan/BACKUP/lagramses_qa
job_root="$qa/hdf5_restart_jobs/${SLURM_JOB_ID:?}"
logs="$job_root/logs"
snapshot="$job_root/snapshot"
continuous="$job_root/continuous"
restart="$job_root/restart"
mkdir -p "$logs" "$snapshot" "$continuous" "$restart"

tree_state=$(git -C "$code" status --porcelain --untracked-files=normal)
if [ -n "$tree_state" ]; then
  echo 'HDF5-RESTART: refusing a dirty or untracked worktree'
  printf '%s\n' "$tree_state"
  exit 1
fi
commit=$(git -C "$code" rev-parse HEAD)
branch=$(git -C "$code" branch --show-current)

snapshot_files=(
  tests/hdf5_restart_smoke.sbatch
  tests/hdf5_restart_smoke.nml.in
  tests/make_cuda_ndgp_ic.py
  tests/cuda_ndgp_matched_ic.sha256
  tests/expand_cpuset.sh
  tests/slurm_rank_bind.sh
  tests/validate_affinity_map.sh
  tests/mpi_affinity_probe.f90
)
for relative in "${snapshot_files[@]}"; do
  cp -p "$code/$relative" "$snapshot/$(basename "$relative")"
done

{
  echo "job_id=$SLURM_JOB_ID"
  echo "host=$(hostname)"
  echo "workdir=$code"
  echo "branch=$branch"
  echo "commit=$commit"
  echo 'scope=two-rank same-order CPU HDF5 restart smoke'
  echo 'build=make -B HDF5=1 USE_CUDA= USE_FFTW='
  echo "started=$(date --iso-8601=seconds)"
  sha256sum "$snapshot"/*
} | tee "$logs/provenance.txt"

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:/opt/ohpc/pub/intel/oneapi/compiler/2025.3/bin:$PATH
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false
export OMP_WAIT_POLICY=passive
export RAMSES_DMO_FINITE_DIAG=1
export I_MPI_PIN=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_LIST I_MPI_PIN_RESPECT_CPUSET
unset KMP_AFFINITY KMP_HW_SUBSET GOMP_CPU_AFFINITY
export I_MPI_FABRICS=shm
unset HDF5 USE_CUDA USE_FFTW USE_CUFFTMP FDMDEBUG FDMTRACE AMR_INDEX_CHECK MORTON128

cpuset_expander="$snapshot/expand_cpuset.sh"
rank_bind="$snapshot/slurm_rank_bind.sh"
affinity_validator="$snapshot/validate_affinity_map.sh"
for helper in /usr/bin/taskset "$cpuset_expander" "$rank_bind" \
  "$affinity_validator"; do
  [ -x "$helper" ] || { echo "HDF5-RESTART: missing helper $helper"; exit 1; }
done

allocation_cpu_mask=$(awk -F: '
  /^Cpus_allowed_list:/ {gsub(/[[:space:]]/,"",$2); print $2; exit}
' /proc/self/status)
if phase4_cpu_list=$("$cpuset_expander" "$allocation_cpu_mask" 2); then
  export PHASE4_CPU_LIST=$phase4_cpu_list
else
  rc=$?
  echo "HDF5-RESTART: invalid allocation CPU mask $allocation_cpu_mask (rc=$rc)"
  exit 1
fi
{
  echo "allocation_cpu_mask=$allocation_cpu_mask"
  echo "rank_cpu_list=$PHASE4_CPU_LIST"
  echo 'binding=taskset_by_SLURM_LOCALID'
} | tee -a "$logs/provenance.txt"

affinity_dir="$job_root/affinity"
affinity_probe="$job_root/mpi_affinity_probe"
mkdir -p "$affinity_dir"
mpiifx -O2 "$snapshot/mpi_affinity_probe.f90" -o "$affinity_probe" \
  >"$logs/affinity_build.log" 2>&1
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 2 -c 1 --cpu-bind=none \
  "$rank_bind" "$affinity_probe" "$affinity_dir" 2 \
  >"$logs/affinity_run.log" 2>&1
find "$affinity_dir" -maxdepth 1 -type f -name 'rank_*.txt' -print0 | \
  sort -z | xargs -0 cat | sort -n >"$logs/affinity_map.txt"
"$affinity_validator" "$PHASE4_CPU_LIST" "$logs/affinity_map.txt" 2
echo 'affinity=PASS ranks=2 unique_single_cpu_masks=2' | tee -a "$logs/provenance.txt"

ic_root="$job_root/ic"
python3 "$snapshot/make_cuda_ndgp_ic.py" "$ic_root" \
  --profile matched-coarse-v3 \
  --verify-manifest "$snapshot/cuda_ndgp_matched_ic.sha256" \
  >"$logs/ic_generation.log" 2>&1
cp -p "$ic_root/MANIFEST.sha256" "$logs/ic_manifest.sha256"
sha256sum "$logs/ic_manifest.sha256" | tee -a "$logs/provenance.txt"

render_nml() {
  nrestart=$1
  target=$2
  sed -e "s|__IC_LEVEL__|$ic_root/level_005|g" \
      -e "s|__IC_LEVEL6__|$ic_root/level_006|g" \
      -e "s|__NRESTART__|$nrestart|g" \
      "$snapshot/hdf5_restart_smoke.nml.in" >"$target"
}
render_nml 0 "$continuous/run.nml"
render_nml 5 "$restart/run.nml"
diff -u \
  <(sed 's/^nrestart=.*/nrestart=__NORMALIZED__/' "$continuous/run.nml") \
  <(sed 's/^nrestart=.*/nrestart=__NORMALIZED__/' "$restart/run.nml") \
  >"$logs/nml_normalized.diff"

make -B -C "$code/bin" HDF5=1 USE_CUDA= USE_FFTW= \
  USE_CUFFTMP= FDMDEBUG= FDMTRACE= AMR_INDEX_CHECK= MORTON128= \
  >"$logs/build.log" 2>&1
grep -Eq -- '-DHDF5([[:space:]]|$)' "$logs/build.log"
if grep -Eq -- '-D(HYDRO_CUDA|USE_FFTW|USE_CUFFTMP|FDMDEBUG|FDMTRACE|AMR_INDEX_CHECK|MORTON128)([[:space:]]|$)|-l(fftw|cudart|cufft)' \
    "$logs/build.log"; then
  echo 'HDF5-RESTART: CUDA/FFTW/debug options entered the CPU HDF5 build'
  exit 1
fi
if [ "$(git -C "$code" rev-parse HEAD)" != "$commit" ] || \
   [ -n "$(git -C "$code" status --porcelain --untracked-files=normal)" ]; then
  echo 'HDF5-RESTART: source changed or became dirty during build'
  exit 1
fi

built="$code/bin/ramses_final3d"
exe="$job_root/ramses_hdf5_cpu"
cp -p "$built" "$exe"
chmod 0555 "$exe"
binary_sha=$(sha256sum "$exe" | awk '{print $1}')
strings "$exe" >"$logs/binary_strings.txt"
grep -Fq "$commit" "$logs/binary_strings.txt"
{
  echo "binary_sha256=$binary_sha"
  echo "binary_mtime=$(stat -c %y "$exe")"
} | tee -a "$logs/provenance.txt"

run_case() {
  run_dir=$1
  if (
    cd "$run_dir"
    srun --kill-on-bad-exit=1 --mpi=pmi2 -n 2 -c 1 --cpu-bind=none \
      "$rank_bind" "$exe" run.nml >run.log 2>&1
  ); then
    rc=0
  else
    rc=$?
  fi
  echo "$rc" >"$run_dir/rc.txt"
  [ "$rc" -eq 0 ] || return "$rc"
}

run_case "$continuous"
cont_outputs=$(find "$continuous" -mindepth 1 -maxdepth 1 -type d \
  -name 'output_*' -printf '%f\n' | sort)
[ "$cont_outputs" = $'output_00001\noutput_00002\noutput_00003\noutput_00004\noutput_00005\noutput_00006' ]
for number in 00001 00002 00003 00004 00005 00006; do
  output="$continuous/output_$number"
  [ -f "$output/COMPLETE" ]
  [ -r "$output/data_$number.h5" ]
done
[ -f "$continuous/output_00005/POISSON_PHI_VALID" ]
[ -f "$continuous/output_00006/POISSON_PHI_VALID" ]
[ "$(grep -c 'Run completed' "$continuous/run.log")" -eq 1 ]

h5diff=/home/kjhan/miniconda3/bin/h5diff
h5dump=/home/kjhan/miniconda3/bin/h5dump
[ -x "$h5diff" ] && [ -x "$h5dump" ]
cp -a "$continuous/output_00005" "$restart/output_00005"
for file in data_00005.h5 COMPLETE POISSON_PHI_VALID; do
  cmp "$continuous/output_00005/$file" "$restart/output_00005/$file"
done
sha256sum "$continuous/output_00005"/{data_00005.h5,COMPLETE,POISSON_PHI_VALID} \
  "$restart/output_00005"/{data_00005.h5,COMPLETE,POISSON_PHI_VALID} \
  >"$logs/checkpoint_copy.sha256"

run_case "$restart"
restart_outputs=$(find "$restart" -mindepth 1 -maxdepth 1 -type d \
  -name 'output_*' -printf '%f\n' | sort)
[ "$restart_outputs" = $'output_00005\noutput_00006' ]
[ -f "$restart/output_00006/COMPLETE" ]
[ -f "$restart/output_00006/POISSON_PHI_VALID" ]
[ -r "$restart/output_00006/data_00006.h5" ]
[ "$(grep -c 'Run completed' "$restart/run.log")" -eq 1 ]
for file in data_00005.h5 COMPLETE POISSON_PHI_VALID; do
  cmp "$continuous/output_00005/$file" "$restart/output_00005/$file"
done

legacy_warning='WARNING: IC header carries no omega_b (legacy grafic); using namelist omega_b'
[ "$(grep -F -c "$legacy_warning" "$continuous/run.log" || true)" -eq 2 ]
[ "$(grep -F -c "$legacy_warning" "$restart/run.log" || true)" -eq 0 ]
for log in "$continuous/run.log" "$restart/run.log"; do
  if grep -iE '(^|[^[:alnum:]_])WARN(ING)?:|FATAL:|ERROR:|HDF5-DIAG|increase ngridmax|increase npartmax|MPI_ABORT|forrtl: severe|segmentation fault|error stop|out of memory|NOT converged|failed to converge|Some grid are outside initial conditions sub-volume' \
      "$log" | grep -Fv "$legacy_warning" >/dev/null; then
    echo "HDF5-RESTART: fatal/nonconvergence marker in $log"
    exit 1
  fi
  if grep -iE '(^|[^[:alnum:]_])OOM([^[:alnum:]_]|$)' "$log" >/dev/null; then
    echo "HDF5-RESTART: OOM marker in $log"
    exit 1
  fi
  if grep -E 'NaN_CHK.*(uold= *[1-9]|f= *[1-9]|d0= *[1-9]|rho= *[1-9]|phi= *[1-9]|scalar= *[1-9])' "$log" >/dev/null; then
    echo "HDF5-RESTART: nonzero NaN census in $log"
    exit 1
  fi
  if grep -Eq '\[CUDA_(MG|NGR|PM)|CUDA pool:' "$log"; then
    echo "HDF5-RESTART: CUDA path entered in CPU smoke: $log"
    exit 1
  fi
done

restore_markers=(
  'HDF5 AMR restore done.'
  'HDF5 particle restore done.'
  'HDF5 poisson restore done.'
  'Modified-gravity scalar field restored from HDF5 checkpoint'
  'Poisson restart: valid phi marker accepted; warm start enabled'
)
for marker in "${restore_markers[@]}"; do
  [ "$(grep -F -c "$marker" "$restart/run.log" || true)" -eq 1 ] || {
    echo "HDF5-RESTART: missing or duplicate restore marker: $marker"
    exit 1
  }
done
if grep -iEq 'VARIABLE NCPU / CROSS-ORDERING RESTART|using variable-ncpu path|Forcing load_balance after variable-ncpu restart|using predictor|Scalar field absent' "$restart/run.log"; then
  echo 'HDF5-RESTART: unexpected variable-ncpu/predictor/scalar fallback'
  exit 1
fi

"$h5dump" -H "$continuous/output_00005/data_00005.h5" \
  >"$logs/continuous_output5.header.txt"
"$h5dump" -H "$continuous/output_00006/data_00006.h5" \
  >"$logs/continuous_output6.header.txt"
"$h5dump" -H "$restart/output_00006/data_00006.h5" \
  >"$logs/restart_output6.header.txt"

attr_scalar() {
  file=$1
  attr=$2
  "$h5dump" -a "/header/$attr" "$file" | awk '
    /\(0\):/ {gsub(/,/,"",$2); print $2; exit}
  '
}
[ "$(attr_scalar "$continuous/output_00005/data_00005.h5" ncpu)" = 2 ]
[ "$(attr_scalar "$continuous/output_00005/data_00005.h5" nstep_coarse)" = 4 ]
[ "$(attr_scalar "$continuous/output_00006/data_00006.h5" ncpu)" = 2 ]
[ "$(attr_scalar "$continuous/output_00006/data_00006.h5" nstep_coarse)" = 5 ]
[ "$(attr_scalar "$restart/output_00006/data_00006.h5" ncpu)" = 2 ]
[ "$(attr_scalar "$restart/output_00006/data_00006.h5" nstep_coarse)" = 5 ]

"$h5diff" -r "$continuous/output_00006/data_00006.h5" \
  "$restart/output_00006/data_00006.h5" >"$logs/h5diff_output6.txt" 2>&1
for attr in time aexp; do
  left=$(attr_scalar "$continuous/output_00006/data_00006.h5" "$attr")
  right=$(attr_scalar "$restart/output_00006/data_00006.h5" "$attr")
  [ "$left" = "$right" ]
  echo "$attr=$left" >>"$logs/final_header.txt"
done

{
  echo 'HDF5_RESTART_SMOKE: PASS'
  echo "job_id=$SLURM_JOB_ID"
  echo "commit=$commit"
  echo "binary_sha256=$binary_sha"
  echo 'scope=two-rank same-order CPU HDF5 AMR+particle+gravity+nDGP roundtrip'
  echo "finished=$(date --iso-8601=seconds)"
} | tee "$job_root/PASS" | tee -a "$logs/provenance.txt"
