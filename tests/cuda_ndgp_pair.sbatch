#!/bin/bash
# [RESIZABLE] First-stage characterization run for the small CUDA/nGR gate.
#SBATCH --job-name=rz_cuda_ndgp
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
#SBATCH --output=/home/kjhan/BACKUP/lagramses_qa/cuda_ndgp_jobs/cal_%j.slurm

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
qa=/home/kjhan/BACKUP/lagramses_qa
cpu_probe=${CUDA_NDGP_CPU_PROBE:-0}
case "$cpu_probe" in
  0|1) ;;
  *) echo 'CUDA-NDGP: CUDA_NDGP_CPU_PROBE must be 0 or 1'; exit 2 ;;
esac
run_mode=characterize
if [ "$cpu_probe" -eq 1 ]; then run_mode=cpu_structure_probe; fi
job_root="$qa/cuda_ndgp_jobs/${SLURM_JOB_ID:?}"
logs="$job_root/logs"
snapshot="$job_root/snapshot"
reports="$job_root/reports"
mkdir -p "$logs" "$snapshot" "$reports"

tree_state=$(git -C "$code" status --porcelain --untracked-files=normal)
if [ -n "$tree_state" ]; then
  echo 'CUDA-NDGP: refusing a dirty or untracked worktree' | tee "$logs/provenance.txt"
  printf '%s\n' "$tree_state" | tee -a "$logs/provenance.txt"
  exit 1
fi
commit=$(git -C "$code" rev-parse HEAD)
branch=$(git -C "$code" branch --show-current)

snapshot_files=(
  tests/cuda_ndgp_pair.sbatch
  tests/cuda_ndgp_pair.nml.in
  tests/make_cuda_ndgp_ic.py
  tests/cuda_ndgp_ic.sha256
  tests/cuda_device_probe.cu
  tests/cuda_ndgp_gate.sh
  tests/validate_cuda_ndgp_logs.py
  tests/test_validate_cuda_ndgp_logs.py
  tests/preflight.sh
  tests/compare_cuda_ndgp_outputs.py
  tests/test_compare_cuda_ndgp_outputs.py
  tests/compare_amr_canonical.py
  tests/test_compare_amr_canonical.py
  tests/cuda_block_index_check.sh
  tests/cuda_block_index_census.py
  tests/cuda_block_index_unit.cpp
  tests/cuda_block_index_device.cu
  tests/run_rank_logged.sh
  tests/expand_cpuset.sh
  tests/slurm_rank_bind.sh
  tests/validate_affinity_map.sh
  tests/mpi_affinity_probe.f90
  patch/cuRamses/force_fine.kjhan.f90
  patch/cuRamses/diag_eint.f90
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
  echo "mode=$run_mode"
  echo 'build=make -B USE_CUDA=1 HDF5= USE_FFTW= USE_CUFFTMP= FDMDEBUG= FDMTRACE= AMR_INDEX_CHECK= MORTON128='
  echo 'cuda_streams_compile_default=16 cuda_streams_runtime=1'
  echo 'dmo_finite_diagnostic=RAMSES_DMO_FINITE_DIAG=1'
  echo "started=$(date --iso-8601=seconds)"
} | tee "$logs/provenance.txt"
sha256sum "$snapshot"/* | tee -a "$logs/provenance.txt"

export PATH=/opt/ohpc/pub/cuda/13.0.2/bin:/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:/opt/ohpc/pub/intel/oneapi/compiler/2025.3/bin:$PATH
export LD_LIBRARY_PATH=/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false
export OMP_WAIT_POLICY=passive
export RAMSES_DMO_FINITE_DIAG=1
export I_MPI_PIN=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_LIST I_MPI_PIN_RESPECT_CPUSET
unset KMP_AFFINITY KMP_HW_SUBSET GOMP_CPU_AFFINITY
export I_MPI_FABRICS=shm
unset HDF5 USE_FFTW USE_CUFFTMP FDMDEBUG FDMTRACE AMR_INDEX_CHECK MORTON128

cpuset_expander="$snapshot/expand_cpuset.sh"
rank_bind="$snapshot/slurm_rank_bind.sh"
affinity_validator="$snapshot/validate_affinity_map.sh"
rank_logger="$snapshot/run_rank_logged.sh"
for helper in /usr/bin/taskset "$cpuset_expander" "$rank_bind" \
  "$affinity_validator" "$rank_logger"; do
  [ -x "$helper" ] || { echo "CUDA-NDGP: helper is not executable: $helper"; exit 1; }
done

allocation_cpu_mask=$(awk -F: '
  /^Cpus_allowed_list:/ {gsub(/[[:space:]]/,"",$2); print $2; exit}
' /proc/self/status)
if phase4_cpu_list=$("$cpuset_expander" "$allocation_cpu_mask" 2); then
  export PHASE4_CPU_LIST=$phase4_cpu_list
else
  parse_rc=$?
  echo "CUDA-NDGP: invalid allocation CPU mask $allocation_cpu_mask (rc=$parse_rc)"
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
if ! mpiifx -O2 "$snapshot/mpi_affinity_probe.f90" -o "$affinity_probe" \
    >"$logs/affinity_build.log" 2>&1; then
  echo 'CUDA-NDGP: affinity probe build failed'
  exit 1
fi
if ! srun --kill-on-bad-exit=1 --mpi=pmi2 -n 2 -c 1 --cpu-bind=none \
    "$rank_bind" "$affinity_probe" "$affinity_dir" \
    2 \
    >"$logs/affinity_run.log" 2>&1; then
  echo 'CUDA-NDGP: affinity probe failed'
  exit 1
fi
find "$affinity_dir" -maxdepth 1 -type f -name 'rank_*.txt' -print0 | \
  sort -z | xargs -0 cat | sort -n >"$logs/affinity_map.txt"
if ! "$affinity_validator" "$PHASE4_CPU_LIST" "$logs/affinity_map.txt" 2; then
  echo 'CUDA-NDGP: affinity mapping failed'
  exit 1
fi
echo 'affinity=PASS ranks=2 unique_single_cpu_masks=2' | tee -a "$logs/provenance.txt"

device_probe="$job_root/cuda_device_probe"
device_dir="$job_root/device"
mkdir -p "$device_dir"
if ! nvcc -O2 -gencode arch=compute_80,code=sm_80 \
    "$snapshot/cuda_device_probe.cu" -o "$device_probe" \
    >"$logs/device_build.log" 2>&1; then
  echo 'CUDA-NDGP: CUDA device probe build failed'
  exit 1
fi
if ! srun --kill-on-bad-exit=1 --mpi=pmi2 -n 2 -c 1 --cpu-bind=none \
    "$rank_bind" "$device_probe" "$device_dir" \
    >"$logs/device_run.log" 2>&1; then
  echo 'CUDA-NDGP: CUDA device probe failed'
  exit 1
fi
cat "$device_dir"/rank_*.txt | sort -n >"$logs/device_map.txt"
if ! awk '
  NF<7 || $1 !~ /^[01]$/ || $2 !~ /^[01]$/ || $3!=1 ||
    $4 !~ /^[0-9a-f]+$/ || length($4)!=32 ||
    $5<8589934592 || $6<8589934592 {bad=1}
  {
    rank[$1]++; localid[$2]++; uuid[$4]++; total[$6]++
    if(index($0,"NVIDIA A10")==0) bad=1
  }
  END {
    nuuid=0; for(x in uuid) nuuid++
    ntotal=0; for(x in total) ntotal++
    if(NR!=2 || rank[0]!=1 || rank[1]!=1 || localid[0]!=1 || localid[1]!=1 ||
       nuuid!=1 || ntotal!=1) bad=1
    exit(bad ? 1 : 0)
  }
' "$logs/device_map.txt"; then
  echo 'CUDA-NDGP: device mapping/free-memory validation failed'
  cat "$logs/device_map.txt"
  exit 1
fi
echo 'cuda_device=PASS ranks=2 visible_device_count=1 shared_uuid=1 free_each_ge_8GiB' | \
  tee -a "$logs/provenance.txt"
cat "$logs/device_map.txt" | tee -a "$logs/provenance.txt"

ic_root="$job_root/ic"
if ! python3 "$snapshot/make_cuda_ndgp_ic.py" "$ic_root" \
    --verify-manifest "$snapshot/cuda_ndgp_ic.sha256" \
    >"$logs/ic_generation.log" 2>&1; then
  echo 'CUDA-NDGP: deterministic IC generation or manifest validation failed'
  exit 1
fi
cp -p "$ic_root/MANIFEST.sha256" "$logs/ic_manifest.sha256"
cat "$logs/ic_manifest.sha256" | tee -a "$logs/provenance.txt"

render_nml() {
  variant=$1
  target=$2
  if [ "$variant" = cpu ]; then
    poisson=.false.; scalar=.false.; particle=.false.
  else
    poisson=.true.; scalar=.true.; particle=.true.
  fi
  sed -e "s|__IC_LEVEL__|$ic_root/level_005|g" \
      -e "s|__IC_LEVEL6__|$ic_root/level_006|g" \
      -e "s|__GPU_POISSON__|$poisson|g" \
      -e "s|__GPU_SCALAR__|$scalar|g" \
      -e "s|__GPU_PARTICLE__|$particle|g" \
      "$snapshot/cuda_ndgp_pair.nml.in" >"$target"
  if [ "$cpu_probe" -eq 1 ] && [ "$variant" = cpu ]; then
    sed -i 's/^nstepmax=4$/nstepmax=1/' "$target"
  fi
}

variants=(cpu gpu)
if [ "$cpu_probe" -eq 1 ]; then variants=(cpu); fi
for variant in "${variants[@]}"; do
  run_dir="$job_root/$variant"
  mkdir -p "$run_dir"
  render_nml "$variant" "$run_dir/run.nml"
done

if ! (
  cd "$snapshot"
  python3 -m unittest -v test_compare_amr_canonical.py \
    test_compare_cuda_ndgp_outputs.py test_validate_cuda_ndgp_logs.py
) >"$logs/unit_tests.log" 2>&1; then
  echo 'CUDA-NDGP: comparator/validator unit tests failed'
  exit 1
fi

if ! RUN_DEVICE=1 "$code/tests/cuda_block_index_check.sh" \
    >"$logs/cuda_static_build.log" 2>&1; then
  echo 'CUDA-NDGP: CUDA static/device block-index gate failed'
  exit 1
fi
if ! make -B -C "$code/bin" USE_CUDA=1 HDF5= USE_FFTW= \
    USE_CUFFTMP= FDMDEBUG= FDMTRACE= AMR_INDEX_CHECK= MORTON128= \
    >"$logs/build.log" 2>&1; then
  echo 'CUDA-NDGP: full CUDA-only build failed'
  exit 1
fi
if grep -Eq -- '-D(USE_FFTW|HDF5|USE_CUFFTMP|FDMDEBUG|FDMTRACE|AMR_INDEX_CHECK|MORTON128)([[:space:]]|$)|-l(fftw|hdf5)' \
    "$logs/build.log"; then
  echo 'CUDA-NDGP: forbidden FFTW/HDF5 flags or libraries entered the build'
  exit 1
fi
if [ "$(git -C "$code" rev-parse HEAD)" != "$commit" ] || \
   [ -n "$(git -C "$code" status --porcelain --untracked-files=normal)" ]; then
  echo 'CUDA-NDGP: source HEAD or clean state changed during build'
  exit 1
fi

built="$code/bin/ramses_final3d"
[ -x "$built" ] || { echo 'CUDA-NDGP: built binary is absent'; exit 1; }
strings "$built" >"$logs/binary_strings.txt"
grep -Fq "$commit" "$logs/binary_strings.txt" || {
  echo 'CUDA-NDGP: binary does not embed the source commit'; exit 1;
}
exe="$job_root/ramses_cuda_ndgp_${commit:0:12}"
cp -p "$built" "$exe"
sha256sum "$exe" | tee -a "$logs/provenance.txt"

for variant in "${variants[@]}"; do
  run_dir="$job_root/$variant"
  preflight_args=(--exe "$exe" --commit "$commit" --nml "$run_dir/run.nml")
  if [ "$variant" = gpu ]; then
    preflight_args+=(--need-switch gpu_poisson=.true. \
      --need-switch gpu_scalar=.true. --need-switch gpu_particle=.true.)
  fi
  if ! "$code/tests/preflight.sh" "${preflight_args[@]}" \
      >"$logs/preflight_${variant}.log" 2>&1; then
    echo "CUDA-NDGP: $variant preflight failed"
    exit 1
  fi
  if (
    cd "$run_dir" || exit 2
    srun --kill-on-bad-exit=1 --mpi=pmi2 -n 2 -c 1 --cpu-bind=none \
      "$rank_bind" "$rank_logger" "$run_dir/rank" "$exe" run.nml
  ); then
    rc=0
  else
    rc=$?
  fi
  cat "$run_dir/rank_0.log" "$run_dir/rank_1.log" >"$run_dir/run.log"
  outputs=$(find "$run_dir" -mindepth 1 -maxdepth 1 -type d -name 'output_*' | wc -l)
  completed=$(grep -c 'Run completed' "$run_dir/run.log" 2>/dev/null || true)
  echo "$variant: rc=$rc outputs=$outputs completed=$completed" | \
    tee -a "$logs/provenance.txt"
  if [ "$rc" -ne 0 ]; then
    echo "CUDA-NDGP: $variant runtime failed"
    exit 1
  fi
done

if [ "$cpu_probe" -eq 1 ]; then
  probe_log="$job_root/cpu/run.log"
  [ "$completed" -eq 1 ] || {
    echo "CUDA-NDGP: CPU probe Run completed count is $completed"; exit 1;
  }
  if grep -Fq 'Some grid are outside initial conditions sub-volume' "$probe_log"; then
    echo 'CUDA-NDGP: CPU probe refinement exceeds the pinned IC volume'
    exit 1
  fi
  for level in 5 6; do
    grep -Eq "Level +$level has +4096 grids" "$probe_log" || {
      echo "CUDA-NDGP: CPU probe initial level-$level grid count is not 4096"
      exit 1
    }
    grep -Eq "nDGP level +$level converged.*res=" "$probe_log" || {
      echo "CUDA-NDGP: CPU probe nDGP level $level did not converge"
      exit 1
    }
  done
  touch "$job_root/CPU_PROBE_RUNTIME_COMPLETE"
  echo "finished=$(date --iso-8601=seconds) rc=0 status=CPU_PROBE_RUNTIME_COMPLETE" | \
    tee -a "$logs/provenance.txt"
  exit 0
fi

if ! "$snapshot/cuda_ndgp_gate.sh" "$job_root/cpu" "$job_root/gpu" \
    "$snapshot" characterize '' "$reports" | tee "$logs/gate.log"; then
  echo "finished=$(date --iso-8601=seconds) rc=1" | tee -a "$logs/provenance.txt"
  exit 1
fi

touch "$job_root/CALIBRATION_COMPLETE"
echo "finished=$(date --iso-8601=seconds) rc=0 status=CALIBRATION_COMPLETE" | \
  tee -a "$logs/provenance.txt"
