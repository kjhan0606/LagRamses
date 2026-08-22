#!/bin/bash
# Current-source legacy-equivalent versus B=64 block grid-major ABBA timing.
# The shared-node result is a directionality characterization, not a dedicated
# node publication benchmark.
#SBATCH --job-name=rz_perf_layout
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
#SBATCH --output=/home/kjhan/BACKUP/lagramses_qa/perf_layout_jobs/perf_%j.slurm

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
qa=/home/kjhan/BACKUP/lagramses_qa
requested_commit=${PERF_LAYOUT_COMMIT:?submit with PERF_LAYOUT_COMMIT=<audited HEAD>}
expected_commit=$(git -C "$code" rev-parse --verify "$requested_commit^{commit}")
job_root="$qa/perf_layout_jobs/${SLURM_JOB_ID:?}"
snapshot="$job_root/snapshot"
logs="$job_root/logs"
reports="$job_root/reports"
runs="$job_root/runs"
mkdir -p "$snapshot" "$logs" "$reports" "$runs" "$job_root/bin"

snapshot_files=(
  tests/perf_layout.sbatch
  tests/perf_layout.nml
  tests/perf_layout_ic.sha256
  tests/perf_layout_legacy.patch
  tests/perf_layout_report.py
  tests/phase2_layout_gate.sh
  tests/compare_amr_canonical.py
  tests/expand_cpuset.sh
  tests/slurm_rank_bind.sh
  tests/validate_affinity_map.sh
  tests/mpi_affinity_probe.f90
)
for relative in "${snapshot_files[@]}"; do
  git -C "$code" show "$expected_commit:$relative" >"$snapshot/$(basename "$relative")"
done
chmod +x "$snapshot"/*.sh "$snapshot"/*.py

{
  echo "job_id=$SLURM_JOB_ID"
  echo "job_name=${SLURM_JOB_NAME:-}"
  echo "workdir=$code"
  echo "source_commit=$expected_commit"
  echo "host=$(hostname -s)"
  echo "started=$(date --iso-8601=seconds)"
  echo 'comparison=current_same_source_legacy_equivalent_vs_block64'
  echo 'order=block_1,legacy_1,legacy_2,block_2'
  echo 'primary_metric=program_Total_elapsed_time'
  echo 'scope=shared_node_directionality_characterization'
  echo 'build=make -B CPU-only O3; HDF5/CUDA/FFTW disabled'
} | tee "$logs/provenance.txt"
sha256sum "$snapshot"/* | tee -a "$logs/provenance.txt"

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:/opt/ohpc/pub/intel/oneapi/compiler/2025.3/bin:$PATH
export LD_LIBRARY_PATH=/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false
export OMP_WAIT_POLICY=passive
export I_MPI_PIN=0
export I_MPI_FABRICS=shm
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_LIST I_MPI_PIN_RESPECT_CPUSET
unset KMP_AFFINITY KMP_HW_SUBSET GOMP_CPU_AFFINITY
unset HDF5 USE_CUDA USE_FFTW USE_CUFFTMP FDMDEBUG FDMTRACE AMR_INDEX_CHECK MORTON128

ic_dir=/home/kjhan/BACKUP/pbh_prod/icz256
(cd "$ic_dir" && sha256sum -c "$snapshot/perf_layout_ic.sha256") | tee "$logs/ic_check.log"
echo "ic_manifest_sha=$(sha256sum "$snapshot/perf_layout_ic.sha256" | awk '{print $1}')" \
  | tee -a "$logs/provenance.txt"

allocation_cpu_mask=$(awk -F: '
  /^Cpus_allowed_list:/ {gsub(/[[:space:]]/,"",$2); print $2; exit}
' /proc/self/status)
if phase4_cpu_list=$("$snapshot/expand_cpuset.sh" "$allocation_cpu_mask" 4); then
  export PHASE4_CPU_LIST=$phase4_cpu_list
else
  parse_rc=$?
  echo "PERF-LAYOUT: invalid allocation CPU mask $allocation_cpu_mask (rc=$parse_rc)"
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
srun --exclusive --kill-on-bad-exit=1 --mpi=pmi2 -n 4 -c 1 --cpu-bind=none \
  "$snapshot/slurm_rank_bind.sh" "$affinity_probe" "$affinity_dir" 4 \
  >"$logs/affinity_run.log" 2>&1
find "$affinity_dir" -maxdepth 1 -type f -name 'rank_*.txt' -print0 | \
  sort -z | xargs -0 cat | sort -n >"$logs/affinity_map.txt"
"$snapshot/validate_affinity_map.sh" "$PHASE4_CPU_LIST" \
  "$logs/affinity_map.txt" 4
echo 'affinity=PASS ranks=4 unique_single_cpu_masks=4' | tee -a "$logs/provenance.txt"

src_block="$job_root/src_block"
src_legacy="$job_root/src_legacy"
git clone --quiet --shared --no-checkout "$code" "$src_block"
git clone --quiet --shared --no-checkout "$code" "$src_legacy"
git -C "$src_block" checkout --quiet --detach "$expected_commit"
git -C "$src_legacy" checkout --quiet --detach "$expected_commit"
test -z "$(git -C "$src_block" status --porcelain --untracked-files=normal)"
test -z "$(git -C "$src_legacy" status --porcelain --untracked-files=normal)"

git -C "$src_legacy" apply --check "$snapshot/perf_layout_legacy.patch"
git -C "$src_legacy" apply "$snapshot/perf_layout_legacy.patch"
legacy_names=$(git -C "$src_legacy" diff --name-only)
expected_names=$'patch/lagRamses/init_amr.f90\npatch/lagRamses/read_params.jaehyun.f90'
test "$legacy_names" = "$expected_names"
test "$(git -C "$src_legacy" diff --numstat | awk '{a+=$1; d+=$2} END {print a, d}')" = '2 2'
test "$(grep -hc '^  amr_block_size=ngridmax$' \
  "$src_legacy/patch/lagRamses/init_amr.f90" \
  "$src_legacy/patch/lagRamses/read_params.jaehyun.f90" | awk '{s+=$1} END{print s}')" = 2
git -C "$src_legacy" diff --check
git -C "$src_legacy" diff >"$logs/legacy_control.diff"

build_one() {
  local variant=$1
  local source=$2
  local build_log="$logs/build_${variant}.log"
  make -B -C "$source/bin" HDF5= USE_CUDA= USE_FFTW= USE_CUFFTMP= \
    FDMDEBUG= FDMTRACE= AMR_INDEX_CHECK= MORTON128= \
    >"$build_log" 2>&1
  if grep -Eq -- '-D(HDF5|USE_CUDA|USE_FFTW|USE_CUFFTMP|FDMDEBUG|FDMTRACE|AMR_INDEX_CHECK|MORTON128)' \
      "$build_log"; then
    echo "PERF-LAYOUT: forbidden build macro in $variant" >&2
    exit 1
  fi
  test -x "$source/bin/ramses_final3d"
  cp -p "$source/bin/ramses_final3d" "$job_root/bin/$variant"
  sha256sum "$job_root/bin/$variant" | tee -a "$logs/provenance.txt"
}

build_one block "$src_block"
build_one legacy "$src_legacy"

capture_co_tenants() {
  local output=$1
  squeue -h -t R -w "$(hostname -s)" -o '%i|%C' | \
    awk -F'|' -v own="$SLURM_JOB_ID" '$1!=own' | sort -u >"$output"
}

run_one() {
  local variant=$1
  local repetition=$2
  local run_dir="$runs/${variant}_${repetition}"
  local executable="$job_root/bin/$variant"
  local start_ns end_ns run_rc
  mkdir -p "$run_dir"
  cp -p "$snapshot/perf_layout.nml" "$run_dir/run.nml"
  {
    echo "started=$(date --iso-8601=seconds)"
    echo "host=$(hostname -s)"
    cat /proc/loadavg
    free -b
  } >"$run_dir/health_before.txt"
  capture_co_tenants "$run_dir/co_tenants.before"
  start_ns=$(date +%s%N)
  if (cd "$run_dir" && srun --exclusive --kill-on-bad-exit=1 --mpi=pmi2 \
      -n 4 -c 1 --cpu-bind=none "$snapshot/slurm_rank_bind.sh" \
      "$executable" run.nml >run.log 2>&1); then
    run_rc=0
  else
    run_rc=$?
  fi
  end_ns=$(date +%s%N)
  echo "$run_rc" >"$run_dir/run.rc"
  echo $((end_ns-start_ns)) >"$run_dir/external_ns.txt"
  capture_co_tenants "$run_dir/co_tenants.after"
  if cmp -s "$run_dir/co_tenants.before" "$run_dir/co_tenants.after"; then
    cp "$run_dir/co_tenants.before" "$run_dir/co_tenants.txt"
  else
    {
      echo CO_TENANTS_CHANGED
      echo BEFORE
      cat "$run_dir/co_tenants.before"
      echo AFTER
      cat "$run_dir/co_tenants.after"
    } >"$run_dir/co_tenants.txt"
  fi
  {
    echo "finished=$(date --iso-8601=seconds)"
    cat /proc/loadavg
    free -b
  } >"$run_dir/health_after.txt"
  if [ "$run_rc" -ne 0 ]; then
    echo "PERF-LAYOUT: ${variant}_${repetition} rc=$run_rc" >&2
    exit "$run_rc"
  fi
}

# ABBA order reduces first-run cache bias and linear shared-node drift.
run_one block 1
run_one legacy 1
run_one legacy 2
run_one block 2

validate_run() {
  local name=$1
  local run_dir="$runs/$name"
  local log="$run_dir/run.log"
  local embedded_commit=$expected_commit
  if [[ $name == legacy_* ]]; then embedded_commit=${expected_commit}-dirty; fi
  test "$(cat "$run_dir/run.rc")" = 0
  test "$(grep -c 'Run completed' "$log")" = 1
  test "$(grep -c 'Total elapsed time:' "$log")" = 1
  test "$(grep -c 'PERF_TIMER_REPORT step=8 interval=8' "$log")" = 1
  test "$(grep 'Main step=' "$log" | sed -E 's/.*Main step=[[:space:]]*([0-9]+).*/\1/' | paste -sd, -)" = '1,2,3,4,5,6,7,8'
  test "$(grep -c 'WARNING: IC header carries no omega_b (legacy grafic); using namelist omega_b' "$log")" = 1
  test -z "$(grep -Ei 'warn' "$log" | grep -Fv 'WARNING: IC header carries no omega_b (legacy grafic); using namelist omega_b' || true)"
  test "$(grep -Ec "^[[:space:]]*last commit[[:space:]]*=[[:space:]]*${embedded_commit}$" "$log")" = 1
  if grep -Eiq 'MPI_ABORT|forrtl|segmentation|ERROR STOP|FATAL|allocation failed|CUDA.*error|failed to converge|NOT converged|increase ngridmax|increase npartmax|out of memory|Some grid are outside' "$log"; then
    echo "PERF-LAYOUT: fatal marker in $name" >&2
    return 1
  fi
  test -n "$(grep 'NaN_CHK' "$log")"
  if grep 'NaN_CHK' "$log" | grep -Eq '=[[:space:]]*[-+]?[1-9][0-9]*'; then
    echo "PERF-LAYOUT: nonzero NaN census in $name" >&2
    return 1
  fi
  test "$(grep -Ec '\[RESIZABLE\].*(GRID_GROW|PARTICLE_GROW)' "$log" || true)" = 0
  test "$(find "$run_dir" -maxdepth 1 -type d -name 'output_*' -printf '%f\n' | sort)" = output_00001
  test -f "$run_dir/output_00001/COMPLETE"
  grep -Eq '^ncpu[[:space:]]*=[[:space:]]*4$' "$run_dir/output_00001/info_00001.txt"
  grep -Eq '^ngridmax[[:space:]]*=[[:space:]]*800000$' "$run_dir/output_00001/info_00001.txt"
  cmp -s "$run_dir/run.nml" "$snapshot/perf_layout.nml"
}

for name in block_1 legacy_1 legacy_2 block_2; do
  validate_run "$name"
done

"$snapshot/phase2_layout_gate.sh" "$runs/legacy_1" "$runs/block_1" \
  | tee "$reports/physics_pair1.log"
"$snapshot/phase2_layout_gate.sh" "$runs/legacy_2" "$runs/block_2" \
  | tee "$reports/physics_pair2.log"
mkdir -p "$job_root/canonical_scratch"
"$snapshot/compare_amr_canonical.py" \
  "$runs/legacy_1/output_00001" "$runs/block_1/output_00001" \
  --scratch "$job_root/canonical_scratch" \
  --json "$reports/canonical_output1.json" \
  | tee "$reports/canonical_output1.log"

python3 "$snapshot/perf_layout_report.py" "$job_root" \
  --json "$reports/performance.json" \
  --markdown "$reports/performance.md" \
  | tee "$reports/performance.log"

{
  echo "finished=$(date --iso-8601=seconds)"
  echo 'correctness=PASS'
  echo 'performance=CHARACTERIZATION_COMPLETE'
} | tee -a "$logs/provenance.txt"
touch "$job_root/PERF_LAYOUT_CHARACTERIZATION_COMPLETE"
