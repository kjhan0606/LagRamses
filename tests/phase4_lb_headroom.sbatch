#!/bin/bash
# [RESIZABLE] Production-default LB headroom/grid-growth paired gate.
#SBATCH --job-name=rz_p4_lbhead
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --chdir=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
#SBATCH --output=/home/kjhan/BACKUP/lagramses_qa/gridgrow_jobs/lbhead_%j.slurm

set -uo pipefail

code=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
qa=/home/kjhan/BACKUP/lagramses_qa
pbh_table=/home/kjhan/BACKUP/lagRamses-PBH/papers/pbh_working/tables/pbh_evap_table_ast_s1.dat
pbh_table_sha=204888b8d5154785b837db280f9b06d2d55caf9a0583fbc0877292da521e906c
job_root="$qa/gridgrow_jobs/${SLURM_JOB_ID:?}"
logs="$job_root/logs"
mkdir -p "$logs"

tree_state=$(git -C "$code" status --porcelain --untracked-files=normal)
if [ -n "$tree_state" ]; then
  echo 'PHASE4-LB: refusing to build a dirty or untracked worktree' | \
    tee "$logs/provenance.txt"
  printf '%s\n' "$tree_state" | tee -a "$logs/provenance.txt"
  exit 1
fi

snapshot="$job_root/snapshot"
mkdir -p "$snapshot"
cp -p "$code/tests/phase4_lb_control.nml" "$snapshot/control.nml"
cp -p "$code/tests/phase4_lb_candidate.nml" "$snapshot/candidate.nml"
cp -p "$code/tests/phase4_lb_gate.sh" "$snapshot/phase4_lb_gate.sh"
cp -p "$code/tests/phase4_lb_headroom.sbatch" "$snapshot/phase4_lb_headroom.sbatch"
snapshot_control="$snapshot/control.nml"
snapshot_candidate="$snapshot/candidate.nml"
snapshot_gate="$snapshot/phase4_lb_gate.sh"

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:/opt/ohpc/pub/intel/oneapi/compiler/2025.3/bin:$PATH
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=false
export OMP_WAIT_POLICY=passive
export I_MPI_PIN=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_LIST I_MPI_PIN_RESPECT_CPUSET
unset KMP_AFFINITY KMP_HW_SUBSET GOMP_CPU_AFFINITY
export I_MPI_FABRICS=shm

cpuset_expander="$code/tests/expand_cpuset.sh"
rank_bind="$code/tests/slurm_rank_bind.sh"
affinity_validator="$code/tests/validate_affinity_map.sh"
if [ ! -x /usr/bin/taskset ] || [ ! -x "$cpuset_expander" ] || \
   [ ! -x "$rank_bind" ] || [ ! -x "$affinity_validator" ]; then
  echo 'PHASE4-LB: taskset or CPU-binding helper is not executable'
  exit 1
fi

allocation_cpu_mask=$(awk -F: '
  /^Cpus_allowed_list:/ {
    gsub(/[[:space:]]/,"",$2)
    print $2
    exit
  }
' /proc/self/status)
if [ -z "$allocation_cpu_mask" ]; then
  echo 'PHASE4-LB: could not read the batch allocation CPU mask'
  exit 1
fi
if phase4_cpu_list=$("$cpuset_expander" "$allocation_cpu_mask" 32); then
  export PHASE4_CPU_LIST=$phase4_cpu_list
else
  parse_rc=$?
  echo "PHASE4-LB: allocation mask is not 32 canonical unique CPUs: $allocation_cpu_mask (rc=$parse_rc)"
  exit 1
fi

commit=$(git -C "$code" rev-parse HEAD)
branch=$(git -C "$code" branch --show-current)
{
  echo "job_id=$SLURM_JOB_ID"
  echo "host=$(hostname)"
  echo "workdir=$code"
  echo "branch=$branch"
  echo "commit=$commit"
  echo "allocation_cpu_mask=$allocation_cpu_mask"
  echo "phase4_cpu_list=$PHASE4_CPU_LIST"
  echo 'binding=taskset_by_SLURM_LOCALID'
  echo 'gate=production_default_lb_headroom'
  echo "started=$(date --iso-8601=seconds)"
} | tee "$logs/provenance.txt"
sha256sum "$snapshot_control" "$snapshot_candidate" "$snapshot_gate" \
  "$snapshot/phase4_lb_headroom.sbatch" | tee -a "$logs/provenance.txt"

if [ ! -r "$pbh_table" ]; then
  echo "PHASE4-LB: PBH table is not readable: $pbh_table" | \
    tee -a "$logs/provenance.txt"
  exit 1
fi
actual_pbh_sha=$(sha256sum "$pbh_table" | awk '{print $1}')
echo "pbh_table=$pbh_table sha256=$actual_pbh_sha" | tee -a "$logs/provenance.txt"
if [ "$actual_pbh_sha" != "$pbh_table_sha" ]; then
  echo 'PHASE4-LB: PBH table checksum mismatch' | tee -a "$logs/provenance.txt"
  exit 1
fi
for nml in "$snapshot_control" "$snapshot_candidate"; do
  if ! grep -Fq "pbh_table_file='$pbh_table'" "$nml"; then
    echo "PHASE4-LB: $nml does not reference the pinned PBH table" | \
      tee -a "$logs/provenance.txt"
    exit 1
  fi
done

affinity_dir="$job_root/affinity"
affinity_probe="$job_root/mpi_affinity_probe"
mkdir -p "$affinity_dir"
if ! mpiifx -O2 "$code/tests/mpi_affinity_probe.f90" -o "$affinity_probe" \
    >"$logs/affinity_build.log" 2>&1; then
  echo 'PHASE4-LB: affinity probe build failed' | tee -a "$logs/provenance.txt"
  exit 1
fi
if ! srun --kill-on-bad-exit=1 --mpi=pmi2 -n 32 -c 1 \
    --cpu-bind=none "$rank_bind" "$affinity_probe" "$affinity_dir" \
    >"$logs/affinity_run.log" 2>&1; then
  echo 'PHASE4-LB: affinity probe execution failed' | tee -a "$logs/provenance.txt"
  exit 1
fi
affinity_file_count=$(find "$affinity_dir" -maxdepth 1 -type f \
  -name 'rank_*.txt' | wc -l)
if [ "$affinity_file_count" -ne 32 ]; then
  echo "PHASE4-LB: affinity probe wrote $affinity_file_count files, expected 32" | \
    tee -a "$logs/provenance.txt"
  exit 1
fi
find "$affinity_dir" -maxdepth 1 -type f -name 'rank_*.txt' -print0 | \
  sort -z | xargs -0 cat | sort -n >"$logs/affinity_map.txt"
if ! "$affinity_validator" "$PHASE4_CPU_LIST" \
    "$logs/affinity_map.txt" 32; then
  echo 'PHASE4-LB: affinity map does not match expected CPUs' | \
    tee -a "$logs/provenance.txt"
  cat "$logs/affinity_map.txt" | tee -a "$logs/provenance.txt"
  exit 1
fi
echo 'affinity=PASS ranks=32 localids=32 expected_unique_single_cpu_masks=32' | \
  tee -a "$logs/provenance.txt"

syntax_files=(
  patch/cuRamses/timing.f90
  patch/cuRamses/amr_commons.kjhan.f90
  patch/cuRamses/radiation_commons.f90
  patch/cuRamses/refine_utils.f90
  patch/lagRamses/pm_commons.f90
  patch/lagRamses/load_balance.kjhan.f90
  patch/lagRamses/read_params.jaehyun.f90
)
if ! "$code/tests/syntax_check.sh" "$code" "${syntax_files[@]}" \
    >"$logs/syntax.log" 2>&1; then
  echo 'PHASE4-LB: syntax gate failed' | tee -a "$logs/provenance.txt"
  exit 1
fi

if ! make -C "$code/bin" HDF5=1 USE_FFTW=1 >"$logs/build.log" 2>&1; then
  echo 'PHASE4-LB: build failed' | tee -a "$logs/provenance.txt"
  exit 1
fi

built="$code/bin/ramses_final3d"
[ -x "$built" ] || { echo 'PHASE4-LB: missing built binary'; exit 1; }
binary_epoch=$(stat -c %Y "$built")
commit_epoch=$(git -C "$code" log -1 --format=%ct "$commit")
[ "$binary_epoch" -ge "$commit_epoch" ] || \
  { echo 'PHASE4-LB: binary predates commit'; exit 1; }

exe="$job_root/ramses_phase4_lb_${commit:0:12}"
cp -p "$built" "$exe"
sha256sum "$exe" | tee -a "$logs/provenance.txt"

overall=0
for variant in control candidate; do
  run_dir="$job_root/$variant"
  mkdir -p "$run_dir"
  cp "$snapshot/${variant}.nml" "$run_dir/run.nml"

  if ! "$code/tests/preflight.sh" --exe "$exe" --commit "$commit" \
      --nml "$run_dir/run.nml" >"$logs/preflight_${variant}.log" 2>&1; then
    echo "$variant: preflight failed" | tee -a "$logs/provenance.txt"
    overall=1
    continue
  fi

  (
    cd "$run_dir" || exit 2
    srun --kill-on-bad-exit=1 --mpi=pmi2 -n 32 -c 1 \
      --cpu-bind=none "$rank_bind" "$exe" run.nml >run.log 2>&1
  )
  rc=$?
  outputs=$(find "$run_dir" -mindepth 1 -maxdepth 1 -type d \
    -name 'output_*' | wc -l)
  completed=$(grep -c 'Run completed' "$run_dir/run.log" 2>/dev/null || true)
  grows=$(grep -c '\[RESIZABLE\] GRID_GROW' "$run_dir/run.log" \
    2>/dev/null || true)
  lb_done=$(grep -c 'load_balance total:' "$run_dir/run.log" \
    2>/dev/null || true)
  echo "$variant: rc=$rc outputs=$outputs completed=$completed grows=$grows lb_done=$lb_done" | \
    tee -a "$logs/provenance.txt"
  [ "$rc" -eq 0 ] || overall=1
done

if ! "$snapshot_gate" "$job_root/control" \
    "$job_root/candidate" | tee "$logs/gate.log"; then
  overall=1
fi

echo "finished=$(date --iso-8601=seconds) rc=$overall" | \
  tee -a "$logs/provenance.txt"
if [ "$overall" -eq 0 ]; then
  touch "$job_root/PASS"
else
  touch "$job_root/FAIL"
fi
exit "$overall"
