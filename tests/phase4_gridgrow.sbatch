#!/bin/bash
# [RESIZABLE] Phase 4 collective grid-growth validation allocation.
#SBATCH --job-name=rz_p4_gridgrow
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
#SBATCH --output=/home/kjhan/BACKUP/lagramses_qa/gridgrow_jobs/gridgrow_%j.slurm

set -uo pipefail

code=/home/kjhan/BACKUP/lagRamses-RESIZABLE/code
qa=/home/kjhan/BACKUP/lagramses_qa
pbh_table=/home/kjhan/BACKUP/lagRamses-PBH/papers/pbh_working/tables/pbh_evap_table_ast_s1.dat
pbh_table_sha=204888b8d5154785b837db280f9b06d2d55caf9a0583fbc0877292da521e906c
job_root="$qa/gridgrow_jobs/${SLURM_JOB_ID:?}"
logs="$job_root/logs"
mkdir -p "$logs"

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:/opt/ohpc/pub/intel/oneapi/compiler/2025.3/bin:$PATH
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=false
export OMP_WAIT_POLICY=passive
# [RESIZABLE] This site has task/cgroup but not task/affinity.  Slurm owns the
# allocation boundary; slurm_rank_bind.sh is the sole per-rank authority.
# Intel MPI pinning combined with srun collapsed all 32 ranks onto CPU 0 in
# invalid job 321015, so it must stay disabled.
export I_MPI_PIN=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_LIST I_MPI_PIN_RESPECT_CPUSET
unset KMP_AFFINITY KMP_HW_SUBSET GOMP_CPU_AFFINITY
export I_MPI_FABRICS=shm

cpuset_expander="$code/tests/expand_cpuset.sh"
rank_bind="$code/tests/slurm_rank_bind.sh"
affinity_validator="$code/tests/validate_affinity_map.sh"
if [ ! -x /usr/bin/taskset ] || [ ! -x "$cpuset_expander" ] || \
   [ ! -x "$rank_bind" ] || [ ! -x "$affinity_validator" ]; then
  echo 'PHASE4: taskset or CPU-binding helper is not executable'
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
  echo 'PHASE4: could not read the batch allocation CPU mask'
  exit 1
fi
if phase4_cpu_list=$("$cpuset_expander" "$allocation_cpu_mask" 32); then
  export PHASE4_CPU_LIST=$phase4_cpu_list
else
  parse_rc=$?
  echo "PHASE4: allocation mask is not 32 canonical unique CPUs: $allocation_cpu_mask (rc=$parse_rc)"
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
  echo "started=$(date --iso-8601=seconds)"
} | tee "$logs/provenance.txt"

# [RESIZABLE] The split PBH project is a read-only input dependency.  Pin the
# exact clean tracked table so a moved symlink or regenerated table cannot
# silently change this allocator-equivalence test.
if [ ! -r "$pbh_table" ]; then
  echo "PHASE4: PBH table is not readable: $pbh_table" | \
    tee -a "$logs/provenance.txt"
  exit 1
fi
actual_pbh_sha=$(sha256sum "$pbh_table" | awk '{print $1}')
echo "pbh_table=$pbh_table sha256=$actual_pbh_sha" | tee -a "$logs/provenance.txt"
if [ "$actual_pbh_sha" != "$pbh_table_sha" ]; then
  echo "PHASE4: PBH table checksum mismatch" | tee -a "$logs/provenance.txt"
  exit 1
fi
for nml in "$code/tests/phase4_roomy.nml" "$code/tests/phase4_tight.nml"; do
  if ! grep -Fq "pbh_table_file='$pbh_table'" "$nml"; then
    echo "PHASE4: $nml does not reference the pinned PBH table" | \
      tee -a "$logs/provenance.txt"
    exit 1
  fi
done

# [RESIZABLE] Check affinity after MPI_Init, where Intel MPI previously
# narrowed every task mask to CPU 0 despite a valid 32-core Slurm allocation.
affinity_dir="$job_root/affinity"
affinity_probe="$job_root/mpi_affinity_probe"
mkdir -p "$affinity_dir"
if ! mpiifx -O2 "$code/tests/mpi_affinity_probe.f90" -o "$affinity_probe" \
    >"$logs/affinity_build.log" 2>&1; then
  echo 'PHASE4: affinity probe build failed' | tee -a "$logs/provenance.txt"
  exit 1
fi
if ! srun --kill-on-bad-exit=1 --mpi=pmi2 -n 32 -c 1 \
    --cpu-bind=none "$rank_bind" "$affinity_probe" "$affinity_dir" \
    >"$logs/affinity_run.log" 2>&1; then
  echo 'PHASE4: affinity probe execution failed' | tee -a "$logs/provenance.txt"
  exit 1
fi
affinity_file_count=$(find "$affinity_dir" -maxdepth 1 -type f \
  -name 'rank_*.txt' | wc -l)
if [ "$affinity_file_count" -ne 32 ]; then
  echo "PHASE4: affinity probe wrote $affinity_file_count files, expected 32" | \
    tee -a "$logs/provenance.txt"
  exit 1
fi
find "$affinity_dir" -maxdepth 1 -type f -name 'rank_*.txt' -print0 | \
  sort -z | xargs -0 cat | sort -n >"$logs/affinity_map.txt"
if ! "$affinity_validator" "$PHASE4_CPU_LIST" \
    "$logs/affinity_map.txt" 32; then
  echo 'PHASE4: rank/localid affinity does not match the expected CPU map' | \
    tee -a "$logs/provenance.txt"
  cat "$logs/affinity_map.txt" | tee -a "$logs/provenance.txt"
  exit 1
fi
echo 'affinity=PASS ranks=32 localids=32 expected_unique_single_cpu_masks=32' | \
  tee -a "$logs/provenance.txt"
cat "$logs/affinity_map.txt" | tee -a "$logs/provenance.txt"

syntax_files=(
  patch/cuRamses/timing.f90
  patch/cuRamses/amr_commons.kjhan.f90
  patch/cuRamses/radiation_commons.f90
  patch/cuRamses/refine_utils.f90
  patch/lagRamses/pm_commons.f90
  patch/lagRamses/read_params.jaehyun.f90
)
if ! "$code/tests/syntax_check.sh" "$code" "${syntax_files[@]}" >"$logs/syntax.log" 2>&1; then
  echo 'PHASE4: syntax gate failed' | tee -a "$logs/provenance.txt"
  exit 1
fi

if ! make -C "$code/bin" HDF5=1 USE_FFTW=1 >"$logs/build.log" 2>&1; then
  echo 'PHASE4: build failed' | tee -a "$logs/provenance.txt"
  exit 1
fi

built="$code/bin/ramses_final3d"
[ -x "$built" ] || { echo 'PHASE4: missing built binary'; exit 1; }
binary_epoch=$(stat -c %Y "$built")
commit_epoch=$(git -C "$code" log -1 --format=%ct "$commit")
[ "$binary_epoch" -ge "$commit_epoch" ] || { echo 'PHASE4: binary predates commit'; exit 1; }

exe="$job_root/ramses_phase4_${commit:0:12}"
cp -p "$built" "$exe"
sha256sum "$exe" | tee -a "$logs/provenance.txt"
sha256sum "$code/tests/phase4_roomy.nml" "$code/tests/phase4_tight.nml" | \
  tee -a "$logs/provenance.txt"

overall=0
for variant in roomy tight; do
  run_dir="$job_root/$variant"
  mkdir -p "$run_dir"
  if [ "$variant" = roomy ]; then
    cp "$code/tests/phase4_roomy.nml" "$run_dir/run.nml"
  else
    cp "$code/tests/phase4_tight.nml" "$run_dir/run.nml"
  fi

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
  outputs=$(find "$run_dir" -mindepth 1 -maxdepth 1 -type d -name 'output_*' | wc -l)
  completed=$(grep -c 'Run completed' "$run_dir/run.log" 2>/dev/null || true)
  grows=$(grep -c '\[RESIZABLE\] GRID_GROW' "$run_dir/run.log" 2>/dev/null || true)
  echo "$variant: rc=$rc outputs=$outputs completed=$completed grows=$grows" | \
    tee -a "$logs/provenance.txt"
  [ "$rc" -eq 0 ] || overall=1
done

if ! "$code/tests/phase4_gridgrow_gate.sh" "$job_root/roomy" "$job_root/tight" \
    | tee "$logs/gate.log"; then
  overall=1
fi

echo "finished=$(date --iso-8601=seconds) rc=$overall" | tee -a "$logs/provenance.txt"
if [ "$overall" -eq 0 ]; then
  touch "$job_root/PASS"
else
  touch "$job_root/FAIL"
fi
exit "$overall"
