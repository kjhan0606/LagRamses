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
export OMP_WAIT_POLICY=passive
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_FABRICS=shm

commit=$(git -C "$code" rev-parse HEAD)
branch=$(git -C "$code" branch --show-current)
{
  echo "job_id=$SLURM_JOB_ID"
  echo "host=$(hostname)"
  echo "workdir=$code"
  echo "branch=$branch"
  echo "commit=$commit"
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
    srun --kill-on-bad-exit=1 --mpi=pmi2 -n 32 -c 1 --cpu-bind=cores \
      "$exe" run.nml >run.log 2>&1
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
