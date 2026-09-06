#!/usr/bin/env bash
#SBATCH --job-name=snrt-fp27-init
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:2
#SBATCH --gpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/runs/fp2_7_initialized_ramses_smoke/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/runs/fp2_7_initialized_ramses_smoke/slurm_%j.err

set -euo pipefail

# Slurm executes a submitted batch script from a spool copy, so BASH_SOURCE is
# not a reliable project-relative path. Keep the qualification root explicit
# and immutable for this /gpfs test case.
ROOT=/gpfs/kjhan/LRD_JWST
SCRIPT_DIR=/gpfs/kjhan/LRD_JWST/simulation/snrt/runs/fp2_7_initialized_ramses_smoke
BINARY="$ROOT/bin/ramses_final3d"
TEMPLATE="$SCRIPT_DIR/effective.nml"
GROUP_CONTRACT="$ROOT/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml"
SECONDARY_CONTRACT="$ROOT/simulation/snrt/config/snrt_secondary_table_contract_v1.nml"
JOB_ROOT="$SCRIPT_DIR/job_${SLURM_JOB_ID}"

[[ -x "$BINARY" ]] || { echo "D4_FAIL missing binary=$BINARY" >&2; exit 1; }
[[ -f "$TEMPLATE" ]] || { echo "D4_FAIL missing namelist=$TEMPLATE" >&2; exit 1; }
[[ -f "$GROUP_CONTRACT" ]] || { echo "D4_FAIL missing group contract=$GROUP_CONTRACT" >&2; exit 1; }
[[ -f "$SECONDARY_CONTRACT" ]] || { echo "D4_FAIL missing secondary contract=$SECONDARY_CONTRACT" >&2; exit 1; }
[[ ! -e "$JOB_ROOT" ]] || { echo "D4_FAIL refusing to overwrite run_dir=$JOB_ROOT" >&2; exit 1; }
mkdir -p "$JOB_ROOT/baseline" "$JOB_ROOT/injected"

export OMP_NUM_THREADS=1
export CUDA_DEVICE_MAX_CONNECTIONS=1
export I_MPI_FABRICS=shm

write_metadata() {
  local case_dir="$1"
  sha256sum "$BINARY" > "$case_dir/binary.sha256"
  sha256sum "$TEMPLATE" > "$case_dir/template_namelist.sha256"
  sha256sum "$case_dir/effective.nml" > "$case_dir/effective_namelist.sha256"
  {
    printf 'purpose=F-P2.7 D4 initialized RAMSES SNRT integration smoke\n'
    printf 'case_dir=%s\n' "$case_dir"
    printf 'effective_namelist=%s\n' "$case_dir/effective.nml"
    printf 'binary=%s\n' "$BINARY"
    printf 'mpi_ranks=%s\n' "${SLURM_NTASKS:-unknown}"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-unknown}"
    printf 'slurm_partition=%s\n' "${SLURM_JOB_PARTITION:-unknown}"
    printf 'cuda_visible_devices=%s\n' "${CUDA_VISIBLE_DEVICES:-unset}"
    printf 'output_policy=noutput=1 aout=2.0 tout=1.0d30 foutput=1000000 fbackup=1000000\n'
    mpiifx --version > "$case_dir/compiler.version" 2>&1 || true
    mpirun --version > "$case_dir/launcher.version" 2>&1 || true
    printf 'compiler=%s\n' "$(sed -n '1p' "$case_dir/compiler.version")"
    printf 'launcher=%s\n' "$(sed -n '1p' "$case_dir/launcher.version")"
  } > "$case_dir/run_metadata.txt"
}

assert_no_outputs() {
  local case_dir="$1"
  if find "$case_dir" -maxdepth 1 -type d -name 'output_*' -print -quit | grep -q .; then
    echo "D4_FAIL unexpected RAMSES output directory in $case_dir" >&2
    return 1
  fi
}

run_case() {
  local case_name="$1"
  local case_dir="$JOB_ROOT/$case_name"
  local log="$case_dir/ramses.log"
  local status

  cp -- "$TEMPLATE" "$case_dir/effective.nml"
  write_metadata "$case_dir"
  if [[ "$case_name" == baseline ]]; then
    set +e
    (
      cd -- "$case_dir"
      env -u SNRT_RT_TX_DIAGNOSTIC_MODE \
        -u SNRT_RT_TX_TEST_FAIL_STAGE -u SNRT_RT_TX_TEST_FAIL_LEAF \
        SNRT_RT_ENABLE=1 SNRT_REDUCED_C=0.01 SNRT_RT_LEVEL=3 \
        SNRT_GROUP_CONTRACT="$GROUP_CONTRACT" \
        SNRT_ALLOW_REFERENCE_CONTROL=1 \
        SNRT_SECONDARY_TABLE_CONTRACT="$SECONDARY_CONTRACT" \
        timeout 240s mpirun -np "${SLURM_NTASKS:-2}" "$BINARY" \
        effective.nml
    ) > "$log" 2>&1
    status=$?
    set -e
    [[ "$status" -eq 0 ]] || { echo "D4_FAIL baseline return_code=$status" >&2; tail -n 100 "$log" >&2; return 1; }
    grep -Fq 'SNRT_RT_TRANSACTION_COMMIT_PASS' "$log" || {
      echo "D4_FAIL baseline commit marker missing" >&2; tail -n 100 "$log" >&2; return 1;
    }
    grep -Fq 'SNRT_RT_CLOSURE_PASS' "$log" || {
      echo "D4_FAIL baseline closure marker missing" >&2; tail -n 100 "$log" >&2; return 1;
    }
    printf 'D4_CASE baseline status=PASS return_code=%s\n' "$status" | tee "$case_dir/status.txt"
  else
    set +e
    (
      cd -- "$case_dir"
      env SNRT_RT_ENABLE=1 SNRT_REDUCED_C=0.01 SNRT_RT_LEVEL=3 \
        SNRT_GROUP_CONTRACT="$GROUP_CONTRACT" \
        SNRT_ALLOW_REFERENCE_CONTROL=1 \
        SNRT_SECONDARY_TABLE_CONTRACT="$SECONDARY_CONTRACT" \
        SNRT_RT_TX_DIAGNOSTIC_MODE=1 \
        SNRT_RT_TX_TEST_FAIL_STAGE=receiver \
        SNRT_RT_TX_TEST_FAIL_LEAF=1 \
        timeout 240s mpirun -np "${SLURM_NTASKS:-2}" "$BINARY" \
        effective.nml
    ) > "$log" 2>&1
    status=$?
    set -e
    # Diagnostic injection is fail-closed inside the runtime and may return
    # normally after recording the rollback.  The markers, not a nonzero
    # launcher status, are the authoritative injected-case result.
    [[ "$status" -eq 0 ]] || { echo "D4_FAIL injected process return_code=$status" >&2; tail -n 100 "$log" >&2; return 1; }
    grep -Fq 'SNRT_RT_DIAGNOSTIC_FAIL_CLOSED class=receiver' "$log" || {
      echo "D4_FAIL injected rollback marker missing" >&2; tail -n 100 "$log" >&2; return 1;
    }
    grep -Fq 'SNRT RT transaction rollback: class=receiver' "$log" || {
      echo "D4_FAIL injected transaction rollback marker missing" >&2; tail -n 100 "$log" >&2; return 1;
    }
    printf 'D4_CASE injected status=PASS return_code=%s expected_fail_closed_markers=1\n' "$status" | tee "$case_dir/status.txt"
  fi
  assert_no_outputs "$case_dir"
}

{
  printf 'D4_BEGIN job_id=%s partition=%s ntasks=%s gres=gpu:2\n' \
    "${SLURM_JOB_ID:-unknown}" "${SLURM_JOB_PARTITION:-unknown}" "${SLURM_NTASKS:-unknown}"
  printf 'D4_BINARY %s\n' "$BINARY"
  printf 'D4_EFFECTIVE_NAMELIST_TEMPLATE %s\n' "$TEMPLATE"
  printf 'D4_RUN_ROOT %s\n' "$JOB_ROOT"
} > "$JOB_ROOT/launch_record.txt"

run_case baseline
run_case injected

printf 'D4_PASS job_id=%s run_root=%s\n' "${SLURM_JOB_ID:-unknown}" "$JOB_ROOT" | tee "$JOB_ROOT/summary.txt"
