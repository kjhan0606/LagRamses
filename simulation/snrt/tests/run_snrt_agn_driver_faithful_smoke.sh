#!/usr/bin/env bash
# F-P2.8: execute the production amr_step -> SNRT AGN source path.
set -euo pipefail

ROOT=/gpfs/kjhan/LRD_JWST
SCRIPT_DIR=$ROOT/simulation/snrt/tests
CONFIG=$ROOT/simulation/snrt/config/snrt_agn_driver_faithful_smoke.nml
GROUP_CONTRACT=$ROOT/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml
SECONDARY_CONTRACT=$ROOT/simulation/snrt/config/snrt_secondary_table_contract_v1.nml
BINARY=${SNRT_BINARY:-$ROOT/bin/ramses_final3d}
STAMP=$(date +%Y%m%dT%H%M%S)
RUN_ROOT=$ROOT/simulation/snrt/runs/fp2_8_agn_driver_faithful_smoke/job_${STAMP}_$$

[[ -x "$BINARY" ]] || { echo "F-P2.8 FAIL: missing executable $BINARY" >&2; exit 1; }
[[ -f "$CONFIG" ]] || { echo "F-P2.8 FAIL: missing namelist $CONFIG" >&2; exit 1; }
[[ -f "$GROUP_CONTRACT" ]] || { echo "F-P2.8 FAIL: missing group contract" >&2; exit 1; }
[[ -f "$SECONDARY_CONTRACT" ]] || { echo "F-P2.8 FAIL: missing secondary contract" >&2; exit 1; }
[[ ! -e "$RUN_ROOT" ]] || { echo "F-P2.8 FAIL: refusing existing run root" >&2; exit 1; }
mkdir -p "$RUN_ROOT/baseline" "$RUN_ROOT/injected"

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false
export OMP_WAIT_POLICY=passive
export CUDA_DEVICE_MAX_CONNECTIONS=1
export I_MPI_FABRICS=shm
export I_MPI_PIN=0

fatal_log_pattern='forrtl: severe|SIGSEGV|AGN reference event arrays have inconsistent size|Invalid AGN reference|Invalid AGN feedback budget|SNRT RT preflight failed closed|Problem in merge_sink'

run_case() {
  local name=$1
  local case_dir=$RUN_ROOT/$name
  local log=$case_dir/ramses.log
  local status

  cp -- "$CONFIG" "$case_dir/effective.nml"
  sha256sum "$BINARY" > "$case_dir/binary.sha256"
  sha256sum "$CONFIG" "$GROUP_CONTRACT" "$SECONDARY_CONTRACT" > "$case_dir/contracts.sha256"

  set +e
  if [[ "$name" == baseline ]]; then
    (
      cd -- "$case_dir"
      env -u SNRT_RT_TX_TEST_FAIL_STAGE -u SNRT_RT_TX_TEST_FAIL_LEAF \
        SNRT_RT_TX_DIAGNOSTIC_MODE=1 \
        SNRT_RT_ENABLE=1 SNRT_AGN_MODEL=partition_reference_v1 \
        SNRT_DRIVER_TEST_SEED_SOURCE=1 \
        SNRT_REDUCED_C=0.01 SNRT_RT_LEVEL=3 \
        SNRT_GROUP_CONTRACT="$GROUP_CONTRACT" \
        SNRT_ALLOW_REFERENCE_CONTROL=1 \
        SNRT_SECONDARY_TABLE_CONTRACT="$SECONDARY_CONTRACT" \
        timeout 300s mpirun -np 1 "$BINARY" effective.nml
    ) > "$log" 2>&1
  else
    (
      cd -- "$case_dir"
      env SNRT_RT_ENABLE=1 SNRT_AGN_MODEL=partition_reference_v1 \
        SNRT_DRIVER_TEST_SEED_SOURCE=1 \
        SNRT_REDUCED_C=0.01 SNRT_RT_LEVEL=3 \
        SNRT_GROUP_CONTRACT="$GROUP_CONTRACT" \
        SNRT_ALLOW_REFERENCE_CONTROL=1 \
        SNRT_SECONDARY_TABLE_CONTRACT="$SECONDARY_CONTRACT" \
        SNRT_RT_TX_DIAGNOSTIC_MODE=1 \
        SNRT_RT_TX_TEST_FAIL_STAGE=receiver SNRT_RT_TX_TEST_FAIL_LEAF=1 \
        timeout 300s mpirun -np 1 "$BINARY" effective.nml
    ) > "$log" 2>&1
  fi
  status=$?
  set -e

  if grep -Eq "$fatal_log_pattern" "$log"; then
    echo "F-P2.8 FAIL $name contains an unexpected fatal/error marker" >&2
    tail -n 120 "$log" >&2
    return 1
  fi

  if [[ "$name" == baseline ]]; then
    [[ "$status" -eq 0 ]] || { echo "F-P2.8 FAIL baseline return=$status" >&2; tail -n 120 "$log" >&2; return 1; }
    grep -Eq 'active sources:[[:space:]]+[1-9][0-9]*' "$log" || {
      echo "F-P2.8 FAIL baseline did not reach an active AGN source" >&2; tail -n 120 "$log" >&2; return 1;
    }
    grep -Fq 'SNRT_RT_TRANSACTION_COMMIT_PASS' "$log" || {
      echo "F-P2.8 FAIL baseline commit marker missing" >&2; tail -n 120 "$log" >&2; return 1;
    }
    grep -Fq 'SNRT_RT_CLOSURE_PASS' "$log" || {
      echo "F-P2.8 FAIL baseline closure marker missing" >&2; tail -n 120 "$log" >&2; return 1;
    }
    printf 'F-P2.8_CASE baseline PASS return_code=%s\n' "$status" | tee "$case_dir/status.txt"
  else
    # clean_stop may terminate through MPI_ABORT or return normally in a
    # diagnostic build; the rollback markers are authoritative here.
    grep -Fq 'SNRT_RT_DIAGNOSTIC_FAIL_CLOSED class=receiver' "$log" || {
      echo "F-P2.8 FAIL injected diagnostic marker missing" >&2; tail -n 120 "$log" >&2; return 1;
    }
    grep -Fq 'SNRT RT transaction rollback: class=receiver' "$log" || {
      echo "F-P2.8 FAIL injected rollback marker missing" >&2; tail -n 120 "$log" >&2; return 1;
    }
    grep -Fq 'SNRT_DRIVER_TEST_SEED_SOURCE applied: NONPRODUCTION' "$log" || {
      echo "F-P2.8 FAIL injected source seed marker missing" >&2; tail -n 120 "$log" >&2; return 1;
    }
    if grep -Fq 'SNRT_RT_TRANSACTION_COMMIT_PASS' "$log"; then
      echo "F-P2.8 FAIL injected case committed after forced failure" >&2
      return 1
    fi
    printf 'F-P2.8_CASE injected PASS return_code=%s expected_mpi_abort=1\n' "$status" | tee "$case_dir/status.txt"
  fi

  if find "$case_dir" -maxdepth 1 -type d -name 'output_*' -print -quit | grep -q .; then
    echo "F-P2.8 FAIL unexpected output directory in $case_dir" >&2
    return 1
  fi
}

{
  printf 'F-P2.8_BEGIN stamp=%s\n' "$STAMP"
  printf 'binary=%s\n' "$BINARY"
  printf 'config=%s\n' "$CONFIG"
  printf 'run_root=%s\n' "$RUN_ROOT"
  printf 'purpose=driver-faithful production call, explicit nonproduction source seed, RT commit/rollback\n'
} > "$RUN_ROOT/launch_record.txt"

run_case baseline
run_case injected

printf 'F-P2.8_PASS run_root=%s\n' "$RUN_ROOT" | tee "$RUN_ROOT/summary.txt"
