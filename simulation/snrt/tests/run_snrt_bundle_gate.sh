#!/usr/bin/env bash
set -euo pipefail

# One bundle-level native gate for the SNRT production path. Focused runners
# remain useful for debugging; this is the single project gate that owns the
# SNRT-enabled full link, native controls, negative controls, timing, and one
# concise summary.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SUMMARY="${SNRT_BUNDLE_GATE_SUMMARY:-$ROOT/simulation/snrt/build/snrt_bundle_gate_summary.txt}"
RUN_DIR="$(mktemp -d /gpfs/kjhan/LRD_JWST/.snrt-bundle-gate.XXXXXX)"
trap 'rm -rf "$RUN_DIR"' EXIT
# The focused native runners intentionally use a private module directory, but
# Fortran also searches the process working directory for .mod files.  Enter
# the repository before invoking them so a make -C bin caller cannot make a
# stale/corrupt bin/*.mod shadow the freshly generated module set.
cd "$ROOT"
mkdir -p "$(dirname "$SUMMARY")"
: > "$SUMMARY"

fail_gate() {
  echo "SNRT_BUNDLE_GATE_FAIL reason=$*" | tee -a "$SUMMARY" >&2
  exit 1
}

require_command() {
  local command_name="$1"
  command -v "$command_name" >/dev/null 2>&1 || fail_gate "missing_command=$command_name"
}

time_now_ns() {
  date +%s%N
}

elapsed_seconds() {
  awk -v start="$1" -v end="$2" 'BEGIN { printf "%.3f", (end-start)/1e9 }'
}

run_stage() {
  local stage="$1"
  shift
  local log="$RUN_DIR/${stage}.log"
  local start_ns end_ns elapsed
  start_ns="$(time_now_ns)"
  if "$@" >"$log" 2>&1; then
    end_ns="$(time_now_ns)"
    elapsed="$(elapsed_seconds "$start_ns" "$end_ns")"
    printf 'STAGE %s status=PASS elapsed_s=%s\n' "$stage" "$elapsed" | tee -a "$SUMMARY"
  else
    end_ns="$(time_now_ns)"
    elapsed="$(elapsed_seconds "$start_ns" "$end_ns")"
    printf 'STAGE %s status=FAIL elapsed_s=%s\n' "$stage" "$elapsed" | tee -a "$SUMMARY" >&2
    tail -n 80 "$log" >&2 || true
    exit 1
  fi
  LAST_STAGE_LOG="$log"
}

require_command mpiifx
require_command gfortran
require_command mpirun
require_command nm

CUDA_ROOT="${CUDA_ROOT:-/opt/ohpc/pub/cuda/13.0.2}"
NVCC_BIN="${NVCC:-$CUDA_ROOT/bin/nvcc}"
[[ -x "$NVCC_BIN" ]] || fail_gate "missing_command=nvcc path=$NVCC_BIN"

printf 'SNRT_BUNDLE_GATE_BEGIN root=%s\n' "$ROOT" | tee -a "$SUMMARY"
printf 'SNRT_BUNDLE_GATE_COMMIT %s\n' "$(git -C "$ROOT" rev-parse HEAD)" | tee -a "$SUMMARY"
printf 'SNRT_BUNDLE_GATE_TOOLS mpiifx=%s gfortran=%s mpirun=%s nvcc=%s\n' \
  "$(command -v mpiifx)" "$(command -v gfortran)" "$(command -v mpirun)" "$NVCC_BIN" | tee -a "$SUMMARY"

run_stage production_build \
  make -C "$ROOT/bin" -B EXEC="$RUN_DIR/ramses_bundle" SNRT=1 USE_CUDA=1 ramses

BINARY="$RUN_DIR/ramses_bundle3d"
[[ -x "$BINARY" ]] || fail_gate "production_binary_missing=$BINARY"

run_stage agn_partition_reference \
  "$ROOT/simulation/snrt/tests/run_fp15_agn_efficiency.sh"
grep -q 'AGN_REFERENCE_PARTITION_SMOKE_OK' "$LAST_STAGE_LOG" || \
  fail_gate 'agn_reference_partition_marker_missing'
printf 'CONSERVATION_ASSERTION name=agn_release_partition source=agn_feedback_deposition_smoke.f90 threshold=1e-14\n' | tee -a "$SUMMARY"

NATIVE_SYMBOLS="$RUN_DIR/native_symbols.txt"
nm -g --defined-only "$BINARY" > "$NATIVE_SYMBOLS" || fail_gate 'native_symbol_dump_failed'
for symbol in \
  snrt_ramses_driver_mp_snrt_ramses_advance_level_ \
  snrt_transport_step_mp_snrt_transport_absorb_multigroup_prepared_trial_ \
  snrt_transport_step_mp_snrt_transport_absorb_multigroup_prepared_dust_trial_ \
  snrt_rt_transaction_mp_snrt_transaction_reduce_decision_ \
  snrt_thermochemistry_mp_snrt_thermochemistry_advance_cell_
do
  grep -Fq -- "$symbol" "$NATIVE_SYMBOLS" || fail_gate "missing_native_symbol=$symbol"
done
printf 'NATIVE_SYMBOLS_CHECK count=5 status=PASS\n' | tee -a "$SUMMARY"

run_stage dust_ledger_receiver \
  "$ROOT/simulation/snrt/tests/run_snrt_native_dust_transaction.sh"
grep -q 'SNRT_NATIVE_DUST_TRANSACTION_ALL_OK' "$LAST_STAGE_LOG" || \
  fail_gate 'dust_ledger_receiver_marker_missing'
printf 'CONSERVATION_ASSERTION name=dust_ledger_receiver source=snrt_dust_transaction.f90 threshold=64_fp32_epsilon\n' | tee -a "$SUMMARY"

run_stage thermochemistry \
  "$ROOT/simulation/snrt/tests/run_snrt_native_thermochemistry.sh"
thermo_negative_count="$(grep -Ec 'PASS: (unset|missing|malformed|wrong).*rejected' "$LAST_STAGE_LOG" || true)"
[[ "$thermo_negative_count" -ge 4 ]] || fail_gate "thermochemistry_negative_cases=$thermo_negative_count"
printf 'NEGATIVE_CASES family=thermochemistry_loader count=%s\n' "$thermo_negative_count" | tee -a "$SUMMARY"

run_stage spectral_contract \
  "$ROOT/simulation/snrt/tests/run_snrt_native_spectral_contract.sh"
spectral_loader_count="$(grep -c 'SNRT_SPECTRAL_CONTRACT_LOADER_OK' "$LAST_STAGE_LOG" || true)"
[[ "$spectral_loader_count" -ge 10 ]] || fail_gate "spectral_loader_cases=$spectral_loader_count"
printf 'NEGATIVE_CASES family=spectral_contract_loader count=%s\n' "$spectral_loader_count" | tee -a "$SUMMARY"

run_stage transaction_mpi \
  "$ROOT/simulation/snrt/tests/run_snrt_native_rt_transaction.sh"
grep -q 'SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_PASS ranks=2' "$LAST_STAGE_LOG" || \
  fail_gate 'mpi_zero_leaf_two_rank_evidence_missing'
transaction_rejection_count="$(grep -Ec 'SNRT_NATIVE_RT_TRANSACTION_(PARTITION|CHEMISTRY)_ROLLBACK_PASS|SNRT_NATIVE_RT_TRANSACTION_MAX_ITER_LIMIT_REJECT_PASS' "$LAST_STAGE_LOG" || true)"
[[ "$transaction_rejection_count" -ge 3 ]] || fail_gate "transaction_rejection_cases=$transaction_rejection_count"
grep -q 'STATIC_SUPPORTING_CHECK driver_failure_routes' "$LAST_STAGE_LOG" || \
  fail_gate 'static_failure_route_check_missing'
grep -q 'STATIC_SUPPORTING_CHECK driver_hydro_preflight' "$LAST_STAGE_LOG" || \
  fail_gate 'static_hydro_preflight_check_missing'
printf 'NEGATIVE_CASES family=transaction_rollback_and_config count=%s\n' "$transaction_rejection_count" | tee -a "$SUMMARY"
printf 'MPI_COVERAGE required_ranks=2 marker=SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_PASS\n' | tee -a "$SUMMARY"
printf 'STATIC_SUPPORTING_CHECK driver_failure_routes\nSTATIC_SUPPORTING_CHECK driver_hydro_preflight\n' | tee -a "$SUMMARY"

run_stage cuda_multigroup \
  env CUDA_ROOT="$CUDA_ROOT" NVCC="$NVCC_BIN" \
  "$ROOT/simulation/snrt/tests/run_snrt_native_cuda_multigroup.sh"
grep -q 'SNRT_CUDA_MULTIGROUP_OK' "$LAST_STAGE_LOG" || fail_gate 'cuda_photon_budget_marker_missing'
grep -q 'SNRT_CUDA_MULTIGROUP_SPECIES_MIX_OK' "$LAST_STAGE_LOG" || fail_gate 'cuda_species_budget_marker_missing'
grep -q 'SNRT_CUDA_MULTIGROUP_SPECIES_DUST_ZERO_DUST_BITWISE_OK' "$LAST_STAGE_LOG" || \
  fail_gate 'cuda_zero_dust_equivalence_marker_missing'
printf 'CONSERVATION_ASSERTION name=cuda_photon_budget source=snrt_cuda_multigroup_smoke.f90 threshold=2e-6\n' | tee -a "$SUMMARY"

if rg -q 'snrt_partition_absorption' "$ROOT/patch/lagRamses/snrt_ramses_driver.f90"; then
  fail_gate 'host_partition_reintroduced'
fi
grep -q 'ZERO_SCAFFOLD' "$ROOT/patch/lagRamses/snrt_ramses_driver.f90" || \
  fail_gate 'dust_scaffold_mode_missing'
printf 'STATIC_SUPPORTING_CHECK dust_direct_hhe_handoff\nSTATIC_SUPPORTING_CHECK dust_zero_scaffold\n' | tee -a "$SUMMARY"

run_stage production_negative \
  "$ROOT/simulation/snrt/tests/run_p04_production_negative.sh"
grep -q 'P04_PRODUCTION_NEGATIVE_OK' "$LAST_STAGE_LOG" || fail_gate 'production_negative_marker_missing'
printf 'NEGATIVE_CASES family=production_fail_closed count=1_plus_snia_runtime\n' | tee -a "$SUMMARY"

run_stage diff_check git -C "$ROOT" diff --check
printf 'CONSERVATION_ASSERTION name=spectral_group_sum_vs_Lbol source=snrt_spectral_contract_smoke.f90 threshold=1e-13\n' | tee -a "$SUMMARY"

BINARY_SHA256="$(sha256sum "$BINARY" | awk '{print $1}')"
printf 'BINARY_SHA256 %s %s\n' "$BINARY_SHA256" "$BINARY" | tee -a "$SUMMARY"
for source_file in \
  patch/lagRamses/agn_feedback_deposition.f90 \
  patch/lagRamses/sink_particle.kjhan.f90 \
  patch/lagRamses/pm_commons.f90 \
  patch/lagRamses/init_sink.f90 \
  patch/lagRamses/read_params.jaehyun.f90 \
  patch/lagRamses/snrt_agn_efficiency.f90
do
  [[ -f "$ROOT/$source_file" ]] || fail_gate "missing_changed_source=$source_file"
  printf 'SOURCE_SHA256 %s %s\n' "$source_file" \
    "$(sha256sum "$ROOT/$source_file" | awk '{print $1}')" | tee -a "$SUMMARY"
done
printf 'SOURCE_SHA256 snrt_ramses_driver %s\n' \
  "$(sha256sum "$ROOT/patch/lagRamses/snrt_ramses_driver.f90" | awk '{print $1}')" | tee -a "$SUMMARY"
printf 'SOURCE_SHA256 snrt_rt_transaction %s\n' \
  "$(sha256sum "$ROOT/patch/lagRamses/snrt_rt_transaction.f90" | awk '{print $1}')" | tee -a "$SUMMARY"
printf 'SOURCE_SHA256 snrt_transport_step %s\n' \
  "$(sha256sum "$ROOT/patch/lagRamses/snrt_transport_step.f90" | awk '{print $1}')" | tee -a "$SUMMARY"
printf 'SOURCE_SHA256 snrt_thermochemistry %s\n' \
  "$(sha256sum "$ROOT/patch/lagRamses/snrt_thermochemistry.f90" | awk '{print $1}')" | tee -a "$SUMMARY"
printf 'SOURCE_SHA256 snrt_cuda_kernels %s\n' \
  "$(sha256sum "$ROOT/patch/lagRamses/snrt_cuda_kernels.cu" | awk '{print $1}')" | tee -a "$SUMMARY"
printf 'SNRT_BUNDLE_GATE_PASS\n' | tee -a "$SUMMARY"
