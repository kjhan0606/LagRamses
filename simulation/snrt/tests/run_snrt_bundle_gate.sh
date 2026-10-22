#!/usr/bin/env bash
set -euo pipefail

# One bundle-level native gate for the SNRT production path. Focused runners
# remain useful for debugging; this is the single project gate that owns the
# SNRT-enabled full link, native controls, negative controls, timing, and one
# concise summary.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo 'SNRT_BUNDLE_GATE_FAIL reason=requires_Slurm_allocation; do_not_run_MPI_CUDA_gate_on_login_node' >&2
  exit 2
fi
SUMMARY="${SNRT_BUNDLE_GATE_SUMMARY:-$ROOT/simulation/snrt/build/snrt_bundle_gate_summary.txt}"
EVIDENCE_DIR="$ROOT/.snrt-bundle-gate-${SLURM_JOB_ID}.evidence"
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

finalize_gate() {
  local exit_status=$? archive_ok=1 stage_log
  trap - EXIT
  if ! mkdir "$EVIDENCE_DIR" 2>/dev/null; then
    echo "SNRT_BUNDLE_GATE_EVIDENCE_FAIL path_exists_or_unwritable=$EVIDENCE_DIR" >&2
    archive_ok=0
  else
    mkdir "$EVIDENCE_DIR/logs" || archive_ok=0
    shopt -s nullglob
    for stage_log in "$RUN_DIR"/*.log; do
      cp -p -- "$stage_log" "$EVIDENCE_DIR/logs/" || archive_ok=0
    done
    if [[ -f "$RUN_DIR/native_symbols.txt" ]]; then
      cp -p -- "$RUN_DIR/native_symbols.txt" "$EVIDENCE_DIR/" || archive_ok=0
    fi
    if [[ -d "$RUN_DIR/m5-live-dust" ]]; then
      for input_file in effective.nml ic_sink inputs.sha256 binary.sha256; do
        if [[ -f "$RUN_DIR/m5-live-dust/$input_file" ]]; then
          cp -p -- "$RUN_DIR/m5-live-dust/$input_file" "$EVIDENCE_DIR/" || archive_ok=0
        fi
      done
    fi
    cp -p -- "$SUMMARY" "$EVIDENCE_DIR/summary.txt" 2>/dev/null || archive_ok=0
    printf 'exit_status=%s\n' "$exit_status" > "$EVIDENCE_DIR/run-status.txt" || archive_ok=0
  fi
  if (( archive_ok )); then
    rm -rf -- "$RUN_DIR"
    echo "SNRT_BUNDLE_GATE_EVIDENCE path=$EVIDENCE_DIR"
  else
    echo "SNRT_BUNDLE_GATE_TEMPORARIES_RETAINED path=$RUN_DIR" >&2
  fi
  return "$exit_status"
}

require_command() {
  local command_name="$1"
  command -v "$command_name" >/dev/null 2>&1 || fail_gate "missing_command=$command_name"
}

[[ ! -e "$EVIDENCE_DIR" ]] || fail_gate "evidence_dir_already_exists=$EVIDENCE_DIR"
RUN_DIR="$(mktemp -d "$ROOT/.snrt-bundle-gate.XXXXXX")"
trap finalize_gate EXIT

time_now_ns() {
  date +%s%N
}

elapsed_seconds() {
  awk -v start="$1" -v end="$2" 'BEGIN { printf "%.3f", (end-start)/1e9 }'
}

run_m5_live_dust_integration() {
  local binary="$1" case_dir="$RUN_DIR/m5-live-dust" template="$ROOT/simulation/snrt/config/snrt_stellar_agn_dust_reference_smoke.nml"
  mkdir "$case_dir"
  cp -- "$template" "$case_dir/effective.nml"
  test "$(grep -c "radiation_transport='snrt_sn'" "$case_dir/effective.nml")" -eq 1 || \
    fail_gate 'M5_live_dust_template_selector_changed'
  test "$(grep -c 'foutput=1$' "$case_dir/effective.nml")" -eq 1 || \
    fail_gate 'M5_live_dust_template_output_clock_changed'
  sed -i \
    -e "s/radiation_transport='snrt_sn'/radiation_transport='snrt_mn', snrt_moment_order=5/" \
    -e 's/foutput=1$/foutput=1000000/' "$case_dir/effective.nml"
  cp -- "$ROOT/simulation/snrt/config/snrt_agn_accretion_smoke.ic_sink" "$case_dir/ic_sink"
  sha256sum "$binary" > "$case_dir/binary.sha256"
  sha256sum "$template" "$case_dir/effective.nml" "$case_dir/ic_sink" \
    "$ROOT/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml" \
    "$ROOT/simulation/snrt/config/snrt_secondary_table_contract_v1.nml" \
    "$ROOT/simulation/snrt/config/snrt_stellar_sed_reference_control_v1.nml" \
    "$ROOT/simulation/snrt/config/dust_dl01_bulk_030_scattering_exchange_reference_v4.nml" \
    /gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_agn_driver_faithful_smoke_yields.dat \
    "$ROOT/simulation/snrt/tests/run_snrt_bundle_gate.sh" \
    > "$case_dir/inputs.sha256"
  grep -Fq "radiation_transport='snrt_mn', snrt_moment_order=5" "$case_dir/effective.nml" || \
    fail_gate 'M5_live_dust_effective_selector_missing'
  grep -Fq 'nstepmax=3' "$case_dir/effective.nml" || fail_gate 'M5_live_dust_step_limit_changed'
  printf 'M5_LIVE_DUST_RUN effective_nml=%s noutput=1 aout=2.0(future) tout=1e30 foutput=fbackup=1000000 expected_dumps=0 nstepmax=3 reference_controls=NONPRODUCTION\n' \
    "$case_dir/effective.nml"
  (
    cd "$case_dir"
    export OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
    export OMP_STACKSIZE=512M KMP_STACKSIZE=512M I_MPI_FABRICS=shm I_MPI_PIN=0
    export SNRT_BACKEND=auto SNRT_DUST_BACKEND=auto SNRT_AGN_MODEL=partition_reference_v1
    export SNRT_REDUCED_C=0.01 SNRT_RT_LEVEL=3 SNRT_ALLOW_REFERENCE_CONTROL=1
    export SNRT_GROUP_CONTRACT="$ROOT/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml"
    export SNRT_SECONDARY_TABLE_CONTRACT="$ROOT/simulation/snrt/config/snrt_secondary_table_contract_v1.nml"
    export SNRT_STELLAR_SED="$ROOT/simulation/snrt/config/snrt_stellar_sed_reference_control_v1.nml"
    export SNRT_DUST_CONTRACT="$ROOT/simulation/snrt/config/dust_dl01_bulk_030_scattering_exchange_reference_v4.nml"
    unset SNRT_DRIVER_TEST_SEED_SOURCE SNRT_RT_TX_TEST_FAIL_STAGE SNRT_RT_TX_TEST_FAIL_LEAF
    timeout 2400s srun --mpi=pmi2 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
      --cpu-bind=none "$binary" effective.nml
  )
  shopt -s nullglob
  local raw_outputs=("$case_dir"/output_[0-9][0-9][0-9][0-9][0-9])
  (( ${#raw_outputs[@]} == 0 )) || fail_gate 'M5_live_dust_unexpected_raw_output'
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
source_head="$(git -C "$ROOT" rev-parse HEAD)"
source_branch="$(git -C "$ROOT" rev-parse --abbrev-ref HEAD)"
source_hash="$source_head"
if ! git -C "$ROOT" diff-index --quiet HEAD --; then source_hash+="-dirty"; fi
source_remote="$(git -C "$ROOT" config --get remote.origin.url)"
printf 'SNRT_BUNDLE_GATE_SOURCE head=%s branch=%s remote=%s tracked_tree=%s\n' \
  "$source_head" "$source_branch" "$source_remote" "$source_hash" | tee -a "$SUMMARY"
printf 'SNRT_BUNDLE_GATE_TOOLS mpiifx=%s gfortran=%s mpirun=%s nvcc=%s\n' \
  "$(command -v mpiifx)" "$(command -v gfortran)" "$(command -v mpirun)" "$NVCC_BIN" | tee -a "$SUMMARY"

BUILD_DIR="$RUN_DIR/build"
mkdir -p "$BUILD_DIR/bin"
for source_dir in patch hydro aton pm poisson amr utils simulation tests mhd; do
  ln -s "$ROOT/$source_dir" "$BUILD_DIR/$source_dir"
done
cp "$ROOT/bin/Makefile" "$BUILD_DIR/bin/Makefile"

# Keep the production Makefile and its ordered VPATH intact, but isolate every
# generated .o/.mod and the linked binary from the shared ROOT/bin directory.
run_stage production_build \
  make -C "$BUILD_DIR/bin" -B GITBRANCH="$source_branch" GITHASH="$source_hash" \
  GITREPO="$source_remote" EXEC="$BUILD_DIR/bin/ramses_bundle" \
  SNRT=1 DUST_LIVE=1 HDF5=1 USE_CUDA=1 USE_FFTW=0 ramses

BINARY="$BUILD_DIR/bin/ramses_bundle3d"
[[ -x "$BINARY" ]] || fail_gate "production_binary_missing=$BINARY"

run_stage m5_live_dust_actual_sources run_m5_live_dust_integration "$BINARY"
grep -q 'Run completed' "$LAST_STAGE_LOG" || fail_gate 'M5_live_dust_run_incomplete'
grep -q 'Radiation transport selector: snrt_mn' "$LAST_STAGE_LOG" || fail_gate 'M5_selector_not_admitted'
grep -Eq 'SNRT_MN_COMMIT level/order/substeps=[[:space:]]*3[[:space:]]+5' "$LAST_STAGE_LOG" || \
  fail_gate 'M5_order5_commit_missing'
mn_commit_count="$(grep -Ec 'SNRT_MN_COMMIT level/order/substeps=[[:space:]]*3[[:space:]]+5' "$LAST_STAGE_LOG")"
(( mn_commit_count >= 2 )) || fail_gate "M5_actual_sources_commit_count=$mn_commit_count"
grep -q 'SNRT_MN_LEDGER source/absorbed/escaped/projection/residual/gas_heat/dust_heat/N=' \
  "$LAST_STAGE_LOG" || fail_gate 'M5_photon_material_ledger_missing'
awk '/SNRT_MN_LEDGER/ { if (($3+0)>0 && ($8+0)>0 && ($9+0)>0) exchanged=1 }
     END { exit !exchanged }' "$LAST_STAGE_LOG" || \
  fail_gate 'M5_actual_sources_gas_and_dust_exchange_missing'
awk '/SNRT_MN_ENERGY/ { seen++; if (sqrt(($6+0)*($6+0))>1.0e-12) bad=1 }
     END { exit !(seen>0 && !bad) }' "$LAST_STAGE_LOG" || \
  fail_gate 'M5_photon_energy_ledger_residual_exceeded'
awk '/SNRT_MN_STATE/ { if (($5+0)>0 && ($6+0)>0) dust_ir_state=1 }
     END { exit !dust_ir_state }' "$LAST_STAGE_LOG" || \
  fail_gate 'M5_live_dust_IR_state_not_updated'
awk '/SNRT_MN_LEDGER/ { if (($3+0)>0 && ($4+0)>0) sourced=1 }
     END { exit !sourced }' "$LAST_STAGE_LOG" || \
  fail_gate 'M5_actual_source_ledger_missing'
if rg -qi 'ERROR:|Aborting|forrtl:|failed to converge|NaN_CHK.*=[[:space:]]*[1-9]' "$LAST_STAGE_LOG"; then
  fail_gate 'M5_live_dust_run_error_marker'
fi
printf 'INTEGRATION_ASSERTION M5_order=5 actual_sources=nonzero gas_dust_IR_energy_ledgers=closed dust_IR_state=updated raw_outputs=0\n' | tee -a "$SUMMARY"

run_stage agn_partition_reference \
  env AMR_MODULE_DIR="$BUILD_DIR/bin" \
  AMR_PARAMETERS_OBJECT="$BUILD_DIR/bin/amr_parameters.jaehyun.o" \
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
  env P04_BINARY="$BINARY" P04_EVIDENCE_DIR="$RUN_DIR/p04_fail_closed" \
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
