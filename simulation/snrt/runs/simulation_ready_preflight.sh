#!/usr/bin/env bash
set -euo pipefail

# Read-only launch preflight for the high-level RT/feedback project. This
# script does not submit a job and does not create a run directory. It makes
# the allowed simulation class explicit before a caller supplies sbatch/srun.

PROJECT_ROOT=/gpfs/kjhan/LRD_JWST
DEFAULT_GROUP_CONTRACT="$PROJECT_ROOT/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml"
DEFAULT_SECONDARY_CONTRACT="$PROJECT_ROOT/simulation/snrt/config/snrt_secondary_table_contract_v1.nml"
READINESS_MANIFEST="$PROJECT_ROOT/manifests/production_readiness_manifest_v1.json"
SNRT_DRIVER="$PROJECT_ROOT/patch/lagRamses/snrt_ramses_driver.f90"

MODE=
BINARY_PATH=
NAMELIST_PATH=
GROUP_CONTRACT="$DEFAULT_GROUP_CONTRACT"
SECONDARY_CONTRACT="$DEFAULT_SECONDARY_CONTRACT"

usage() {
  cat >&2 <<'EOF'
Usage:
  simulation_ready_preflight.sh --mode MODE --binary PATH --namelist PATH [options]

MODE:
  reference-control    bounded nine-group SNRT + hydro wiring qualification
  legacy-comparison    legacy stellar/feedback comparison binary, no SNRT
  physical-production  physical production request; currently fail-closed

Options:
  --group-contract PATH       reference-control spectral contract
  --secondary-contract PATH   reference-control secondary-ionization contract
EOF
}

fail_preflight() {
  echo "SIMULATION_READY_BLOCKED reason=$*" >&2
  exit 2
}

repo_path() {
  local candidate="$1"
  if [[ "$candidate" == /* ]]; then
    printf '%s\n' "$candidate"
  else
    printf '%s/%s\n' "$PROJECT_ROOT" "$candidate"
  fi
}

require_file() {
  local label="$1"
  local path="$2"
  [[ -f "$path" ]] || fail_preflight "missing_${label}=$path"
}

require_nml_value() {
  local key="$1"
  local expected="$2"
  local actual
  actual="$(sed -nE "s/^[[:space:]]*${key}[[:space:]]*=[[:space:]]*([^!,[:space:]]+).*/\1/p" "$NAMELIST_PATH" | head -n 1 | tr '[:upper:]' '[:lower:]')"
  [[ "$actual" == "$expected" ]] || \
    fail_preflight "namelist_${key}_expected=${expected}_actual=${actual:-missing}"
}

require_nml_pattern() {
  local label="$1"
  local pattern="$2"
  grep -Eiq -- "$pattern" "$NAMELIST_PATH" || \
    fail_preflight "namelist_${label}_missing"
}

has_snrt_symbols() {
  nm -g --defined-only "$BINARY_PATH" 2>/dev/null \
    | grep -Ei 'snrt_(ramses_driver|rt_transaction|transport_step)' >/dev/null
}

while (($# > 0)); do
  case "$1" in
    --mode)
      (($# >= 2)) || { usage; exit 64; }
      MODE="$2"
      shift 2
      ;;
    --binary)
      (($# >= 2)) || { usage; exit 64; }
      BINARY_PATH="$(repo_path "$2")"
      shift 2
      ;;
    --namelist)
      (($# >= 2)) || { usage; exit 64; }
      NAMELIST_PATH="$(repo_path "$2")"
      shift 2
      ;;
    --group-contract)
      (($# >= 2)) || { usage; exit 64; }
      GROUP_CONTRACT="$(repo_path "$2")"
      shift 2
      ;;
    --secondary-contract)
      (($# >= 2)) || { usage; exit 64; }
      SECONDARY_CONTRACT="$(repo_path "$2")"
      shift 2
      ;;
    -h|--help)
      usage >&1
      exit 0
      ;;
    *)
      usage
      exit 64
      ;;
  esac
done

[[ "$MODE" == reference-control || "$MODE" == legacy-comparison || \
   "$MODE" == physical-production ]] || {
  usage
  exit 64
}
[[ -n "$BINARY_PATH" ]] || { usage; exit 64; }
[[ -n "$NAMELIST_PATH" ]] || { usage; exit 64; }
[[ -x "$BINARY_PATH" ]] || fail_preflight "binary_not_executable=$BINARY_PATH"
require_file namelist "$NAMELIST_PATH"
command -v nm >/dev/null 2>&1 || fail_preflight missing_command=nm

if [[ "$MODE" == physical-production ]]; then
  # This is deliberately an explicit refusal. A future physical release
  # must change the manifest, source contracts, and this admission boundary
  # together; a caller must not bypass it by selecting a different binary.
  if grep -Fq '"production_status": "blocked_until_all_required_assets_are_approved"' \
       "$READINESS_MANIFEST"; then
    fail_preflight 'physical_production_manifest_blocked'
  fi
  if grep -Fq 'ZERO_SCAFFOLD' "$SNRT_DRIVER"; then
    fail_preflight 'physical_production_live_dust_is_zero_scaffold'
  fi
  fail_preflight physical_production_requires_explicit_release_record
fi

require_nml_pattern hydro '^[[:space:]]*hydro[[:space:]]*=[[:space:]]*\.true\.'
require_nml_pattern feedback_mode "feedback_mode[[:space:]]*=[[:space:]]*'legacy'"

case "$MODE" in
  reference-control)
    has_snrt_symbols || fail_preflight "snrt_symbols_missing=$BINARY_PATH"
    require_file group_contract "$GROUP_CONTRACT"
    require_file secondary_contract "$SECONDARY_CONTRACT"
    grep -Fq "contract_status='reference_control'" "$GROUP_CONTRACT" || \
      fail_preflight group_contract_is_not_reference_control
    grep -Fq 'nenergy_contract=258' "$SECONDARY_CONTRACT" || \
      fail_preflight secondary_contract_energy_axis_mismatch
    grep -Fq 'nxi_contract=14' "$SECONDARY_CONTRACT" || \
      fail_preflight secondary_contract_xi_axis_mismatch
    require_nml_value cosmo .false.
    require_nml_value pic .false.
    require_nml_value poisson .false.
    require_nml_value sink .false.
    require_nml_value noutput 1
    require_nml_value aout 2.0d0
    require_nml_value tout 1.0d30
    require_nml_value foutput 1000000
    require_nml_value fbackup 1000000
    ;;
  legacy-comparison)
    has_snrt_symbols && fail_preflight "legacy_binary_contains_snrt_symbols=$BINARY_PATH"
    ;;
esac

printf 'SIMULATION_READY_BEGIN mode=%s\n' "$MODE"
printf 'SIMULATION_READY_BINARY %s\n' "$BINARY_PATH"
printf 'SIMULATION_READY_BINARY_SHA256 %s\n' "$(sha256sum "$BINARY_PATH" | awk '{print $1}')"
printf 'SIMULATION_READY_NAMELIST %s\n' "$NAMELIST_PATH"
printf 'SIMULATION_READY_NAMELIST_SHA256 %s\n' "$(sha256sum "$NAMELIST_PATH" | awk '{print $1}')"
printf 'SIMULATION_READY_GROUP_CONTRACT %s\n' "$GROUP_CONTRACT"
printf 'SIMULATION_READY_SECONDARY_CONTRACT %s\n' "$SECONDARY_CONTRACT"
printf 'SIMULATION_READY_REPOSITORY_HEAD %s\n' "$(git -C "$PROJECT_ROOT" rev-parse HEAD)"
printf 'SIMULATION_READY_WORKTREE_ENTRIES %s\n' "$(git -C "$PROJECT_ROOT" status --porcelain | wc -l)"
if [[ "$MODE" == reference-control ]]; then
  printf 'SIMULATION_READY_OUTPUT_POLICY noutput=1 aout=2.0d0 tout=1.0d30 foutput=1000000 fbackup=1000000\n'
else
  printf 'SIMULATION_READY_OUTPUT_POLICY_NOTE caller_must_audit_noutput_aout_tout_foutput_fbackup\n'
fi
printf 'SIMULATION_READY_PHYSICAL_STATUS blocked_until_required_assets_and_live_receivers\n'
printf 'SIMULATION_READY_PASS mode=%s\n' "$MODE"
