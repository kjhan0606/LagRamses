#!/usr/bin/env bash
set -euo pipefail

# DUST-10 native admission smoke.  The tracked fixtures are contract tests,
# not approved opacity or thermal science assets.
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d /gpfs/kjhan/LRD_JWST/.snrt-dust-contract.XXXXXX)"
trap 'rm -rf "$build_dir"' EXIT
module_dir="$build_dir/modules"
mkdir -p "$module_dir"
valid="$repo_root/simulation/snrt/config/dust_native_contract_test.nml"
invalid="$repo_root/simulation/snrt/config/dust_native_contract_invalid_status.nml"

check_reference_modes() {
  local binary="$1" setting expected reference version
  for version in 2 3; do
  reference="$repo_root/simulation/snrt/config/dust_native_reference_control_v${version}.nml"
  for setting in 0 1 invalid 11111111111111111111111111111111; do
    expected=0
    [[ "$setting" != 1 ]] || expected=1
    SNRT_ALLOW_REFERENCE_CONTROL="$setting" SNRT_DUST_CONTRACT="$valid" \
      "$binary" "$valid" "$invalid" "$reference" "$expected"
  done
  env -u SNRT_ALLOW_REFERENCE_CONTROL SNRT_DUST_CONTRACT="$valid" \
    "$binary" "$valid" "$invalid" "$reference" 0
  done
}

if command -v ifx >/dev/null 2>&1; then
  ifx -fpp -DWITHOUTMPI -module "$module_dir" -I"$module_dir" \
    -c "$repo_root/patch/lagRamses/snrt_dust_contract.f90" \
    -o "$build_dir/snrt_dust_contract_ifx.o"
  ifx -fpp -DWITHOUTMPI -module "$module_dir" -I"$module_dir" \
    "$repo_root/patch/lagRamses/snrt_dust_ir.f90" \
    "$repo_root/patch/lagRamses/snrt_dust_contract_smoke.f90" \
    "$build_dir/snrt_dust_contract_ifx.o" -o "$build_dir/snrt_dust_contract_ifx"
  SNRT_DUST_CONTRACT="$valid" "$build_dir/snrt_dust_contract_ifx" "$valid" "$invalid"
  check_reference_modes "$build_dir/snrt_dust_contract_ifx"
  echo SNRT_NATIVE_DUST_CONTRACT_IFX_PASS
else
  echo SNRT_NATIVE_DUST_CONTRACT_IFX_SKIP: ifx_unavailable
fi

if command -v gfortran >/dev/null 2>&1; then
  gnu_module_dir="$build_dir/gnu-modules"
  mkdir -p "$gnu_module_dir"
  gfortran -cpp -ffree-line-length-none -DWITHOUTMPI \
    -J"$gnu_module_dir" -I"$gnu_module_dir" \
    -c "$repo_root/patch/lagRamses/snrt_dust_contract.f90" \
    -o "$build_dir/snrt_dust_contract_gnu.o"
  gfortran -cpp -ffree-line-length-none -DWITHOUTMPI \
    -J"$gnu_module_dir" -I"$gnu_module_dir" \
    "$repo_root/patch/lagRamses/snrt_dust_ir.f90" \
    "$repo_root/patch/lagRamses/snrt_dust_contract_smoke.f90" \
    "$build_dir/snrt_dust_contract_gnu.o" -o "$build_dir/snrt_dust_contract_gnu"
  SNRT_DUST_CONTRACT="$valid" "$build_dir/snrt_dust_contract_gnu" "$valid" "$invalid"
  check_reference_modes "$build_dir/snrt_dust_contract_gnu"
  echo SNRT_NATIVE_DUST_CONTRACT_GNU_PASS
else
  echo SNRT_NATIVE_DUST_CONTRACT_GNU_SKIP: gfortran_unavailable
fi

echo SNRT_NATIVE_DUST_CONTRACT_RUN_PASS
