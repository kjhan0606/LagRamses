#!/usr/bin/env bash
set -euo pipefail

# DUST-9 native boundary smoke.  This intentionally does not launch RAMSES or
# Python/JAX: the live RAMSES driver has no dedicated dust state yet.
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d /gpfs/kjhan/LRD_JWST/.snrt-dust-mapping-receiver.XXXXXX)"
trap 'rm -rf "$build_dir"' EXIT
module_dir="$build_dir/modules"
mkdir -p "$module_dir"

if command -v ifx >/dev/null 2>&1; then
  ifx -fpp -DWITHOUTMPI -module "$module_dir" -I"$module_dir" \
    -c "$repo_root/patch/lagRamses/snrt_dust_receiver.f90" \
    -o "$build_dir/snrt_dust_receiver_ifx.o"
  ifx -fpp -DWITHOUTMPI -module "$module_dir" -I"$module_dir" \
    "$repo_root/patch/lagRamses/snrt_dust_receiver_smoke.f90" \
    "$build_dir/snrt_dust_receiver_ifx.o" -o "$build_dir/snrt_dust_receiver_ifx"
  "$build_dir/snrt_dust_receiver_ifx"
  echo SNRT_NATIVE_DUST_MAPPING_RECEIVER_IFX_PASS
else
  echo SNRT_NATIVE_DUST_MAPPING_RECEIVER_IFX_SKIP: ifx_unavailable
fi

if command -v gfortran >/dev/null 2>&1; then
  gnu_module_dir="$build_dir/gnu-modules"
  mkdir -p "$gnu_module_dir"
  gfortran -cpp -ffree-line-length-none -DWITHOUTMPI \
    -J"$gnu_module_dir" -I"$gnu_module_dir" \
    -c "$repo_root/patch/lagRamses/snrt_dust_receiver.f90" \
    -o "$build_dir/snrt_dust_receiver_gnu.o"
  gfortran -cpp -ffree-line-length-none -DWITHOUTMPI \
    -J"$gnu_module_dir" -I"$gnu_module_dir" \
    "$repo_root/patch/lagRamses/snrt_dust_receiver_smoke.f90" \
    "$build_dir/snrt_dust_receiver_gnu.o" -o "$build_dir/snrt_dust_receiver_gnu"
  "$build_dir/snrt_dust_receiver_gnu"
  echo SNRT_NATIVE_DUST_MAPPING_RECEIVER_GNU_PASS
else
  echo SNRT_NATIVE_DUST_MAPPING_RECEIVER_GNU_SKIP: gfortran_unavailable
fi

echo SNRT_NATIVE_DUST_MAPPING_RECEIVER_RUN_PASS
