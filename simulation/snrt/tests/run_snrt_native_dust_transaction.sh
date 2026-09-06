#!/usr/bin/env bash
set -euo pipefail

# Focused native FP64 receiver smoke.  It does not launch RAMSES, CUDA, or
# the Python/JAX reference solver; the actual CUDA boundary is covered by the
# existing A10 multigroup smoke.
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d /gpfs/kjhan/LRD_JWST/.snrt-dust-transaction.XXXXXX)"
trap 'rm -rf "$build_dir"' EXIT
module_dir="$build_dir/modules"
mkdir -p "$module_dir"

fc="${FC:-gfortran}"
fc_name="${fc##*/}"
if [[ "$fc_name" == "ifx" || "$fc_name" == "mpiifx" ]]; then
  flags=(-fpp -DWITHOUTMPI -module "$module_dir" -I"$module_dir")
else
  flags=(-cpp -ffree-line-length-none -DWITHOUTMPI -J"$module_dir" -I"$module_dir")
fi

"$fc" "${flags[@]}" -c \
  "$repo_root/simulation/snrt/tests/fortran/amr_parameters_transaction_smoke_stub.f90" \
  -o "$build_dir/amr_parameters.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_dust_transaction.f90" \
  -o "$build_dir/snrt_dust_transaction.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_dust_transaction_smoke.f90" \
  -o "$build_dir/snrt_dust_transaction_smoke.o"
"$fc" "$build_dir/amr_parameters.o" "$build_dir/snrt_dust_transaction.o" \
  "$build_dir/snrt_dust_transaction_smoke.o" -o "$build_dir/snrt_dust_transaction_smoke"

"$build_dir/snrt_dust_transaction_smoke"
echo SNRT_NATIVE_DUST_TRANSACTION_ALL_OK
