#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/.." && pwd)
source_dir="$repo_dir/patch/lagRamses"
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT
compiler=${F90:-ifx}

sources=(
  stellar_enrichment_config.f90
  stellar_enrichment_contract.f90
  stellar_snia_physical_contract.f90
  stellar_native_units.f90
  stellar_snia_cell_deposition.f90
  stellar_cell_deposition.f90
  stellar_ramses_bridge.f90
  stellar_ramses_field_map.f90
  stellar_ramses_mapped_bridge.f90
)
objects=()
for source in "${sources[@]}"; do
  "$compiler" -O0 -g -traceback -warn all -check all \
    -module "$build_dir" -I"$build_dir" -c "$source_dir/$source" \
    -o "$build_dir/${source%.f90}.o"
  objects+=("$build_dir/${source%.f90}.o")
done

"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" \
  "$repo_dir/tests/test_stellar_residual_deposition.f90" \
  "${objects[@]}" -o "$build_dir/test_stellar_residual_deposition"

"$build_dir/test_stellar_residual_deposition"
