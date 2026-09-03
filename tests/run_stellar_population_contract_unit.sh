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
  stellar_yield_tables.f90
  stellar_yield_interpolation.f90
  stellar_yield_provider.f90
  stellar_ssp_sources.f90
  stellar_source_increment.f90
  stellar_population_ledger.f90
  stellar_enrichment_driver.f90
)
objects=()

for source in "${sources[@]}"; do
  "$compiler" -O0 -g -traceback -warn all -check all \
    -module "$build_dir" -I"$build_dir" -c "$source_dir/$source" \
    -o "$build_dir/${source%.f90}.o"
  objects+=("$build_dir/${source%.f90}.o")
done

"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/simulation/snrt/native/phase0/g2_population_ledger_test.f90" \
  -o "$build_dir/g2_population_ledger_test.o"

"$compiler" -O0 -g -traceback -check all \
  "$build_dir/g2_population_ledger_test.o" \
  "${objects[@]}" \
  -o "$build_dir/test_stellar_population_contract"

"$build_dir/test_stellar_population_contract"
echo "STELLAR_POPULATION_CONTRACT_UNIT_OK"
