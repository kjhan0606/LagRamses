#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/.." && pwd)
source_dir="$repo_dir/patch/lagRamses"
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT

compiler=${F90:-ifx}
for source in stellar_enrichment_config.f90 stellar_yield_tables.f90 \
  stellar_yield_audit.f90
do
  "$compiler" -O0 -g -traceback -warn all -check all \
    -module "$build_dir" -I"$build_dir" -c "$source_dir/$source" \
    -o "$build_dir/${source%.f90}.o"
done

"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" \
  "$repo_dir/tests/test_stellar_yield_audit_contract.f90" \
  "$build_dir/stellar_enrichment_config.o" \
  "$build_dir/stellar_yield_tables.o" \
  "$build_dir/stellar_yield_audit.o" \
  -o "$build_dir/test_stellar_yield_audit_contract"

"$build_dir/test_stellar_yield_audit_contract"
