#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/.." && pwd)
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT

compiler=${F90:-ifx}
for source in \
  stellar_enrichment_config.f90 \
  stellar_yield_tables.f90 \
  stellar_yield_interpolation.f90
do
  object=${source%.f90}.o
  "$compiler" -module "$build_dir" -I"$build_dir" -c \
    "$repo_dir/patch/lagRamses/$source" -o "$build_dir/$object"
done
"$compiler" -module "$build_dir" -I"$build_dir" \
  "$repo_dir/tests/test_stellar_yield_fail_closed.f90" \
  "$build_dir/stellar_enrichment_config.o" \
  "$build_dir/stellar_yield_tables.o" \
  "$build_dir/stellar_yield_interpolation.o" \
  -o "$build_dir/test_stellar_yield_fail_closed"
"$build_dir/test_stellar_yield_fail_closed"
