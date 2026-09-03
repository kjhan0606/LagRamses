#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/.." && pwd)
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT
compiler=${F90:-ifx}

"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/patch/lagRamses/stellar_enrichment_config.f90" \
  -o "$build_dir/stellar_enrichment_config.o"
"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/patch/lagRamses/stellar_enrichment_contract.f90" \
  -o "$build_dir/stellar_enrichment_contract.o"
"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/patch/lagRamses/stellar_yield_tables.f90" \
  -o "$build_dir/stellar_yield_tables.o"
"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/patch/lagRamses/stellar_yield_interpolation.f90" \
  -o "$build_dir/stellar_yield_interpolation.o"
"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/patch/lagRamses/stellar_yield_provider.f90" \
  -o "$build_dir/stellar_yield_provider.o"
"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/patch/lagRamses/stellar_ssp_sources.f90" \
  -o "$build_dir/stellar_ssp_sources.o"
"$compiler" -O0 -g -traceback -warn all -check all \
  -module "$build_dir" -I"$build_dir" \
  "$repo_dir/tests/stellar_imf_reference.f90" \
  "$build_dir/stellar_enrichment_config.o" \
  "$build_dir/stellar_enrichment_contract.o" \
  "$build_dir/stellar_yield_tables.o" \
  "$build_dir/stellar_yield_interpolation.o" \
  "$build_dir/stellar_yield_provider.o" \
  "$build_dir/stellar_ssp_sources.o" \
  -o "$build_dir/stellar_imf_reference"

"$build_dir/stellar_imf_reference" > "$build_dir/fortran_imf.txt"
"$repo_dir/simulation/snrt/.venv/bin/python" \
  "$repo_dir/simulation/snrt/tests/fp1_imf_jax_differential.py" \
  "$build_dir/fortran_imf.txt" \
  "$repo_dir/simulation/snrt/data/fp1_imf_jax_differential.json"
