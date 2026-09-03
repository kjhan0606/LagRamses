#!/usr/bin/env bash
set -euo pipefail

# Reproducible G1 contract runner.  It builds only the RAMSES-independent
# native mirror in an isolated /gpfs build directory, then runs the native
# and JAX differential tests.  No external checkout is modified.

ROOT=/gpfs/kjhan/LRD_JWST
SOURCE_DIR="$ROOT/simulation/snrt/native/phase0"
BUILD_DIR="$ROOT/build/g1_native"
DATA_DIR="$ROOT/simulation/snrt/data"
FC=mpiifx

mkdir -p "$BUILD_DIR" "$DATA_DIR"

# Do not let a native-mirror PASS masquerade as a production-linked PASS. The
# normal path requires fresh production build/linkage/smoke evidence. A
# diagnostic-only path keeps the native differential oracle runnable without
# claiming production closure; it emits a distinct terminal marker.
if [[ "${P0_DIAGNOSTIC:-0}" == 1 ]]; then
  python3 "$ROOT/simulation/snrt/tools/validate_stellar_source_parity.py"
else
  python3 "$ROOT/simulation/snrt/tools/validate_stellar_source_parity.py" --require-pass
fi

sources=(
  stellar_enrichment_config.f90
  stellar_enrichment_contract.f90
  stellar_snia_physical_contract.f90
  stellar_native_units.f90
  stellar_progress_contract.f90
  stellar_yield_tables.f90
  stellar_yield_interpolation.f90
  stellar_yield_provider.f90
  stellar_ssp_sources.f90
  stellar_source_increment.f90
  stellar_population_ledger.f90
  stellar_enrichment_driver.f90
  stellar_yield_audit.f90
  stellar_ramses_field_map.f90
  g1_contract_test.f90
)

for source in "${sources[@]}"; do
  "$FC" -O2 -g -traceback -warn all -check all -fpp \
    -module "$BUILD_DIR" -c "$SOURCE_DIR/$source" \
    -o "$BUILD_DIR/${source%.f90}.o"
done

objects=(
  g1_contract_test.o
  stellar_enrichment_config.o
  stellar_enrichment_contract.o
  stellar_snia_physical_contract.o
  stellar_native_units.o
  stellar_progress_contract.o
  stellar_yield_tables.o
  stellar_yield_interpolation.o
  stellar_yield_provider.o
  stellar_ssp_sources.o
  stellar_source_increment.o
  stellar_population_ledger.o
  stellar_enrichment_driver.o
  stellar_yield_audit.o
  stellar_ramses_field_map.o
)

"$FC" -O2 -g -traceback -check all \
  "${objects[@]/#/$BUILD_DIR/}" -o "$BUILD_DIR/g1_contract_test"

G1_TABLE_PATH="$BUILD_DIR/g1_synthetic_yields.dat" \
G1_NATIVE_RESULT="$BUILD_DIR/g1_native_interpolation.txt" \
  "$BUILD_DIR/g1_contract_test"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/g1_native_jax_differential.py" \
  --table "$BUILD_DIR/g1_synthetic_yields.dat" \
  --native-result "$BUILD_DIR/g1_native_interpolation.txt" \
  --evidence "$DATA_DIR/g1_native_jax_differential.json"

if [[ "${P0_DIAGNOSTIC:-0}" == 1 ]]; then
  echo "G1_NATIVE_DIAGNOSTIC_ONLY"
else
  echo "G1_NATIVE_CONTRACT_RUN_OK"
fi
