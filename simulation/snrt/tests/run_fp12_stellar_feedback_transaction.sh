#!/usr/bin/env bash
set -euo pipefail

# Focused F-P1.2 native bridge/transaction evidence.  It uses only the
# RAMSES-independent patch modules and never activates a production run.
ROOT=/gpfs/kjhan/LRD_JWST
SOURCE_DIR="$ROOT/patch/lagRamses"
TEST_SOURCE="$ROOT/simulation/snrt/native/phase0/fp12_stellar_feedback_transaction_test.f90"
BUILD_DIR="$ROOT/build/fp12_stellar_feedback_transaction"
FC=ifx

# The synchronization evidence must exercise more than one OpenMP worker.
# Pin the test configuration so the reported full-row lock case is
# reproducible across login-node defaults.
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=FALSE

mkdir -p "$BUILD_DIR"

sources=(
  stellar_enrichment_config.f90
  stellar_enrichment_contract.f90
  stellar_snia_physical_contract.f90
  stellar_native_units.f90
  stellar_snia_cell_deposition.f90
  stellar_ramses_field_map.f90
  stellar_ramses_bridge.f90
)

for source in "${sources[@]}"; do
  "$FC" -O0 -g -traceback -warn all -check all -qopenmp \
    -module "$BUILD_DIR" -c "$SOURCE_DIR/$source" \
    -o "$BUILD_DIR/${source%.f90}.o"
done

"$FC" -O0 -g -traceback -warn all -check all -qopenmp \
  -module "$BUILD_DIR" -I"$BUILD_DIR" -c "$TEST_SOURCE" \
  -o "$BUILD_DIR/fp12_stellar_feedback_transaction_test.o"

objects=(
  fp12_stellar_feedback_transaction_test.o
  stellar_ramses_bridge.o
  stellar_ramses_field_map.o
  stellar_snia_cell_deposition.o
  stellar_snia_physical_contract.o
  stellar_native_units.o
  stellar_enrichment_contract.o
  stellar_enrichment_config.o
)
"$FC" -qopenmp "${objects[@]/#/$BUILD_DIR/}" \
  -o "$BUILD_DIR/fp12_stellar_feedback_transaction_test"

"$BUILD_DIR/fp12_stellar_feedback_transaction_test"
echo "FP12_NATIVE_TRANSACTION_RUN_OK"
