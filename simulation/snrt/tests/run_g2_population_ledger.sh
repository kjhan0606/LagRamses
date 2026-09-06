#!/usr/bin/env bash
set -euo pipefail

# G2 population ledger contract.  This is a RAMSES-independent native test;
# it checks channel closure and proves that unsupported event channels cannot
# fall through to ordinary IMF integration.

ROOT=/gpfs/kjhan/LRD_JWST
SOURCE_DIR="$ROOT/simulation/snrt/native/phase0"
BUILD_DIR="$ROOT/build/g2_population_ledger"
FC=mpiifx

mkdir -p "$BUILD_DIR"

sources=(
  stellar_enrichment_config.f90
  stellar_enrichment_contract.f90
  stellar_snia_physical_contract.f90
  stellar_yield_tables.f90
  stellar_yield_interpolation.f90
  stellar_yield_provider.f90
  stellar_ssp_sources.f90
  stellar_source_increment.f90
  stellar_population_ledger.f90
  stellar_enrichment_driver.f90
  g2_population_ledger_test.f90
)

for source in "${sources[@]}"; do
  "$FC" -O2 -g -traceback -warn all -check all -fpp \
    -module "$BUILD_DIR" -c "$SOURCE_DIR/$source" \
    -o "$BUILD_DIR/${source%.f90}.o"
done

objects=(
  g2_population_ledger_test.o
  stellar_enrichment_config.o
  stellar_enrichment_contract.o
  stellar_snia_physical_contract.o
  stellar_yield_tables.o
  stellar_yield_interpolation.o
  stellar_yield_provider.o
  stellar_ssp_sources.o
  stellar_source_increment.o
  stellar_population_ledger.o
  stellar_enrichment_driver.o
)

"$FC" -O2 -g -traceback -check all \
  "${objects[@]/#/$BUILD_DIR/}" -o "$BUILD_DIR/g2_population_ledger_test"

(cd "$BUILD_DIR" && ./g2_population_ledger_test)

echo "G2_POPULATION_LEDGER_RUN_OK"
