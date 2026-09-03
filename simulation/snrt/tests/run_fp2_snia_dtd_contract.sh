#!/usr/bin/env bash
set -euo pipefail

# F-P2 contract and source-admission runner.  This proves interval-integrated
# DTD behavior and the approved physical baseline while SNIa runtime activation
# remains blocked until its AMR/MPI caller is connected.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "$SCRIPT_DIR/../../.." && pwd)"
SOURCE_DIR="$ROOT/simulation/snrt/native/phase0"
BUILD_DIR="$ROOT/build/fp2_snia_dtd"
PRODUCTION_SOURCE_DIR="$ROOT/patch/lagRamses"
PRODUCTION_BUILD_DIR="$ROOT/build/fp2_snia_dtd_production"
FC=mpiifx

mkdir -p "$BUILD_DIR"
for source in stellar_enrichment_config.f90 stellar_enrichment_contract.f90 \
    stellar_snia_dtd.f90 stellar_snia_population_contract.f90 \
    stellar_snia_physical_contract.f90 \
    stellar_snia_event_ledger.f90 stellar_native_units.f90 \
    stellar_snia_cell_deposition.f90 fp2_snia_dtd_test.f90 \
    fp2_snia_population_contract_test.f90 fp2_snia_physical_contract_test.f90 \
    fp2_snia_cell_deposition_test.f90 \
    fp2_snia_event_ledger_test.f90; do
  "$FC" -O2 -g -traceback -warn all -check all -fpp \
    -module "$BUILD_DIR" -c "$SOURCE_DIR/$source" \
    -o "$BUILD_DIR/${source%.f90}.o"
done

"$FC" -O2 -g -traceback -check all \
  "$BUILD_DIR/fp2_snia_dtd_test.o" \
  "$BUILD_DIR/stellar_snia_dtd.o" \
  "$BUILD_DIR/stellar_enrichment_config.o" \
  -o "$BUILD_DIR/fp2_snia_dtd_test"

"$FC" -O2 -g -traceback -check all \
  "$BUILD_DIR/fp2_snia_physical_contract_test.o" \
  "$BUILD_DIR/stellar_snia_physical_contract.o" \
  "$BUILD_DIR/stellar_enrichment_config.o" \
  -o "$BUILD_DIR/fp2_snia_physical_contract_test"

"$FC" -O2 -g -traceback -check all \
  "$BUILD_DIR/fp2_snia_population_contract_test.o" \
  "$BUILD_DIR/stellar_snia_population_contract.o" \
  "$BUILD_DIR/stellar_snia_dtd.o" \
  "$BUILD_DIR/stellar_enrichment_config.o" \
  -o "$BUILD_DIR/fp2_snia_population_contract_test"

"$FC" -O2 -g -traceback -check all \
  "$BUILD_DIR/fp2_snia_cell_deposition_test.o" \
  "$BUILD_DIR/stellar_snia_cell_deposition.o" \
  "$BUILD_DIR/stellar_snia_physical_contract.o" \
  "$BUILD_DIR/stellar_native_units.o" \
  "$BUILD_DIR/stellar_enrichment_contract.o" \
  "$BUILD_DIR/stellar_enrichment_config.o" \
  -o "$BUILD_DIR/fp2_snia_cell_deposition_test"

"$FC" -O2 -g -traceback -check all \
  "$BUILD_DIR/fp2_snia_event_ledger_test.o" \
  "$BUILD_DIR/stellar_snia_event_ledger.o" \
  "$BUILD_DIR/stellar_enrichment_contract.o" \
  "$BUILD_DIR/stellar_enrichment_config.o" \
  -o "$BUILD_DIR/fp2_snia_event_ledger_test"

(cd "$BUILD_DIR" && ./fp2_snia_dtd_test)
(cd "$BUILD_DIR" && ./fp2_snia_population_contract_test)
(cd "$BUILD_DIR" && ./fp2_snia_physical_contract_test)
(cd "$BUILD_DIR" && ./fp2_snia_cell_deposition_test)
(cd "$BUILD_DIR" && ./fp2_snia_event_ledger_test)

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/fp2_snia_event_yield_converter.py"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tools/audit_fp2_snia_event_yield_asset.py" \
  --json-out "$ROOT/simulation/snrt/data/fp2_snia_event_yield_asset_audit.json" \
  > "$BUILD_DIR/fp2_snia_event_yield_asset_audit.stdout"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tools/audit_fp2_snia_keegans_format.py" \
  --json-out "$ROOT/simulation/snrt/data/fp2_snia_keegans_format_audit.json" \
  > "$BUILD_DIR/fp2_snia_keegans_format_audit.stdout"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tools/audit_fp2_snia_hesma.py" \
  --json-out "$ROOT/simulation/snrt/data/fp2_snia_hesma_source_audit.json" \
  > "$BUILD_DIR/fp2_snia_hesma_source_audit.stdout"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/fp2_snia_hesma_source.py"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/fp2_snia_hesma_adapter.py"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/promote_hesma_snia_source.py"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tools/adapt_hesma_snia_source.py" \
  --model n100 \
  --json-out "$ROOT/simulation/snrt/data/fp2_snia_hesma_n100_review_normalized.json" \
  > "$BUILD_DIR/fp2_snia_hesma_n100_review_normalized.stdout"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/fp2_snia_hesma_model_comparison.py"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tools/build_hesma_snia_model_comparison.py" \
  --json-out "$ROOT/simulation/snrt/data/fp2_snia_hesma_model_comparison.json" \
  > "$BUILD_DIR/fp2_snia_hesma_model_comparison.stdout"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/fp2_snia_hesma_profile_estimators.py"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tools/audit_hesma_profile_estimators.py" \
  --json-out "$ROOT/simulation/snrt/data/fp2_snia_hesma_profile_estimator_comparison.json" \
  > "$BUILD_DIR/fp2_snia_hesma_profile_estimators.stdout"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/fp2_snia_hesma_selection_packet.py"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tools/build_hesma_snia_selection_packet.py" \
  --json-out "$BUILD_DIR/fp2_snia_hesma_source_selection_packet_review.json" \
  > "$BUILD_DIR/fp2_snia_hesma_selection_packet.stdout"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tools/audit_fp2_snia_event_source_admission.py" \
  --json-out "$ROOT/simulation/snrt/data/fp2_snia_event_source_admission_audit.json" \
  > "$BUILD_DIR/fp2_snia_event_source_admission.stdout"

"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/fp2_snia_event_source_admission.py"

# Compile the exact production mirror sources as a source-order smoke test.
# This intentionally does not link or activate SNIa in the runtime binary.
mkdir -p "$PRODUCTION_BUILD_DIR"
for source in stellar_enrichment_config.f90 stellar_enrichment_contract.f90 \
    stellar_snia_dtd.f90 stellar_snia_population_contract.f90 \
    stellar_snia_physical_contract.f90 \
    stellar_snia_event_ledger.f90 stellar_native_units.f90 \
    stellar_snia_cell_deposition.f90; do
  "$FC" -O2 -g -traceback -warn all -check all -fpp \
    -module "$PRODUCTION_BUILD_DIR" -I "$PRODUCTION_BUILD_DIR" \
    -c "$PRODUCTION_SOURCE_DIR/$source" \
    -o "$PRODUCTION_BUILD_DIR/${source%.f90}.o"
done
echo "FP2_SNIa_PRODUCTION_MIRROR_COMPILE_OK"

"$ROOT/simulation/snrt/.venv/bin/python" \
"$ROOT/simulation/snrt/tools/audit_fp2_snia_dtd_contract.py" \
  --json-out "$ROOT/simulation/snrt/data/fp2_snia_dtd_contract_audit.json" \
  > "$BUILD_DIR/fp2_snia_dtd_contract_audit.stdout"
"$ROOT/simulation/snrt/.venv/bin/python" \
  "$ROOT/simulation/snrt/tests/fp2_snia_dtd_contract.py"
echo "FP2_SNIa_DTD_CONTRACT_RUN_OK"
