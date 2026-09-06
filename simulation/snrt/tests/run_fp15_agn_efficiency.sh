#!/usr/bin/env bash
set -euo pipefail

# RAMSES-independent native contract for the shared AGN efficiency helper and
# accepted-fuel/energy photon budgets and native deposition. It uses a
# scratch directory and never enables SNRT runtime execution.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
FC=${FC:-mpiifx}
BUILD=$(mktemp -d "$ROOT/build/snrt-agn-efficiency.XXXXXX")
trap 'rm -rf "$BUILD"' EXIT

if [[ ! -f "$ROOT/bin/amr_parameters.mod" || ! -f "$ROOT/bin/amr_parameters.jaehyun.o" ]]; then
  echo "missing production amr_parameters module/object under $ROOT/bin" >&2
  exit 2
fi

FFLAGS=(-qopenmp -fpp -O0 -g -check all -check noarg_temp_created -I"$BUILD" -I"$ROOT/bin" -I"$ROOT/patch/lagRamses" -module "$BUILD")
"$FC" "${FFLAGS[@]}" -c "$ROOT/patch/lagRamses/snrt_agn_efficiency.f90" -o "$BUILD/snrt_agn_efficiency.o"
"$FC" "${FFLAGS[@]}" "$ROOT/patch/lagRamses/snrt_agn_efficiency_smoke.f90" \
  "$BUILD/snrt_agn_efficiency.o" "$ROOT/bin/amr_parameters.jaehyun.o" -o "$BUILD/efficiency_smoke"
"$BUILD/efficiency_smoke"

"$FC" "${FFLAGS[@]}" -c "$ROOT/patch/lagRamses/snrt_agn_source.f90" -o "$BUILD/snrt_agn_source.o"
"$FC" "${FFLAGS[@]}" "$ROOT/patch/lagRamses/snrt_agn_source_smoke.f90" \
  "$BUILD/snrt_agn_source.o" "$ROOT/bin/amr_parameters.jaehyun.o" -o "$BUILD/source_smoke"
"$BUILD/source_smoke"

"$FC" "${FFLAGS[@]}" -c "$ROOT/patch/lagRamses/agn_feedback_deposition.f90" -o "$BUILD/agn_feedback_deposition.o"
"$FC" "${FFLAGS[@]}" "$ROOT/patch/lagRamses/agn_feedback_deposition_smoke.f90" \
  "$BUILD/agn_feedback_deposition.o" "$ROOT/bin/amr_parameters.jaehyun.o" -o "$BUILD/deposition_smoke"
"$BUILD/deposition_smoke"

echo "SNRT_AGN_EFFICIENCY_NATIVE_TEST_OK helper=compiled source_api=compiled runtime=disabled"
