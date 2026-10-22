#!/usr/bin/env bash
set -euo pipefail

# RAMSES-independent native contract for the shared AGN efficiency helper and
# accepted-fuel/energy photon budgets and native deposition. It uses a
# scratch directory and never enables SNRT runtime execution.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
FC=${FC:-mpiifx}
BUILD=$(mktemp -d "$ROOT/.snrt-agn-efficiency.XXXXXX")
trap 'rm -rf "$BUILD"' EXIT
AMR_MODULE_DIR=${AMR_MODULE_DIR:-$ROOT/bin}
AMR_PARAMETERS_OBJECT=${AMR_PARAMETERS_OBJECT:-$AMR_MODULE_DIR/amr_parameters.jaehyun.o}

if [[ ! -f "$AMR_MODULE_DIR/amr_parameters.mod" || ! -f "$AMR_PARAMETERS_OBJECT" ]]; then
  echo "missing production amr_parameters module/object under $AMR_MODULE_DIR" >&2
  exit 2
fi

FFLAGS=(-qopenmp -fpp -O0 -g -check all -check noarg_temp_created -I"$BUILD" -I"$AMR_MODULE_DIR" -I"$ROOT/patch/lagRamses" -module "$BUILD")
"$FC" "${FFLAGS[@]}" -c "$ROOT/patch/lagRamses/snrt_agn_efficiency.f90" -o "$BUILD/snrt_agn_efficiency.o"
"$FC" "${FFLAGS[@]}" "$ROOT/patch/lagRamses/snrt_agn_efficiency_smoke.f90" \
  "$BUILD/snrt_agn_efficiency.o" "$AMR_PARAMETERS_OBJECT" -o "$BUILD/efficiency_smoke"
"$BUILD/efficiency_smoke"

"$FC" "${FFLAGS[@]}" -c "$ROOT/patch/lagRamses/snrt_agn_source.f90" -o "$BUILD/snrt_agn_source.o"
"$FC" "${FFLAGS[@]}" -DWITHOUTMPI -c "$ROOT/patch/lagRamses/snrt_rt_transaction.f90" \
  -o "$BUILD/snrt_rt_transaction.o"
"$FC" "${FFLAGS[@]}" "$ROOT/patch/lagRamses/snrt_agn_source_smoke.f90" \
  "$BUILD/snrt_agn_source.o" "$BUILD/snrt_rt_transaction.o" \
  "$AMR_PARAMETERS_OBJECT" -o "$BUILD/source_smoke"
"$BUILD/source_smoke"

"$FC" "${FFLAGS[@]}" -c "$ROOT/patch/lagRamses/agn_feedback_deposition.f90" -o "$BUILD/agn_feedback_deposition.o"
"$FC" "${FFLAGS[@]}" "$ROOT/patch/lagRamses/agn_feedback_deposition_smoke.f90" \
  "$BUILD/agn_feedback_deposition.o" "$AMR_PARAMETERS_OBJECT" -o "$BUILD/deposition_smoke"
"$BUILD/deposition_smoke"

echo "SNRT_AGN_EFFICIENCY_NATIVE_TEST_OK helper=compiled source_api=compiled runtime=disabled"
