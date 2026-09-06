#!/usr/bin/env bash
set -euo pipefail

# Supporting wiring check for the live SNRT hydro path.  This is intentionally
# a source-order assertion: the native initialized-RAMSES gate owns runtime
# behavior, while this check prevents a later edit from moving the required
# post-SNRT restriction back before the receiver.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE="$ROOT/patch/lagRamses/amr_step.jaehyun.f90"

[[ -f "$SOURCE" ]] || { echo "SNRT_HYDRO_RESTRICTION_FAIL missing_source=$SOURCE" >&2; exit 1; }

snrt_line="$(rg -n '^ *call snrt_ramses_advance_level\(ilevel\)' "$SOURCE" | cut -d: -f1)"
post_upload_line="$(awk -v start="$snrt_line" 'NR > start && /if\(hydro \.and\. snrt_agn_rt_requested\(\)\) call upload_fine\(ilevel\)/ { print NR; exit }' "$SOURCE")"
diagnose_line="$(rg -n '^ *call snrt_ramses_diagnose_level\(ilevel\)' "$SOURCE" | cut -d: -f1)"

[[ -n "$snrt_line" ]] || { echo "SNRT_HYDRO_RESTRICTION_FAIL advance_call_missing" >&2; exit 1; }
[[ -n "$post_upload_line" ]] || { echo "SNRT_HYDRO_RESTRICTION_FAIL post_receiver_restriction_missing" >&2; exit 1; }
[[ -n "$diagnose_line" ]] || { echo "SNRT_HYDRO_RESTRICTION_FAIL diagnostic_call_missing" >&2; exit 1; }
(( post_upload_line > snrt_line )) || { echo "SNRT_HYDRO_RESTRICTION_FAIL restriction_precedes_receiver" >&2; exit 1; }
(( diagnose_line > post_upload_line )) || { echo "SNRT_HYDRO_RESTRICTION_FAIL diagnostic_precedes_restriction" >&2; exit 1; }

awk -v advance="$snrt_line" '
  NR <= advance && /^#ifdef SNRT$/ { opened = 1 }
  NR <= advance && opened && /^#endif$/ { closed = 1 }
  END { exit(opened && closed ? 0 : 1) }
' "$SOURCE" || { echo "SNRT_HYDRO_RESTRICTION_FAIL missing_SNRT_guard" >&2; exit 1; }

rg -q '^  use snrt_agn_efficiency, only: snrt_agn_rt_requested$' "$SOURCE" || {
  echo "SNRT_HYDRO_RESTRICTION_FAIL request_latch_import_missing" >&2
  exit 1
}

git -C "$ROOT" diff --check
printf 'SNRT_HYDRO_RESTRICTION_WIRING_OK receiver_line=%s post_restriction_line=%s diagnostic_line=%s\n' \
  "$snrt_line" "$post_upload_line" "$diagnose_line"
