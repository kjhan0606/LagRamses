#!/bin/bash
# Layout-change gate. The amr_* files store father/nbor/son as CELL INDICES,
# which a block layout deliberately renumbers, so comparing them bitwise can
# never pass and says nothing about physics. This gate compares the physical
# payloads (hydro, part, grav) bitwise and requires the per-step energy lines
# to match exactly. Use it for any chunk that changes amr_block_size; use
# phase0_bitwise_gate.sh only when the layout is unchanged.
#   phase2_layout_gate.sh <ref_run_dir> <cand_run_dir>
set -uo pipefail
ref=$1; cand=$2; fail=0
for odir in "$ref"/output_*; do
  o=$(basename "$odir"); [ -d "$cand/$o" ] || { echo "MISSING $o"; fail=1; continue; }
  for f in "$odir"/hydro_*.out* "$odir"/part_*.out* "$odir"/grav_*.out*; do
    [ -e "$f" ] || continue
    b="$cand/$o/$(basename "$f")"
    cmp -s "$f" "$b" || { echo "DIFFER $o/$(basename "$f")"; fail=1; }
  done
done
filt(){ grep -E "Main step=|SGS_DT" "$1/run.log" | sed "s/[[:space:]]\+/ /g"; }
diff -q <(filt "$ref") <(filt "$cand") >/dev/null || {
  echo "DIFFER energy/timestep lines:"; diff <(filt "$ref") <(filt "$cand") | head -10; fail=1; }
[ $fail -eq 0 ] && echo "GATE: PASS (physics bitwise identical; amr indices renumbered by design)" || echo "GATE: FAIL"
exit $fail
