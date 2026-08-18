#!/bin/bash
# Bitwise acceptance gate for the Phase-1 index-API conversion
# (BLOCK_GRID_PHASE01_BRIEF.md). Compares two run directories produced by the
# reference and candidate binaries from identical namelists and rank counts.
#
#   tests/phase0_bitwise_gate.sh <ref_run_dir> <cand_run_dir>
#
# PASS requires every output payload byte-identical and the physics lines of
# the logs identical. Timing lines are excluded; nothing else is.
set -uo pipefail
ref=${1:?ref run dir}; cand=${2:?candidate run dir}
fail=0

# 1. payload files: every amr_/hydro_/part_/info_ file in every output_*
for rdir in "$ref"/output_*; do
  odir=$(basename "$rdir")
  cdir="$cand/$odir"
  if [[ ! -d $cdir ]]; then echo "MISSING $odir in candidate"; fail=1; continue; fi
  while IFS= read -r f; do
    rel=${f#"$rdir"/}
    # build provenance differs between builds by design, not physics
    case $(basename "$f") in compilation.txt|makefile.txt|patches.txt) continue;; esac
    if ! cmp -s "$f" "$cdir/$rel"; then
      echo "DIFFER  $odir/$rel"
      fail=1
    fi
  done < <(find "$rdir" -type f | sort)
done
# and no extra outputs on the candidate side
for cdir in "$cand"/output_*; do
  odir=$(basename "$cdir")
  [[ -d $ref/$odir ]] || { echo "EXTRA $odir in candidate"; fail=1; }
done

# 2. physics lines of the log (timings vary run to run, physics must not)
filt() {
  grep -E "PBHDIAG|PBHCACHE|SGS_DT|Main step=|Fine step=|aexp=| Error=" "$1" \
    | sed 's/[[:space:]]\+/ /g'
}
if ! diff -q <(filt "$ref/run.log") <(filt "$cand/run.log") >/dev/null; then
  echo "DIFFER  run.log physics lines:"
  diff <(filt "$ref/run.log") <(filt "$cand/run.log") | head -20
  fail=1
fi

if [[ $fail -eq 0 ]]; then echo "GATE: PASS (bitwise identical)"; else echo "GATE: FAIL"; fi
exit $fail
