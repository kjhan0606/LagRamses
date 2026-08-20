#!/bin/bash
# [RESIZABLE] Phase 4 grid-capacity growth gate.  Capacity metadata and legacy-output AMR
# indices legitimately differ when roomy and tight runs end with different
# ngridmax values, so this gate compares physical payloads and physics logs.
set -uo pipefail

ref=${1:?roomy run directory}
cand=${2:?tight run directory}
fail=0

bad() {
  printf 'FAIL: %s\n' "$*"
  fail=1
}

for item in "$ref" "$cand" "$ref/run.log" "$cand/run.log"; do
  [ -e "$item" ] || bad "missing $item"
done

if [ -f "$ref/run.log" ] && [ -f "$cand/run.log" ]; then
  for run in "$ref" "$cand"; do
    grep -q 'Run completed' "$run/run.log" || bad "$(basename "$run") did not complete"
    if grep -qiE 'Increase ngridmax|MPI_ABORT|forrtl: severe|Segmentation fault|error stop' "$run/run.log"; then
      bad "$(basename "$run") contains a fatal runtime marker"
    fi
    if grep -E 'NaN_CHK.*(uold= *[1-9]|f= *[1-9]|d0= *[1-9])' "$run/run.log" >/dev/null; then
      bad "$(basename "$run") reports nonzero NaN counters"
    fi
  done

  awk '
    /\[RESIZABLE\] GRID_GROW/ {
      old=""; new=""
      for (i=1; i<=NF; i++) {
        if ($i ~ /^old=/) { split($i,a,"="); old=a[2]+0 }
        if ($i ~ /^new=/) { split($i,a,"="); new=a[2]+0 }
      }
      if (old=="" || new=="" || new<=old) bad=1
      count++
    }
    END { exit(count>0 && !bad ? 0 : 1) }
  ' "$cand/run.log" || bad 'tight run has no valid GRID_GROW old/new event'
fi

ref_outputs=$(find "$ref" -mindepth 1 -maxdepth 1 -type d -name 'output_*' -printf '%f\n' 2>/dev/null | sort)
cand_outputs=$(find "$cand" -mindepth 1 -maxdepth 1 -type d -name 'output_*' -printf '%f\n' 2>/dev/null | sort)
[ -n "$ref_outputs" ] || bad 'roomy run has no outputs'
[ "$ref_outputs" = "$cand_outputs" ] || bad 'roomy/tight output directory sets differ'

while IFS= read -r output; do
  [ -n "$output" ] || continue
  [ -f "$ref/$output/COMPLETE" ] || bad "roomy $output lacks COMPLETE"
  [ -f "$cand/$output/COMPLETE" ] || bad "tight $output lacks COMPLETE"

  ref_files=$(find "$ref/$output" -maxdepth 1 -type f \
    \( -name 'hydro_*.out*' -o -name 'part_*.out*' -o -name 'grav_*.out*' \) \
    -printf '%f\n' 2>/dev/null | sort)
  cand_files=$(find "$cand/$output" -maxdepth 1 -type f \
    \( -name 'hydro_*.out*' -o -name 'part_*.out*' -o -name 'grav_*.out*' \) \
    -printf '%f\n' 2>/dev/null | sort)
  [ "$ref_files" = "$cand_files" ] || bad "$output physical payload sets differ"

  while IFS= read -r payload; do
    [ -n "$payload" ] || continue
    cmp -s "$ref/$output/$payload" "$cand/$output/$payload" || \
      bad "$output/$payload differs"
  done <<< "$ref_files"
done <<< "$ref_outputs"

physics_lines() {
  grep -E 'PBHDIAG|PBHCACHE|SGS_DT|Main step=|Fine step=|aexp=| Error=' "$1" \
    | sed 's/mem=[0-9.]*% *[0-9.]*%//' \
    | sed 's/[[:space:]]\+/ /g'
}

if [ -f "$ref/run.log" ] && [ -f "$cand/run.log" ]; then
  if ! diff -q <(physics_lines "$ref/run.log") <(physics_lines "$cand/run.log") >/dev/null; then
    bad 'roomy/tight physics logs differ'
    diff <(physics_lines "$ref/run.log") <(physics_lines "$cand/run.log") | head -20
  fi
fi

if [ "$fail" -eq 0 ]; then
  echo 'GATE: PASS (grid growth fired; physical payloads and physics logs are bitwise identical)'
else
  echo 'GATE: FAIL'
fi
exit "$fail"
