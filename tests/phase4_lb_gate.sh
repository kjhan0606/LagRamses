#!/bin/bash
# [RESIZABLE] Production-default load-balance/grid-growth equivalence gate.
# The paired runs intentionally differ only in initial ngridtot.  Candidate
# growth must make the final capacity and all legacy output bytes identical.

set -uo pipefail

control=${1:?control run directory}
candidate=${2:?candidate run directory}
fail=0

bad() {
  printf 'FAIL: %s\n' "$*"
  fail=1
}

for item in "$control" "$candidate" "$control/run.nml" \
  "$candidate/run.nml" "$control/run.log" "$candidate/run.log"; do
  [ -e "$item" ] || bad "missing $item"
done

normalize_nml() {
  sed -e '/^! \[RESIZABLE\] Phase 4 production-LB headroom/d' \
      -e '/^ngridtot=/d' "$1"
}

if [ -f "$control/run.nml" ] && [ -f "$candidate/run.nml" ]; then
  diff -q <(normalize_nml "$control/run.nml") \
    <(normalize_nml "$candidate/run.nml") >/dev/null || \
    bad 'paired namelists differ beyond identity comment and ngridtot'
  grep -qx 'ngridtot=3629056' "$control/run.nml" || \
    bad 'control ngridtot is not 3629056'
  grep -qx 'ngridtot=2900000' "$candidate/run.nml" || \
    bad 'candidate ngridtot is not 2900000'
  for nml in "$control/run.nml" "$candidate/run.nml"; do
    grep -qx 'lb_grid_headroom=0.85' "$nml" || \
      bad "$(basename "$(dirname "$nml")") does not use default LB headroom"
  done
fi

for run in "$control" "$candidate"; do
  [ -f "$run/run.log" ] || continue
  grep -q 'Run completed' "$run/run.log" || \
    bad "$(basename "$run") did not complete"
  if grep -qiE 'Increase ngridmax|MPI_ABORT|forrtl: severe|Segmentation fault|error stop' \
      "$run/run.log"; then
    bad "$(basename "$run") contains a fatal runtime marker"
  fi
  if grep -E 'NaN_CHK.*(uold= *[1-9]|f= *[1-9]|d0= *[1-9])' \
      "$run/run.log" >/dev/null; then
    bad "$(basename "$run") reports nonzero NaN counters"
  fi
  grep -q 'LB grid usage:' "$run/run.log" || \
    bad "$(basename "$run") never exercised the LB grid guard"
  grep -q 'load_balance total:' "$run/run.log" || \
    bad "$(basename "$run") never completed a real load balance"
  if grep -q 'Bounded remap: no safe progress' "$run/run.log"; then
    bad "$(basename "$run") used the no-op LB path"
  fi
done

if [ -f "$control/run.log" ]; then
  control_grows=$(grep -c '\[RESIZABLE\] GRID_GROW' \
    "$control/run.log" 2>/dev/null || true)
  [ "$control_grows" -eq 0 ] || bad "control has $control_grows GRID_GROW events"
fi

if [ -f "$candidate/run.log" ]; then
  if ! awk '
    /\[RESIZABLE\] GRID_GROW/ {
      old=""; new=""; rank=""
      for(i=1;i<=NF;i++) {
        if($i~/^rank=/) {split($i,a,"="); rank=a[2]+0}
        if($i~/^old=/)  {split($i,a,"="); old=a[2]+0}
        if($i~/^new=/)  {split($i,a,"="); new=a[2]+0}
      }
      count++
      if(rank!=1 || old!=90688 || new!=113408) bad=1
    }
    END {exit(count==1 && !bad ? 0 : 1)}
  ' "$candidate/run.log"; then
    bad 'candidate growth is not exactly rank=1 old=90688 new=113408 once'
  fi

  grow_line=$(grep -n -m1 '\[RESIZABLE\] GRID_GROW' \
    "$candidate/run.log" | cut -d: -f1)
  first_lb_line=$(grep -n -m1 'Load balancing AMR grid' \
    "$candidate/run.log" | cut -d: -f1)
  if [ -z "$grow_line" ] || [ -z "$first_lb_line" ] || \
     [ "$grow_line" -ge "$first_lb_line" ]; then
    bad 'candidate did not grow before the first load-balance mutation'
  fi
fi

control_outputs=$(find "$control" -mindepth 1 -maxdepth 1 -type d \
  -name 'output_*' -printf '%f\n' 2>/dev/null | sort)
candidate_outputs=$(find "$candidate" -mindepth 1 -maxdepth 1 -type d \
  -name 'output_*' -printf '%f\n' 2>/dev/null | sort)
[ -n "$control_outputs" ] || bad 'control has no outputs'
[ "$control_outputs" = "$candidate_outputs" ] || \
  bad 'control/candidate output directory sets differ'

amr_files=0
info_files=0
while IFS= read -r output; do
  [ -n "$output" ] || continue
  [ -f "$control/$output/COMPLETE" ] || bad "control $output lacks COMPLETE"
  [ -f "$candidate/$output/COMPLETE" ] || bad "candidate $output lacks COMPLETE"

  control_files=$(find "$control/$output" -type f -printf '%P\n' \
    2>/dev/null | sort)
  candidate_files=$(find "$candidate/$output" -type f -printf '%P\n' \
    2>/dev/null | sort)
  [ "$control_files" = "$candidate_files" ] || \
    bad "$output file sets differ"

  while IFS= read -r relative; do
    [ -n "$relative" ] || continue
    case $(basename "$relative") in
      compilation.txt|makefile.txt|patches.txt)
        continue
        ;;
      # The archived source namelist intentionally records different initial
      # ngridtot.  Its normalized equivalence is checked above.
      namelist.txt)
        continue
        ;;
    esac
    case $(basename "$relative") in
      amr_*.out*) amr_files=$((amr_files+1));;
      info_*.txt) info_files=$((info_files+1));;
    esac
    cmp -s "$control/$output/$relative" \
      "$candidate/$output/$relative" || bad "$output/$relative differs"
  done <<< "$control_files"
done <<< "$control_outputs"

[ "$amr_files" -gt 0 ] || bad 'no AMR binary was compared'
[ "$info_files" -gt 0 ] || bad 'no info file was compared'

physics_lines() {
  grep -E 'PBHDIAG|PBHCACHE|SGS_DT|Main step=|Fine step=|aexp=| Error=' "$1" \
    | sed 's/mem=[0-9.]*% *[0-9.]*%//' \
    | sed 's/[[:space:]]\+/ /g'
}

if [ -f "$control/run.log" ] && [ -f "$candidate/run.log" ]; then
  if ! diff -q <(physics_lines "$control/run.log") \
      <(physics_lines "$candidate/run.log") >/dev/null; then
    bad 'control/candidate physics logs differ'
    diff <(physics_lines "$control/run.log") \
      <(physics_lines "$candidate/run.log") | head -20
  fi
fi

if [ "$fail" -eq 0 ]; then
  echo "GATE: PASS (production LB growth; $amr_files AMR and $info_files info files included in full byte comparison)"
else
  echo 'GATE: FAIL'
fi
exit "$fail"
