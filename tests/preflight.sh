#!/bin/bash
# Pre-flight for a lagRamses QA run. Every check here corresponds to a test
# that was actually wasted by skipping it. Sourced or called at the TOP of a
# test script; exits nonzero and prints why, so the test never runs on a
# wrong premise.
#
#   tests/preflight.sh --exe <binary> [--commit <sha>] [--nml <namelist>]
#                      [--need-switch name=value ...] [--gpu-free-gb N]
set -uo pipefail
exe=""; commit=""; nml=""; gpu_need=""; declare -a switches=()
while [ $# -gt 0 ]; do
  case "$1" in
    --exe) exe=$2; shift 2;;
    --commit) commit=$2; shift 2;;
    --nml) nml=$2; shift 2;;
    --need-switch) switches+=("$2"); shift 2;;
    --gpu-free-gb) gpu_need=$2; shift 2;;
    *) echo "preflight: unknown arg $1" >&2; exit 2;;
  esac
done
fail=0
say(){ printf 'preflight: %s\n' "$*"; }
bad(){ printf 'PREFLIGHT FAIL: %s\n' "$*" >&2; fail=1; }

# 1. binary must exist and be NEWER than the commit under test.
#    A growth test once ran on a binary built 75 min before the allocator commit.
if [ -n "$exe" ]; then
  [ -x "$exe" ] || bad "binary not executable: $exe"
  if [ -x "$exe" ] && [ -n "$commit" ]; then
    bt=$(stat -c %Y "$exe")
    ct=$(git -C "$(dirname "$0")/.." log -1 --format=%ct "$commit" 2>/dev/null)
    if [ -n "$ct" ]; then
      if [ "$bt" -lt "$ct" ]; then
        bad "binary predates $commit ($(date -d @$bt '+%m-%d %H:%M') < $(date -d @$ct '+%m-%d %H:%M')) — rebuild"
      else
        say "binary $(basename "$exe") built $(date -d @$bt '+%m-%d %H:%M'), after $commit. ok"
      fi
    fi
  fi
fi

# 2. namelist sanity: levelmin must match the IC resolution.
#    Lowering levelmin once aborted a run: the grafic IC is 256^3 (level 8).
if [ -n "$nml" ] && [ -f "$nml" ]; then
  lmin=$(grep -oE '^levelmin=[0-9]+' "$nml" | cut -d= -f2)
  icdir=$(grep -oE "^initfile\(1\)='[^']+'" "$nml" | cut -d"'" -f2)
  if [ -n "$lmin" ] && [ -n "$icdir" ]; then
    # Read the grafic header's n1 (first int of the first record) rather than
    # guessing from the directory name: icz256 and level_010 both carry digits
    # that are not the level.
    hdr=""
    for c in "$icdir/ic_deltab" "$icdir/level_$(printf '%03d' "$lmin")/ic_deltab"; do
      [ -f "$c" ] && { hdr=$c; break; }
    done
    if [ -n "$hdr" ]; then
      n1=$(python3 - "$hdr" <<'PY2'
import struct,sys
with open(sys.argv[1],'rb') as f:
    f.read(4)
    print(struct.unpack('<i', f.read(4))[0])
PY2
)
      want=$((1 << lmin))
      if [ "$n1" != "$want" ]; then
        bad "levelmin=$lmin expects a ${want}^3 IC but $hdr is ${n1}^3 — the run aborts on grid mismatch"
      else
        say "levelmin=$lmin matches the ${n1}^3 IC. ok"
      fi
    else
      say "levelmin=$lmin; no ic_deltab found under $icdir to verify against (skipped)"
    fi
  fi
  # 3. switches the test depends on must actually be set.
  #    A CUDA test aimed at gpu_poisson never ran it: the default is .false.
  for s in "${switches[@]:-}"; do
    [ -z "$s" ] && continue
    k=${s%%=*}; v=${s#*=}
    got=$(grep -oE "^${k}=[^ ]+" "$nml" | tail -1 | cut -d= -f2)
    if [ -z "$got" ]; then
      bad "$k not set in $nml but the test depends on it (defaults are NOT what you assume; check amr_parameters)"
    elif [ "$got" != "$v" ]; then
      bad "$k=$got in $nml, test needs $v"
    else
      say "$k=$got. ok"
    fi
  done
fi

# 4. GPU must actually be free.
#    A CUDA test read as a clean pass while falling back to CPU with
#    "out of memory" because 31.4 of 32.8 GB were already held.
if [ -n "$gpu_need" ]; then
  if command -v nvidia-smi >/dev/null 2>&1; then
    used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    tot=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
    free_gb=$(( (tot - used) / 1024 ))
    if [ "$free_gb" -lt "$gpu_need" ]; then
      bad "GPU has ${free_gb}GB free, test needs ${gpu_need}GB — kernels would silently fall back to CPU"
    else
      say "GPU ${free_gb}GB free. ok"
    fi
  else
    bad "--gpu-free-gb requested but nvidia-smi is absent"
  fi
fi

[ $fail -eq 0 ] && say "all checks passed" || printf 'preflight: ABORTING\n' >&2
exit $fail
