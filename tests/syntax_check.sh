#!/bin/bash
# Per-file compile check: catches undeclared names, missing USE imports, wrong
# rank -- the errors a gate run should never be the first to discover.
#
#   tests/syntax_check.sh <worktree> [file ...]
#
# No file list => Fortran files modified vs HEAD.  Modified modules are compiled
# first into a scratch dir so dependants see the NEW interface; everything else
# resolves against REFBIN, a bin/ from a complete build.
set -u
export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:/opt/ohpc/pub/intel/oneapi/compiler/2025.3/bin:$PATH
REFBIN="${REFBIN:-/home/kjhan/BACKUP/lagRamses/bin}"

WT="${1:?usage: syntax_check.sh <worktree> [file ...]}"; shift
cd "$WT" || exit 2

if [ $# -gt 0 ]; then FILES="$*"; else FILES=$(git diff --name-only HEAD -- '*.f90' '*.F90'); fi
[ -z "$FILES" ] && { echo "SYNTAX: nothing to check"; exit 0; }

command -v mpiifx >/dev/null || { echo "SYNTAX: SKIP (no mpiifx)"; exit 3; }
ls "$REFBIN"/*.mod >/dev/null 2>&1 || { echo "SYNTAX: SKIP (REFBIN $REFBIN has no .mod)"; exit 3; }

D=$(mktemp -d); trap 'rm -rf "$D"' EXIT
FF="-fpp -syntax-only -I patch/lagRamses -I $D -I $REFBIN -module $D
    -DNPRE=8 -DNDIM=3 -DNVAR=6 -DNENER=0 -DSOLVER=hydro
    -DLONGINT -DQUADHILBERT -DOUTPUT_PARTICLE_POTENTIAL"

# pass 1 -- publish fresh .mod for modified modules (repeat to settle ordering)
MODFILES=""
for f in $FILES; do
    [ -f "$f" ] && grep -qiE '^[[:space:]]*module[[:space:]]+[a-z_]' "$f" && MODFILES="$MODFILES $f"
done
for pass in 1 2 3; do
    for f in $MODFILES; do mpiifx $FF "$f" >/dev/null 2>&1; done
done

# pass 2 -- report
rc=0; n=0
for f in $FILES; do
    [ -f "$f" ] || continue
    n=$((n+1))
    out=$(mpiifx $FF "$f" 2>&1)
    if echo "$out" | grep -qiE '^.*error #'; then
        echo "--- $f"; echo "$out" | grep -iE 'error #' | head -10; rc=1
    fi
done
if [ $rc -eq 0 ]; then echo "SYNTAX: PASS ($n file(s))"; else echo "SYNTAX: FAIL"; fi
exit $rc
