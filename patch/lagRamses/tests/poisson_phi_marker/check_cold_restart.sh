#!/bin/bash
# Run after fresh.nml has generated output_00002. Preserve all existing dumps.
set -euo pipefail
if (( $# != 3 )); then
  echo "Usage: $0 BINARY CHECKPOINT_00002 NEW_RUN_DIRECTORY" >&2
  exit 2
fi
binary=$(realpath "$1")
checkpoint=$(realpath "$2")
run=$3
test -x "$binary"
test -f "$checkpoint/COMPLETE"
test ! -e "$run"
mkdir "$run"
run=$(realpath "$run")
here=$(cd "$(dirname "$0")" && pwd)
cp "$here/restart.nml" "$run/restart.nml"
ln -s "$checkpoint" "$run/output_00002"
export OMP_NUM_THREADS=1
cd "$run"
mpirun -np "${TEST_NCPU:-2}" "$binary" restart.nml > restart.log 2>&1
grep -q 'restart_phi_warm_start is deprecated and ignored' restart.log
grep -q 'Restart phi policy = cold predictor' restart.log
grep -q 'Run completed' restart.log
if grep -q 'Poisson warm start from restored phi' restart.log; then
  echo 'Unsafe saved-phi warm initialization was used' >&2
  exit 1
fi
if grep -q 'WARN: Fine multigrid Poisson failed to converge' restart.log; then
  echo 'MG did not converge' >&2
  exit 1
fi
# A fresh output must still include the gravity payload.
find "$run" -mindepth 2 -type f -name 'grav_*.out*' -print -quit | grep -q .
echo 'PASS: legacy warm request uses cold restart and preserves gravity output'
