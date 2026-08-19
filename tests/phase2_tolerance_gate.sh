#!/bin/bash
# Conservation-tolerance gate for a block-layout chunk. The cell values are the
# same as legacy but stored in a different memory order, so full-cell reduction
# loops reorder FP sums and results drift at the rounding level (~1e-8). A real
# layout bug corrupts the mesh -> order-unity difference or NaN/crash. Threshold
# 1e-6 on the final-step physical energies cleanly separates the two.
#   phase2_tolerance_gate.sh <ref_run_dir> <cand_run_dir>
set -uo pipefail
ref=$1; cand=$2
python3 - "$ref" "$cand" <<'PY'
import sys, re
ref, cand = sys.argv[1], sys.argv[2]
def last(path):
    v=None
    for l in open(f"{path}/run.log", errors="replace"):
        m=re.search(r"Main step=\s*(\d+)\s+mcons=\s*(\S+)\s+econs=\s*(\S+)\s+epot=\s*(\S+)\s+ekin=\s*(\S+)\s+eint=\s*(\S+)", l)
        if m: v=m
    return v
a=last(ref); b=last(cand)
if a is None or b is None:
    print("GATE: FAIL (no Main step lines)"); sys.exit(1)
import os
na=len([d for d in os.listdir(ref) if d.startswith("output_")])
nc=len([d for d in os.listdir(cand) if d.startswith("output_")])
if na!=nc:
    print(f"GATE: FAIL (output count {na} vs {nc})"); sys.exit(1)
ok=True
print(f"{'qty':>6} {'ref':>15} {'cand':>15} {'rel':>10}")
for i,name in ((4,'epot'),(5,'ekin'),(6,'eint')):
    ra=float(a.group(i)); ca=float(b.group(i))
    rel=abs(ca/ra-1) if ra!=0 else abs(ca)
    print(f"{name:>6} {ra:15.6e} {ca:15.6e} {rel:10.2e}")
    if rel>1e-6: ok=False
print("GATE: PASS (within 1e-6, FP-reorder level)" if ok else "GATE: FAIL (exceeds 1e-6)")
sys.exit(0 if ok else 1)
PY
