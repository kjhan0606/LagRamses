#!/usr/bin/env python3
import re
import sys
from pathlib import Path


LEVEL_RE = re.compile(r"^\s*Level\s+(\d+)\s+has\s+(\d+)\s+grids", re.MULTILINE)


def last_mesh(log_path):
    text = Path(log_path).read_text()
    meshes = text.split("Mesh structure")
    matches = LEVEL_RE.findall(meshes[-1] if len(meshes) > 1 else text)
    return {int(level): int(count) for level, count in matches}, text


if len(sys.argv) != 5:
    raise SystemExit("usage: check_logs.py OFF.log ON1.log ON4.log INVALID.log")

off, _ = last_mesh(sys.argv[1])
on1, on1_text = last_mesh(sys.argv[2])
on4, on4_text = last_mesh(sys.argv[3])
invalid, invalid_text = last_mesh(sys.argv[4])

assert max(off) == 3, f"disabled run refined unexpectedly: {off}"
assert max(on1) == 5 and on1[5] > 0, f"serial floor missing: {on1}"
assert max(on4) == 5 and on4[5] > 0, f"MPI floor missing: {on4}"
assert on1 == on4, f"MPI decomposition changed mesh counts: {on1} != {on4}"
assert "Void refinement floor enabled at level" in on1_text
assert "Void refinement floor enabled at level" in on4_text
assert "void refinement requires r_refine>0 at level" in invalid_text
assert "Aborting" in invalid_text
assert not invalid, f"invalid input advanced to a mesh: {invalid}"

print(f"PASSED: off={off}, on-1rank={on1}, on-4rank={on4}, invalid input rejected")
