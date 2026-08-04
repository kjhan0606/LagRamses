#!/usr/bin/env python3
import sys
from pathlib import Path

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "utils" / "py"))
import ramses_io  # noqa: E402


if len(sys.argv) != 2:
    raise SystemExit("usage: check_spatial.py RUN_DIRECTORY")

run_directory = Path(sys.argv[1]).resolve()
outputs = sorted(run_directory.glob("output_*"))
if not outputs:
    raise SystemExit(f"no output directory in {run_directory}")
nout = int(outputs[-1].name.split("_")[-1])

old_directory = Path.cwd()
try:
    import os

    os.chdir(run_directory)
    amr = ramses_io.rd_amr(nout)
finally:
    os.chdir(old_directory)

boxlen = amr[0].boxlen
offset = np.array(
    [
        [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5],
        [-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5],
        [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5],
    ]
)

coordinates = []
levels = []
for level_number, level in enumerate(amr, start=1):
    dx = boxlen / 2**level_number
    for child in range(8):
        leaf = ~level.refined[child]
        if np.any(leaf):
            coordinates.append(
                level.xg[:, leaf] + offset[:, child, None] * dx
            )
            levels.append(np.full(np.count_nonzero(leaf), level_number))

xyz = np.concatenate(coordinates, axis=1)
level = np.concatenate(levels)
radius = np.sqrt(np.sum((xyz - 0.5) ** 2, axis=0))
core = level[radius < 0.15]
far = level[radius > 0.65]
counts = {
    int(key): int(value)
    for key, value in zip(*np.unique(level, return_counts=True))
}

assert core.size > 0 and core.min() == core.max() == 5
assert far.size > 0 and far.min() == far.max() == 3
print(
    "PASSED: "
    f"leaf counts={counts}, "
    f"core cells={core.size} at level 5, "
    f"far cells={far.size} at level 3"
)
