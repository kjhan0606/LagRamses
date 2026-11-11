#!/usr/bin/env python3
"""Inspect raw GRAFIC gas contrast, passive field and refinement masks."""
import json
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

root, output = map(Path, sys.argv[1:])


def read_field(path):
    with FortranFile(path, "r") as f:
        header = f.read_record(np.uint8).tobytes()
        if len(header) != 48:
            raise ValueError(f"invalid header: {path}")
        nx, ny, nz, *_ = struct.unpack("<3i9f", header)
        if not all(0 < x <= 1024 for x in (nx, ny, nz)):
            raise ValueError(f"invalid dimensions: {path}")
        data = np.empty((nz, ny, nx), dtype=np.float32)
        for z in range(nz):
            data[z] = f.read_reals(np.float32).reshape(ny, nx)
    if not np.isfinite(data).all():
        raise ValueError(f"nonfinite values: {path}")
    return header, data


report = {"scope": "raw IC arrays and same-level refmap overlap; not an AMR leaf inventory", "levels": {}}
for level in (8, 9, 10, 11):
    directory = root / "ics" / f"level_{level:03d}"
    headers, arrays = zip(*(read_field(directory / name) for name in
                            ("ic_deltab", "ic_pvar_00001", "ic_refmap")))
    if len(set(headers)) != 1:
        raise ValueError(f"mismatched headers at level {level}")
    delta, tracer, refmap = arrays
    if not np.all((refmap == 0) | (refmap == 1)):
        raise ValueError(f"nonbinary refinement mask at level {level}")
    flagged = refmap > 0.5
    nonpositive = delta <= -1.0
    report["levels"][str(level)] = {
        "cells": int(delta.size),
        "delta_min": float(delta.min()),
        "delta_le_minus_one_count": int(np.count_nonzero(nonpositive)),
        "refmap_one_count": int(np.count_nonzero(flagged)),
        "delta_le_minus_one_in_refmap_one": int(np.count_nonzero(nonpositive & flagged)),
        "tracer_min": float(tracer.min()),
        "tracer_max": float(tracer.max()),
        "tracer_nonzero_count": int(np.count_nonzero(tracer)),
        "tracer_above_half_count": int(np.count_nonzero(tracer > 0.5)),
        "tracer_above_half_in_refmap_one": int(np.count_nonzero((tracer > 0.5) & flagged)),
    }
    del headers, arrays, delta, tracer, refmap, flagged, nonpositive
output.write_text(json.dumps(report, indent=2) + "\n")
