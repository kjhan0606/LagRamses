#!/usr/bin/env python3
"""Measure signed tracer-bound excursions in stored AMR leaves, without modifying them."""
import json
import os
from pathlib import Path

import h5py
import numpy as np

results = {}
for name in ("EARLY", "LATE"):
    result = {"levels": {}}
    with h5py.File(os.environ[f"RAR_{name}_SNAPSHOT"], "r") as h:
        result["scale_factor"] = float(h["header"].attrs["aexp"])
        for level in range(8, 12):
            flag = h[f"amr/level_{level}/son_flag"]
            rho = h[f"hydro/level_{level}/uold_1"]
            tracer = h[f"hydro/level_{level}/uold_18"]
            if not len(flag) == len(rho) == len(tracer):
                raise ValueError("dataset length mismatch")
            row = {"negative_excess_mass": 0., "upper_excess_mass": 0.,
                   "total_tracer_mass": 0., "affected_gas_mass": 0.,
                   "out_of_bounds_cells": 0, "beyond_1e-3_cells": 0}
            volume = 2. ** (-3 * level)
            for lo in range(0, len(flag), 500_000):
                hi = min(lo + 500_000, len(flag))
                f = np.asarray(flag[lo:hi])
                if not np.isin(f, [0, 1]).all():
                    raise ValueError("unexpected son flag")
                leaf = f == 0
                d = np.asarray(rho[lo:hi])[leaf]
                q = np.asarray(tracer[lo:hi])[leaf]
                if not np.isfinite(d).all() or not np.isfinite(q).all() or np.any(d <= 0):
                    raise ValueError("invalid data")
                bad = (q < 0) | (q > d)
                row["negative_excess_mass"] += float(np.sum(np.maximum(-q, 0), dtype=np.float64)) * volume
                row["upper_excess_mass"] += float(np.sum(np.maximum(q-d, 0), dtype=np.float64)) * volume
                row["total_tracer_mass"] += float(np.sum(q, dtype=np.float64)) * volume
                row["affected_gas_mass"] += float(np.sum(d[bad], dtype=np.float64)) * volume
                row["out_of_bounds_cells"] += int(np.count_nonzero(bad))
                row["beyond_1e-3_cells"] += int(np.count_nonzero((q < -0.001*d) | (q > 1.001*d)))
            result["levels"][str(level)] = row
    totals = {key: sum(row[key] for row in result["levels"].values()) for key in row}
    for key in ("negative_excess_mass", "upper_excess_mass", "affected_gas_mass"):
        totals[key + "_over_tracer_mass"] = totals[key] / totals["total_tracer_mass"]
    result["totals"] = totals
    results[name.lower()] = result
results["limits"] = ("Excess masses describe hypothetical clipping magnitudes only; no clipping is performed. "
                     "Gas in affected cells is not itself an error mass. Separate runs have different output and timestep schedules. "
                     "This is not a force, metallicity, halo membership, or production validation.")
Path(os.environ["RAR_OUTPUT"], "tracer_bounds.json").write_text(json.dumps(results, indent=2) + "\n")
