#!/usr/bin/env python3
"""Measure gas/tracer integrals and out-of-range tracer in paired AMR snapshots."""
import json
import os
from pathlib import Path

import h5py
import numpy as np

results = {}
for name, expected_a in (("EARLY", 0.0101), ("LATE", 0.02)):
    result = {"levels": {}}
    with h5py.File(os.environ[f"RAR_{name}_SNAPSHOT"], "r") as h:
        a = np.asarray(h["header"].attrs["aexp"]).item()
        if abs(a - expected_a) > 1e-6:
            raise ValueError(f"{name}: expected a={expected_a}, found {a}")
        result["scale_factor"] = a
        for level in range(8, 12):
            flag = h[f"amr/level_{level}/son_flag"]
            gas = h[f"hydro/level_{level}/uold_1"]
            tracer = h[f"hydro/level_{level}/uold_18"]
            if not len(flag) == len(gas) == len(tracer):
                raise ValueError("AMR and hydro dataset size mismatch")
            row = {"leaf_cells": 0, "gas_mass_code": 0., "tracer_mass_code": 0.,
                   "negative_excess_mass": 0., "upper_excess_mass": 0.,
                   "affected_gas_mass": 0., "out_of_bounds_cells": 0,
                   "beyond_1e-3_cells": 0}
            volume = 2. ** (-3 * level)
            for lo in range(0, len(flag), 500_000):
                hi = min(lo + 500_000, len(flag))
                f = np.asarray(flag[lo:hi])
                if not np.isin(f, [0, 1]).all():
                    raise ValueError("unexpected son flag")
                leaf = f == 0
                d = np.asarray(gas[lo:hi])[leaf]
                q = np.asarray(tracer[lo:hi])[leaf]
                if not np.isfinite(d).all() or not np.isfinite(q).all() or np.any(d <= 0):
                    raise ValueError("invalid leaf gas density or tracer")
                bad = (q < 0) | (q > d)
                row["leaf_cells"] += len(d)
                row["gas_mass_code"] += float(np.sum(d, dtype=np.float64)) * volume
                row["tracer_mass_code"] += float(np.sum(q, dtype=np.float64)) * volume
                row["negative_excess_mass"] += float(np.sum(np.maximum(-q, 0), dtype=np.float64)) * volume
                row["upper_excess_mass"] += float(np.sum(np.maximum(q-d, 0), dtype=np.float64)) * volume
                row["affected_gas_mass"] += float(np.sum(d[bad], dtype=np.float64)) * volume
                row["out_of_bounds_cells"] += int(np.count_nonzero(bad))
                row["beyond_1e-3_cells"] += int(np.count_nonzero((q < -0.001*d) | (q > 1.001*d)))
            result["levels"][str(level)] = row
    result["total_leaf_gas_mass_code"] = sum(r["gas_mass_code"] for r in result["levels"].values())
    result["total_leaf_tracer_mass_code"] = sum(r["tracer_mass_code"] for r in result["levels"].values())
    result["total_negative_excess_mass"] = sum(r["negative_excess_mass"] for r in result["levels"].values())
    result["total_upper_excess_mass"] = sum(r["upper_excess_mass"] for r in result["levels"].values())
    result["negative_excess_over_tracer_mass"] = result["total_negative_excess_mass"] / result["total_leaf_tracer_mass_code"]
    result["upper_excess_over_tracer_mass"] = result["total_upper_excess_mass"] / result["total_leaf_tracer_mass_code"]
    result["out_of_bounds_cells"] = sum(r["out_of_bounds_cells"] for r in result["levels"].values())
    result["beyond_1e-3_cells"] = sum(r["beyond_1e-3_cells"] for r in result["levels"].values())
    results[name.lower()] = result

early, late = results["early"], results["late"]
results["same_trajectory_comparison"] = {
    "gas_mass_relative_change": late["total_leaf_gas_mass_code"] / early["total_leaf_gas_mass_code"] - 1.,
    "tracer_mass_relative_change": late["total_leaf_tracer_mass_code"] / early["total_leaf_tracer_mass_code"] - 1.,
    "limits": "Same integration trajectory, but AMR leaf meshes differ between epochs. Global integrals test mass accounting only; they do not establish local mixing accuracy, correct AMR footprint, or force convergence."
}
results["limits"] = (
    "Excess masses are hypothetical clipping magnitudes; no clipping is performed. "
    "Gas mass in affected cells is not itself an error mass. This is not a physical metallicity, "
    "target-halo inventory, or full-physics validation.")
Path(os.environ["RAR_OUTPUT"], "paired_tracer_bounds.json").write_text(json.dumps(results, indent=2) + "\n")
