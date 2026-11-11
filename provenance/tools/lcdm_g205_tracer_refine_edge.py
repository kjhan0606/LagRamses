#!/usr/bin/env python3
"""Measure gas-tracer fractions around the configured AMR mask threshold."""
import json
import os
from pathlib import Path

import h5py
import numpy as np

threshold = 0.5
bands = ((0.45, 0.50), (0.50, 0.55), (0.45, 0.55))
results = {}
for epoch in ("EARLY", "LATE"):
    snapshot = Path(os.environ[f"RAR_{epoch}_SNAPSHOT"])
    result = {"levels": {}}
    with h5py.File(snapshot, "r") as h:
        result["scale_factor"] = float(h["header"].attrs["aexp"])
        for level in range(8, 12):
            flags = h[f"amr/level_{level}/son_flag"]
            gas = h[f"hydro/level_{level}/uold_1"]
            tracer = h[f"hydro/level_{level}/uold_18"]
            if not len(flags) == len(gas) == len(tracer):
                raise ValueError("AMR and hydro dataset size mismatch")
            volume = 2.0 ** (-3 * level)
            row = {scope: {"cells": 0, "ratio_gt_0_5_cells": 0,
                           "ratio_gt_0_5_gas_mass_code": 0.0,
                           "gas_mass_code": 0.0,
                           "bands": {f"{lo:.2f}_{hi:.2f}":
                                     {"cells": 0, "gas_mass_code": 0.0}
                                     for lo, hi in bands}}
                   for scope in ("all", "leaf")}
            for lo in range(0, len(flags), 500_000):
                hi = min(lo + 500_000, len(flags))
                f = np.asarray(flags[lo:hi])
                d = np.asarray(gas[lo:hi])
                q = np.asarray(tracer[lo:hi])
                if not np.isin(f, [0, 1]).all() or not np.isfinite(d).all() \
                        or not np.isfinite(q).all() or np.any(d <= 0):
                    raise ValueError("invalid flag, gas density, or tracer data")
                ratio = q / d
                for scope, mask in (("all", np.ones(len(f), dtype=bool)),
                                    ("leaf", f == 0)):
                    counts = row[scope]
                    r = ratio[mask]
                    rho = d[mask]
                    counts["cells"] += int(mask.sum())
                    counts["gas_mass_code"] += float(np.sum(rho, dtype=np.float64)) * volume
                    selected = r > threshold
                    counts["ratio_gt_0_5_cells"] += int(selected.sum())
                    counts["ratio_gt_0_5_gas_mass_code"] += \
                        float(np.sum(rho[selected], dtype=np.float64)) * volume
                    for band_lo, band_hi in bands:
                        band = (r >= band_lo) & (r < band_hi)
                        label = f"{band_lo:.2f}_{band_hi:.2f}"
                        counts["bands"][label]["cells"] += int(band.sum())
                        counts["bands"][label]["gas_mass_code"] += \
                            float(np.sum(rho[band], dtype=np.float64)) * volume
            for scope, counts in row.items():
                counts["fraction_cells_ratio_gt_0_5"] = \
                    counts["ratio_gt_0_5_cells"] / counts["cells"]
                counts["fraction_gas_mass_ratio_gt_0_5"] = \
                    counts["ratio_gt_0_5_gas_mass_code"] / counts["gas_mass_code"]
                for band_row in counts["bands"].values():
                    band_row["fraction_gas_mass"] = \
                        band_row["gas_mass_code"] / counts["gas_mass_code"]
            result["levels"][str(level)] = row
    results[epoch.lower()] = result
results["criterion"] = {
    "ivar_refine": 18,
    "var_cut_refine": threshold,
    "mask": "uold_18/uold_1 > var_cut_refine",
    "bands": [[lo, hi] for lo, hi in bands],
}
results["limits"] = (
    "Per-level fractions summarize the saved mesh, not one-to-one Lagrangian cells. "
    "The actual density map also applies the particle/CIC density threshold m_refine=8; "
    "these tracer-only statistics do not reconstruct the final AMR decision map. "
    "Separate snapshots have different time and AMR histories, so differences are descriptive.")
Path(os.environ["RAR_OUTPUT"], "tracer_refine_edge.json").write_text(
    json.dumps(results, indent=2) + "\n")
