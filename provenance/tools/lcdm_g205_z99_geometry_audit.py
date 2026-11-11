#!/usr/bin/env python3
"""Audit Group 205 multi-fluid GRAFIC2 IC geometry through level 11."""

import json
import re
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

root, parent, trace, output = map(Path, sys.argv[1:5])
levels = (8, 9, 10, 11)
fields = (
    "ic_refmap",
    "ic_deltab",
    "ic_velbx",
    "ic_velby",
    "ic_velbz",
    "ic_poscx",
    "ic_poscy",
    "ic_poscz",
    "ic_velcx",
    "ic_velcy",
    "ic_velcz",
    "ic_pvar_00001",
)
expected_cosmology = np.array([0.01, 0.315, 0.685, 67.36, 0.0493017])
cosmology_atol = np.array([1.0e-6, 2.0e-6, 2.0e-6, 1.0e-3, 2.0e-6])


def read_grafic(path, expected=expected_cosmology):
    with FortranFile(path, "r") as stream:
        raw = stream.read_record(np.uint8).tobytes()
        if len(raw) != 48:
            raise ValueError(f"Unexpected GRAFIC header size in {path}")
        header = struct.unpack("<3i9f", raw)
        nx, ny, nz = header[:3]
        if not (0 < nx <= 1024 and 0 < ny <= 1024 and 0 < nz <= 1024):
            raise ValueError(f"Unexpected GRAFIC dimensions in {path}")
        values = np.empty((nz, ny, nx), dtype=np.float32)
        for k in range(nz):
            values[k] = stream.read_reals(np.float32).reshape(ny, nx)
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite GRAFIC values in {path}")
    cosmology = np.asarray(header[7:12], dtype=np.float64)
    if not np.allclose(cosmology, expected, rtol=0.0, atol=cosmology_atol):
        raise ValueError(f"Unexpected cosmology header in {path}: {cosmology}")
    return header, values


summary = json.loads((trace / "summary.json").read_text())
candidate = next(row for row in summary["candidates"] if row["group_id"] == 205)
trace_record = next(
    row for row in candidate["traces"] if row["radius_multiple_R200c"] == 3.0
)
if trace_record["unique_parent_IC_cells"] != 4536:
    raise ValueError("Unexpected Group 205 3R200c trace size")
with np.load(trace / "candidate_traces.npz") as archive:
    cells = np.asarray(archive["g205_r3_cell_xyz"])
if cells.shape != (4536, 3) or cells.dtype != np.uint8:
    raise ValueError("Unexpected Group 205 Lagrangian cell array")
if np.unique(cells, axis=0).shape[0] != len(cells):
    raise ValueError("Duplicate cells in Group 205 Lagrangian trace")

points = np.loadtxt(root / "points.txt")
if points.shape != (len(cells), 3) or not np.isfinite(points).all():
    raise ValueError("Unexpected point-list shape or values")
point_cells = np.floor(points * 256).astype(np.int64)
if not np.array_equal(point_cells, cells.astype(np.int64)):
    raise ValueError("Point list does not map back to selected IC cells")

log = (root / "music.log").read_text()
match = re.search(
    r"Domain will be shifted by \(\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\)",
    log,
)
if match is None:
    raise RuntimeError("MUSIC log lacks an unambiguous integer domain shift")
shift = np.array([int(value) for value in match.groups()], dtype=np.int64)

parent_header, parent_delta = read_grafic(
    parent / "ics_lcdm/level_008/ic_deltab",
    expected=np.array([0.02, 0.315, 0.685, 67.36, 0.0493017])
)
if parent_header[:3] != (256, 256, 256):
    raise ValueError("Unexpected parent base-grid dimensions")

mask_stats = {}
field_stats = {}
for level in levels:
    directory = root / f"ics/level_{level:03d}"
    headers = {}
    for field in fields:
        path = directory / field
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing GRAFIC field {path}")
        header, values = read_grafic(path)
        headers[field] = header[:3]
        if field == "ic_refmap":
            if not np.all((values == 0) | (values == 1)):
                raise ValueError(f"Refinement mask is not binary at level {level}")
            active = int(np.count_nonzero(values > 0.5))
        if field in ("ic_deltab", "ic_velbx", "ic_velby", "ic_velbz"):
            values64 = values.astype(np.float64)
            field_stats[f"L{level}_{field}"] = {
                "minimum": float(values64.min()),
                "maximum": float(values64.max()),
                "mean": float(values64.mean()),
                "rms": float(np.sqrt(np.mean(values64**2))),
            }
        del values
    if len(set(headers.values())) != 1:
        raise ValueError(f"GRAFIC field dimensions differ at level {level}")
    mask_stats[str(level)] = {
        "shape_xyz": list(headers["ic_refmap"]),
        "active_cells": active,
    }
    if active <= 0:
        raise ValueError(f"Empty refinement mask at level {level}")

if mask_stats["8"]["shape_xyz"] != [256, 256, 256]:
    raise ValueError("Unexpected base refinement-mask dimensions")
active_counts = [mask_stats[str(level)]["active_cells"] for level in levels]
if any(next_count <= previous for previous, next_count in zip(active_counts, active_counts[1:])):
    raise ValueError("Nested refinement masks do not grow with level")

base_mask_header, base_mask = read_grafic(root / "ics/level_008/ic_refmap")
child_header, child_delta = read_grafic(root / "ics/level_008/ic_deltab")
if base_mask_header[:3] != parent_header[:3] or child_header[:3] != parent_header[:3]:
    raise ValueError("Base-grid dimensions differ from the parent IC")
index = np.floor(points * 256).astype(np.int64)
registered_index = (index + shift) % 256
covered = base_mask[
    registered_index[:, 2], registered_index[:, 1], registered_index[:, 0]
] > 0.5
if not covered.all():
    raise ValueError("A 3R200c trace point falls outside the base refinement mask")

registered_parent = np.roll(
    parent_delta,
    tuple(int(value) for value in shift[::-1]),
    axis=(0, 1, 2),
)
difference = child_delta.astype(np.float64) - registered_parent.astype(np.float64)
parent_rms = float(np.sqrt(np.mean(registered_parent.astype(np.float64) ** 2)))
child_rms = float(np.sqrt(np.mean(child_delta.astype(np.float64) ** 2)))
difference_rms = float(np.sqrt(np.mean(difference**2)))
correlation = float(
    np.corrcoef(child_delta.ravel(), registered_parent.ravel())[0, 1]
)

result = {
    "status": "NESTED_L11_LCDM_Z99_TWOFLUID_IC_GEOMETRY_PASS_NO_EVOLUTION",
    "candidate_group_id": 205,
    "candidate_selection": "provisional z=0 mass-screen candidate; not production approval",
    "model": "two-fluid Lambda-CDM IC geometry diagnostic; not validated for production evolution",
    "trace_radius_multiple_R200c": 3.0,
    "trace_unique_parent_cells": int(len(cells)),
    "provisional_parent_M200c_msun": float(candidate["provisional_M200c_msun"]),
    "domain_shift_base_cells_xyz": shift.tolist(),
    "GRAFIC_cosmology_header_expected": expected_cosmology.tolist(),
    "trace_points_covered_by_level8_mask": int(covered.sum()),
    "trace_point_coverage_fraction": float(covered.mean()),
    "refinement_masks": mask_stats,
    "field_statistics": field_stats,
    "base_density_parent_rms": parent_rms,
    "base_density_zoom_rms": child_rms,
    "base_density_difference_rms": difference_rms,
    "base_density_correlation_after_integer_registration": correlation,
    "limits": [
        "The parent density is at z=49 and the zoom density at z=99; their amplitude and difference are cross-epoch diagnostics, not a matched-time comparison.",
        "Geometry and finite-field integrity do not validate transfer-function normalization or velocity power spectra.",
        "The GRAFIC header and IC reader require separate RAMSES validation.",
        "No evolved particle mass count, contamination, relaxation, galaxy, feedback, or RAR result was measured.",
        "A million-particle target is a resolution estimate until the evolved R200c CDM count is measured.",
    ],
}
(output / "geometry_audit.json").write_text(
    json.dumps(result, indent=2, allow_nan=False) + "\n"
)
print(json.dumps(result, indent=2))
