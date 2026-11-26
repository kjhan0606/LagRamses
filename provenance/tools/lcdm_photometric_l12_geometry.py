#!/usr/bin/env python3
"""Check traced-region coverage and the generated GRAFIC fields."""

import configparser
import json
import re
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.io import FortranFile


root = Path(sys.argv[1])
group = int(sys.argv[2])
expected_points = int(sys.argv[3])
candidate_dir = root / f"g{group}"
config = configparser.ConfigParser(interpolation=None)
if not config.read(candidate_dir / "music.conf"):
    raise ValueError("Missing effective MUSIC configuration")
levelmin = config.getint("setup", "levelmin")
levelmax = config.getint("setup", "levelmax")
if levelmin != 8 or levelmax < levelmin:
    raise ValueError("Unexpected level range")
points = np.loadtxt(candidate_dir / "points.txt")
if points.shape != (expected_points, 3) or not np.isfinite(points).all():
    raise ValueError("Unexpected point array")
if np.any(points < 0) or np.any(points >= 1):
    raise ValueError("Trace points must lie inside the periodic parent box")
cells = np.floor(points * 256).astype(np.int64)
if len(np.unique(cells, axis=0)) != len(cells):
    raise ValueError("Duplicate parent-grid trace cells")

log = (root / "music.log").read_text(errors="replace")
matches = re.findall(
    r"Domain will be shifted by \(\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\)", log
)
if len(matches) != 1:
    raise ValueError("MUSIC must report exactly one integer domain shift")
shift = np.asarray(matches[0], dtype=np.int64)
registered = (cells + shift) % 256

fields = (
    "ic_refmap", "ic_deltab", "ic_velbx", "ic_velby", "ic_velbz",
    "ic_poscx", "ic_poscy", "ic_poscz", "ic_velcx", "ic_velcy",
    "ic_velcz", "ic_pvar_00001",
)
expected_cosmology = np.array([0.02, 0.315, 0.685, 67.36, 0.0493017])
cosmology_atol = np.array([1.0e-6, 2.0e-6, 2.0e-6, 1.0e-3, 2.0e-6])


def read_header(stream, path):
    raw = stream.read_record(np.uint8).tobytes()
    if len(raw) != 48:
        raise ValueError(f"Unexpected GRAFIC header length in {path}")
    header = struct.unpack("<3i9f", raw)
    if any(dimension <= 0 or dimension > 4096 for dimension in header[:3]):
        raise ValueError(f"Unexpected GRAFIC dimensions in {path}")
    if not np.allclose(header[7:12], expected_cosmology, rtol=0, atol=cosmology_atol):
        raise ValueError(f"Unexpected cosmology in {path}")
    return header


def read_plane(stream, nx, ny, path):
    plane = stream.read_reals(np.float32)
    if plane.size != nx * ny or not np.isfinite(plane).all():
        raise ValueError(f"Invalid data plane in {path}")
    return plane.reshape(ny, nx)


levels = {}
for level in range(levelmin, levelmax + 1):
    directory = root / "ics" / f"level_{level:03d}"
    headers = {}
    stats = {}
    for field in fields:
        path = directory / field
        with FortranFile(path, "r") as stream:
            header = read_header(stream, path)
            headers[field] = header
            nx, ny, nz = header[:3]
            minimum, maximum = float("inf"), float("-inf")
            nonpositive_density = 0
            active = 0
            covered = 0
            for z in range(nz):
                plane = read_plane(stream, nx, ny, path)
                minimum = min(minimum, float(plane.min()))
                maximum = max(maximum, float(plane.max()))
                if field == "ic_refmap":
                    if not np.all((plane == 0) | (plane == 1)):
                        raise ValueError(f"Nonbinary refinement mask in {path}")
                    active += int(np.count_nonzero(plane))
                    if level == 8:
                        in_plane = registered[registered[:, 2] == z]
                        covered += int(np.count_nonzero(plane[in_plane[:, 1], in_plane[:, 0]]))
                elif field == "ic_deltab":
                    nonpositive_density += int(np.count_nonzero(plane <= -1))
            stats[field] = {"min": minimum, "max": maximum}
            if field == "ic_refmap":
                stats[field]["one_count"] = active
                if level == 8:
                    stats[field]["trace_cells_covered"] = covered
            if field == "ic_deltab":
                stats[field]["delta_le_minus_one_count"] = nonpositive_density
    if len(set(headers.values())) != 1:
        raise ValueError(f"Field headers differ at level {level}")
    if stats["ic_refmap"]["one_count"] == 0:
        raise ValueError(f"Empty refinement mask at level {level}")
    if level == 8:
        if headers["ic_refmap"][:3] != (256, 256, 256):
            raise ValueError("Unexpected base grid")
        if stats["ic_refmap"]["trace_cells_covered"] != len(cells):
            raise ValueError("Base refinement mask does not cover all trace points")
    levels[str(level)] = {"shape_xyz": list(headers["ic_refmap"][:3]), "fields": stats}

report = {
    "status": "GRAFIC_GEOMETRY_AND_FINITE_FIELDS_CHECKED",
    "candidate_group": group,
    "domain_shift_base_cells_xyz": shift.tolist(),
    "traced_parent_cells": len(cells),
    "levels": levels,
    "scope": (
        "The density counts cover raw arrays, including possible buffer or covered cells. "
        "They are not an AMR leaf inventory or proof of a physical hydro state. "
        "Power normalization, evolved halo membership and contamination require additional checks."
    ),
}
(root / "geometry_audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
print(json.dumps(report, indent=2, allow_nan=False))
