#!/usr/bin/env python3
"""Cross-check finest-mass particle levels against z=49 level-11 AMR boxes."""
import argparse
import json
import os
from pathlib import Path

import h5py
import numpy as np

N_SIDE = 1024  # each level-11 AMR grid spans a 2^3-cell box at level 11
MASS_FINE = (1.0 - 0.0493017 / 0.315) * 2.0**-33
EXPECTED_A = 0.02


def cell_keys(xyz):
    xyz = np.asarray(xyz, dtype=np.int64)
    return (xyz[:, 0] * N_SIDE + xyz[:, 1]) * N_SIDE + xyz[:, 2]


def analyze(label, snapshot, snapshot_sha256):
    with h5py.File(snapshot, "r") as handle:
        aexp = float(np.asarray(handle["header"].attrs["aexp"]).item())
        if abs(aexp - EXPECTED_A) > 1e-6:
            raise ValueError(f"{label}: expected a={EXPECTED_A}, found {aexp}")

        grid = handle["amr/level_11"]
        centres = np.column_stack([grid[f"xg_{axis}"][:] for axis in (1, 2, 3)])
        lattice = centres * N_SIDE - 0.5
        if np.max(np.abs(lattice - np.rint(lattice)), initial=0.0) > 1e-7:
            raise ValueError(f"{label}: unexpected level-11 grid-centre convention")
        grid_xyz = np.rint(lattice).astype(np.int64) % N_SIDE
        grid_keys = np.sort(cell_keys(grid_xyz))
        if len(np.unique(grid_keys)) != len(grid_keys):
            raise ValueError(f"{label}: duplicate level-11 grid boxes")

        particles = handle["particles"]
        counts = {
            str(level): {"total": 0, "inside_level11_box": 0, "outside_level11_box": 0}
            for level in (10, 11)
        }
        counts["other"] = {
            "total": 0, "inside_level11_box": 0, "outside_level11_box": 0
        }
        total_particles = len(particles["mass"])
        for start in range(0, total_particles, 500_000):
            end = min(start + 500_000, total_particles)
            mass = np.asarray(particles["mass"][start:end])
            fine = np.isclose(mass, MASS_FINE, rtol=1e-5, atol=0.0)
            if not fine.any():
                continue

            level = np.asarray(particles["levelp"][start:end])[fine]
            xyz = np.column_stack(
                [particles[f"x_{axis}"][start:end] for axis in (1, 2, 3)]
            )[fine]
            if not np.isfinite(xyz).all() or np.any((xyz < 0.0) | (xyz >= 1.0)):
                raise ValueError(f"{label}: invalid periodic particle coordinates")
            particle_keys = cell_keys(np.floor(xyz * N_SIDE).astype(np.int64))
            locations = np.searchsorted(grid_keys, particle_keys)
            valid = locations < len(grid_keys)
            inside = np.zeros(len(locations), dtype=bool)
            inside[valid] = grid_keys[locations[valid]] == particle_keys[valid]

            for lev in (10, 11):
                selected = level == lev
                row = counts[str(lev)]
                row["total"] += int(np.count_nonzero(selected))
                row["inside_level11_box"] += int(np.count_nonzero(selected & inside))
                row["outside_level11_box"] += int(np.count_nonzero(selected & ~inside))
            other = (level != 10) & (level != 11)
            row = counts["other"]
            row["total"] += int(np.count_nonzero(other))
            row["inside_level11_box"] += int(np.count_nonzero(other & inside))
            row["outside_level11_box"] += int(np.count_nonzero(other & ~inside))

    fine_total = sum(row["total"] for row in counts.values())
    if fine_total == 0:
        raise ValueError(f"{label}: no particles at the expected finest mass")
    for level, row in counts.items():
        if row["total"] != row["inside_level11_box"] + row["outside_level11_box"]:
            raise ValueError(f"{label}: inconsistent particle cross-tabulation for level {level}")
        row["inside_fraction"] = (
            row["inside_level11_box"] / row["total"] if row["total"] else None
        )

    return {
        "label": label,
        "snapshot": str(snapshot),
        "snapshot_sha256": snapshot_sha256,
        "scale_factor": aexp,
        "level11_grid_boxes": int(len(grid_keys)),
        "finest_mass_code": MASS_FINE,
        "finest_mass_particles": fine_total,
        "counts_by_particle_level": counts,
        "level11_labelled_outside_fraction": counts["11"]["outside_level11_box"]
        / counts["11"]["total"] if counts["11"]["total"] else None,
        "limits": (
            "Level-11-box containment is a geometric mesh check. It does not test the "
            "full CIC force stencil, force error, target-halo membership, contamination, "
            "or late-time galaxy resolution. The adaptive and fixed-map snapshots have "
            "different meshes and timestep histories."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--adaptive", type=Path, required=True)
    parser.add_argument("--fixed-map", type=Path, required=True)
    parser.add_argument("--adaptive-sha256", required=True)
    parser.add_argument("--fixed-map-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result = {
        "status": "Z49_LEVEL11_PARTICLE_BOX_COVERAGE_REVIEW_REQUIRED",
        "cases": [
            analyze("adaptive", args.adaptive, args.adaptive_sha256),
            analyze("fixed_map", args.fixed_map, args.fixed_map_sha256),
        ],
        "comparison_limit": (
            "This tests the actual saved a=0.02 meshes against particle positions; it is "
            "not a matched force-convergence experiment and does not establish the "
            "million-particle count inside the z=0 target halo."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
