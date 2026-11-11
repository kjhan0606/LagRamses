#!/usr/bin/env python3
"""Bound distances of z=49 level-11-labelled particles outside level-11 boxes."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree

N_SIDE = 1024
MASS_FINE = (1.0 - 0.0493017 / 0.315) * 2.0**-33
EXPECTED_A = 0.02


def cell_keys(xyz):
    return (xyz[:, 0] * N_SIDE + xyz[:, 1]) * N_SIDE + xyz[:, 2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    parser.add_argument("--expected-outside", type=int, default=8326)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with h5py.File(args.snapshot, "r") as handle:
        aexp = float(np.asarray(handle["header"].attrs["aexp"]).item())
        if abs(aexp - EXPECTED_A) > 1e-6:
            raise ValueError(f"expected a={EXPECTED_A}, found {aexp}")

        grid = handle["amr/level_11"]
        centres = np.column_stack([grid[f"xg_{axis}"][:] for axis in (1, 2, 3)])
        lattice = centres * N_SIDE - 0.5
        if np.max(np.abs(lattice - np.rint(lattice)), initial=0.0) > 1e-7:
            raise ValueError("unexpected level-11 grid-centre convention")
        grid_xyz = np.rint(lattice).astype(np.int64) % N_SIDE
        grid_keys = np.sort(cell_keys(grid_xyz))
        if len(np.unique(grid_keys)) != len(grid_keys):
            raise ValueError("duplicate level-11 grid boxes")

        particles = handle["particles"]
        outside = []
        total_particles = len(particles["mass"])
        for start in range(0, total_particles, 500_000):
            end = min(start + 500_000, total_particles)
            mass = np.asarray(particles["mass"][start:end])
            level = np.asarray(particles["levelp"][start:end])
            selected = np.isclose(mass, MASS_FINE, rtol=1e-5, atol=0.0) & (level == 11)
            if not selected.any():
                continue

            xyz = np.column_stack(
                [particles[f"x_{axis}"][start:end] for axis in (1, 2, 3)]
            )[selected]
            if not np.isfinite(xyz).all() or np.any((xyz < 0.0) | (xyz >= 1.0)):
                raise ValueError("invalid periodic particle coordinates")
            particle_keys = cell_keys(np.floor(xyz * N_SIDE).astype(np.int64))
            locations = np.searchsorted(grid_keys, particle_keys)
            valid = locations < len(grid_keys)
            inside = np.zeros(len(locations), dtype=bool)
            inside[valid] = grid_keys[locations[valid]] == particle_keys[valid]
            outside.append(xyz[~inside])

    outside = np.concatenate(outside, axis=0) if outside else np.empty((0, 3))
    if len(outside) != args.expected_outside:
        raise ValueError(
            f"expected {args.expected_outside} level-11 particles outside boxes, "
            f"found {len(outside)}"
        )
    if not len(outside):
        raise ValueError("no outside particles to analyze")

    tree = cKDTree(centres, boxsize=1.0)
    _, candidates = tree.query(outside, k=8, workers=2)
    delta = np.abs(outside[:, None, :] - centres[candidates])
    delta = np.minimum(delta, 1.0 - delta)
    distance = (
        np.sqrt(np.sum(np.maximum(delta - 0.5 / N_SIDE, 0.0) ** 2, axis=2)).min(axis=1)
        * 2048.0
    )
    thresholds = (0.01, 0.1, 0.5, 1.0, 2.0)
    result = {
        "status": "Z49_LEVEL11_PARTICLE_BOX_DISTANCE_REVIEW_REQUIRED",
        "snapshot": str(args.snapshot),
        "snapshot_sha256": args.snapshot_sha256,
        "scale_factor": aexp,
        "mismatched_level11_particles": int(len(outside)),
        "distance_units": "level-11 cell widths (25/2048 comoving Mpc/h)",
        "quantiles_upper_bound": {
            str(q): float(np.quantile(distance, q)) for q in (0.0, 0.5, 0.9, 0.99, 1.0)
        },
        "within_upper_bound_counts": {
            str(t): int(np.count_nonzero(distance <= t)) for t in thresholds
        },
        "method": (
            "Periodic distance to the closest among eight nearest level-11 grid centres, "
            "subtracting the grid-box half-width componentwise. These are upper bounds "
            "on distance to any level-11 box."
        ),
        "limits": (
            "This localizes geometric level/box discrepancies only. It cannot identify "
            "their cause, prove force-stencil coverage, measure gravity error, or establish "
            "the final target-halo particle count."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
