#!/usr/bin/env python3
"""Map finest-mass particle positions to the deepest actual AMR grid level."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

LEVELS = tuple(range(8, 12))
CHILD_BITS = np.asarray(
    [[(child >> axis) & 1 for axis in range(3)] for child in range(8)],
    dtype=np.int64,
)
MASS_FINE = (1.0 - 0.0493017 / 0.315) * 2.0**-33
EXPECTED_A = 0.02


def cell_keys(xyz, nside):
    xyz = np.asarray(xyz, dtype=np.int64)
    return (xyz[:, 0] * nside + xyz[:, 1]) * nside + xyz[:, 2]


def grid_parent_coordinates(grid_centres, level):
    nside = 2 ** (level - 1)
    lattice = grid_centres * nside - 0.5
    rounded = np.rint(lattice)
    if np.max(np.abs(lattice - rounded), initial=0.0) > 1e-7:
        raise ValueError(f"level {level}: unexpected grid-centre convention")
    return rounded.astype(np.int64) % nside


def read_grid_keys(handle, level):
    group = handle[f"amr/level_{level}"]
    centres = np.column_stack([group[f"xg_{axis}"][:] for axis in (1, 2, 3)])
    parent = grid_parent_coordinates(centres, level)
    keys = cell_keys(parent, 2 ** (level - 1))
    if len(np.unique(keys)) != len(keys):
        raise ValueError(f"level {level}: duplicate grid-parent coordinates")
    return np.sort(keys), centres


def verify_nested_mesh(handle, level):
    """Verify every level-l grid has a refined parent cell at level l-1."""
    parent_level = level - 1
    group = handle[f"amr/level_{parent_level}"]
    centres = np.column_stack([group[f"xg_{axis}"][:] for axis in (1, 2, 3)])
    parent = grid_parent_coordinates(centres, parent_level)
    nside = 2**parent_level
    cells = (2 * parent[:, None, :] + CHILD_BITS[None, :, :]).reshape(-1, 3) % nside
    cell_key_array = cell_keys(cells, nside)
    son = np.asarray(group["son_flag"][:]).reshape(-1)
    if len(son) != len(cell_key_array) or not np.isin(son, (0, 1)).all():
        raise ValueError(f"level {parent_level}: malformed son_flag data")
    order = np.argsort(cell_key_array, kind="mergesort")
    sorted_cell_keys = cell_key_array[order]
    if len(sorted_cell_keys) > 1 and np.any(sorted_cell_keys[1:] == sorted_cell_keys[:-1]):
        raise ValueError(f"level {parent_level}: duplicate AMR cell coordinates")
    child_parent_keys, _ = read_grid_keys(handle, level)
    positions = np.searchsorted(sorted_cell_keys, child_parent_keys)
    found = positions < len(sorted_cell_keys)
    inside = np.zeros(len(positions), dtype=bool)
    inside[found] = sorted_cell_keys[positions[found]] == child_parent_keys[found]
    refined = np.zeros(len(positions), dtype=bool)
    refined[inside] = son[order[positions[inside]]] == 1
    if not np.all(refined):
        raise ValueError(
            f"level {level}: {np.count_nonzero(~refined)} grids lack a refined parent cell"
        )
    return int(len(child_parent_keys))


def analyze(label, snapshot, snapshot_sha256):
    with h5py.File(snapshot, "r") as handle:
        aexp = float(np.asarray(handle["header"].attrs["aexp"]).item())
        if abs(aexp - EXPECTED_A) > 1e-6:
            raise ValueError(f"{label}: expected a={EXPECTED_A}, found {aexp}")

        grid_keys = {}
        grid_counts = {}
        for level in LEVELS:
            grid_keys[level], _ = read_grid_keys(handle, level)
            grid_counts[level] = int(len(grid_keys[level]))
        base_nside = 2 ** (LEVELS[0] - 1)
        if grid_counts[LEVELS[0]] != base_nside**3:
            raise ValueError(f"{label}: level-8 grid blocks do not cover the full base mesh")
        nested_counts = {
            str(level): verify_nested_mesh(handle, level) for level in LEVELS[1:]
        }

        highest_son = np.asarray(handle[f"amr/level_{LEVELS[-1]}/son_flag"][:]).reshape(-1)
        if not np.isin(highest_son, (0, 1)).all() or np.any(highest_son != 0):
            raise ValueError(f"{label}: level 11 is not the leaf-only maximum level")

        particles = handle["particles"]
        if not all(name in particles for name in ("mass", "levelp", "x_1", "x_2", "x_3")):
            raise ValueError(f"{label}: particle positions, masses, or levelp are missing")
        counts = np.zeros((len(LEVELS), len(LEVELS)), dtype=np.int64)
        spatial_counts = np.zeros(len(LEVELS), dtype=np.int64)
        integration_counts = np.zeros(len(LEVELS), dtype=np.int64)
        finest_particles = 0
        total_particles = len(particles["mass"])
        for start in range(0, total_particles, 500_000):
            end = min(start + 500_000, total_particles)
            mass = np.asarray(particles["mass"][start:end])
            integration_level = np.asarray(particles["levelp"][start:end])
            fine = np.isclose(mass, MASS_FINE, rtol=1e-5, atol=0.0)
            if not fine.any():
                continue

            xyz = np.column_stack(
                [particles[f"x_{axis}"][start:end] for axis in (1, 2, 3)]
            )[fine]
            integration_level = integration_level[fine]
            if not np.isfinite(xyz).all() or np.any((xyz < 0.0) | (xyz >= 1.0)):
                raise ValueError(f"{label}: invalid periodic particle coordinates")
            if np.any((integration_level < LEVELS[0]) | (integration_level > LEVELS[-1])):
                raise ValueError(f"{label}: finest-mass integration levels outside 8--11")

            spatial_level = np.zeros(len(xyz), dtype=np.int8)
            for level in LEVELS:
                nside = 2 ** (level - 1)
                parent_cell = np.floor(xyz * nside).astype(np.int64)
                query = cell_keys(parent_cell, nside)
                keys = grid_keys[level]
                locations = np.searchsorted(keys, query)
                valid = locations < len(keys)
                inside = np.zeros(len(locations), dtype=bool)
                inside[valid] = keys[locations[valid]] == query[valid]
                spatial_level[inside] = level

            if np.any(spatial_level < LEVELS[0]):
                raise ValueError(f"{label}: particles fall outside the level-8 base mesh")
            finest_particles += len(xyz)
            for idx, level in enumerate(LEVELS):
                spatial_counts[idx] += np.count_nonzero(spatial_level == level)
                integration_counts[idx] += np.count_nonzero(integration_level == level)
            for integration in LEVELS:
                selected = integration_level == integration
                if selected.any():
                    hist = np.bincount(
                        spatial_level[selected] - LEVELS[0], minlength=len(LEVELS)
                    )
                    counts[integration - LEVELS[0], :] += hist

    if finest_particles != 6_272_696:
        raise ValueError(f"{label}: expected 6,272,696 finest-mass particles, found {finest_particles}")
    if int(counts.sum()) != finest_particles:
        raise ValueError(f"{label}: incomplete integration/spatial-level cross-tabulation")

    spatial_level_counts = {str(level): int(spatial_counts[i]) for i, level in enumerate(LEVELS)}
    integration_level_counts = {
        str(level): int(integration_counts[i]) for i, level in enumerate(LEVELS)
    }
    cross_tab = {
        str(integration): {
            str(spatial): int(counts[integration - LEVELS[0], spatial - LEVELS[0]])
            for spatial in LEVELS
        }
        for integration in LEVELS
    }
    return {
        "label": label,
        "snapshot": str(snapshot),
        "snapshot_sha256": snapshot_sha256,
        "scale_factor": aexp,
        "finest_mass_particles": finest_particles,
        "level_grid_block_counts": {str(k): v for k, v in grid_counts.items()},
        "nested_child_grid_counts_by_level": nested_counts,
        "saved_integration_level_counts": integration_level_counts,
        "spatial_refinement_level_counts_from_grid_geometry": spatial_level_counts,
        "integration_level_by_spatial_level": cross_tab,
        "fraction_of_finest_mass_particles_in_level11_spatial_mesh": (
            spatial_counts[-1] / finest_particles
        ),
        "interpretation": (
            "Spatial level is the deepest AMR grid whose parent cell contains the particle. "
            "The saved levelp field is a particle integration/tree level and is cross-tabulated "
            "separately; it is not used to assign spatial resolution."
        ),
        "limits": (
            "This is a whole-domain mesh-resolution inventory. It does not identify final-halo "
            "membership, low-resolution contamination inside R200c, force error, or the final "
            "one-million-CDM-particle target. The adaptive and fixed-map runs have different "
            "mesh histories, so their counts are descriptive rather than a matched convergence test."
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
        "status": "Z49_FINEST_MASS_SPATIAL_LEVEL_REVIEW_REQUIRED",
        "cases": [
            analyze("adaptive", args.adaptive, args.adaptive_sha256),
            analyze("fixed_map", args.fixed_map, args.fixed_map_sha256),
        ],
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
