#!/usr/bin/env python3
"""Locate late-time level-11 tracer excursions relative to coarse/fine faces."""
import json
import os
from pathlib import Path

import h5py
import numpy as np

CHILD_BITS = np.asarray(
    [[(child >> axis) & 1 for axis in range(3)] for child in range(8)],
    dtype=np.int64,
)


def cell_keys(xyz, nside):
    xyz = np.asarray(xyz, dtype=np.int64)
    return (xyz[:, 0] * nside + xyz[:, 1]) * nside + xyz[:, 2]


def coordinate_indices(grid_centers, grid_level):
    """Return integer indices of parent cells whose child grids sit at grid_level."""
    parent_level = grid_level - 1
    nside = 2**parent_level
    lattice = grid_centers * nside - 0.5
    rounded = np.rint(lattice)
    if np.max(np.abs(lattice - rounded), initial=0.0) > 1e-7:
        raise ValueError(f"level {grid_level}: unexpected GRAFIC/RAMSES grid-centre lattice")
    return rounded.astype(np.int64) % nside


def level10_cell_lookup(handle):
    """Map all stored level-10 cells to their child-grid flags."""
    level = 10
    nside = 2**level
    group = handle[f"amr/level_{level}"]
    grids = np.column_stack([group[f"xg_{axis}"][:] for axis in (1, 2, 3)])
    parent = coordinate_indices(grids, level)
    son = np.asarray(group["son_flag"][:])
    if son.size != len(grids) * 8 or not np.isin(son, [0, 1]).all():
        raise ValueError("level 10: malformed son_flag array")
    cells = (2 * parent[:, None, :] + CHILD_BITS[None, :, :]).reshape(-1, 3) % nside
    keys = cell_keys(cells, nside)
    order = np.argsort(keys, kind="mergesort")
    sorted_keys = keys[order]
    if len(sorted_keys) > 1 and np.any(sorted_keys[1:] == sorted_keys[:-1]):
        raise ValueError("level 10: duplicate cell coordinates")
    return sorted_keys, son[order]


def lookup_status(sorted_keys, sorted_son, query_keys):
    pos = np.searchsorted(sorted_keys, query_keys)
    found = pos < len(sorted_keys)
    safe = np.minimum(pos, max(len(sorted_keys) - 1, 0))
    if len(sorted_keys):
        found &= sorted_keys[safe] == query_keys
    else:
        found[:] = False
    status = np.full(len(query_keys), -1, dtype=np.int8)
    status[found] = sorted_son[safe[found]]
    return status, found


def summarize_mask(mask, leaf, boundary_same_level, boundary_missing):
    selected = mask & leaf
    any_boundary = boundary_same_level | boundary_missing
    n = int(np.count_nonzero(selected))
    n_same = int(np.count_nonzero(selected & boundary_same_level))
    n_missing = int(np.count_nonzero(selected & boundary_missing))
    n_any = int(np.count_nonzero(selected & any_boundary))
    return {
        "leaf_cells": n,
        "touch_level10_leaf_face": n_same,
        "touch_missing_level10_face": n_missing,
        "touch_any_coarser_face": n_any,
        "interior_to_level10_refinement_boundary": n - n_any,
        "boundary_fraction": (n_any / n if n else None),
    }


def analyze(snapshot, bound_json, label):
    with h5py.File(snapshot, "r") as handle:
        aexp = float(np.asarray(handle["header"].attrs["aexp"]).item())
        if abs(aexp - 0.02) > 1e-6:
            raise ValueError(f"{label}: expected a=0.02, found {aexp}")

        parent_keys, parent_son = level10_cell_lookup(handle)
        group11 = handle["amr/level_11"]
        grids11 = np.column_stack([group11[f"xg_{axis}"][:] for axis in (1, 2, 3)])
        parents11 = coordinate_indices(grids11, 11)
        parent11_keys = cell_keys(parents11, 2**10)
        if len(np.unique(parent11_keys)) != len(parent11_keys):
            raise ValueError(f"{label}: duplicate level-11 parent cells")
        parent_status, parent_found = lookup_status(parent_keys, parent_son, parent11_keys)
        if not np.all(parent_found) or np.any(parent_status != 1):
            raise ValueError(f"{label}: level-11 grids do not map to refined level-10 cells")

        # Find which faces of each level-11 parent cell border an unrefined or
        # absent level-10 cell. Absent cells are coarser than level 10 here.
        coarse_face = np.zeros((len(parents11), 3, 2), dtype=bool)
        missing_face = np.zeros_like(coarse_face)
        for axis in range(3):
            for side, delta in enumerate((-1, 1)):
                neighbor = parents11.copy()
                neighbor[:, axis] = (neighbor[:, axis] + delta) % 1024
                neighbor_keys = cell_keys(neighbor, 1024)
                status, found = lookup_status(parent_keys, parent_son, neighbor_keys)
                coarse_face[:, axis, side] = found & (status == 0)
                missing_face[:, axis, side] = ~found

        touch_same = np.zeros((len(parents11), 8), dtype=bool)
        touch_missing = np.zeros_like(touch_same)
        for axis in range(3):
            touch_same |= coarse_face[:, axis, 0, None] & (CHILD_BITS[None, :, axis] == 0)
            touch_same |= coarse_face[:, axis, 1, None] & (CHILD_BITS[None, :, axis] == 1)
            touch_missing |= missing_face[:, axis, 0, None] & (CHILD_BITS[None, :, axis] == 0)
            touch_missing |= missing_face[:, axis, 1, None] & (CHILD_BITS[None, :, axis] == 1)
        boundary_same = touch_same.reshape(-1)
        boundary_missing = touch_missing.reshape(-1)

        son11 = np.asarray(group11["son_flag"][:])
        gas = np.asarray(handle["hydro/level_11/uold_1"][:], dtype=np.float64)
        tracer = np.asarray(handle["hydro/level_11/uold_18"][:], dtype=np.float64)
        if not (len(son11) == len(gas) == len(tracer) == len(grids11) * 8):
            raise ValueError(f"{label}: level-11 AMR/hydro dimensions disagree")
        if not np.isin(son11, [0, 1]).all() or np.any(son11 != 0):
            raise ValueError(f"{label}: level 11 is not the expected leaf-only maximum level")
        if not np.isfinite(gas).all() or not np.isfinite(tracer).all() or np.any(gas <= 0):
            raise ValueError(f"{label}: invalid level-11 gas/tracer values")

        ratio = tracer / gas
        leaf = son11 == 0
        any_outside = (ratio < 0.0) | (ratio > 1.0)
        severe = (ratio < -0.001) | (ratio > 1.001)
        threshold_band = (ratio >= 0.45) & (ratio <= 0.55)
        result = {
            "snapshot": str(snapshot),
            "scale_factor": aexp,
            "level11_grid_count": int(len(grids11)),
            "level11_leaf_cells": int(np.count_nonzero(leaf)),
            "level10_cells_with_child_count": int(np.count_nonzero(parent_son == 1)),
            "level11_parent_mapping": "Every level-11 grid centre maps to one level-10 son_flag=1 cell; periodic face neighbors are checked at level 10.",
            "coarser_face_definition": "A level-11 child cell touches a coarser face when its child bit is on that face and the neighboring level-10 cell is unrefined (son_flag=0) or absent from the stored level-10 mesh.",
            "all_level11_leaf_cells": summarize_mask(leaf, leaf, boundary_same, boundary_missing),
            "all_out_of_bounds": summarize_mask(any_outside, leaf, boundary_same, boundary_missing),
            "beyond_0.001": summarize_mask(severe, leaf, boundary_same, boundary_missing),
            "ratio_0.45_to_0.55": summarize_mask(threshold_band, leaf, boundary_same, boundary_missing),
        }

    expected = json.loads(Path(bound_json).read_text())["late"]["levels"]["11"]
    if result["all_out_of_bounds"]["leaf_cells"] != expected["out_of_bounds_cells"]:
        raise ValueError(f"{label}: out-of-bounds count disagrees with saved bounds audit")
    if result["beyond_0.001"]["leaf_cells"] != expected["beyond_1e-3_cells"]:
        raise ValueError(f"{label}: severe-excursion count disagrees with saved bounds audit")
    return result


adaptive = analyze(
    os.environ["RAR_ADAPTIVE_SNAPSHOT"],
    os.environ["RAR_ADAPTIVE_BOUNDS_JSON"],
    "adaptive",
)
fixed = analyze(
    os.environ["RAR_FIXED_SNAPSHOT"],
    os.environ["RAR_FIXED_BOUNDS_JSON"],
    "fixed_map",
)
result = {
    "adaptive": adaptive,
    "fixed_map": fixed,
    "limits": (
        "This is geometric face adjacency from saved AMR cells, not a force or transport-error measurement. "
        "The fixed-IC-map run still evolves AMR support. Spatial co-location can identify a boundary association, "
        "but cannot by itself distinguish hydrodynamic advection from restriction/prolongation."
    ),
}
Path(os.environ["RAR_OUTPUT"], "refinement_boundary_excursions.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
