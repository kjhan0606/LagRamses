"""Conservatively resample RAMSES HDF5 AMR leaves onto the fixed P4 mesh."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys

import h5py
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.snapshot import GridSpec, neutral_primordial_input, write_static_rt_input
from snrt_core.source_ledger import read_photon_source_ledger_csv
from snrt_core.thermal_atlas import read_thermal_atlas


_INFO_VALUE = re.compile(r"^\s*([A-Za-z_]+)\s*=\s*([^!#]+)")
_CHILD_BITS = np.asarray(
    tuple((child & 1, (child >> 1) & 1, (child >> 2) & 1) for child in range(8)), dtype=np.float64
)


def _read_info(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text().splitlines():
        match = _INFO_VALUE.match(line)
        if match is not None:
            try:
                values[match.group(1)] = float(match.group(2).strip().replace("D", "E"))
            except ValueError:
                continue
    required = {"aexp", "unit_l", "unit_d"}
    missing = required - set(values)
    if missing:
        raise ValueError(f"RAMSES info file missing values: {sorted(missing)}")
    return values


def _read_grid_blocks(dataset: h5py.Dataset, grid_indices: np.ndarray) -> np.ndarray:
    """Read eight RAMSES cell values per selected grid without a full field read."""

    result = np.empty((len(grid_indices), 8), dtype=np.float64)
    starts = np.r_[0, np.flatnonzero(np.diff(grid_indices) != 1) + 1]
    stops = np.r_[starts[1:], len(grid_indices)]
    for begin, end in zip(starts, stops, strict=True):
        first_grid = int(grid_indices[begin])
        last_grid = int(grid_indices[end - 1]) + 1
        result[begin:end] = np.asarray(dataset[first_grid * 8 : last_grid * 8]).reshape((-1, 8))
    return result


def _candidate_grid_batches(
    amr: h5py.Group,
    *,
    cell_width_code: float,
    left_edge_code: np.ndarray,
    right_edge_code: np.ndarray,
    chunk_grids: int,
):
    """Yield local AMR grids whose eight-cell boxes intersect the P4 cube."""

    count = len(amr["xg_1"])
    for start in range(0, count, chunk_grids):
        stop = min(start + chunk_grids, count)
        centers = np.column_stack(
            (
                np.asarray(amr["xg_1"][start:stop]),
                np.asarray(amr["xg_2"][start:stop]),
                np.asarray(amr["xg_3"][start:stop]),
            )
        )
        intersects = np.all(centers + cell_width_code > left_edge_code, axis=1)
        intersects &= np.all(centers - cell_width_code < right_edge_code, axis=1)
        local = np.flatnonzero(intersects)
        if len(local):
            yield start + local, centers[local]


def _deposit_leaf(
    density_sum: np.ndarray,
    coverage: np.ndarray,
    *,
    center_code: np.ndarray,
    density_code: float,
    leaf_width_code: float,
    left_edge_code: np.ndarray,
    analysis_cell_width_code: float,
) -> None:
    """Volume-average one AMR leaf cell into every overlapping analysis cell."""

    shape = np.asarray(density_sum.shape, dtype=np.int64)
    cell_left = center_code - 0.5 * leaf_width_code
    cell_right = center_code + 0.5 * leaf_width_code
    start = np.floor((cell_left - left_edge_code) / analysis_cell_width_code + 1.0e-10).astype(np.int64)
    stop = np.ceil((cell_right - left_edge_code) / analysis_cell_width_code - 1.0e-10).astype(np.int64)
    start = np.maximum(start, 0)
    stop = np.minimum(stop, shape)
    if np.any(stop <= start):
        return
    for index_x in range(start[0], stop[0]):
        for index_y in range(start[1], stop[1]):
            for index_z in range(start[2], stop[2]):
                index = np.asarray((index_x, index_y, index_z))
                analysis_left = left_edge_code + index * analysis_cell_width_code
                analysis_right = analysis_left + analysis_cell_width_code
                overlap = np.minimum(cell_right, analysis_right) - np.maximum(cell_left, analysis_left)
                if np.all(overlap > 0.0):
                    volume_fraction = float(np.prod(overlap / analysis_cell_width_code))
                    density_sum[index_x, index_y, index_z] += density_code * volume_fraction
                    coverage[index_x, index_y, index_z] += volume_fraction


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", required=True)
    parser.add_argument("--info", required=True)
    parser.add_argument("--zoom-manifest", required=True)
    parser.add_argument("--thermal-atlas", required=True)
    parser.add_argument("--source-ledger", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata-output", required=True)
    parser.add_argument("--analysis-level", type=int, default=15)
    parser.add_argument("--scan-chunk-grids", type=int, default=1_000_000)
    parser.add_argument("--metallicity-solar", type=float, default=1.0e-6)
    args = parser.parse_args()
    if args.scan_chunk_grids < 1:
        raise ValueError("scan-chunk-grids must be positive")
    if args.metallicity_solar <= 0.0:
        raise ValueError("metallicity-solar must be positive")

    manifest = json.loads(Path(args.zoom_manifest).read_text())
    final = manifest["final"]
    left = np.asarray(final["left_edge_code"], dtype=np.float64)
    shape = np.asarray(final["shape"], dtype=np.int64)
    width = float(final["width_code"])
    if left.shape != (3,) or shape.shape != (3,) or np.any(shape <= 0) or width <= 0.0:
        raise ValueError("invalid P4 final-cube geometry")
    analysis_cell_width = 2.0 ** (-args.analysis_level)
    if not np.allclose(width, analysis_cell_width * shape[0], rtol=0.0, atol=analysis_cell_width * 1.0e-8):
        raise ValueError("P4 cube does not align with the requested RAMSES analysis level")
    right = left + width
    density_sum = np.zeros(tuple(shape), dtype=np.float64)
    coverage = np.zeros(tuple(shape), dtype=np.float64)
    level_summary = []

    with h5py.File(Path(args.snapshot), "r") as handle:
        levels = sorted(
            int(name.removeprefix("level_"))
            for name in handle["amr"].keys()
            if name.startswith("level_") and f"hydro/{name}" in handle
        )
        if not levels:
            raise ValueError("HDF5 snapshot has no AMR/hydro levels")
        for level in levels:
            amr = handle[f"amr/level_{level}"]
            hydro_density = handle[f"hydro/level_{level}/uold_1"]
            leaf_width = 2.0 ** (-level)
            candidate_grids = 0
            deposited_leaves = 0
            refined_cells = 0
            for grid_indices, grid_centers in _candidate_grid_batches(
                amr,
                cell_width_code=leaf_width,
                left_edge_code=left,
                right_edge_code=right,
                chunk_grids=args.scan_chunk_grids,
            ):
                candidate_grids += len(grid_indices)
                density_blocks = _read_grid_blocks(hydro_density, grid_indices)
                son_blocks = _read_grid_blocks(amr["son_flag"], grid_indices)
                for grid_center, density_values, son_flags in zip(
                    grid_centers, density_blocks, son_blocks, strict=True
                ):
                    child_centers = grid_center + (_CHILD_BITS - 0.5) * leaf_width
                    for child_center, density_code, son_flag in zip(child_centers, density_values, son_flags, strict=True):
                        child_left = child_center - 0.5 * leaf_width
                        child_right = child_center + 0.5 * leaf_width
                        if np.any(child_right <= left) or np.any(child_left >= right):
                            continue
                        if son_flag != 0.0:
                            refined_cells += 1
                            continue
                        _deposit_leaf(
                            density_sum,
                            coverage,
                            center_code=child_center,
                            density_code=float(density_code),
                            leaf_width_code=leaf_width,
                            left_edge_code=left,
                            analysis_cell_width_code=analysis_cell_width,
                        )
                        deposited_leaves += 1
            level_summary.append(
                {
                    "level": level,
                    "candidate_grids": candidate_grids,
                    "deposited_leaf_cells": deposited_leaves,
                    "refined_cells_in_cube": refined_cells,
                }
            )
            print(
                f"P4_HDF5_SCAN level={level} candidates={candidate_grids} "
                f"leaves={deposited_leaves} refined={refined_cells}",
                flush=True,
            )

    if not np.allclose(coverage, 1.0, rtol=0.0, atol=1.0e-8):
        raise ValueError(
            "adaptive AMR leaves do not cover the P4 mesh exactly: "
            f"min={coverage.min():.16g} max={coverage.max():.16g}"
        )
    density_code = density_sum / coverage
    if not np.isfinite(density_code).all() or np.any(density_code <= 0.0):
        raise ValueError("resampled AMR gas density is non-finite or non-positive")

    info = _read_info(Path(args.info))
    atlas = read_thermal_atlas(args.thermal_atlas)
    density_cgs = density_code * info["unit_d"]
    hydrogen_number_density = density_cgs * 0.76 / 1.67262192369e-24
    temperature = atlas.equilibrium_temperature(info["aexp"], hydrogen_number_density, args.metallicity_solar)
    source_ledger = read_photon_source_ledger_csv(args.source_ledger)
    sources = source_ledger.source_catalog_in_cube(left, width, tuple(shape))
    if sources is None:
        raise ValueError("photon ledger has no source inside the P4 cube")
    snapshot = neutral_primordial_input(
        GridSpec(cell_width_cm=info["unit_l"] * analysis_cell_width, left_edge_cm=left * info["unit_l"]),
        density_cgs,
        temperature,
        sources=sources,
    )
    write_static_rt_input(args.output, snapshot)

    metadata = {
        "snapshot": str(Path(args.snapshot).resolve()),
        "info": str(Path(args.info).resolve()),
        "aexp": info["aexp"],
        "analysis_level": args.analysis_level,
        "analysis_rule": "volume-conservative resampling of all intersecting RAMSES AMR leaf cells",
        "populated_amr_levels": levels,
        "level_summary": level_summary,
        "left_edge_code": left.tolist(),
        "shape": shape.tolist(),
        "cell_width_code": analysis_cell_width,
        "cell_width_cm": info["unit_l"] * analysis_cell_width,
        "density_unit_g_cm3": info["unit_d"],
        "mean_n_h_cm3": float(hydrogen_number_density.mean()),
        "max_n_h_cm3": float(hydrogen_number_density.max()),
        "coverage_min": float(coverage.min()),
        "coverage_max": float(coverage.max()),
        "thermal_initialization": "offline Grackle thermal-atlas equilibrium fallback",
        "thermal_atlas": str(Path(args.thermal_atlas).resolve()),
        "metallicity_solar": args.metallicity_solar,
        "source_ledger": str(Path(args.source_ledger).resolve()),
        "photon_source_count": int(len(sources.cell_index)),
    }
    metadata_path = Path(args.metadata_output)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(
        "P4_HDF5_STAGE_OK "
        f"shape={tuple(shape)} sources={len(sources.cell_index)} mean_n_h={hydrogen_number_density.mean():.6g}",
        flush=True,
    )


if __name__ == "__main__":
    main()
