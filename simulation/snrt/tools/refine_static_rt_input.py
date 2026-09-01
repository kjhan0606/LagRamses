#!/usr/bin/env python3
"""Build a conservative factor-of-two static-grid refinement control input."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import h5py
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.snapshot import GridSpec, SourceCatalog, StaticRTInput, read_static_rt_input, write_static_rt_input


def _refine_field(field: np.ndarray, factor: int, axes: tuple[int, ...]) -> np.ndarray:
    refined = np.asarray(field)
    for axis in axes:
        refined = np.repeat(refined, factor, axis=axis)
    return refined


def refine_static_rt_input(snapshot: StaticRTInput, factor: int = 2) -> StaticRTInput:
    """Refine every uniform-grid field by repetition and split source luminosity.

    This is a controlled resolution experiment, not an AMR reconstruction. Gas
    fields are piecewise constant in each coarse cell. A source at a coarse
    cell centre is represented by equal luminosity in all ``factor**3`` child
    cells, preserving its centre of luminosity and total photon rate.
    """

    if factor < 2 or int(factor) != factor:
        raise ValueError("factor must be an integer greater than or equal to two")
    factor = int(factor)
    scalar_axes = (0, 1, 2)
    refined_sources = None
    if snapshot.sources is not None:
        child_offsets = np.asarray(list(np.ndindex((factor, factor, factor))), dtype=np.int64)
        coarse_index = snapshot.sources.cell_index
        refined_index = (
            factor * coarse_index[:, None, :] + child_offsets[None, :, :]
        ).reshape((-1, 3))
        refined_luminosity = np.repeat(
            snapshot.sources.photon_luminosity_s,
            factor**3,
            axis=0,
        ) / float(factor**3)
        refined_sources = SourceCatalog(refined_index, refined_luminosity)

    return StaticRTInput(
        grid=GridSpec(snapshot.grid.cell_width_cm / factor, snapshot.grid.left_edge_cm.copy()),
        hydrogen_number_density_cm3=_refine_field(snapshot.hydrogen_number_density_cm3, factor, scalar_axes),
        helium_number_density_cm3=_refine_field(snapshot.helium_number_density_cm3, factor, scalar_axes),
        temperature_k=_refine_field(snapshot.temperature_k, factor, scalar_axes),
        dust_relative_abundance=_refine_field(snapshot.dust_relative_abundance, factor, scalar_axes),
        x_hii=_refine_field(snapshot.x_hii, factor, scalar_axes),
        x_heii=_refine_field(snapshot.x_heii, factor, scalar_axes),
        x_heiii=_refine_field(snapshot.x_heiii, factor, scalar_axes),
        sources=refined_sources,
        velocity_cm_s=None
        if snapshot.velocity_cm_s is None
        else _refine_field(snapshot.velocity_cm_s, factor, (1, 2, 3)),
        metallicity_solar=None
        if snapshot.metallicity_solar is None
        else _refine_field(snapshot.metallicity_solar, factor, scalar_axes),
        dust_to_metal=None
        if snapshot.dust_to_metal is None
        else _refine_field(snapshot.dust_to_metal, factor, scalar_axes),
        x_h2=None if snapshot.x_h2 is None else _refine_field(snapshot.x_h2, factor, scalar_axes),
        cell_level=None
        if snapshot.cell_level is None
        else _refine_field(snapshot.cell_level, factor, scalar_axes),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--factor", type=int, default=2)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite existing refined input: {args.output}")
    snapshot = read_static_rt_input(args.input)
    refined = refine_static_rt_input(snapshot, args.factor)
    write_static_rt_input(args.output, refined)
    with h5py.File(args.output, "a") as handle:
        handle.attrs["synthetic_refinement_factor"] = args.factor
        handle.attrs["refinement_parent_path"] = str(args.input.resolve())
        handle.attrs["refinement_gas_interpolation"] = "piecewise_constant"
        handle.attrs["refinement_source_mapping"] = "equal_luminosity_to_all_child_cells"
    coarse_total = (
        np.zeros(snapshot.sources.photon_luminosity_s.shape[1], dtype=np.float64)
        if snapshot.sources is None
        else snapshot.sources.photon_luminosity_s.sum(axis=0, dtype=np.float64)
    )
    refined_total = (
        np.zeros_like(coarse_total)
        if refined.sources is None
        else refined.sources.photon_luminosity_s.sum(axis=0, dtype=np.float64)
    )
    if not np.allclose(refined_total, coarse_total, rtol=1.0e-13, atol=0.0):
        raise RuntimeError("source luminosity changed during refinement")
    print(
        "REFINE_STATIC_RT_INPUT_OK "
        f"coarse_shape={snapshot.shape} refined_shape={refined.shape} factor={args.factor} "
        f"coarse_sources={0 if snapshot.sources is None else len(snapshot.sources.cell_index)} "
        f"refined_sources={0 if refined.sources is None else len(refined.sources.cell_index)} "
        f"groups={len(coarse_total)} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
