"""Source deposition for static-grid S_N transport."""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp
import numpy as np


class PointSources(NamedTuple):
    """Cell-centered point sources with luminosity ordered as [source, group]."""

    cell_index: jnp.ndarray
    luminosity: jnp.ndarray


def deposit_point_sources(
    shape: tuple[int, int, int],
    cell_width: tuple[float, float, float],
    sources: PointSources,
    dtype=jnp.float32,
    deposition_mode: str = "point",
) -> jnp.ndarray:
    """Deposit source luminosity as a group emissivity field.

    Luminosity and cell widths must use a mutually consistent unit system.
    The returned array has shape ``[group, nx, ny, nz]``.

    ``point`` preserves the original cell-centred deposition.  ``compact3``
    is an opt-in numerical control that distributes each source over a
    3-by-3-by-3 compact kernel with one-dimensional weights ``[1/4, 1/2,
    1/4]``.  The kernel is renormalized at domain boundaries, so each source
    conserves its group luminosity exactly in host-side double precision.
    ``compact3`` is a spatial-regularization control, not a physical source
    size model.
    """

    number_of_groups = sources.luminosity.shape[1]
    source_luminosity = np.asarray(sources.luminosity, dtype=np.float64)
    cell_volume = np.prod(np.asarray(cell_width, dtype=np.float64))
    if deposition_mode == "point":
        emissivity = jnp.zeros((number_of_groups, *shape), dtype=dtype)
        # This setup operation retains stellar/AGN Q in host float64 until
        # after division by a parsec-scale volume. The resulting source field
        # is finite in float32 and is the only quantity required by the JAX
        # transport kernel.
        source_emissivity = jnp.asarray(source_luminosity.T / cell_volume, dtype=dtype)
        ix, iy, iz = (sources.cell_index[:, 0], sources.cell_index[:, 1], sources.cell_index[:, 2])
        return emissivity.at[:, ix, iy, iz].add(source_emissivity)
    if deposition_mode != "compact3":
        raise ValueError(f"unsupported source deposition mode: {deposition_mode!r}")

    cell_index = np.asarray(sources.cell_index, dtype=np.int64)
    if cell_index.ndim != 2 or cell_index.shape[1] != 3:
        raise ValueError("source cell_index must have shape [source, 3]")
    if source_luminosity.ndim != 2 or source_luminosity.shape[0] != cell_index.shape[0]:
        raise ValueError("source luminosity must have shape [source, group]")
    if np.any(cell_index < 0) or np.any(cell_index >= np.asarray(shape, dtype=np.int64)[None, :]):
        raise ValueError("source cell_index contains an out-of-bounds cell")

    # Construct the compact kernel on the host so that duplicate sources and
    # boundary renormalization are handled deterministically before conversion
    # to the JAX transport dtype.
    compact_weights = np.asarray((0.25, 0.5, 0.25), dtype=np.float64)
    normalizer = np.ones(cell_index.shape[0], dtype=np.float64)
    for axis, extent in enumerate(shape):
        axis_normalizer = np.zeros(cell_index.shape[0], dtype=np.float64)
        for offset, weight in zip((-1, 0, 1), compact_weights, strict=True):
            target = cell_index[:, axis] + offset
            valid = (target >= 0) & (target < extent)
            axis_normalizer += np.where(valid, weight, 0.0)
        normalizer *= axis_normalizer

    deposited_luminosity = np.zeros((number_of_groups, *shape), dtype=np.float64)
    for offsets in np.ndindex((3, 3, 3)):
        delta = tuple(index - 1 for index in offsets)
        target = cell_index + np.asarray(delta, dtype=np.int64)[None, :]
        valid = np.all((target >= 0) & (target < np.asarray(shape, dtype=np.int64)[None, :]), axis=1)
        if not np.any(valid):
            continue
        weight = float(np.prod(compact_weights[list(offsets)], dtype=np.float64))
        factor = weight / normalizer[valid]
        for group in range(number_of_groups):
            np.add.at(
                deposited_luminosity[group],
                (target[valid, 0], target[valid, 1], target[valid, 2]),
                source_luminosity[valid, group] * factor,
            )
    return jnp.asarray(deposited_luminosity / cell_volume, dtype=dtype)
