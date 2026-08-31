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
) -> jnp.ndarray:
    """Deposit source luminosity as a group emissivity field.

    Luminosity and cell widths must use a mutually consistent unit system.
    The returned array has shape ``[group, nx, ny, nz]``.
    """

    number_of_groups = sources.luminosity.shape[1]
    emissivity = jnp.zeros((number_of_groups, *shape), dtype=dtype)
    # This setup operation retains stellar/AGN Q in host float64 until after
    # division by a parsec-scale volume. The resulting source field is finite
    # in float32 and is the only quantity required by the JAX transport kernel.
    source_luminosity = np.asarray(sources.luminosity, dtype=np.float64)
    cell_volume = np.prod(np.asarray(cell_width, dtype=np.float64))
    source_emissivity = jnp.asarray(source_luminosity.T / cell_volume, dtype=dtype)
    ix, iy, iz = (sources.cell_index[:, 0], sources.cell_index[:, 1], sources.cell_index[:, 2])
    return emissivity.at[:, ix, iy, iz].add(source_emissivity)
