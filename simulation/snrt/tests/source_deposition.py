#!/usr/bin/env python3
"""Validate source-deposition modes and group-luminosity conservation."""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np

from snrt_core.sources import PointSources, deposit_point_sources


jax.config.update("jax_enable_x64", True)


def main() -> int:
    shape = (5, 4, 3)
    cell_width = (2.0, 3.0, 5.0)
    cell_volume = float(np.prod(cell_width))
    sources = PointSources(
        jnp.asarray(((0, 0, 0), (2, 2, 1), (2, 2, 1)), dtype=jnp.int32),
        jnp.asarray(((3.0e49, 7.0e48), (2.0e49, 5.0e48), (1.0e49, 2.0e48)), dtype=jnp.float64),
    )
    expected = np.asarray(sources.luminosity, dtype=np.float64).sum(axis=0)

    point = np.asarray(deposit_point_sources(shape, cell_width, sources, dtype=jnp.float64))
    compact = np.asarray(
        deposit_point_sources(shape, cell_width, sources, dtype=jnp.float64, deposition_mode="compact3")
    )
    for mode, field in (("point", point), ("compact3", compact)):
        recovered = field.sum(axis=(1, 2, 3), dtype=np.float64) * cell_volume
        assert np.allclose(recovered, expected, rtol=1.0e-13, atol=1.0e34), (mode, recovered, expected)
        assert np.isfinite(field).all()
        assert np.min(field) >= 0.0

    # The boundary source must not wrap around the opposite edge. The compact
    # kernel has support only in the first two cells along each boundary axis.
    boundary_field = np.asarray(
        deposit_point_sources(
            shape,
            cell_width,
            PointSources(sources.cell_index[:1], sources.luminosity[:1]),
            dtype=jnp.float64,
            deposition_mode="compact3",
        )
    )
    boundary_only = boundary_field[:, :2, :2, :2]
    assert np.allclose(boundary_only.sum(axis=(1, 2, 3)) * cell_volume, sources.luminosity[0], rtol=1.0e-13, atol=1.0e34)
    assert np.all(boundary_field[:, -1, :, :] == 0.0)
    assert np.all(boundary_field[:, :, -1, :] == 0.0)
    assert np.all(boundary_field[:, :, :, -1] == 0.0)
    assert np.count_nonzero(compact) > np.count_nonzero(point)

    try:
        deposit_point_sources(shape, cell_width, sources, deposition_mode="unknown")
    except ValueError:
        pass
    else:
        raise AssertionError("unsupported deposition mode was accepted")
    print("SOURCE_DEPOSITION_TEST_OK modes=point,compact3 luminosity_conserved=true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
