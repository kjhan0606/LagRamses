"""Static-opacity shadow benchmark for angular-transport verification."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .quadrature import level_symmetric_quadrature
from .sources import PointSources, deposit_point_sources
from .transport import TransportConfig, advance_explicit


class ShadowProblem(NamedTuple):
    """One-group point-source shadow through a fixed, opaque spherical clump."""

    config: TransportConfig
    directions: jnp.ndarray
    weights: jnp.ndarray
    emissivity: jnp.ndarray
    absorption: jnp.ndarray
    shadow_mask: jnp.ndarray
    control_mask: jnp.ndarray


def make_opaque_clump_problem(
    shape: tuple[int, int, int] = (64, 64, 64),
    order: int = 4,
    courant: float = 0.2,
    clump_radius_cells: float = 8.0,
    clump_absorption: float = 64.0,
    source_luminosity: float = 1.0,
    dtype: jnp.dtype = jnp.float32,
) -> ShadowProblem:
    """Build B03: a point source and an optically thick, non-scattering clump.

    Code units set the cube side and reduced light speed to one. The explicit
    time step meets both the directional transport CFL condition and the clump
    absorption condition for the supplied default parameters.
    """
    if min(shape) < 32:
        raise ValueError("The shadow benchmark requires at least 32 cells per dimension.")

    directions, weights = level_symmetric_quadrature(order, dtype)
    cell_width = jnp.asarray(tuple(1.0 / size for size in shape), dtype=dtype)
    directional_rate = jnp.max(jnp.sum(jnp.abs(directions) / cell_width[None, :], axis=1))
    dt = courant / directional_rate
    config = TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=1.0)

    source_cell = jnp.asarray([[max(2, shape[0] // 8), shape[1] // 2, shape[2] // 2]], dtype=jnp.int32)
    luminosity = jnp.asarray([[source_luminosity]], dtype=dtype)
    emissivity = deposit_point_sources(shape, cell_width, PointSources(source_cell, luminosity))

    x, y, z = jnp.meshgrid(
        jnp.arange(shape[0], dtype=dtype) + 0.5,
        jnp.arange(shape[1], dtype=dtype) + 0.5,
        jnp.arange(shape[2], dtype=dtype) + 0.5,
        indexing="ij",
    )
    center_x = 0.50 * shape[0]
    center_y = 0.50 * shape[1]
    center_z = 0.50 * shape[2]
    transverse_radius = jnp.sqrt((y - center_y) ** 2 + (z - center_z) ** 2)
    clump = (x - center_x) ** 2 + (y - center_y) ** 2 + (z - center_z) ** 2 <= clump_radius_cells**2
    absorption = jnp.where(clump, clump_absorption, 0.0)[None, ...].astype(dtype)

    downstream = x > center_x + clump_radius_cells
    shadow_mask = downstream & (transverse_radius < 0.75 * clump_radius_cells)
    control_mask = downstream & (transverse_radius > 1.5 * clump_radius_cells) & (transverse_radius < 2.25 * clump_radius_cells)
    return ShadowProblem(config, directions, weights, emissivity, absorption, shadow_mask, control_mask)


def build_shadow_runner(problem: ShadowProblem, n_steps: int):
    """Create a compiled explicit B03 evolution with the clump held fixed."""
    if n_steps < 1:
        raise ValueError("n_steps must be positive.")

    initial_intensity = jnp.zeros(
        (problem.emissivity.shape[0], problem.directions.shape[0], *problem.emissivity.shape[1:]),
        dtype=problem.emissivity.dtype,
    )

    def run() -> jnp.ndarray:
        def step(_: int, intensity: jnp.ndarray) -> jnp.ndarray:
            return advance_explicit(
                problem.config,
                problem.directions,
                intensity,
                problem.emissivity,
                problem.absorption,
            )

        return jax.lax.fori_loop(0, n_steps, step, initial_intensity)

    return jax.jit(run)


def shadow_contrast(number_density: jnp.ndarray, problem: ShadowProblem) -> jnp.ndarray:
    """Return mean shadow intensity divided by a nearby unblocked control region."""
    field = number_density[0]
    return jnp.mean(field[problem.shadow_mask]) / jnp.mean(field[problem.control_mask])
