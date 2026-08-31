"""Two-source crossing-beam benchmark independent of any AGN application."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .quadrature import level_symmetric_quadrature
from .sources import PointSources, deposit_point_sources
from .transport import TransportConfig, advance_explicit


class CrossingBeamProblem(NamedTuple):
    config: TransportConfig
    directions: jnp.ndarray
    weights: jnp.ndarray
    emissivity: jnp.ndarray
    absorption: jnp.ndarray
    midpoint: tuple[int, int, int]


def make_crossing_beam_problem(
    shape: tuple[int, int, int] = (64, 64, 64),
    order: int = 4,
    courant: float = 0.2,
    source_luminosity: float = 1.0,
    dtype: jnp.dtype = jnp.float32,
) -> CrossingBeamProblem:
    """Build B04: two equal, transparent-medium sources whose fields overlap."""
    directions, weights = level_symmetric_quadrature(order, dtype)
    cell_width = jnp.asarray(tuple(1.0 / size for size in shape), dtype=dtype)
    directional_rate = jnp.max(jnp.sum(jnp.abs(directions) / cell_width[None, :], axis=1))
    config = TransportConfig(cell_width=cell_width, dt=courant / directional_rate, reduced_light_speed=1.0)

    midpoint = (shape[0] // 2, shape[1] // 2, shape[2] // 2)
    source_cells = jnp.asarray(
        [
            [shape[0] // 4, midpoint[1], midpoint[2]],
            [3 * shape[0] // 4, midpoint[1], midpoint[2]],
        ],
        dtype=jnp.int32,
    )
    luminosity = jnp.full((2, 1), source_luminosity, dtype=dtype)
    emissivity = deposit_point_sources(shape, cell_width, PointSources(source_cells, luminosity))
    absorption = jnp.zeros_like(emissivity)
    return CrossingBeamProblem(config, directions, weights, emissivity, absorption, midpoint)


def build_crossing_beam_runner(problem: CrossingBeamProblem, n_steps: int):
    """Return the compiled transparent-medium B04 evolution."""
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


def midpoint_flux_factor(intensity: jnp.ndarray, problem: CrossingBeamProblem) -> jnp.ndarray:
    """Return |F|/(cN) at the symmetry point; equal beams should approach zero."""
    number_density = jnp.einsum("d,gdxyz->gxyz", problem.weights, intensity)
    flux_factor = jnp.einsum("d,di,gdxyz->gixyz", problem.weights, problem.directions, intensity)
    x, y, z = problem.midpoint
    return jnp.linalg.norm(flux_factor[0, :, x, y, z]) / number_density[0, x, y, z]
