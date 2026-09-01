"""Deterministic problem constructors for the first S_N RT benchmarks."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .coupling import build_hydrogen_radiation_step
from .quadrature import s4_quadrature
from .sources import PointSources, deposit_point_sources
from .transport import TransportConfig, initial_intensity


PARSEC_CM = 3.0856775814913673e18
EV_ERG = 1.602176634e-12


class StromgrenProblem(NamedTuple):
    """All static inputs for the B01 isothermal Stromgren-sphere benchmark."""

    config: TransportConfig
    directions: jnp.ndarray
    weights: jnp.ndarray
    intensity: jnp.ndarray
    emissivity: jnp.ndarray
    n_h: jnp.ndarray
    x_hii: jnp.ndarray
    temperature: jnp.ndarray
    group_cross_section: jnp.ndarray
    group_excess_energy: jnp.ndarray
    analytic_radius_cm: float


class StromgrenState(NamedTuple):
    """Dynamical state advanced by the B01 runner."""

    intensity: jnp.ndarray
    x_hii: jnp.ndarray


def make_stromgren_problem(
    shape: tuple[int, int, int] = (64, 64, 64),
    source_photon_rate: float = 1.0e49,
    hydrogen_number_density: float = 1.0,
    reduced_light_speed: float = 2.99792458e10,
    courant: float = 0.2,
    dtype=jnp.float32,
) -> StromgrenProblem:
    """Construct B01 with a single HI-ionizing photon group.

    The gas is held isothermal at 1e4 K so the case-B coefficient in the
    chemistry core matches the analytic Stromgren radius. The source is placed
    in the central cell and all outer boundaries are vacuum boundaries.
    """

    alpha_b = 2.59e-13
    analytic_radius_cm = (3.0 * source_photon_rate / (4.0 * jnp.pi * alpha_b * hydrogen_number_density**2)) ** (1.0 / 3.0)
    domain_length = 6.0 * analytic_radius_cm
    cell_width = tuple(float(domain_length / cells) for cells in shape)
    directions, weights = s4_quadrature(dtype)
    inverse_width = jnp.asarray([1.0 / value for value in cell_width], dtype=dtype)
    maximum_directional_sum = float(jnp.max(jnp.sum(jnp.abs(directions) * inverse_width[None, :], axis=1)))
    dt = courant / (reduced_light_speed * maximum_directional_sum)
    config = TransportConfig(cell_width, dt, reduced_light_speed)

    central_cell = jnp.asarray([[shape[0] // 2, shape[1] // 2, shape[2] // 2]], dtype=jnp.int32)
    sources = PointSources(central_cell, jnp.asarray([[source_photon_rate]], dtype=dtype))
    emissivity = deposit_point_sources(shape, cell_width, sources, dtype)
    intensity = initial_intensity(1, directions.shape[0], shape, dtype)
    n_h = jnp.full(shape, hydrogen_number_density, dtype=dtype)
    x_hii = jnp.zeros(shape, dtype=dtype)
    temperature = jnp.full(shape, 1.0e4, dtype=dtype)
    group_cross_section = jnp.asarray([6.30e-18], dtype=dtype)
    group_excess_energy = jnp.asarray([0.0 * EV_ERG], dtype=dtype)
    return StromgrenProblem(
        config,
        directions,
        weights,
        intensity,
        emissivity,
        n_h,
        x_hii,
        temperature,
        group_cross_section,
        group_excess_energy,
        float(analytic_radius_cm),
    )


def build_stromgren_runner(problem: StromgrenProblem, number_of_steps: int):
    """Build a fixed-step JIT B01 evolution function without executing it."""

    step = build_hydrogen_radiation_step(
        problem.config,
        problem.directions,
        problem.weights,
        problem.group_cross_section,
        problem.group_excess_energy,
    )

    @jax.jit
    def run(state: StromgrenState) -> StromgrenState:
        def body(_: int, current: StromgrenState) -> StromgrenState:
            next_intensity, chemistry = step(
                current.intensity,
                problem.emissivity,
                problem.n_h,
                current.x_hii,
                problem.temperature,
            )
            return StromgrenState(next_intensity, chemistry.x_hii)

        return jax.lax.fori_loop(0, number_of_steps, body, state)

    return run
