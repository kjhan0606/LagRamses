"""Physical-unit primordial Stromgren benchmark for the coupled S_N solver."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpy as np

from .photon_coupling import build_primordial_radiation_step
from .primordial import PrimordialState, hui_gnedin_case_b_hydrogen, primordial_cross_sections
from .quadrature import level_symmetric_quadrature
from .sources import PointSources, deposit_point_sources
from .transport import TransportConfig


PARSEC_CM = 3.0856775814913673e18
LIGHT_SPEED_CM_S = 2.99792458e10


class PrimordialStromgrenProblem(NamedTuple):
    config: TransportConfig
    directions: jnp.ndarray
    weights: jnp.ndarray
    group_energy_ev: jnp.ndarray
    emissivity: jnp.ndarray
    initial_state: PrimordialState
    temperature_k: jnp.ndarray
    stromgren_radius_cm: float


class PrimordialStromgrenState(NamedTuple):
    intensity: jnp.ndarray
    chemistry: PrimordialState


def make_primordial_stromgren_problem(
    shape: tuple[int, int, int] = (64, 64, 64),
    cell_size_parsec: float = 4.0,
    ionizing_photon_rate: float = 1.0e49,
    hydrogen_density_cm3: float = 1.0,
    helium_number_ratio: float = 0.079,
    temperature_k: float = 1.0e4,
    reduced_light_speed_fraction: float = 1.0e-2,
    courant: float = 0.2,
    order: int = 4,
    group_energy_ev: float = 20.0,
    dtype: jnp.dtype = jnp.float32,
) -> PrimordialStromgrenProblem:
    """Build B01 in cgs units with one hydrogen-ionizing photon group.

    The reduced speed of light is a numerical parameter whose validity must be
    established against a fixed-duration analytic or higher-c reference; no
    universal front-speed margin is implied by the default. The analytic
    reference is hydrogen-only. At the default 20 eV, the He I cross section is
    zero, so the structurally present helium remains neutral and does not
    contribute to the measured radius discrepancy.
    """
    if not 0.0 < reduced_light_speed_fraction <= 1.0:
        raise ValueError("reduced_light_speed_fraction must lie in (0, 1].")

    directions, weights = level_symmetric_quadrature(order, dtype)
    cell_width = jnp.full((3,), cell_size_parsec * PARSEC_CM, dtype=dtype)
    reduced_light_speed = reduced_light_speed_fraction * LIGHT_SPEED_CM_S
    group_energies = jnp.asarray([group_energy_ev], dtype=dtype)
    directional_rate = jnp.max(jnp.sum(jnp.abs(directions) / cell_width[None, :], axis=1))
    dt = float(courant / (reduced_light_speed * directional_rate))
    config = TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_speed)

    source_cell = jnp.asarray([[shape[0] // 2, shape[1] // 2, shape[2] // 2]], dtype=jnp.int32)
    luminosity = np.asarray([[ionizing_photon_rate]], dtype=np.float64)
    emissivity = deposit_point_sources(shape, cell_width, PointSources(source_cell, luminosity))
    state = PrimordialState(
        n_hydrogen=jnp.full(shape, hydrogen_density_cm3, dtype=dtype),
        n_helium=jnp.full(shape, helium_number_ratio * hydrogen_density_cm3, dtype=dtype),
        x_hydrogen_ii=jnp.zeros(shape, dtype=dtype),
        x_helium_ii=jnp.zeros(shape, dtype=dtype),
        x_helium_iii=jnp.zeros(shape, dtype=dtype),
    )
    temperature = jnp.full(shape, temperature_k, dtype=dtype)
    alpha_b = float(hui_gnedin_case_b_hydrogen(jnp.asarray(temperature_k)))
    stromgren_radius = (3.0 * ionizing_photon_rate / (4.0 * jnp.pi * alpha_b * hydrogen_density_cm3**2)) ** (1.0 / 3.0)
    return PrimordialStromgrenProblem(
        config=config,
        directions=directions,
        weights=weights,
        group_energy_ev=group_energies,
        emissivity=emissivity,
        initial_state=state,
        temperature_k=temperature,
        stromgren_radius_cm=float(stromgren_radius),
    )


def build_primordial_stromgren_runner(problem: PrimordialStromgrenProblem, n_steps: int):
    """Return a compiled B01 run; execution remains caller-controlled."""
    if n_steps < 1:
        raise ValueError("n_steps must be positive.")

    cross_sections = primordial_cross_sections(problem.group_energy_ev)
    step = build_primordial_radiation_step(
        problem.directions,
        problem.weights,
        problem.config,
        cross_sections,
        problem.group_energy_ev,
    )
    initial_intensity = jnp.zeros(
        (problem.emissivity.shape[0], problem.directions.shape[0], *problem.emissivity.shape[1:]),
        dtype=problem.emissivity.dtype,
    )

    def run() -> PrimordialStromgrenState:
        def advance(_: int, carry: PrimordialStromgrenState) -> PrimordialStromgrenState:
            result = step(carry.intensity, problem.emissivity, carry.chemistry, problem.temperature_k)
            return PrimordialStromgrenState(result.intensity, result.state)

        initial = PrimordialStromgrenState(initial_intensity, problem.initial_state)
        return jax.lax.fori_loop(0, n_steps, advance, initial)

    return jax.jit(run)


def ionized_volume_radius(state: PrimordialState, cell_width_cm: jnp.ndarray) -> jnp.ndarray:
    """Infer an equivalent H II radius from cells with x_HII greater than 0.5."""
    cell_volume_root = jnp.prod(jnp.asarray(cell_width_cm) ** (1.0 / 3.0))
    occupied_cells = jnp.sum(state.x_hydrogen_ii > 0.5)
    return (3.0 * occupied_cells / (4.0 * jnp.pi)) ** (1.0 / 3.0) * cell_volume_root
