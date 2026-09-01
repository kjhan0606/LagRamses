"""Fixed-iteration local implicit closure for primordial recombination."""

from __future__ import annotations

import jax
import jax.numpy as jnp

from .primordial import PrimordialState, case_b_helium_recombination, hui_gnedin_case_b_hydrogen
from .primordial_cooling import collisional_ionization_coefficients


def hydrogen_photoionization_relaxation(
    x_hydrogen_ii: jnp.ndarray,
    photoionization_rate: jnp.ndarray,
    electron_density: jnp.ndarray,
    temperature_k: jnp.ndarray,
    dt: float,
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return exact H I/H II relaxation and its time-averaged H II fraction.

    The photoionization and electron densities are held fixed over ``dt``. This
    is the local analytic closure used by time-averaged-opacity iterations:
    the mean neutral fraction used by transport is consistent with the H rate
    solution rather than with a capped one-shot photon transfer.
    """

    next_fraction, mean_fraction, _, _ = hydrogen_neutral_relaxation(
        1.0 - x_hydrogen_ii,
        photoionization_rate,
        electron_density,
        temperature_k,
        dt,
    )
    return next_fraction, mean_fraction


def hydrogen_neutral_relaxation(
    x_hydrogen_i: jnp.ndarray,
    photoionization_rate: jnp.ndarray,
    electron_density: jnp.ndarray,
    temperature_k: jnp.ndarray,
    dt: float,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Solve H I/H II relaxation while retaining the small neutral fraction.

    The neutral fraction is the opacity carrier in an ionized source cell.
    Evolving it directly avoids subtracting two nearly equal float32 H II
    fractions before each transport/chemistry fixed-point iteration.
    """

    alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
    recombination_rate = alpha_hii * electron_density
    total_rate = photoionization_rate + recombination_rate
    safe_rate = jnp.maximum(total_rate, jnp.finfo(total_rate.dtype).tiny)
    equilibrium_neutral = recombination_rate / safe_rate
    rate_dt = total_rate * dt
    decay = jnp.exp(-rate_dt)
    initial_neutral = jnp.clip(x_hydrogen_i, 0.0, 1.0)
    next_neutral = equilibrium_neutral + (initial_neutral - equilibrium_neutral) * decay
    mean_factor = jnp.where(
        rate_dt > 1.0e-4,
        -jnp.expm1(-rate_dt) / jnp.maximum(rate_dt, jnp.finfo(rate_dt.dtype).tiny),
        1.0 - 0.5 * rate_dt + rate_dt**2 / 6.0,
    )
    mean_neutral = equilibrium_neutral + (initial_neutral - equilibrium_neutral) * mean_factor
    next_neutral = jnp.clip(next_neutral, 0.0, 1.0)
    mean_neutral = jnp.clip(mean_neutral, 0.0, 1.0)
    return 1.0 - next_neutral, 1.0 - mean_neutral, next_neutral, mean_neutral


def helium_photoionization_backward_euler(
    x_helium_ii: jnp.ndarray,
    x_helium_iii: jnp.ndarray,
    photoionization_hei_rate: jnp.ndarray,
    photoionization_heii_rate: jnp.ndarray,
    electron_density: jnp.ndarray,
    temperature_k: jnp.ndarray,
    dt: float,
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Solve the three-state He network with a local backward-Euler update.

    The H I, He I, and He II photoionization rates and the electron density are
    fixed by one time-averaged-opacity iteration. The resulting 3-by-3 linear
    system has an analytic scalar elimination, avoiding a per-cell matrix solve
    while preserving positivity and He-number conservation.
    """

    alpha_heii, alpha_heiii = case_b_helium_recombination(temperature_k)
    recombination_heii = alpha_heii * electron_density
    recombination_heiii = alpha_heiii * electron_density
    heiii_factor = 1.0 + dt * recombination_heiii
    heiii_constant = x_helium_iii / heiii_factor
    heiii_from_heii = dt * photoionization_heii_rate / heiii_factor
    numerator = x_helium_ii + dt * (
        photoionization_hei_rate * (1.0 - heiii_constant) + recombination_heiii * heiii_constant
    )
    denominator = 1.0 + dt * (
        photoionization_hei_rate * (1.0 + heiii_from_heii)
        + photoionization_heii_rate
        + recombination_heii
        - recombination_heiii * heiii_from_heii
    )
    next_heii = numerator / jnp.maximum(denominator, jnp.finfo(denominator.dtype).tiny)
    next_heiii = heiii_constant + heiii_from_heii * next_heii
    next_heii = jnp.clip(next_heii, 0.0, 1.0)
    next_heiii = jnp.clip(next_heiii, 0.0, 1.0 - next_heii)
    return next_heii, next_heiii


def implicit_case_b_recombination(
    state: PrimordialState,
    temperature_k: jnp.ndarray,
    dt: float,
    iterations: int = 24,
) -> PrimordialState:
    """Apply coupled backward-Euler H/He recombination by electron-density bisection."""

    return implicit_case_b_recombination_with_recombinations(
        state,
        temperature_k,
        dt,
        iterations,
    )[0]


def implicit_case_b_recombination_with_recombinations(
    state: PrimordialState,
    temperature_k: jnp.ndarray,
    dt: float,
    iterations: int = 24,
) -> tuple[PrimordialState, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Apply coupled backward-Euler H/He recombination by electron-density bisection.

    The input is the state after photon-driven ionizations. At fixed electron
    density the backward-Euler fractions are analytic, leaving one monotonic
    scalar equation per cell. Fixed-count bisection keeps the kernel static,
    positive, and robust at large alpha*n_e*dt. The returned recombination
    densities are ordered H II -> H I, He II -> He I, and He III -> He II.
    """
    if iterations < 1:
        raise ValueError("iterations must be positive.")

    alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
    alpha_heii, alpha_heiii = case_b_helium_recombination(temperature_k)
    x_hii_initial = state.x_hydrogen_ii
    x_heii_initial = state.x_helium_ii
    x_heiii_initial = state.x_helium_iii

    def fractions_from_electron_density(electron_density: jnp.ndarray) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
        h_factor = dt * alpha_hii * electron_density
        heii_factor = dt * alpha_heii * electron_density
        heiii_factor = dt * alpha_heiii * electron_density
        next_hii = x_hii_initial / (1.0 + h_factor)
        next_heiii = x_heiii_initial / (1.0 + heiii_factor)
        next_heii = (x_heii_initial + heiii_factor * next_heiii) / (1.0 + heii_factor)
        return next_hii, next_heii, next_heiii

    initial_electron_density = state.n_hydrogen * x_hii_initial + state.n_helium * (x_heii_initial + 2.0 * x_heiii_initial)

    def iterate(_: int, bounds: tuple[jnp.ndarray, jnp.ndarray]) -> tuple[jnp.ndarray, jnp.ndarray]:
        lower, upper = bounds
        midpoint = 0.5 * (lower + upper)
        x_hii, x_heii, x_heiii = fractions_from_electron_density(midpoint)
        residual = midpoint - (state.n_hydrogen * x_hii + state.n_helium * (x_heii + 2.0 * x_heiii))
        return jnp.where(residual > 0.0, lower, midpoint), jnp.where(residual > 0.0, midpoint, upper)

    lower, upper = jax.lax.fori_loop(
        0,
        iterations,
        iterate,
        (jnp.zeros_like(initial_electron_density), initial_electron_density),
    )
    solved_electron_density = 0.5 * (lower + upper)
    x_hii, x_heii, x_heiii = fractions_from_electron_density(solved_electron_density)
    next_state = PrimordialState(
        n_hydrogen=state.n_hydrogen,
        n_helium=state.n_helium,
        x_hydrogen_ii=jnp.clip(x_hii, 0.0, 1.0),
        x_helium_ii=jnp.clip(x_heii, 0.0, 1.0),
        x_helium_iii=jnp.clip(x_heiii, 0.0, 1.0),
    )
    recombination_hydrogen = (
        state.n_hydrogen * dt * alpha_hii * solved_electron_density * next_state.x_hydrogen_ii
    )
    recombination_helium_ii = (
        state.n_helium * dt * alpha_heii * solved_electron_density * next_state.x_helium_ii
    )
    recombination_helium_iii = (
        state.n_helium * dt * alpha_heiii * solved_electron_density * next_state.x_helium_iii
    )
    return (
        next_state,
        recombination_hydrogen,
        recombination_helium_ii,
        recombination_helium_iii,
    )


def implicit_atomic_chemistry_with_transitions(
    state: PrimordialState,
    temperature_k: jnp.ndarray,
    dt: float,
    iterations: int = 24,
) -> tuple[
    PrimordialState,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
]:
    """Apply coupled collisional ionization and case-B recombination.

    At fixed electron density the backward-Euler H and three-state He systems
    are analytic. A scalar bisection closes the electron density against the
    resulting ion fractions. Returned transition densities are recombinations
    H II, He II, He III followed by collisional ionizations H I, He I, He II.
    """

    if iterations < 1:
        raise ValueError("iterations must be positive")
    if dt < 0.0:
        raise ValueError("dt must be non-negative")

    alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
    alpha_heii, alpha_heiii = case_b_helium_recombination(temperature_k)
    collisional = collisional_ionization_coefficients(temperature_k)

    def fractions_from_electron_density(
        electron_density: jnp.ndarray,
    ) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
        hydrogen_ionization = collisional.hydrogen_i * electron_density
        hydrogen_recombination = alpha_hii * electron_density
        next_hii = (
            state.x_hydrogen_ii + dt * hydrogen_ionization
        ) / (1.0 + dt * (hydrogen_ionization + hydrogen_recombination))

        helium_i_ionization = collisional.helium_i * electron_density
        helium_ii_ionization = collisional.helium_ii * electron_density
        helium_ii_recombination = alpha_heii * electron_density
        helium_iii_recombination = alpha_heiii * electron_density
        heiii_factor = 1.0 + dt * helium_iii_recombination
        heiii_constant = state.x_helium_iii / heiii_factor
        heiii_from_heii = dt * helium_ii_ionization / heiii_factor
        numerator = state.x_helium_ii + dt * (
            helium_i_ionization * (1.0 - heiii_constant)
            + helium_iii_recombination * heiii_constant
        )
        denominator = 1.0 + dt * (
            helium_i_ionization * (1.0 + heiii_from_heii)
            + helium_ii_ionization
            + helium_ii_recombination
            - helium_iii_recombination * heiii_from_heii
        )
        next_heii = numerator / jnp.maximum(denominator, jnp.finfo(denominator.dtype).tiny)
        next_heiii = heiii_constant + heiii_from_heii * next_heii
        return next_hii, next_heii, next_heiii

    maximum_electron_density = state.n_hydrogen + 2.0 * state.n_helium

    def iterate(_: int, bounds: tuple[jnp.ndarray, jnp.ndarray]) -> tuple[jnp.ndarray, jnp.ndarray]:
        lower, upper = bounds
        midpoint = 0.5 * (lower + upper)
        x_hii, x_heii, x_heiii = fractions_from_electron_density(midpoint)
        implied = state.n_hydrogen * x_hii + state.n_helium * (x_heii + 2.0 * x_heiii)
        residual = midpoint - implied
        return jnp.where(residual > 0.0, lower, midpoint), jnp.where(residual > 0.0, midpoint, upper)

    lower, upper = jax.lax.fori_loop(
        0,
        iterations,
        iterate,
        (jnp.zeros_like(maximum_electron_density), maximum_electron_density),
    )
    electron_density = 0.5 * (lower + upper)
    x_hii, x_heii, x_heiii = fractions_from_electron_density(electron_density)
    next_state = PrimordialState(
        state.n_hydrogen,
        state.n_helium,
        jnp.clip(x_hii, 0.0, 1.0),
        jnp.clip(x_heii, 0.0, 1.0),
        jnp.clip(x_heiii, 0.0, 1.0 - jnp.clip(x_heii, 0.0, 1.0)),
    )
    next_hei = 1.0 - next_state.x_helium_ii - next_state.x_helium_iii
    recombination_hydrogen = state.n_hydrogen * dt * alpha_hii * electron_density * next_state.x_hydrogen_ii
    recombination_helium_ii = state.n_helium * dt * alpha_heii * electron_density * next_state.x_helium_ii
    recombination_helium_iii = state.n_helium * dt * alpha_heiii * electron_density * next_state.x_helium_iii
    ionization_hydrogen_i = state.n_hydrogen * dt * collisional.hydrogen_i * electron_density * (1.0 - next_state.x_hydrogen_ii)
    ionization_helium_i = state.n_helium * dt * collisional.helium_i * electron_density * next_hei
    ionization_helium_ii = state.n_helium * dt * collisional.helium_ii * electron_density * next_state.x_helium_ii
    return (
        next_state,
        recombination_hydrogen,
        recombination_helium_ii,
        recombination_helium_iii,
        ionization_hydrogen_i,
        ionization_helium_i,
        ionization_helium_ii,
    )
