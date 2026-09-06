"""Static excess-IR energy transport with local equilibrium reprocessing.

Energy fields are erg/cm^3 per normalized angular direction. Primary dust
heating is erg/cm^3/s. Only configured IR bands are transported; the spectral
complement is an explicit escape sink. No grain heat capacity or gas coupling.
"""

from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpy as np

from .dust import DustThermalClosure, EV_ERG
from .transport import TransportConfig, advance_with_absorption, angular_integral


class ExcessTable(NamedTuple):
    power: jnp.ndarray
    log_temperature: jnp.ndarray
    band_power: jnp.ndarray
    band_photons: jnp.ndarray
    outside_power: jnp.ndarray
    background_power: jnp.ndarray


class ExcessEmission(NamedTuple):
    temperature: jnp.ndarray
    energy_rate: jnp.ndarray
    photon_rate: jnp.ndarray
    outside_rate: jnp.ndarray
    invalid: jnp.ndarray


def prepare_excess_table(closure: DustThermalClosure, background_k: float) -> ExcessTable:
    """Validate a piecewise-linear emission curve parameterized by total power.

    All channels use the same interpolation parameter, so their differential
    powers sum to the total increment even next to the background floor.
    """
    if not jax.config.x64_enabled:
        raise ValueError("IR reference transport requires JAX float64")
    temperature = closure.temperature_k
    if not np.isfinite(background_k) or not temperature[0] <= background_k <= temperature[-1]:
        raise ValueError("IR thermal table does not cover CMB temperature")
    indices = closure.ir_group_indices
    if np.any(closure.group_edges_ev[indices + 1] > 1.0):
        raise ValueError("IR transport only admits complete groups below 1 eV")
    power = closure.emitted_power_per_h_erg_s
    if not np.isfinite(power).all() or np.any(power <= 0) or np.any(np.diff(power) <= 0):
        raise ValueError("IR total power must be finite, positive and strictly increasing")
    band = power[:, None] * closure.ir_energy_fraction
    photons = band / (EV_ERG * closure.ir_mean_photon_energy_ev)
    outside = power * closure.untracked_energy_fraction
    for array in (band, photons, outside):
        if not np.isfinite(array).all() or np.any(np.diff(array, axis=0) < 0):
            raise ValueError("differential IR channels must be finite and monotone")
    if not np.allclose(band.sum(axis=1) + outside, power, rtol=1e-10, atol=0):
        raise ValueError("differential IR channels do not close")
    background = np.interp(np.log(background_k), np.log(temperature), power)
    return ExcessTable(*(jnp.asarray(x, dtype=jnp.float64) for x in
                         (power, np.log(temperature), band, photons, outside, background)))


def excess_emission(table: ExcessTable, heating: jnp.ndarray,
                    dust_density: jnp.ndarray) -> ExcessEmission:
    """Integrate channel slopes over a power increment without bath subtraction."""
    active = (dust_density > 0) & (heating > 0)
    safe_density = jnp.where(dust_density > 0, dust_density, 1.0)
    increment = jnp.where(active, heating / safe_density, 0.0)
    # Express both ends relative to the bath, not bath+tiny_increment. This
    # preserves increments much smaller than machine epsilon times P_CMB.
    starts = jnp.maximum(table.power[:-1] - table.background_power, 0.0)
    ends = jnp.maximum(table.power[1:] - table.background_power, 0.0)
    widths = jnp.maximum(jnp.minimum(increment[..., None], ends) - starts, 0.0)
    delta_power = jnp.diff(table.power)

    def integrate(values):
        slopes = jnp.diff(values, axis=0) / (
            delta_power[:, None] if values.ndim == 2 else delta_power)
        result = jnp.einsum("...t,tg->g...", widths, slopes) if values.ndim == 2 else widths @ slopes
        return result * dust_density

    log_t = jnp.interp(table.background_power + increment, table.power, table.log_temperature)
    invalid = (~jnp.isfinite(heating) | ~jnp.isfinite(dust_density)
               | (heating < 0) | (dust_density < 0)
               | ((heating > 0) & (dust_density == 0))
               | (increment > table.power[-1] - table.background_power))
    return ExcessEmission(jnp.where(active, jnp.exp(log_t), 0.0),
                          integrate(table.band_power), integrate(table.band_photons),
                          integrate(table.outside_power), invalid)


def outward_energy(config, directions, weights, energy):
    """Outgoing face-flux integral in erg for the explicit transport stage."""
    total = jnp.asarray(0.0, dtype=energy.dtype)
    for axis in range(3):
        area = np.prod([config.cell_width[k] for k in range(3) if k != axis])
        for side, sign in ((0, -1), (-1, 1)):
            face = jnp.take(energy, side, axis=axis + 2)
            velocity = jnp.maximum(sign * directions[:, axis], 0.0)
            total += config.dt * config.reduced_light_speed * area * jnp.sum(
                face * (weights * velocity)[None, :, None, None])
    return total


class IRStep(NamedTuple):
    energy: jnp.ndarray
    temperature: jnp.ndarray
    emitted_photons: jnp.ndarray
    outside_energy: jnp.ndarray
    absorbed_energy: jnp.ndarray
    escaped_energy: jnp.ndarray
    balance_residual: jnp.ndarray
    balance_relative: jnp.ndarray
    local_relative: jnp.ndarray
    iterations: jnp.ndarray
    thermal_invalid: jnp.ndarray
    valid: jnp.ndarray


def build_ir_step(config: TransportConfig, directions, weights, table: ExcessTable,
                  dust_density, absorption, *, tolerance=1e-9, max_iterations=128):
    """Build a bounded implicit re-emission step for frozen band opacity.

    Each iterate calls the transport kernel on exactly the same old field;
    only its constant local source changes. Reabsorbed energy is recycled
    within the step, rather than injected as new primary heating.
    """
    dirs, w = np.asarray(directions), np.asarray(weights)
    density, opacity = np.asarray(dust_density), np.asarray(absorption)
    if (dirs.ndim != 2 or dirs.shape[1] != 3 or w.shape != (len(dirs),)
            or not np.isfinite(dirs).all() or not np.isfinite(w).all()
            or np.any(w <= 0) or not np.isclose(w.sum(), 1, rtol=0, atol=1e-12)
            or not np.allclose(np.linalg.norm(dirs, axis=1), 1, atol=1e-6)):
        raise ValueError("IR angular directions/weights are invalid")
    if (len(config.cell_width) != 3 or not np.isfinite(config.cell_width).all()
            or min(config.cell_width) <= 0 or not np.isfinite(config.dt) or config.dt <= 0
            or not np.isfinite(config.reduced_light_speed) or config.reduced_light_speed <= 0
            or not np.isfinite(tolerance) or not 0 < tolerance < 1 or max_iterations < 1):
        raise ValueError("invalid IR timestep or iteration controls")
    cfl = config.dt * config.reduced_light_speed * np.max(
        (np.abs(dirs) / np.asarray(config.cell_width)).sum(axis=1))
    if cfl > 1 + 1e-12:
        raise ValueError("IR transport CFL exceeds one")
    if (density.ndim != 3 or min(density.shape) < 1
            or opacity.shape != (table.band_power.shape[1], *density.shape)
            or not np.isfinite(density).all() or not np.isfinite(opacity).all()
            or np.any(density < 0) or np.any(opacity < 0)
            or np.any(opacity[:, density == 0] != 0)):
        raise ValueError("IR density or band absorption is invalid")
    dust_density, absorption = jnp.asarray(density), jnp.asarray(opacity)
    directions, weights = jnp.asarray(dirs), jnp.asarray(w)
    volume = float(np.prod(config.cell_width))

    @jax.jit
    def step(energy, primary):
        if energy.shape != (opacity.shape[0], len(dirs), *density.shape) or primary.shape != density.shape:
            raise ValueError("IR state/source shapes are inconsistent")
        escaped = outward_energy(config, directions, weights, energy)
        injected = jnp.sum(primary) * config.dt * volume
        old_total = jnp.sum(angular_integral(energy, weights)) * volume
        scale = jnp.maximum(injected, jnp.finfo(energy.dtype).tiny)
        # No primary source: use the existing radiation inventory as scale.
        scale = jnp.where(injected > 0, scale, jnp.maximum(old_total, scale))

        def evaluate(guess):
            emission = excess_emission(table, primary + guess / config.dt, dust_density)
            next_energy, absorbed_directional = advance_with_absorption(
                config, directions, energy, emission.energy_rate / weights.sum(), absorption)
            absorbed = angular_integral(absorbed_directional, weights).sum(axis=0)
            outside = emission.outside_rate * config.dt
            remaining = jnp.sum(angular_integral(next_energy, weights)) * volume
            balance = remaining - old_total + escaped + jnp.sum(outside) * volume - injected
            local_scale = jnp.maximum(primary * config.dt + absorbed, jnp.finfo(energy.dtype).tiny)
            local = jnp.max(jnp.abs(absorbed - guess) / local_scale)
            valid = (jnp.all(~emission.invalid) & jnp.all(jnp.isfinite(next_energy))
                     & jnp.all(next_energy >= 0) & jnp.all(jnp.isfinite(primary))
                     & jnp.all(primary >= 0) & jnp.all(jnp.isfinite(energy)) & jnp.all(energy >= 0))
            return IRStep(next_energy, emission.temperature, emission.photon_rate * config.dt,
                          outside, absorbed, escaped, balance, jnp.abs(balance) / scale,
                          local, jnp.asarray(0), jnp.any(emission.invalid), valid)

        initial_guess = jnp.zeros_like(primary)
        initial = evaluate(initial_guess)

        def cond(carry):
            iteration, _, result = carry
            return ((iteration < max_iterations) & result.valid
                    & ((result.local_relative > tolerance) | (result.balance_relative > tolerance)))

        def body(carry):
            iteration, guess, result = carry
            guess = 0.5 * (guess + result.absorbed_energy)
            return iteration + 1, guess, evaluate(guess)

        count, _, result = jax.lax.while_loop(cond, body, (jnp.asarray(0), initial_guess, initial))
        return result._replace(iterations=count, valid=(result.valid
            & (result.local_relative <= tolerance) & (result.balance_relative <= tolerance)))

    return step
