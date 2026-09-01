"""JAX/XLA interpolation of a provenance-validated metal-only atlas."""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp
import numpy as np


class JaxThermalAtlas(NamedTuple):
    """Static 3-D solar-metallicity tables for a compiled S_N kernel.

    ``net_rate_erg_s_cm3`` contains only the UVB-free metal contribution and
    uses the SNRT convention: heating is positive and cooling is negative.
    Non-equilibrium primordial rates are evaluated from the live chemistry.
    """

    scale_factor: jnp.ndarray
    log_hydrogen_number_density_cm3: jnp.ndarray
    log_temperature_k: jnp.ndarray
    net_rate_erg_s_cm3: jnp.ndarray
    mean_molecular_weight: jnp.ndarray


def from_numpy_atlas(atlas, dtype=jnp.float32) -> JaxThermalAtlas:
    """Transfer an already validated host-side thermal atlas to JAX arrays."""

    return JaxThermalAtlas(
        scale_factor=jnp.asarray(np.asarray(atlas.scale_factor), dtype=dtype),
        log_hydrogen_number_density_cm3=jnp.asarray(np.asarray(atlas.log_hydrogen_number_density_cm3), dtype=dtype),
        log_temperature_k=jnp.asarray(np.asarray(atlas.log_temperature_k), dtype=dtype),
        net_rate_erg_s_cm3=jnp.asarray(np.asarray(atlas.net_rate_erg_s_cm3), dtype=dtype),
        mean_molecular_weight=jnp.asarray(np.asarray(atlas.mean_molecular_weight), dtype=dtype),
    )


def _indices(axis: jnp.ndarray, values: jnp.ndarray) -> tuple[jnp.ndarray, jnp.ndarray]:
    clipped = jnp.clip(values, axis[0], axis[-1])
    index = jnp.clip(jnp.searchsorted(axis, clipped, side="right") - 1, 0, len(axis) - 2)
    weight = (clipped - axis[index]) / (axis[index + 1] - axis[index])
    return index, weight


def _interpolate(
    atlas: JaxThermalAtlas,
    values: jnp.ndarray,
    scale_factor: float | jnp.ndarray,
    temperature_k: jnp.ndarray,
    n_hydrogen_cm3: jnp.ndarray,
) -> jnp.ndarray:
    """Interpolate a solar-metallicity (a, n_H, T) table with edge clamping."""

    temperature = jnp.maximum(jnp.asarray(temperature_k), jnp.finfo(jnp.asarray(temperature_k).dtype).tiny)
    density = jnp.maximum(jnp.asarray(n_hydrogen_cm3), jnp.finfo(jnp.asarray(n_hydrogen_cm3).dtype).tiny)
    time_index, time_weight = _indices(atlas.scale_factor, jnp.asarray(scale_factor, dtype=temperature.dtype))
    density_index, density_weight = _indices(atlas.log_hydrogen_number_density_cm3, jnp.log10(density))
    temperature_index, temperature_weight = _indices(atlas.log_temperature_k, jnp.log10(temperature))
    result = jnp.zeros_like(temperature)
    for time_offset in (0, 1):
        time_factor = time_weight if time_offset else 1.0 - time_weight
        for density_offset in (0, 1):
            density_factor = density_weight if density_offset else 1.0 - density_weight
            for temperature_offset in (0, 1):
                temperature_factor = temperature_weight if temperature_offset else 1.0 - temperature_weight
                result = result + (
                    time_factor
                    * density_factor
                    * temperature_factor
                    * values[
                        time_index + time_offset,
                        density_index + density_offset,
                        temperature_index + temperature_offset,
                    ]
                )
    return result


def net_rate(
    atlas: JaxThermalAtlas,
    scale_factor: float | jnp.ndarray,
    temperature_k: jnp.ndarray,
    n_hydrogen_cm3: jnp.ndarray,
    metallicity_solar: float | jnp.ndarray,
) -> jnp.ndarray:
    """Return signed UVB-free metal heating/cooling in erg cm^-3 s^-1."""

    temperature = jnp.asarray(temperature_k)
    temperature, density, metallicity = jnp.broadcast_arrays(
        temperature,
        jnp.asarray(n_hydrogen_cm3, dtype=temperature.dtype),
        jnp.asarray(metallicity_solar, dtype=temperature.dtype),
    )
    metallicity = jnp.where(jnp.isfinite(metallicity) & (metallicity >= 0.0), metallicity, jnp.nan)
    solar_rate = _interpolate(atlas, atlas.net_rate_erg_s_cm3, scale_factor, temperature, density)
    return solar_rate * metallicity


def mean_mu(
    atlas: JaxThermalAtlas,
    scale_factor: float | jnp.ndarray,
    temperature_k: jnp.ndarray,
    n_hydrogen_cm3: jnp.ndarray,
    metallicity_solar: float | jnp.ndarray,
) -> jnp.ndarray:
    """Return the table mean molecular weight with the same interpolation."""

    temperature = jnp.asarray(temperature_k)
    temperature, density, metallicity = jnp.broadcast_arrays(
        temperature,
        jnp.asarray(n_hydrogen_cm3, dtype=temperature.dtype),
        jnp.asarray(metallicity_solar, dtype=temperature.dtype),
    )
    validity = jnp.where(jnp.isfinite(metallicity) & (metallicity >= 0.0), 1.0, jnp.nan)
    return _interpolate(atlas, atlas.mean_molecular_weight, scale_factor, temperature, density) * validity
