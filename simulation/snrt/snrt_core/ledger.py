"""Photon and energy ledgers for finite-volume S_N transport."""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp

from .transport import TransportConfig, advance_with_absorption, angular_integral, radiation_moments


class PhotonLedger(NamedTuple):
    """Per-group photon accounting over a single explicit transport step."""

    initial: jnp.ndarray
    emitted: jnp.ndarray
    absorbed: jnp.ndarray
    escaped: jnp.ndarray
    final: jnp.ndarray
    residual: jnp.ndarray


def _inventory(number_density: jnp.ndarray, cell_volume: float) -> jnp.ndarray:
    return jnp.sum(number_density, axis=(1, 2, 3)) * cell_volume


def _face_escape(
    face_intensity: jnp.ndarray,
    direction_cosine: jnp.ndarray,
    weights: jnp.ndarray,
    area: float,
    light_speed: float,
    outward_sign: float,
) -> jnp.ndarray:
    """Return outward photon rate through one vacuum boundary face per group."""

    cosine = direction_cosine[None, :, None, None]
    angular_weight = weights[None, :, None, None]
    outward_cosine = jnp.maximum(outward_sign * cosine, 0.0)
    return light_speed * area * jnp.sum(face_intensity * angular_weight * outward_cosine, axis=(1, 2, 3))


def boundary_escape_rate(
    intensity: jnp.ndarray,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    config: TransportConfig,
) -> jnp.ndarray:
    """Return the total outward photon rate through all six vacuum faces."""

    dx, dy, dz = config.cell_width
    escaped = jnp.zeros((intensity.shape[0],), dtype=intensity.dtype)
    escaped = escaped + _face_escape(intensity[:, :, 0, :, :], directions[:, 0], weights, dy * dz, config.reduced_light_speed, -1.0)
    escaped = escaped + _face_escape(intensity[:, :, -1, :, :], directions[:, 0], weights, dy * dz, config.reduced_light_speed, 1.0)
    escaped = escaped + _face_escape(intensity[:, :, :, 0, :], directions[:, 1], weights, dx * dz, config.reduced_light_speed, -1.0)
    escaped = escaped + _face_escape(intensity[:, :, :, -1, :], directions[:, 1], weights, dx * dz, config.reduced_light_speed, 1.0)
    escaped = escaped + _face_escape(intensity[:, :, :, :, 0], directions[:, 2], weights, dx * dy, config.reduced_light_speed, -1.0)
    escaped = escaped + _face_escape(intensity[:, :, :, :, -1], directions[:, 2], weights, dx * dy, config.reduced_light_speed, 1.0)
    return escaped


def photon_ledger(
    config: TransportConfig,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    intensity_before: jnp.ndarray,
    intensity_after: jnp.ndarray,
    emissivity: jnp.ndarray,
    absorption: jnp.ndarray,
) -> PhotonLedger:
    """Account for sources, absorption, boundary escape, and finite-volume residual.

    The absorption term is reconstructed by the exact same
    ``advance_with_absorption`` operation used by the transport step.  In
    particular, it includes photons emitted during this step before the local
    exponential attenuation is applied.  Boundary flux is evaluated from the
    beginning-of-step upwind state, which is the finite-volume flux used by the
    explicit spatial operator.  If ``intensity_after`` does not match the
    supplied operator inputs, that mismatch is exposed by ``residual``.
    """

    _, absorbed_directional = advance_with_absorption(
        config,
        directions,
        intensity_before,
        emissivity,
        absorption,
    )
    absorbed_photons = angular_integral(absorbed_directional, weights)
    return photon_ledger_from_absorbed(
        config,
        directions,
        weights,
        intensity_before,
        intensity_after,
        emissivity,
        absorbed_photons,
    )


def photon_ledger_from_absorbed(
    config: TransportConfig,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    intensity_before: jnp.ndarray,
    intensity_after: jnp.ndarray,
    emissivity: jnp.ndarray,
    absorbed_photons: jnp.ndarray,
) -> PhotonLedger:
    """Account for a supplied local absorption result and boundary escape.

    ``absorbed_photons`` must be the angularly integrated photon loss actually
    used to produce ``intensity_after``.  This variant is required when the
    opacity is time-averaged or iterated inside a coupled chemistry kernel and
    therefore cannot be reconstructed from the caller's old state alone.
    """

    cell_volume = config.cell_width[0] * config.cell_width[1] * config.cell_width[2]
    number_before, _ = radiation_moments(intensity_before, directions, weights, config.reduced_light_speed)
    number_after, _ = radiation_moments(intensity_after, directions, weights, config.reduced_light_speed)
    initial = _inventory(number_before, cell_volume)
    final = _inventory(number_after, cell_volume)
    absorbed = _inventory(absorbed_photons, cell_volume)
    source_directional = jnp.broadcast_to(emissivity[:, None, :, :, :], intensity_before.shape)
    emitted = _inventory(
        config.dt * angular_integral(source_directional, weights),
        cell_volume,
    )
    escaped = config.dt * boundary_escape_rate(intensity_before, directions, weights, config)
    residual = final - initial - emitted + absorbed + escaped
    return PhotonLedger(initial, emitted, absorbed, escaped, final, residual)


def energy_from_photons(photon_count: jnp.ndarray, group_photon_energy: jnp.ndarray) -> jnp.ndarray:
    """Convert a per-group photon inventory or ledger entry to energy."""

    return photon_count * group_photon_energy
