"""Dust absorption primitives and audited group-opacity loading."""

from __future__ import annotations

import json
from pathlib import Path
from typing import NamedTuple

import jax.numpy as jnp
import numpy as np


EV_ERG = 1.602176634e-12
LIGHT_SPEED_CM_S = 2.99792458e10


class DustModel(NamedTuple):
    """Absorbing dust scaled to a supplied reference cross section per H.

    ``absorption_cross_section_per_h`` is [group] in cm^2 per H nucleus for
    the reference dust mixture. ``relative_abundance`` is a non-negative
    cell field relative to that mixture; it may encode metallicity and a
    dust-to-metal prescription outside the transport kernel.
    ``absorption_weighted_energy_ev`` is the per-group photon energy used for
    dust heating and absorption momentum; it is optional for compatibility
    with legacy callers, which fall back to the supplied transport group
    energy.
    """

    absorption_cross_section_per_h: jnp.ndarray
    relative_abundance: jnp.ndarray
    absorption_weighted_energy_ev: jnp.ndarray | None = None


class DustOpacityClosure(NamedTuple):
    """Validated, group-averaged dust opacity metadata.

    The cross section is per H nucleus for the reference dust mixture.  The
    energy is weighted by the same dust absorption opacity and is therefore
    the appropriate energy per absorbed dust photon for local heating.
    """

    group_edges_ev: np.ndarray
    absorption_cross_section_per_h_cm2: np.ndarray
    absorption_weighted_energy_ev: np.ndarray


def read_dust_opacity_metadata(
    path: str | Path,
    *,
    expected_group_edges_ev: np.ndarray | None = None,
) -> DustOpacityClosure:
    """Read and validate a source-SED-dependent dust opacity closure.

    The JSON schema is ``snrt_dust_opacity_v1`` and requires
    ``group_edges_ev``, ``absorption_cross_section_per_h_cm2``, and
    ``absorption_weighted_energy_ev`` plus non-empty ``reference_mixture``,
    ``opacity_source``, and ``spectral_weighting`` provenance strings.  No
    opacity normalization or group ordering is inferred.
    ``expected_group_edges_ev`` should be supplied by the photon-ledger
    metadata when the closure is attached to a run.
    """

    opacity_path = Path(path)
    with opacity_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    if metadata.get("schema") != "snrt_dust_opacity_v1":
        raise ValueError(f"{opacity_path}: unsupported dust opacity schema")
    required = (
        "group_edges_ev",
        "absorption_cross_section_per_h_cm2",
        "absorption_weighted_energy_ev",
    )
    missing = [name for name in required if name not in metadata]
    if missing:
        raise ValueError(f"{opacity_path}: missing dust opacity fields {missing}")
    for name in ("reference_mixture", "opacity_source", "spectral_weighting"):
        if not isinstance(metadata.get(name), str) or not metadata[name].strip():
            raise ValueError(f"{opacity_path}: {name} must be a non-empty provenance string")

    edges = np.asarray(metadata["group_edges_ev"], dtype=np.float64)
    cross_section = np.asarray(metadata["absorption_cross_section_per_h_cm2"], dtype=np.float64)
    weighted_energy = np.asarray(metadata["absorption_weighted_energy_ev"], dtype=np.float64)
    if (
        edges.ndim != 1
        or edges.size < 2
        or not np.isfinite(edges).all()
        or np.any(edges <= 0.0)
        or np.any(np.diff(edges) <= 0.0)
    ):
        raise ValueError(f"{opacity_path}: group edges must be finite, positive, and increasing")
    number_of_groups = edges.size - 1
    if cross_section.shape != (number_of_groups,) or weighted_energy.shape != (number_of_groups,):
        raise ValueError(f"{opacity_path}: dust arrays do not match the number of groups")
    if (
        not np.isfinite(cross_section).all()
        or np.any(cross_section < 0.0)
        or not np.isfinite(weighted_energy).all()
        or np.any(weighted_energy <= 0.0)
    ):
        raise ValueError(f"{opacity_path}: dust cross sections/energies are invalid")
    tolerance = 1.0e-12 * np.maximum(1.0, edges[1:])
    if np.any(weighted_energy < edges[:-1] - tolerance) or np.any(weighted_energy > edges[1:] + tolerance):
        raise ValueError(f"{opacity_path}: absorption-weighted energies lie outside their groups")
    if expected_group_edges_ev is not None:
        expected = np.asarray(expected_group_edges_ev, dtype=np.float64)
        if expected.shape != edges.shape or not np.allclose(expected, edges, rtol=0.0, atol=1.0e-12):
            raise ValueError(f"{opacity_path}: dust groups do not match the photon-ledger groups")
    return DustOpacityClosure(edges, cross_section, weighted_energy)


def dust_model_from_metadata(
    path: str | Path,
    relative_abundance: jnp.ndarray,
    *,
    dtype: jnp.dtype = jnp.float32,
    expected_group_edges_ev: np.ndarray | None = None,
) -> DustModel:
    """Build a JAX dust model from a validated opacity sidecar."""

    closure = read_dust_opacity_metadata(path, expected_group_edges_ev=expected_group_edges_ev)
    abundance = jnp.asarray(relative_abundance, dtype=dtype)
    if not np.isfinite(np.asarray(abundance)).all() or np.any(np.asarray(abundance) < 0.0):
        raise ValueError("relative dust abundance must be finite and non-negative")
    return DustModel(
        absorption_cross_section_per_h=jnp.asarray(closure.absorption_cross_section_per_h_cm2, dtype=dtype),
        relative_abundance=abundance,
        absorption_weighted_energy_ev=jnp.asarray(closure.absorption_weighted_energy_ev, dtype=dtype),
    )


def absorption_coefficient(n_hydrogen: jnp.ndarray, dust: DustModel) -> jnp.ndarray:
    """Return dust absorption coefficient [group, cell] in cm^-1."""
    extra_axes = (1,) * n_hydrogen.ndim
    cross_section = dust.absorption_cross_section_per_h.reshape((-1,) + extra_axes)
    return cross_section * n_hydrogen[None, ...] * jnp.maximum(dust.relative_abundance[None, ...], 0.0)


def absorbed_dust_momentum_rate(
    absorbed_intensity: jnp.ndarray,
    dust_fraction: jnp.ndarray,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    absorption_weighted_energy_ev: jnp.ndarray,
    dt: float,
) -> jnp.ndarray:
    """Return dust absorption momentum deposition per volume and time.

    ``absorbed_intensity`` is the directional photon-number loss over one
    transport step, ordered ``[group, direction, x, y, z]``.  The returned
    vector is in dyn cm^-3 and uses the physical speed of light for photon
    momentum, even when transport uses a reduced light speed.  This is an
    absorption-only force diagnostic; scattering and IR re-emission are not
    included.
    """

    absorbed = jnp.asarray(absorbed_intensity)
    fraction = jnp.asarray(dust_fraction)
    direction = jnp.asarray(directions)
    angular_weight = jnp.asarray(weights)
    energy = jnp.asarray(absorption_weighted_energy_ev)
    if absorbed.ndim != 5 or fraction.shape != (absorbed.shape[0], *absorbed.shape[2:]):
        raise ValueError("absorbed intensity/dust fraction must have shapes [group,direction,x,y,z] and [group,x,y,z]")
    if direction.shape != (absorbed.shape[1], 3) or angular_weight.shape != (absorbed.shape[1],):
        raise ValueError("directions and weights must match the directional intensity axis")
    if energy.shape != (absorbed.shape[0],) or not np.isfinite(np.asarray(energy)).all() or np.any(np.asarray(energy) <= 0.0):
        raise ValueError("absorption-weighted photon energy must have shape (group,) and be positive")
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be positive and finite")
    energy_shape = (-1, 1, 1, 1, 1)
    dust_absorbed_energy = absorbed * fraction[:, None, ...] * energy.reshape(energy_shape) * EV_ERG
    return jnp.einsum(
        "d,gdxyz,di->igxyz",
        angular_weight,
        dust_absorbed_energy / (LIGHT_SPEED_CM_S * dt),
        direction,
    ).sum(axis=1)


def zero_dust(number_of_groups: int, shape: tuple[int, ...], dtype: jnp.dtype = jnp.float32) -> DustModel:
    """Return a no-dust model compatible with a static transport shape."""
    return DustModel(
        absorption_cross_section_per_h=jnp.zeros((number_of_groups,), dtype=dtype),
        relative_abundance=jnp.zeros(shape, dtype=dtype),
    )
