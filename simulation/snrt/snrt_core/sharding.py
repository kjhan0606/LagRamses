"""Spatially sharded JAX/XLA execution for static-grid S_N transport."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import jax
import jax.numpy as jnp
import numpy as np
from jax.sharding import Mesh, NamedSharding, PartitionSpec

from .transport import TransportConfig, advance_explicit


@dataclass(frozen=True)
class XShardings:
    """One-dimensional x-domain layouts for intensity and cell fields."""

    mesh: Mesh
    intensity: NamedSharding
    group_field: NamedSharding
    scalar_field: NamedSharding


def make_x_shardings(
    devices: Sequence[jax.Device] | None = None,
    axis_name: str = "x",
) -> XShardings:
    """Create a 1-D spatial mesh; groups and directions remain replicated."""
    selected_devices = tuple(jax.devices() if devices is None else devices)
    if not selected_devices:
        raise ValueError("At least one JAX device is required.")
    mesh = Mesh(np.asarray(selected_devices), (axis_name,))
    intensity = NamedSharding(mesh, PartitionSpec(None, None, axis_name, None, None))
    group_field = NamedSharding(mesh, PartitionSpec(None, axis_name, None, None))
    scalar_field = NamedSharding(mesh, PartitionSpec(axis_name, None, None))
    return XShardings(mesh, intensity, group_field, scalar_field)


def validate_x_partition(shape: tuple[int, int, int], shardings: XShardings) -> None:
    """Reject grids whose x dimension cannot be evenly split across devices."""
    device_count = shardings.mesh.shape[next(iter(shardings.mesh.axis_names))]
    if shape[0] % device_count:
        raise ValueError(f"x dimension {shape[0]} is not divisible by {device_count} devices.")


def place_transport_fields(
    intensity: jnp.ndarray,
    emissivity: jnp.ndarray,
    absorption: jnp.ndarray,
    shardings: XShardings,
) -> tuple[jax.Array, jax.Array, jax.Array]:
    """Place transport fields using the x-domain partition layouts."""
    validate_x_partition(tuple(intensity.shape[2:]), shardings)
    return (
        jax.device_put(intensity, shardings.intensity),
        jax.device_put(emissivity, shardings.group_field),
        jax.device_put(absorption, shardings.group_field),
    )


def build_x_sharded_transport_step(
    config: TransportConfig,
    directions: jnp.ndarray,
    shardings: XShardings,
):
    """Build a compiler-partitioned transport step with x-boundary halo exchange.

    The global-array expression stays identical to the CPU reference. XLA SPMD
    partitions the x dimension and inserts the required neighbor communication
    for the upwind stencil, avoiding direction-dependent Python control flow.
    """

    def step(intensity: jnp.ndarray, emissivity: jnp.ndarray, absorption: jnp.ndarray) -> jnp.ndarray:
        return advance_explicit(config, directions, intensity, emissivity, absorption)

    return jax.jit(
        step,
        in_shardings=(shardings.intensity, shardings.group_field, shardings.group_field),
        out_shardings=shardings.intensity,
    )
