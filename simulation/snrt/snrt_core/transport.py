"""JAX/XLA-friendly explicit multigroup discrete-ordinates transport."""

from __future__ import annotations

from dataclasses import dataclass

import jax
import jax.numpy as jnp
from jax import lax


@dataclass(frozen=True)
class TransportConfig:
    """Static transport parameters in one internally consistent unit system."""

    cell_width: tuple[float, float, float]
    dt: float
    reduced_light_speed: float


def initial_intensity(
    number_of_groups: int,
    number_of_directions: int,
    shape: tuple[int, int, int],
    dtype=jnp.float32,
) -> jnp.ndarray:
    """Create the zero-radiation initial state [group, direction, x, y, z]."""

    return jnp.zeros((number_of_groups, number_of_directions, *shape), dtype=dtype)


def _neighbors_with_vacuum(field: jnp.ndarray, axis: int) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return lower and upper neighbors with zero-intensity outer boundaries."""

    if axis == 2:
        lower = jnp.concatenate((jnp.zeros_like(field[:, :, :1]), field[:, :, :-1]), axis=2)
        upper = jnp.concatenate((field[:, :, 1:], jnp.zeros_like(field[:, :, :1])), axis=2)
    elif axis == 3:
        lower = jnp.concatenate((jnp.zeros_like(field[:, :, :, :1]), field[:, :, :, :-1]), axis=3)
        upper = jnp.concatenate((field[:, :, :, 1:], jnp.zeros_like(field[:, :, :, :1])), axis=3)
    elif axis == 4:
        lower = jnp.concatenate((jnp.zeros_like(field[:, :, :, :, :1]), field[:, :, :, :, :-1]), axis=4)
        upper = jnp.concatenate((field[:, :, :, :, 1:], jnp.zeros_like(field[:, :, :, :, :1])), axis=4)
    else:
        raise ValueError(f"unsupported spatial axis: {axis}")
    return lower, upper


def _transport_rhs(
    intensity: jnp.ndarray,
    directions: jnp.ndarray,
    cell_width: tuple[float, float, float],
    reduced_light_speed: float,
) -> jnp.ndarray:
    rhs = jnp.zeros_like(intensity)
    for direction_axis, field_axis in enumerate((2, 3, 4)):
        lower, upper = _neighbors_with_vacuum(intensity, field_axis)
        backward = (intensity - lower) / cell_width[direction_axis]
        forward = (upper - intensity) / cell_width[direction_axis]
        velocity = directions[:, direction_axis][None, :, None, None, None]
        gradient = jnp.where(velocity >= 0.0, backward, forward)
        rhs = rhs - reduced_light_speed * velocity * gradient
    return rhs


def cfl_number(config: TransportConfig, directions: jnp.ndarray) -> jnp.ndarray:
    """Return the maximum directional Courant number for one explicit step."""

    inverse_width = jnp.asarray([1.0 / value for value in config.cell_width], dtype=directions.dtype)
    directional_sum = jnp.sum(jnp.abs(directions) * inverse_width[None, :], axis=1)
    return config.reduced_light_speed * config.dt * jnp.max(directional_sum)


def angular_integral(
    directional_field: jnp.ndarray,
    weights: jnp.ndarray,
    *,
    use_tensor_core_reduction: bool = False,
) -> jnp.ndarray:
    """Integrate a ``[group, direction, x, y, z]`` field over direction.

    The default is the direct reduction used by the reference solver.  The
    Tensor Core path partitions directions into 16 weighted angular bins and
    evaluates the contraction as a static GEMM.  Summing those bins restores
    the identical discrete angular integral while making the contraction shape
    suitable for TF32 Tensor Cores on NVIDIA GPUs.
    """

    if not use_tensor_core_reduction:
        return jnp.einsum("d,gdxyz->gxyz", weights, directional_field)

    number_of_directions = directional_field.shape[1]
    number_of_bins = 16
    bin_index = jnp.arange(number_of_directions) * number_of_bins // number_of_directions
    projection = jax.nn.one_hot(bin_index, number_of_bins, dtype=directional_field.dtype)
    projection = projection * weights[:, None]
    flattened = jnp.moveaxis(directional_field, 1, -1).reshape((-1, number_of_directions))
    binned = jnp.matmul(flattened, projection, precision=lax.Precision.DEFAULT)
    return jnp.sum(binned, axis=-1).reshape((directional_field.shape[0], *directional_field.shape[2:]))


def advance_with_absorption(
    config: TransportConfig,
    directions: jnp.ndarray,
    intensity: jnp.ndarray,
    emissivity: jnp.ndarray,
    absorption: jnp.ndarray,
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Advance one S_N step and return the exact local absorption loss.

    ``intensity`` is ordered [group, direction, x, y, z]. ``emissivity`` and
    ``absorption`` are ordered [group, x, y, z]. Transport and emission use an
    explicit step; the local absorption operator is applied exactly. Positivity
    therefore requires only the directional transport CFL condition.
    """

    transport = _transport_rhs(
        intensity,
        directions,
        config.cell_width,
        config.reduced_light_speed,
    )
    source = emissivity[:, None, :, :, :]
    pre_absorption = intensity + config.dt * (transport + source)
    optical_depth = config.reduced_light_speed * config.dt * absorption[:, None, :, :, :]
    next_intensity = pre_absorption * jnp.exp(-optical_depth)
    return next_intensity, pre_absorption - next_intensity


def advance_explicit(
    config: TransportConfig,
    directions: jnp.ndarray,
    intensity: jnp.ndarray,
    emissivity: jnp.ndarray,
    absorption: jnp.ndarray,
) -> jnp.ndarray:
    """Advance an S_N system with explicit transport and exact local absorption."""
    next_intensity, _ = advance_with_absorption(config, directions, intensity, emissivity, absorption)
    return next_intensity


def build_explicit_step(config: TransportConfig, directions: jnp.ndarray):
    """Build a statically shaped JIT-compiled transport step for TPU execution."""

    directions = jnp.asarray(directions)

    @jax.jit
    def step(intensity: jnp.ndarray, emissivity: jnp.ndarray, absorption: jnp.ndarray) -> jnp.ndarray:
        return advance_explicit(config, directions, intensity, emissivity, absorption)

    return step


def radiation_moments(
    intensity: jnp.ndarray,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    reduced_light_speed: float,
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return angular number density and flux from an S_N intensity field."""

    angular_weight = weights[None, :, None, None, None]
    number_density = jnp.sum(intensity * angular_weight, axis=1)
    flux_components = []
    for axis in range(3):
        direction = directions[:, axis][None, :, None, None, None]
        flux_components.append(reduced_light_speed * jnp.sum(intensity * angular_weight * direction, axis=1))
    return number_density, jnp.stack(flux_components, axis=1)
