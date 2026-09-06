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
    ``absorption`` are ordered [group, x, y, z]. Transport is explicit. After
    that transport stage, the local constant-source/absorption problem is
    integrated analytically over the step. This keeps source-cell emission
    time-centred with respect to local absorption and returns the corresponding
    local photon loss. Positivity therefore requires only the directional
    transport CFL condition.
    """

    transport = _transport_rhs(
        intensity,
        directions,
        config.cell_width,
        config.reduced_light_speed,
    )
    source = emissivity[:, None, :, :, :]
    transported = intensity + config.dt * transport
    optical_depth = config.reduced_light_speed * config.dt * absorption[:, None, :, :, :]
    transmission = jnp.exp(-optical_depth)
    one_minus_transmission = -jnp.expm1(-optical_depth)

    # phi(tau) = (1 - exp(-tau)) / tau is the exact source response. The
    # A second-order polynomial branch avoids losing the source contribution
    # to cancellation when a cell is optically thin during one step. Keep the
    # polynomial quadratic: XLA may reassociate powers of a very large dt, and
    # a cubic term can overflow in float32 even when the resulting tau is zero.
    safe_optical_depth = jnp.maximum(optical_depth, jnp.finfo(optical_depth.dtype).tiny)
    source_response = one_minus_transmission / safe_optical_depth
    source_response_series = (
        1.0
        - 0.5 * optical_depth
        + optical_depth * optical_depth / 6.0
    )
    source_response = jnp.where(optical_depth < 1.0e-4, source_response_series, source_response)

    next_intensity = transported * transmission + config.dt * source * source_response
    absorbed_intensity = (
        transported * one_minus_transmission
        + config.dt * source * (1.0 - source_response)
    )
    return next_intensity, absorbed_intensity


def _exponential_integral(rate: jnp.ndarray, dt: float) -> jnp.ndarray:
    """Return ``(1-exp(-rate*dt))/rate`` with a thin-cell branch."""

    safe_rate = jnp.maximum(rate, jnp.finfo(rate.dtype).tiny)
    optical_depth = rate * dt
    exact = -jnp.expm1(-optical_depth) / safe_rate
    series = dt * (1.0 - 0.5 * optical_depth + optical_depth * optical_depth / 6.0)
    return jnp.where(optical_depth < 1.0e-4, series, exact)


def _exponential_double_integral(rate: jnp.ndarray, dt: float) -> jnp.ndarray:
    """Return ``(dt-F(rate))/rate`` for a constant local source."""

    safe_rate = jnp.maximum(rate, jnp.finfo(rate.dtype).tiny)
    optical_depth = rate * dt
    first = _exponential_integral(rate, dt)
    exact = (dt - first) / safe_rate
    # The quadratic branch avoids an unnecessary high power in float32 while
    # retaining the correct limit dt**2/2 as rate -> 0.
    series = dt * dt * (0.5 - optical_depth / 6.0 + optical_depth * optical_depth / 24.0)
    return jnp.where(optical_depth < 1.0e-4, series, exact)


def advance_with_isotropic_scattering(
    config: TransportConfig,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    intensity: jnp.ndarray,
    emissivity: jnp.ndarray,
    absorption: jnp.ndarray,
    scattering: jnp.ndarray,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Advance transport with exact local isotropic within-group scattering.

    The explicit upwind transport stage is unchanged.  In the subsequent
    local constant-coefficient solve, ``absorption`` and ``scattering`` are
    coefficients in cm^-1 and ``emissivity`` is an isotropic source.  The
    returned arrays are, respectively, the next intensity, absorption loss,
    incoming scattering events, and isotropic outgoing scattering events.  All
    three event arrays are directional photon-number integrals over the step;
    ``absorbed`` is intentionally absorption-only for the H/He ledgers.

    For the isotropic closure, the angular mean obeys
    ``dJ/dt = -c_hat*kappa_abs*J + source``.  Thus scattering redistributes
    directions but cannot remove photons or add dust heating.
    """

    transport = _transport_rhs(
        intensity,
        directions,
        config.cell_width,
        config.reduced_light_speed,
    )
    transported = intensity + config.dt * transport
    source = emissivity[:, None, :, :, :]
    absorption_rate = config.reduced_light_speed * absorption[:, None, :, :, :]
    scattering_rate = config.reduced_light_speed * scattering[:, None, :, :, :]
    total_rate = absorption_rate + scattering_rate
    absorption_integral = _exponential_integral(absorption_rate, config.dt)
    total_integral = _exponential_integral(total_rate, config.dt)
    absorption_double_integral = _exponential_double_integral(absorption_rate, config.dt)
    transmission_absorption = jnp.exp(-absorption_rate * config.dt)
    transmission_total = jnp.exp(-total_rate * config.dt)

    mean_initial = jnp.einsum("d,gdxyz->gxyz", weights, transported)
    mean_initial_directional = mean_initial[:, None, :, :, :]
    mean_integral = (
        mean_initial * absorption_integral[:, 0, :, :, :]
        + emissivity * absorption_double_integral[:, 0, :, :, :]
    )
    directional_integral = (
        transported * total_integral
        + mean_initial_directional * (absorption_integral - total_integral)
        + source * absorption_double_integral
    )
    next_intensity = (
        transported * transmission_total
        + mean_initial_directional * (transmission_absorption - transmission_total)
        + source * absorption_integral
    )
    absorbed_intensity = absorption_rate * directional_integral
    scattered_incoming = scattering_rate * directional_integral
    scattered_outgoing = jnp.broadcast_to(
        scattering_rate * mean_integral[:, None, :, :, :], directional_integral.shape
    )
    return next_intensity, absorbed_intensity, scattered_incoming, scattered_outgoing


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
