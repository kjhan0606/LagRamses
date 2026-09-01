"""Furlanetto--Stoever fast-electron energy deposition."""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import jax.numpy as jnp
import numpy as np


HYDROGEN_I_IONIZATION_ENERGY_EV = 13.60
HELIUM_I_IONIZATION_ENERGY_EV = 24.59
HELIUM_II_IONIZATION_ENERGY_EV = 54.42

_IONIZED_FRACTION_GRID = np.asarray(
    [
        1.0e-4,
        2.318e-4,
        4.677e-4,
        1.0e-3,
        2.318e-3,
        4.677e-3,
        1.0e-2,
        2.318e-2,
        4.677e-2,
        1.0e-1,
        0.5,
        0.9,
        0.99,
        0.999,
    ],
    dtype=np.float64,
)
_TABLE_FILENAMES = (
    "log_xi_-4.0.dat",
    "log_xi_-3.6.dat",
    "log_xi_-3.3.dat",
    "log_xi_-3.0.dat",
    "log_xi_-2.6.dat",
    "log_xi_-2.3.dat",
    "log_xi_-2.0.dat",
    "log_xi_-1.6.dat",
    "log_xi_-1.3.dat",
    "log_xi_-1.0.dat",
    "xi_0.500.dat",
    "xi_0.900.dat",
    "xi_0.990.dat",
    "xi_0.999.dat",
)
_TABLE_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "data" / "furlanetto_stoever_2010"
)


class SecondaryFractions(NamedTuple):
    """Fractions of one photoelectron's kinetic energy sent to each channel."""

    heating: jnp.ndarray
    hydrogen_i_ionization: jnp.ndarray
    helium_i_ionization: jnp.ndarray
    helium_ii_ionization: jnp.ndarray
    excitation: jnp.ndarray


class _FurlanettoStoeverTables(NamedTuple):
    energy_ev: np.ndarray
    total_ionization: np.ndarray
    heating: np.ndarray
    excitation: np.ndarray
    hydrogen_i_ionizations: np.ndarray
    helium_i_ionizations: np.ndarray
    helium_ii_ionizations: np.ndarray


def _load_furlanetto_stoever_tables() -> _FurlanettoStoeverTables:
    """Load and validate the 21cmFAST distribution of the FS2010 tables."""

    rows = []
    for filename in _TABLE_FILENAMES:
        path = _TABLE_DIRECTORY / filename
        values = np.loadtxt(path, skiprows=3, dtype=np.float64)
        if values.shape != (258, 9):
            raise RuntimeError(f"invalid Furlanetto--Stoever table shape in {path}: {values.shape}")
        if not np.isfinite(values).all():
            raise RuntimeError(f"non-finite Furlanetto--Stoever table value in {path}")
        rows.append(values)

    tables = np.stack(rows, axis=0)
    energy_ev = tables[0, :, 0]
    if not np.all(np.diff(energy_ev) > 0.0):
        raise RuntimeError("Furlanetto--Stoever energy grid is not strictly increasing")
    if not all(np.array_equal(table[:, 0], energy_ev) for table in tables):
        raise RuntimeError("Furlanetto--Stoever files do not share one energy grid")

    result = _FurlanettoStoeverTables(
        energy_ev=energy_ev,
        total_ionization=tables[:, :, 1],
        heating=tables[:, :, 2],
        excitation=tables[:, :, 3],
        hydrogen_i_ionizations=tables[:, :, 5],
        helium_i_ionizations=tables[:, :, 6],
        helium_ii_ionizations=tables[:, :, 7],
    )
    for value in result:
        value.setflags(write=False)
    return result


_TABLES = _load_furlanetto_stoever_tables()


def _bilinear_interpolate(
    values: np.ndarray,
    energy_ev: jnp.ndarray,
    ionized_fraction: jnp.ndarray,
) -> jnp.ndarray:
    """Linearly interpolate a table in electron energy and ionized fraction."""

    dtype = energy_ev.dtype
    energy_grid = jnp.asarray(_TABLES.energy_ev, dtype=dtype)
    ionized_grid = jnp.asarray(_IONIZED_FRACTION_GRID, dtype=dtype)
    table = jnp.asarray(values, dtype=dtype)

    energy_high = jnp.clip(
        jnp.searchsorted(energy_grid, energy_ev, side="right"),
        1,
        energy_grid.size - 1,
    )
    energy_low = energy_high - 1
    ionized_high = jnp.clip(
        jnp.searchsorted(ionized_grid, ionized_fraction, side="right"),
        1,
        ionized_grid.size - 1,
    )
    ionized_low = ionized_high - 1

    energy_weight = (energy_ev - energy_grid[energy_low]) / (
        energy_grid[energy_high] - energy_grid[energy_low]
    )
    ionized_weight = (ionized_fraction - ionized_grid[ionized_low]) / (
        ionized_grid[ionized_high] - ionized_grid[ionized_low]
    )
    low_ionized_value = (
        table[ionized_low, energy_low] * (1.0 - energy_weight)
        + table[ionized_low, energy_high] * energy_weight
    )
    high_ionized_value = (
        table[ionized_high, energy_low] * (1.0 - energy_weight)
        + table[ionized_high, energy_high] * energy_weight
    )
    return low_ionized_value * (1.0 - ionized_weight) + high_ionized_value * ionized_weight


def furlanetto_stoever_2010(
    electron_energy_ev: jnp.ndarray,
    hydrogen_ionized_fraction: jnp.ndarray,
) -> SecondaryFractions:
    """Interpolate FS2010 energy fractions for a primordial H/He gas.

    The tables assume equal H II and He II fractions, negligible He III, and a
    primordial abundance ratio. Their tabulated parameter is therefore the H II
    fraction rather than the electron fraction. Values are bilinearly
    interpolated in energy and ionized fraction, matching the authors' public
    21cmFAST implementation. Below 10 eV, all energy becomes heat. Above the
    9937.21 eV table limit, the asymptotic energy fractions are held fixed.

    The table's total ionization fraction is split among H I, He I, and He II
    using its tabulated ionization counts and threshold energies. All five
    channels are normalized together to enforce exact local energy closure
    despite Monte Carlo/table-rounding residuals.
    """

    dtype = jnp.result_type(electron_energy_ev, hydrogen_ionized_fraction, jnp.float32)
    energy = jnp.asarray(electron_energy_ev, dtype=dtype)
    ionized = jnp.asarray(hydrogen_ionized_fraction, dtype=dtype)
    while energy.ndim < ionized.ndim + 1:
        energy = energy[..., None]
    ionized = ionized[None, ...]

    minimum_energy = jnp.asarray(_TABLES.energy_ev[0], dtype=dtype)
    maximum_energy = jnp.asarray(_TABLES.energy_ev[-1], dtype=dtype)
    query_energy = jnp.clip(energy, minimum_energy, maximum_energy)
    query_ionized = jnp.clip(
        ionized,
        jnp.asarray(_IONIZED_FRACTION_GRID[0], dtype=dtype),
        jnp.asarray(_IONIZED_FRACTION_GRID[-1], dtype=dtype),
    )

    total_ionization = jnp.maximum(
        _bilinear_interpolate(_TABLES.total_ionization, query_energy, query_ionized),
        0.0,
    )
    heating = jnp.maximum(
        _bilinear_interpolate(_TABLES.heating, query_energy, query_ionized),
        0.0,
    )
    excitation = jnp.maximum(
        _bilinear_interpolate(_TABLES.excitation, query_energy, query_ionized),
        0.0,
    )
    hydrogen_count = jnp.maximum(
        _bilinear_interpolate(_TABLES.hydrogen_i_ionizations, query_energy, query_ionized),
        0.0,
    )
    helium_i_count = jnp.maximum(
        _bilinear_interpolate(_TABLES.helium_i_ionizations, query_energy, query_ionized),
        0.0,
    )
    helium_ii_count = jnp.maximum(
        _bilinear_interpolate(_TABLES.helium_ii_ionizations, query_energy, query_ionized),
        0.0,
    )

    hydrogen_weight = hydrogen_count * HYDROGEN_I_IONIZATION_ENERGY_EV
    helium_i_weight = helium_i_count * HELIUM_I_IONIZATION_ENERGY_EV
    helium_ii_weight = helium_ii_count * HELIUM_II_IONIZATION_ENERGY_EV
    ionization_weight = hydrogen_weight + helium_i_weight + helium_ii_weight
    safe_ionization_weight = jnp.maximum(
        ionization_weight,
        jnp.finfo(dtype).tiny,
    )
    hydrogen_ionization = total_ionization * hydrogen_weight / safe_ionization_weight
    helium_i_ionization = total_ionization * helium_i_weight / safe_ionization_weight
    helium_ii_ionization = total_ionization * helium_ii_weight / safe_ionization_weight

    total_fraction = (
        heating
        + hydrogen_ionization
        + helium_i_ionization
        + helium_ii_ionization
        + excitation
    )
    safe_total_fraction = jnp.maximum(total_fraction, jnp.finfo(dtype).tiny)
    active = energy >= minimum_energy
    return SecondaryFractions(
        heating=jnp.where(active, heating / safe_total_fraction, 1.0),
        hydrogen_i_ionization=jnp.where(
            active,
            hydrogen_ionization / safe_total_fraction,
            0.0,
        ),
        helium_i_ionization=jnp.where(
            active,
            helium_i_ionization / safe_total_fraction,
            0.0,
        ),
        helium_ii_ionization=jnp.where(
            active,
            helium_ii_ionization / safe_total_fraction,
            0.0,
        ),
        excitation=jnp.where(active, excitation / safe_total_fraction, 0.0),
    )
