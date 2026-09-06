"""Reader and thermal inversion for the local Grackle equilibrium-grid format."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


HYDROGEN_MASS_FRACTION = 0.76
PROTON_MASS_G = 1.67262192369e-24
BOLTZMANN_ERG_K = 1.380649e-16


@dataclass(frozen=True)
class GrackleEquilibriumTable:
    """NLTE equilibrium values on (log n_H, log Z/Zsun, log T) axes.

    The source table supplies thermal equilibrium and cooling rates, but not
    the H/He ion fractions needed by the S_N chemistry state.  It therefore
    initializes temperature only; radiation chemistry remains in SNRT.
    """

    redshift: float
    log_hydrogen_number_density_cm3: np.ndarray
    log_metallicity_solar: np.ndarray
    log_temperature_k: np.ndarray
    net_rate_erg_s_cm3: np.ndarray
    mean_molecular_weight: np.ndarray

    def __post_init__(self) -> None:
        axes = {
            "log_hydrogen_number_density_cm3": self.log_hydrogen_number_density_cm3,
            "log_metallicity_solar": self.log_metallicity_solar,
            "log_temperature_k": self.log_temperature_k,
        }
        normalized_axes = {name: np.asarray(value, dtype=np.float64) for name, value in axes.items()}
        if any(axis.ndim != 1 or len(axis) < 2 or np.any(np.diff(axis) <= 0.0) for axis in normalized_axes.values()):
            raise ValueError("Grackle table axes must be strictly increasing one-dimensional arrays")
        expected_shape = tuple(len(axis) for axis in normalized_axes.values())
        rates = np.asarray(self.net_rate_erg_s_cm3, dtype=np.float64)
        mu = np.asarray(self.mean_molecular_weight, dtype=np.float64)
        if rates.shape != expected_shape or mu.shape != expected_shape:
            raise ValueError("Grackle data must have shape (n_density, n_metallicity, n_temperature)")
        if not np.isfinite(rates).all() or not np.isfinite(mu).all() or np.any(mu <= 0.0):
            raise ValueError("Grackle table contains invalid rates or mean molecular weights")
        for name, axis in normalized_axes.items():
            object.__setattr__(self, name, axis)
        object.__setattr__(self, "net_rate_erg_s_cm3", rates)
        object.__setattr__(self, "mean_molecular_weight", mu)

    @staticmethod
    def _axis_indices(axis: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        clipped = np.clip(values, axis[0], axis[-1])
        index = np.clip(np.searchsorted(axis, clipped, side="right") - 1, 0, len(axis) - 2)
        weight = (clipped - axis[index]) / (axis[index + 1] - axis[index])
        return index, weight

    def _interpolate(
        self,
        values: np.ndarray,
        temperature_k: np.ndarray,
        n_h_cm3: np.ndarray,
        metallicity_solar: np.ndarray | float,
    ) -> np.ndarray:
        temperature, n_h, metallicity = np.broadcast_arrays(
            np.asarray(temperature_k, dtype=np.float64),
            np.asarray(n_h_cm3, dtype=np.float64),
            np.asarray(metallicity_solar, dtype=np.float64),
        )
        if np.any(temperature <= 0.0) or np.any(n_h <= 0.0) or np.any(metallicity <= 0.0):
            raise ValueError("Grackle interpolation requires positive temperature, n_H, and metallicity")
        density_index, density_weight = self._axis_indices(
            self.log_hydrogen_number_density_cm3, np.log10(n_h)
        )
        metallicity_index, metallicity_weight = self._axis_indices(
            self.log_metallicity_solar, np.log10(metallicity)
        )
        temperature_index, temperature_weight = self._axis_indices(self.log_temperature_k, np.log10(temperature))
        result = np.zeros_like(temperature, dtype=np.float64)
        for density_offset in (0, 1):
            density_factor = density_weight if density_offset else 1.0 - density_weight
            for metallicity_offset in (0, 1):
                metallicity_factor = metallicity_weight if metallicity_offset else 1.0 - metallicity_weight
                for temperature_offset in (0, 1):
                    temperature_factor = temperature_weight if temperature_offset else 1.0 - temperature_weight
                    result += (
                        density_factor
                        * metallicity_factor
                        * temperature_factor
                        * values[
                            density_index + density_offset,
                            metallicity_index + metallicity_offset,
                            temperature_index + temperature_offset,
                        ]
                    )
        return result

    def mean_mu(
        self, temperature_k: np.ndarray, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray | float
    ) -> np.ndarray:
        return self._interpolate(self.mean_molecular_weight, temperature_k, n_h_cm3, metallicity_solar)

    def temperature_from_pressure(
        self,
        pressure_dyn_cm2: np.ndarray,
        mass_density_g_cm3: np.ndarray,
        metallicity_solar: np.ndarray | float,
        *,
        iterations: int = 32,
        hydrogen_mass_fraction: float = HYDROGEN_MASS_FRACTION,
    ) -> np.ndarray:
        """Solve T=P*mu(T,n_H,Z)*m_p/(rho*k_B) with damped fixed-point updates."""

        if iterations < 1 or not 0.0 < hydrogen_mass_fraction < 1.0:
            raise ValueError("invalid Grackle temperature-inversion configuration")
        pressure, density, metallicity = np.broadcast_arrays(
            np.asarray(pressure_dyn_cm2, dtype=np.float64),
            np.asarray(mass_density_g_cm3, dtype=np.float64),
            np.asarray(metallicity_solar, dtype=np.float64),
        )
        if np.any(pressure <= 0.0) or np.any(density <= 0.0) or np.any(metallicity <= 0.0):
            raise ValueError("pressure, density, and metallicity must be positive")
        n_h = density * hydrogen_mass_fraction / PROTON_MASS_G
        temperature = pressure * 1.2195121951219512 * PROTON_MASS_G / (density * BOLTZMANN_ERG_K)
        for _ in range(iterations):
            mu = self.mean_mu(temperature, n_h, metallicity)
            updated = pressure * mu * PROTON_MASS_G / (density * BOLTZMANN_ERG_K)
            temperature = np.sqrt(temperature * updated)
        return temperature


def read_grackle_equilibrium_table(path: str | Path) -> GrackleEquilibriumTable:
    """Read the binary layout written by Cooling_Grackle/grackle_cooling_grid.c."""

    path = Path(path)
    with path.open("rb") as handle:
        dimensions = np.fromfile(handle, dtype=np.int32, count=3)
        metadata = np.fromfile(handle, dtype=np.float64, count=7)
        if len(dimensions) != 3 or len(metadata) != 7 or np.any(dimensions < 2):
            raise ValueError("invalid Grackle equilibrium-grid header")
        n_density, n_metallicity, n_temperature = (int(value) for value in dimensions)
        log_temperature = np.fromfile(handle, dtype=np.float64, count=n_temperature)
        count = n_density * n_metallicity * n_temperature
        net_rate = np.fromfile(handle, dtype=np.float64, count=count)
        mean_mu = np.fromfile(handle, dtype=np.float64, count=count)
        if len(log_temperature) != n_temperature or len(net_rate) != count or len(mean_mu) != count:
            raise ValueError("truncated Grackle equilibrium-grid file")
        if handle.read(1):
            raise ValueError("unexpected trailing bytes in Grackle equilibrium-grid file")
    redshift, log_n_min, log_n_max, log_z_min, log_z_max, _, _ = metadata
    return GrackleEquilibriumTable(
        redshift=float(redshift),
        log_hydrogen_number_density_cm3=np.linspace(log_n_min, log_n_max, n_density),
        log_metallicity_solar=np.linspace(log_z_min, log_z_max, n_metallicity),
        log_temperature_k=log_temperature,
        net_rate_erg_s_cm3=net_rate.reshape((n_density, n_metallicity, n_temperature)),
        mean_molecular_weight=mean_mu.reshape((n_density, n_metallicity, n_temperature)),
    )
