"""Primordial mean-molecular-weight tables for hydro-to-RT staging."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


HYDROGEN_MASS_FRACTION = 0.76
PROTON_MASS_G = 1.67262192369e-24
BOLTZMANN_ERG_K = 1.380649e-16


def primordial_mean_molecular_weight(
    x_hii: np.ndarray | float,
    x_heii: np.ndarray | float,
    x_heiii: np.ndarray | float,
    *,
    hydrogen_mass_fraction: float = HYDROGEN_MASS_FRACTION,
) -> np.ndarray:
    """Return mu for H/He gas from nuclei-number ionization fractions.

    ``1/mu = X(1+x_HII) + Y/4(1+x_HeII+2*x_HeIII)``.  Electrons are included
    in the particle count; metals are intentionally excluded because the
    current hydro output has no metallicity scalar.
    """

    if not 0.0 < hydrogen_mass_fraction < 1.0:
        raise ValueError("hydrogen_mass_fraction must be between zero and one")
    x_hii, x_heii, x_heiii = np.broadcast_arrays(
        np.asarray(x_hii, dtype=np.float64),
        np.asarray(x_heii, dtype=np.float64),
        np.asarray(x_heiii, dtype=np.float64),
    )
    if np.any((x_hii < 0.0) | (x_hii > 1.0)):
        raise ValueError("x_hii must be in [0, 1]")
    if np.any((x_heii < 0.0) | (x_heii > 1.0) | (x_heiii < 0.0) | (x_heiii > 1.0)):
        raise ValueError("helium ionization fractions must be in [0, 1]")
    if np.any(x_heii + x_heiii > 1.0):
        raise ValueError("x_heii + x_heiii must not exceed one")
    helium_mass_fraction = 1.0 - hydrogen_mass_fraction
    inverse_mu = hydrogen_mass_fraction * (1.0 + x_hii)
    inverse_mu += helium_mass_fraction * (1.0 + x_heii + 2.0 * x_heiii) / 4.0
    return 1.0 / inverse_mu


@dataclass(frozen=True)
class PrimordialMuTable:
    """Bilinearly interpolated primordial thermal and ionization table.

    Table axes are log10(T/K) and log10(n_H/cm^-3).  Values outside the table
    use edge clamping rather than extrapolating unvalidated chemistry.
    """

    log_temperature_k: np.ndarray
    log_hydrogen_number_density_cm3: np.ndarray
    mean_molecular_weight: np.ndarray
    x_hii: np.ndarray
    x_heii: np.ndarray
    x_heiii: np.ndarray

    def __post_init__(self) -> None:
        log_temperature = np.asarray(self.log_temperature_k, dtype=np.float64)
        log_density = np.asarray(self.log_hydrogen_number_density_cm3, dtype=np.float64)
        if log_temperature.ndim != 1 or log_density.ndim != 1:
            raise ValueError("mu table axes must be one-dimensional")
        if len(log_temperature) < 2 or len(log_density) < 2:
            raise ValueError("mu table requires at least two points per axis")
        if np.any(np.diff(log_temperature) <= 0.0) or np.any(np.diff(log_density) <= 0.0):
            raise ValueError("mu table axes must be strictly increasing")
        expected_shape = (len(log_temperature), len(log_density))
        values = {
            "mean_molecular_weight": self.mean_molecular_weight,
            "x_hii": self.x_hii,
            "x_heii": self.x_heii,
            "x_heiii": self.x_heiii,
        }
        normalized = {name: np.asarray(value, dtype=np.float64) for name, value in values.items()}
        if any(value.shape != expected_shape for value in normalized.values()):
            raise ValueError("mu table values must have shape (n_temperature, n_density)")
        if any(not np.isfinite(value).all() for value in normalized.values()):
            raise ValueError("mu table values must be finite")
        if np.any(normalized["mean_molecular_weight"] <= 0.0):
            raise ValueError("mean molecular weight must be positive")
        primordial_mean_molecular_weight(
            normalized["x_hii"], normalized["x_heii"], normalized["x_heiii"]
        )
        object.__setattr__(self, "log_temperature_k", log_temperature)
        object.__setattr__(self, "log_hydrogen_number_density_cm3", log_density)
        for name, value in normalized.items():
            object.__setattr__(self, name, value)

    def _interpolate(self, values: np.ndarray, temperature_k: np.ndarray, n_h_cm3: np.ndarray) -> np.ndarray:
        temperature, density = np.broadcast_arrays(
            np.asarray(temperature_k, dtype=np.float64), np.asarray(n_h_cm3, dtype=np.float64)
        )
        if np.any(temperature <= 0.0) or np.any(density <= 0.0):
            raise ValueError("temperature and hydrogen number density must be positive")
        log_temperature = np.clip(
            np.log10(temperature), self.log_temperature_k[0], self.log_temperature_k[-1]
        )
        log_density = np.clip(
            np.log10(density),
            self.log_hydrogen_number_density_cm3[0],
            self.log_hydrogen_number_density_cm3[-1],
        )
        temperature_index = np.clip(
            np.searchsorted(self.log_temperature_k, log_temperature, side="right") - 1,
            0,
            len(self.log_temperature_k) - 2,
        )
        density_index = np.clip(
            np.searchsorted(self.log_hydrogen_number_density_cm3, log_density, side="right") - 1,
            0,
            len(self.log_hydrogen_number_density_cm3) - 2,
        )
        temperature_weight = (
            (log_temperature - self.log_temperature_k[temperature_index])
            / (self.log_temperature_k[temperature_index + 1] - self.log_temperature_k[temperature_index])
        )
        density_weight = (
            (log_density - self.log_hydrogen_number_density_cm3[density_index])
            / (
                self.log_hydrogen_number_density_cm3[density_index + 1]
                - self.log_hydrogen_number_density_cm3[density_index]
            )
        )
        lower = (1.0 - density_weight) * values[temperature_index, density_index]
        lower += density_weight * values[temperature_index, density_index + 1]
        upper = (1.0 - density_weight) * values[temperature_index + 1, density_index]
        upper += density_weight * values[temperature_index + 1, density_index + 1]
        return (1.0 - temperature_weight) * lower + temperature_weight * upper

    def state(self, temperature_k: np.ndarray, n_h_cm3: np.ndarray) -> tuple[np.ndarray, ...]:
        """Return interpolated (mu, x_HII, x_HeII, x_HeIII)."""

        return (
            self._interpolate(self.mean_molecular_weight, temperature_k, n_h_cm3),
            self._interpolate(self.x_hii, temperature_k, n_h_cm3),
            self._interpolate(self.x_heii, temperature_k, n_h_cm3),
            self._interpolate(self.x_heiii, temperature_k, n_h_cm3),
        )

    def temperature_from_pressure(
        self,
        pressure_dyn_cm2: np.ndarray,
        mass_density_g_cm3: np.ndarray,
        n_h_cm3: np.ndarray,
        *,
        iterations: int = 32,
    ) -> np.ndarray:
        """Solve T = P*mu(T,n_H)*m_p/(rho*k_B) by damped fixed-point iteration."""

        if iterations < 1:
            raise ValueError("iterations must be positive")
        pressure, density, n_h = np.broadcast_arrays(
            np.asarray(pressure_dyn_cm2, dtype=np.float64),
            np.asarray(mass_density_g_cm3, dtype=np.float64),
            np.asarray(n_h_cm3, dtype=np.float64),
        )
        if np.any(pressure <= 0.0) or np.any(density <= 0.0) or np.any(n_h <= 0.0):
            raise ValueError("pressure, density, and n_H must be positive")
        neutral_mu = float(primordial_mean_molecular_weight(0.0, 0.0, 0.0))
        temperature = pressure * neutral_mu * PROTON_MASS_G / (density * BOLTZMANN_ERG_K)
        for _ in range(iterations):
            mu = self._interpolate(self.mean_molecular_weight, temperature, n_h)
            updated = pressure * mu * PROTON_MASS_G / (density * BOLTZMANN_ERG_K)
            temperature = np.sqrt(temperature * updated)
        return temperature


def neutral_primordial_mu_table(
    log_temperature_k: np.ndarray,
    log_hydrogen_number_density_cm3: np.ndarray,
) -> PrimordialMuTable:
    """Create the valid neutral-H/He baseline table, not a fitted ionized gas model."""

    shape = (len(log_temperature_k), len(log_hydrogen_number_density_cm3))
    zeros = np.zeros(shape, dtype=np.float64)
    return PrimordialMuTable(
        log_temperature_k=np.asarray(log_temperature_k),
        log_hydrogen_number_density_cm3=np.asarray(log_hydrogen_number_density_cm3),
        mean_molecular_weight=np.full(shape, primordial_mean_molecular_weight(0.0, 0.0, 0.0)),
        x_hii=zeros,
        x_heii=zeros,
        x_heiii=zeros,
    )


def write_primordial_mu_table(path: str | Path, table: PrimordialMuTable) -> None:
    """Store a portable table without introducing a second binary convention."""

    np.savez_compressed(
        Path(path),
        log_temperature_k=table.log_temperature_k,
        log_hydrogen_number_density_cm3=table.log_hydrogen_number_density_cm3,
        mean_molecular_weight=table.mean_molecular_weight,
        x_hii=table.x_hii,
        x_heii=table.x_heii,
        x_heiii=table.x_heiii,
    )


def read_primordial_mu_table(path: str | Path) -> PrimordialMuTable:
    """Read a table generated by :func:`write_primordial_mu_table`."""

    with np.load(Path(path)) as data:
        return PrimordialMuTable(
            log_temperature_k=data["log_temperature_k"],
            log_hydrogen_number_density_cm3=data["log_hydrogen_number_density_cm3"],
            mean_molecular_weight=data["mean_molecular_weight"],
            x_hii=data["x_hii"],
            x_heii=data["x_heii"],
            x_heiii=data["x_heiii"],
        )
