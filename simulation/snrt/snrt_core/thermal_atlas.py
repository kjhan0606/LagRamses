"""Offline Grackle thermal atlas and runtime interpolation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from snrt_core.grackle import BOLTZMANN_ERG_K, HYDROGEN_MASS_FRACTION, PROTON_MASS_G, GrackleEquilibriumTable


@dataclass(frozen=True)
class ThermalAtlas:
    """Background thermal quantities on (a, log n_H, log Z, log T).

    This is an offline product. Runtime users interpolate it but never invoke
    Grackle or regenerate chemistry tables inside a hydro/RT timestep.
    """

    scale_factor: np.ndarray
    log_hydrogen_number_density_cm3: np.ndarray
    log_metallicity_solar: np.ndarray
    log_temperature_k: np.ndarray
    net_rate_erg_s_cm3: np.ndarray
    mean_molecular_weight: np.ndarray
    equilibrium_log_temperature_k: np.ndarray

    def __post_init__(self) -> None:
        axes = {
            "scale_factor": self.scale_factor,
            "log_hydrogen_number_density_cm3": self.log_hydrogen_number_density_cm3,
            "log_metallicity_solar": self.log_metallicity_solar,
            "log_temperature_k": self.log_temperature_k,
        }
        normalized = {name: np.asarray(value, dtype=np.float64) for name, value in axes.items()}
        if any(axis.ndim != 1 or len(axis) < 2 or np.any(np.diff(axis) <= 0.0) for axis in normalized.values()):
            raise ValueError("thermal atlas axes must be strictly increasing one-dimensional arrays")
        expected_shape = tuple(len(axis) for axis in normalized.values())
        rates = np.asarray(self.net_rate_erg_s_cm3, dtype=np.float64)
        mu = np.asarray(self.mean_molecular_weight, dtype=np.float64)
        equilibrium_temperature = np.asarray(self.equilibrium_log_temperature_k, dtype=np.float64)
        if rates.shape != expected_shape or mu.shape != expected_shape:
            raise ValueError("thermal atlas data must have shape (n_a, n_nH, n_Z, n_T)")
        if equilibrium_temperature.shape != expected_shape[:-1]:
            raise ValueError("equilibrium temperatures must have shape (n_a, n_nH, n_Z)")
        if (
            not np.isfinite(rates).all()
            or not np.isfinite(mu).all()
            or not np.isfinite(equilibrium_temperature).all()
            or np.any(mu <= 0.0)
        ):
            raise ValueError("thermal atlas contains invalid values")
        for name, axis in normalized.items():
            object.__setattr__(self, name, axis)
        object.__setattr__(self, "net_rate_erg_s_cm3", rates)
        object.__setattr__(self, "mean_molecular_weight", mu)
        object.__setattr__(self, "equilibrium_log_temperature_k", equilibrium_temperature)

    @staticmethod
    def _indices(axis: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        clipped = np.clip(values, axis[0], axis[-1])
        index = np.clip(np.searchsorted(axis, clipped, side="right") - 1, 0, len(axis) - 2)
        weight = (clipped - axis[index]) / (axis[index + 1] - axis[index])
        return index, weight

    def _spatial_interpolate(
        self, values: np.ndarray, temperature_k: np.ndarray, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray
    ) -> np.ndarray:
        temperature, density, metallicity = np.broadcast_arrays(
            np.asarray(temperature_k, dtype=np.float64),
            np.asarray(n_h_cm3, dtype=np.float64),
            np.asarray(metallicity_solar, dtype=np.float64),
        )
        if np.any(temperature <= 0.0) or np.any(density <= 0.0) or np.any(metallicity <= 0.0):
            raise ValueError("thermal-atlas interpolation requires positive T, n_H, and Z")
        density_index, density_weight = self._indices(self.log_hydrogen_number_density_cm3, np.log10(density))
        metallicity_index, metallicity_weight = self._indices(self.log_metallicity_solar, np.log10(metallicity))
        temperature_index, temperature_weight = self._indices(self.log_temperature_k, np.log10(temperature))
        result = np.zeros_like(temperature, dtype=np.float64)
        for density_offset in (0, 1):
            for metallicity_offset in (0, 1):
                for temperature_offset in (0, 1):
                    result += (
                        (density_weight if density_offset else 1.0 - density_weight)
                        * (metallicity_weight if metallicity_offset else 1.0 - metallicity_weight)
                        * (temperature_weight if temperature_offset else 1.0 - temperature_weight)
                        * values[
                            density_index + density_offset,
                            metallicity_index + metallicity_offset,
                            temperature_index + temperature_offset,
                        ]
                    )
        return result

    def _interpolate(
        self,
        values: np.ndarray,
        scale_factor: np.ndarray | float,
        temperature_k: np.ndarray,
        n_h_cm3: np.ndarray,
        metallicity_solar: np.ndarray | float,
    ) -> np.ndarray:
        time = np.asarray(scale_factor, dtype=np.float64)
        time_index, time_weight = self._indices(self.scale_factor, time)
        lower = self._spatial_interpolate(values[time_index], temperature_k, n_h_cm3, metallicity_solar)
        upper = self._spatial_interpolate(values[time_index + 1], temperature_k, n_h_cm3, metallicity_solar)
        return (1.0 - time_weight) * lower + time_weight * upper

    def mean_mu(
        self, scale_factor: np.ndarray | float, temperature_k: np.ndarray, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray | float
    ) -> np.ndarray:
        return self._interpolate(self.mean_molecular_weight, scale_factor, temperature_k, n_h_cm3, metallicity_solar)

    def net_rate(
        self, scale_factor: np.ndarray | float, temperature_k: np.ndarray, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray | float
    ) -> np.ndarray:
        """Return signed background rate. Rates interpolate linearly, not logarithmically."""

        return self._interpolate(self.net_rate_erg_s_cm3, scale_factor, temperature_k, n_h_cm3, metallicity_solar)

    def equilibrium_temperature(
        self, scale_factor: np.ndarray | float, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray | float
    ) -> np.ndarray:
        """Interpolate the precomputed net-rate-zero temperature.

        This is a fallback for snapshots that do not store usable pressure or
        temperature. It is not a substitute for the evolved hydro temperature.
        """

        time = np.asarray(scale_factor, dtype=np.float64)
        time_index, time_weight = self._indices(self.scale_factor, time)
        dummy_temperature = np.full(np.broadcast(n_h_cm3, metallicity_solar).shape, 10.0)
        lower = self._spatial_interpolate_2d(
            self.equilibrium_log_temperature_k[time_index], dummy_temperature, n_h_cm3, metallicity_solar
        )
        upper = self._spatial_interpolate_2d(
            self.equilibrium_log_temperature_k[time_index + 1], dummy_temperature, n_h_cm3, metallicity_solar
        )
        return 10.0 ** ((1.0 - time_weight) * lower + time_weight * upper)

    def _spatial_interpolate_2d(
        self, values: np.ndarray, dummy_temperature: np.ndarray, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray
    ) -> np.ndarray:
        _, density, metallicity = np.broadcast_arrays(dummy_temperature, np.asarray(n_h_cm3), np.asarray(metallicity_solar))
        if np.any(density <= 0.0) or np.any(metallicity <= 0.0):
            raise ValueError("equilibrium-temperature lookup requires positive n_H and metallicity")
        density_index, density_weight = self._indices(self.log_hydrogen_number_density_cm3, np.log10(density))
        metallicity_index, metallicity_weight = self._indices(self.log_metallicity_solar, np.log10(metallicity))
        return (
            (1.0 - density_weight) * (1.0 - metallicity_weight) * values[density_index, metallicity_index]
            + density_weight * (1.0 - metallicity_weight) * values[density_index + 1, metallicity_index]
            + (1.0 - density_weight) * metallicity_weight * values[density_index, metallicity_index + 1]
            + density_weight * metallicity_weight * values[density_index + 1, metallicity_index + 1]
        )

    def temperature_from_pressure(
        self,
        scale_factor: float,
        pressure_dyn_cm2: np.ndarray,
        mass_density_g_cm3: np.ndarray,
        metallicity_solar: np.ndarray | float,
        *,
        iterations: int = 32,
        hydrogen_mass_fraction: float = HYDROGEN_MASS_FRACTION,
    ) -> np.ndarray:
        """Invert P/rho using interpolated mu(a, n_H, Z, T)."""

        if iterations < 1 or not 0.0 < hydrogen_mass_fraction < 1.0:
            raise ValueError("invalid thermal-atlas temperature-inversion configuration")
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
            mu = self.mean_mu(scale_factor, temperature, n_h, metallicity)
            temperature = np.sqrt(temperature * pressure * mu * PROTON_MASS_G / (density * BOLTZMANN_ERG_K))
        return temperature


def thermal_atlas_from_grackle(scale_factor: np.ndarray, tables: list[GrackleEquilibriumTable]) -> ThermalAtlas:
    """Combine identically gridded offline Grackle subtables into a runtime atlas."""

    scale_factor = np.asarray(scale_factor, dtype=np.float64)
    if len(tables) != len(scale_factor) or len(tables) < 2:
        raise ValueError("at least two matched Grackle subtables are required")
    first = tables[0]
    for table in tables[1:]:
        if not (
            np.array_equal(table.log_hydrogen_number_density_cm3, first.log_hydrogen_number_density_cm3)
            and np.array_equal(table.log_metallicity_solar, first.log_metallicity_solar)
            and np.array_equal(table.log_temperature_k, first.log_temperature_k)
        ):
            raise ValueError("all Grackle subtables must share identical spatial axes")
    net_rate = np.stack([table.net_rate_erg_s_cm3 for table in tables])
    equilibrium_log_temperature = _equilibrium_log_temperature(first.log_temperature_k, net_rate)
    return ThermalAtlas(
        scale_factor=scale_factor,
        log_hydrogen_number_density_cm3=first.log_hydrogen_number_density_cm3,
        log_metallicity_solar=first.log_metallicity_solar,
        log_temperature_k=first.log_temperature_k,
        net_rate_erg_s_cm3=net_rate,
        mean_molecular_weight=np.stack([table.mean_molecular_weight for table in tables]),
        equilibrium_log_temperature_k=equilibrium_log_temperature,
    )


def _equilibrium_log_temperature(log_temperature: np.ndarray, net_rate: np.ndarray) -> np.ndarray:
    """Choose the first stable net-rate-zero crossing in every table column."""

    result = np.empty(net_rate.shape[:-1], dtype=np.float64)
    for index in np.ndindex(result.shape):
        profile = net_rate[index]
        crossings = np.flatnonzero(profile[:-1] * profile[1:] <= 0.0)
        if len(crossings):
            slopes = (profile[crossings + 1] - profile[crossings]) / (log_temperature[crossings + 1] - log_temperature[crossings])
            stable = crossings[slopes < 0.0]
            crossing = int(stable[0] if len(stable) else crossings[0])
            left, right = profile[crossing], profile[crossing + 1]
            fraction = 0.0 if right == left else -left / (right - left)
            result[index] = log_temperature[crossing] + fraction * (log_temperature[crossing + 1] - log_temperature[crossing])
        else:
            result[index] = log_temperature[int(np.argmin(np.abs(profile)))]
    return result


def write_thermal_atlas(path: str | Path, atlas: ThermalAtlas) -> None:
    import h5py

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as handle:
        handle.attrs["format"] = "snrt_thermal_atlas"
        handle.attrs["format_version"] = 1
        handle.create_dataset("scale_factor", data=atlas.scale_factor)
        handle.create_dataset("log_hydrogen_number_density_cm3", data=atlas.log_hydrogen_number_density_cm3)
        handle.create_dataset("log_metallicity_solar", data=atlas.log_metallicity_solar)
        handle.create_dataset("log_temperature_k", data=atlas.log_temperature_k)
        handle.create_dataset("net_rate_erg_s_cm3", data=atlas.net_rate_erg_s_cm3, compression="gzip")
        handle.create_dataset("mean_molecular_weight", data=atlas.mean_molecular_weight, compression="gzip")
        handle.create_dataset("equilibrium_log_temperature_k", data=atlas.equilibrium_log_temperature_k, compression="gzip")


def read_thermal_atlas(path: str | Path) -> ThermalAtlas:
    import h5py

    with h5py.File(Path(path), "r") as handle:
        if handle.attrs.get("format", "") != "snrt_thermal_atlas" or int(handle.attrs.get("format_version", -1)) != 1:
            raise ValueError("not an SNRT thermal-atlas file")
        return ThermalAtlas(
            scale_factor=np.asarray(handle["scale_factor"]),
            log_hydrogen_number_density_cm3=np.asarray(handle["log_hydrogen_number_density_cm3"]),
            log_metallicity_solar=np.asarray(handle["log_metallicity_solar"]),
            log_temperature_k=np.asarray(handle["log_temperature_k"]),
            net_rate_erg_s_cm3=np.asarray(handle["net_rate_erg_s_cm3"]),
            mean_molecular_weight=np.asarray(handle["mean_molecular_weight"]),
            equilibrium_log_temperature_k=np.asarray(handle["equilibrium_log_temperature_k"]),
        )
