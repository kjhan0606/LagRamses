"""Offline Grackle thermal atlas and runtime interpolation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from types import MappingProxyType
from typing import Mapping

import numpy as np

from snrt_core.grackle import BOLTZMANN_ERG_K, HYDROGEN_MASS_FRACTION, PROTON_MASS_G


THERMAL_ATLAS_FORMAT = "snrt_thermal_atlas"
THERMAL_ATLAS_FORMAT_VERSION = 3
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_REQUIRED_PROVENANCE_VALUES = {
    "thermal_component": "metal_only",
    "primordial_rates_included": "false",
    "uv_background_included": "false",
    "photoheating_included": "false",
    "metallicity_scaling": "linear_z_solar",
    "metallicity_application": "analytic_runtime_multiplier",
    "rate_sign_convention": "heating_positive_cooling_negative",
}
_REQUIRED_PROVENANCE_KEYS = (
    "source_data_name",
    "source_data_sha256",
    "source_cooling_dataset",
    "source_repository",
    "source_repository_revision",
    "generator_name",
    "generator_version",
    "generator_sha256",
    "grackle_version",
    "grackle_revision",
    "generated_utc",
    "cmb_metal_floor",
)


def validate_thermal_atlas_provenance(provenance: Mapping[str, str]) -> dict[str, str]:
    """Normalize and validate the production metal-atlas provenance contract."""

    normalized = {str(key): str(value) for key, value in provenance.items()}
    required = set(_REQUIRED_PROVENANCE_VALUES) | set(_REQUIRED_PROVENANCE_KEYS)
    missing = sorted(required - set(normalized))
    if missing:
        raise ValueError(f"thermal atlas lacks required provenance: {', '.join(missing)}")
    for key, expected in _REQUIRED_PROVENANCE_VALUES.items():
        if normalized[key] != expected:
            raise ValueError(f"thermal atlas provenance {key!r} must be {expected!r}")
    for key in ("source_data_sha256", "generator_sha256"):
        if not _SHA256_PATTERN.fullmatch(normalized[key]):
            raise ValueError(f"thermal atlas provenance {key!r} is not a SHA-256")
    if not re.fullmatch(r"[0-9a-f]{40}", normalized["grackle_revision"]):
        raise ValueError("thermal atlas provenance grackle_revision is not a full git revision")
    if normalized["cmb_metal_floor"] != "continuous_subtraction_at_tcmb":
        raise ValueError("thermal atlas provenance cmb_metal_floor must be continuous_subtraction_at_tcmb")
    if not all(normalized[key].strip() for key in _REQUIRED_PROVENANCE_KEYS):
        raise ValueError("thermal atlas provenance contains an empty required value")
    return normalized


@dataclass(frozen=True)
class ThermalAtlas:
    """Solar-metallicity thermal quantities on (a, log n_H, log T).

    This is an offline product. Runtime users interpolate it but never invoke
    Grackle or regenerate chemistry tables inside a hydro/RT timestep. Atomic
    primordial cooling is deliberately excluded and evaluated from the live
    non-equilibrium H/He state instead.
    """

    scale_factor: np.ndarray
    log_hydrogen_number_density_cm3: np.ndarray
    log_temperature_k: np.ndarray
    net_rate_erg_s_cm3: np.ndarray
    mean_molecular_weight: np.ndarray
    equilibrium_log_temperature_k: np.ndarray
    provenance: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        axes = {
            "scale_factor": self.scale_factor,
            "log_hydrogen_number_density_cm3": self.log_hydrogen_number_density_cm3,
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
            raise ValueError("thermal atlas data must have shape (n_a, n_nH, n_T)")
        if equilibrium_temperature.shape != expected_shape[:-1]:
            raise ValueError("equilibrium temperatures must have shape (n_a, n_nH)")
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
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType({str(key): str(value) for key, value in self.provenance.items()}),
        )

    @staticmethod
    def _indices(axis: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        clipped = np.clip(values, axis[0], axis[-1])
        index = np.clip(np.searchsorted(axis, clipped, side="right") - 1, 0, len(axis) - 2)
        weight = (clipped - axis[index]) / (axis[index + 1] - axis[index])
        return index, weight

    def validate_runtime_domain(
        self,
        scale_factor: float,
        temperature_k: np.ndarray,
        n_h_cm3: np.ndarray,
        *,
        label: str = "thermal runtime",
    ) -> None:
        """Reject a runtime state that would otherwise be edge-clamped.

        The interpolation helpers retain clamping for controlled offline
        studies.  A production-facing caller must prove that its epoch and
        initial gas state lie inside the admitted table domain instead of
        silently using an edge value.
        """

        if not np.isfinite(scale_factor) or not self.scale_factor[0] <= scale_factor <= self.scale_factor[-1]:
            raise ValueError(
                f"{label} scale factor {scale_factor:.17g} is outside the admitted "
                f"range [{self.scale_factor[0]:.17g}, {self.scale_factor[-1]:.17g}]"
            )
        temperature, density = np.broadcast_arrays(
            np.asarray(temperature_k, dtype=np.float64),
            np.asarray(n_h_cm3, dtype=np.float64),
        )
        if not np.isfinite(temperature).all() or np.any(temperature <= 0.0):
            raise ValueError(f"{label} temperature must be finite and positive")
        if not np.isfinite(density).all() or np.any(density <= 0.0):
            raise ValueError(f"{label} n_H must be finite and positive")
        log_density = np.log10(density)
        log_temperature = np.log10(temperature)
        if np.any(
            (log_density < self.log_hydrogen_number_density_cm3[0])
            | (log_density > self.log_hydrogen_number_density_cm3[-1])
        ):
            raise ValueError(
                f"{label} n_H lies outside the admitted range "
                f"[{10.0 ** self.log_hydrogen_number_density_cm3[0]:.17g}, "
                f"{10.0 ** self.log_hydrogen_number_density_cm3[-1]:.17g}] cm^-3"
            )
        if np.any(
            (log_temperature < self.log_temperature_k[0])
            | (log_temperature > self.log_temperature_k[-1])
        ):
            raise ValueError(
                f"{label} temperature lies outside the admitted range "
                f"[{10.0 ** self.log_temperature_k[0]:.17g}, "
                f"{10.0 ** self.log_temperature_k[-1]:.17g}] K"
            )

    def _spatial_interpolate(
        self, values: np.ndarray, temperature_k: np.ndarray, n_h_cm3: np.ndarray
    ) -> np.ndarray:
        temperature, density = np.broadcast_arrays(
            np.asarray(temperature_k, dtype=np.float64),
            np.asarray(n_h_cm3, dtype=np.float64),
        )
        if np.any(temperature <= 0.0) or np.any(density <= 0.0):
            raise ValueError("thermal-atlas interpolation requires positive T and n_H")
        density_index, density_weight = self._indices(self.log_hydrogen_number_density_cm3, np.log10(density))
        temperature_index, temperature_weight = self._indices(self.log_temperature_k, np.log10(temperature))
        result = np.zeros_like(temperature, dtype=np.float64)
        for density_offset in (0, 1):
            for temperature_offset in (0, 1):
                result += (
                    (density_weight if density_offset else 1.0 - density_weight)
                    * (temperature_weight if temperature_offset else 1.0 - temperature_weight)
                    * values[
                        density_index + density_offset,
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
    ) -> np.ndarray:
        time = np.asarray(scale_factor, dtype=np.float64)
        time_index, time_weight = self._indices(self.scale_factor, time)
        lower = self._spatial_interpolate(values[time_index], temperature_k, n_h_cm3)
        upper = self._spatial_interpolate(values[time_index + 1], temperature_k, n_h_cm3)
        return (1.0 - time_weight) * lower + time_weight * upper

    def mean_mu(
        self, scale_factor: np.ndarray | float, temperature_k: np.ndarray, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray | float
    ) -> np.ndarray:
        temperature, density, metallicity = np.broadcast_arrays(
            np.asarray(temperature_k, dtype=np.float64),
            np.asarray(n_h_cm3, dtype=np.float64),
            np.asarray(metallicity_solar, dtype=np.float64),
        )
        if not np.isfinite(metallicity).all() or np.any(metallicity < 0.0):
            raise ValueError("metallicity must be finite and non-negative")
        return self._interpolate(self.mean_molecular_weight, scale_factor, temperature, density)

    def net_rate(
        self, scale_factor: np.ndarray | float, temperature_k: np.ndarray, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray | float
    ) -> np.ndarray:
        """Return the signed solar table multiplied analytically by linear Z."""

        temperature, density, metallicity = np.broadcast_arrays(
            np.asarray(temperature_k, dtype=np.float64),
            np.asarray(n_h_cm3, dtype=np.float64),
            np.asarray(metallicity_solar, dtype=np.float64),
        )
        if not np.isfinite(metallicity).all() or np.any(metallicity < 0.0):
            raise ValueError("metallicity must be finite and non-negative")
        solar_rate = self._interpolate(self.net_rate_erg_s_cm3, scale_factor, temperature, density)
        return solar_rate * metallicity

    def equilibrium_temperature(
        self, scale_factor: np.ndarray | float, n_h_cm3: np.ndarray, metallicity_solar: np.ndarray | float
    ) -> np.ndarray:
        """Interpolate the precomputed net-rate-zero temperature.

        This is a fallback for snapshots that do not store usable pressure or
        temperature. It is not a substitute for the evolved hydro temperature.
        """

        time = np.asarray(scale_factor, dtype=np.float64)
        time_index, time_weight = self._indices(self.scale_factor, time)
        density, metallicity = np.broadcast_arrays(
            np.asarray(n_h_cm3, dtype=np.float64),
            np.asarray(metallicity_solar, dtype=np.float64),
        )
        if not np.isfinite(metallicity).all() or np.any(metallicity <= 0.0):
            raise ValueError("metal-only equilibrium lookup requires positive metallicity")
        lower = self._density_interpolate(self.equilibrium_log_temperature_k[time_index], density)
        upper = self._density_interpolate(self.equilibrium_log_temperature_k[time_index + 1], density)
        return 10.0 ** ((1.0 - time_weight) * lower + time_weight * upper)

    def _density_interpolate(self, values: np.ndarray, n_h_cm3: np.ndarray) -> np.ndarray:
        density = np.asarray(n_h_cm3, dtype=np.float64)
        if np.any(density <= 0.0):
            raise ValueError("equilibrium-temperature lookup requires positive n_H")
        density_index, density_weight = self._indices(self.log_hydrogen_number_density_cm3, np.log10(density))
        return (1.0 - density_weight) * values[density_index] + density_weight * values[density_index + 1]

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
        """Invert P/rho using the atlas mu(a, n_H, T) staging aid."""

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


def thermal_atlas_from_grackle(*_args, **_kwargs) -> ThermalAtlas:
    """Reject the retired equilibrium H/He+metal atlas construction path."""

    raise RuntimeError(
        "equilibrium Grackle atlases cannot satisfy thermal-atlas format v3; "
        "use tools/build_metal_thermal_atlas.py"
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
    provenance = validate_thermal_atlas_provenance(atlas.provenance)
    with h5py.File(path, "x") as handle:
        handle.attrs["format"] = THERMAL_ATLAS_FORMAT
        handle.attrs["format_version"] = THERMAL_ATLAS_FORMAT_VERSION
        provenance_group = handle.create_group("provenance")
        for key, value in sorted(provenance.items()):
            provenance_group.attrs[key] = value
        handle.create_dataset("scale_factor", data=atlas.scale_factor)
        handle.create_dataset("log_hydrogen_number_density_cm3", data=atlas.log_hydrogen_number_density_cm3)
        handle.create_dataset("log_temperature_k", data=atlas.log_temperature_k)
        handle.create_dataset("net_rate_erg_s_cm3", data=atlas.net_rate_erg_s_cm3, compression="gzip")
        handle.create_dataset("mean_molecular_weight", data=atlas.mean_molecular_weight, compression="gzip")
        handle.create_dataset("equilibrium_log_temperature_k", data=atlas.equilibrium_log_temperature_k, compression="gzip")


def read_thermal_atlas(path: str | Path) -> ThermalAtlas:
    import h5py

    with h5py.File(Path(path), "r") as handle:
        if (
            handle.attrs.get("format", "") != THERMAL_ATLAS_FORMAT
            or int(handle.attrs.get("format_version", -1)) != THERMAL_ATLAS_FORMAT_VERSION
        ):
            raise ValueError("thermal atlas is not the provenance-enforced SNRT format v3")
        if "provenance" not in handle:
            raise ValueError("thermal atlas has no provenance group")
        provenance = validate_thermal_atlas_provenance(
            {
                key: value.decode() if isinstance(value, bytes) else str(value)
                for key, value in handle["provenance"].attrs.items()
            }
        )
        return ThermalAtlas(
            scale_factor=np.asarray(handle["scale_factor"]),
            log_hydrogen_number_density_cm3=np.asarray(handle["log_hydrogen_number_density_cm3"]),
            log_temperature_k=np.asarray(handle["log_temperature_k"]),
            net_rate_erg_s_cm3=np.asarray(handle["net_rate_erg_s_cm3"]),
            mean_molecular_weight=np.asarray(handle["mean_molecular_weight"]),
            equilibrium_log_temperature_k=np.asarray(handle["equilibrium_log_temperature_k"]),
            provenance=provenance,
        )
