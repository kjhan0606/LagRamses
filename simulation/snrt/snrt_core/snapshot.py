"""Static-grid input contract for S_N RT runs.

All arrays in this module are host-side NumPy arrays.  The resulting
``StaticRTInput`` is the explicit boundary between a RAMSES snapshot reader
and the JAX transport/chemistry kernels.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import re
import shutil
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np

from snrt_core.grackle import GrackleEquilibriumTable
from snrt_core.mu_table import PrimordialMuTable
from snrt_core.thermal_atlas import ThermalAtlas


HYDROGEN_MASS_FRACTION = 0.76
HELIUM_MASS_FRACTION = 1.0 - HYDROGEN_MASS_FRACTION
PROTON_MASS_G = 1.67262192369e-24
BOLTZMANN_ERG_K = 1.380649e-16
FORMAT_NAME = "snrt_static_rt_input"
FORMAT_VERSION = 1


@dataclass(frozen=True)
class GridSpec:
    """Uniform Cartesian grid geometry in proper cgs units."""

    cell_width_cm: float
    left_edge_cm: np.ndarray

    def __post_init__(self) -> None:
        edge = np.asarray(self.left_edge_cm, dtype=np.float64)
        if edge.shape != (3,):
            raise ValueError("left_edge_cm must contain exactly three coordinates")
        if not np.isfinite(edge).all():
            raise ValueError("left_edge_cm must be finite")
        if not np.isfinite(self.cell_width_cm) or self.cell_width_cm <= 0.0:
            raise ValueError("cell_width_cm must be positive and finite")
        object.__setattr__(self, "left_edge_cm", edge)


@dataclass(frozen=True)
class SourceCatalog:
    """Photon-number luminosities assigned to cells for each RT group."""

    cell_index: np.ndarray
    photon_luminosity_s: np.ndarray

    def __post_init__(self) -> None:
        cell_index = np.asarray(self.cell_index, dtype=np.int64)
        luminosity = np.asarray(self.photon_luminosity_s, dtype=np.float64)
        if cell_index.ndim != 2 or cell_index.shape[1] != 3:
            raise ValueError("source cell_index must have shape (n_source, 3)")
        if luminosity.ndim != 2 or luminosity.shape[0] != cell_index.shape[0]:
            raise ValueError("source photon_luminosity_s must have shape (n_source, n_group)")
        if luminosity.shape[1] == 0 or not np.isfinite(luminosity).all() or np.any(luminosity < 0.0):
            raise ValueError("source photon luminosities must be finite, non-negative, and grouped")
        object.__setattr__(self, "cell_index", cell_index)
        object.__setattr__(self, "photon_luminosity_s", luminosity)


@dataclass(frozen=True)
class StaticRTInput:
    """Gas state, initial ionization state, and optional discrete sources.

    ``hydrogen_number_density_cm3`` and ``helium_number_density_cm3`` are
    nuclei number densities.  Fractions are number fractions in [0, 1], and
    temperature is in kelvin.  ``dust_relative_abundance`` is dimensionless:
    one means the dust cross-section normalization supplied to the RT run.
    """

    grid: GridSpec
    hydrogen_number_density_cm3: np.ndarray
    helium_number_density_cm3: np.ndarray
    temperature_k: np.ndarray
    dust_relative_abundance: np.ndarray
    x_hii: np.ndarray
    x_heii: np.ndarray
    x_heiii: np.ndarray
    sources: SourceCatalog | None = None

    def __post_init__(self) -> None:
        arrays = {
            "hydrogen_number_density_cm3": self.hydrogen_number_density_cm3,
            "helium_number_density_cm3": self.helium_number_density_cm3,
            "temperature_k": self.temperature_k,
            "dust_relative_abundance": self.dust_relative_abundance,
            "x_hii": self.x_hii,
            "x_heii": self.x_heii,
            "x_heiii": self.x_heiii,
        }
        normalized = {name: np.asarray(value, dtype=np.float64) for name, value in arrays.items()}
        shape = normalized["hydrogen_number_density_cm3"].shape
        if len(shape) != 3 or any(length <= 0 for length in shape):
            raise ValueError("gas fields must be non-empty three-dimensional arrays")
        if any(value.shape != shape for value in normalized.values()):
            raise ValueError("all gas fields must share one Cartesian shape")
        if any(not np.isfinite(value).all() for value in normalized.values()):
            raise ValueError("gas fields must contain finite values")
        if np.any(normalized["hydrogen_number_density_cm3"] <= 0.0):
            raise ValueError("hydrogen_number_density_cm3 must be positive")
        if np.any(normalized["helium_number_density_cm3"] < 0.0):
            raise ValueError("helium_number_density_cm3 must be non-negative")
        if np.any(normalized["temperature_k"] <= 0.0):
            raise ValueError("temperature_k must be positive")
        if np.any(normalized["dust_relative_abundance"] < 0.0):
            raise ValueError("dust_relative_abundance must be non-negative")
        for name in ("x_hii", "x_heii", "x_heiii"):
            if np.any((normalized[name] < 0.0) | (normalized[name] > 1.0)):
                raise ValueError(f"{name} must be a fraction in [0, 1]")
        if np.any(normalized["x_heii"] + normalized["x_heiii"] > 1.0):
            raise ValueError("x_heii + x_heiii must not exceed one")
        if self.sources is not None:
            if np.any(self.sources.cell_index < 0) or np.any(self.sources.cell_index >= np.asarray(shape)):
                raise ValueError("source cell_index lies outside the gas grid")
        for name, value in normalized.items():
            object.__setattr__(self, name, value)

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.hydrogen_number_density_cm3.shape


def neutral_primordial_input(
    grid: GridSpec,
    mass_density_g_cm3: np.ndarray,
    temperature_k: np.ndarray,
    dust_relative_abundance: np.ndarray | float = 0.0,
    sources: SourceCatalog | None = None,
    *,
    hydrogen_mass_fraction: float = HYDROGEN_MASS_FRACTION,
) -> StaticRTInput:
    """Build a neutral H/He input from RAMSES gas density in cgs units."""

    if not 0.0 < hydrogen_mass_fraction < 1.0:
        raise ValueError("hydrogen_mass_fraction must be between zero and one")
    density = np.asarray(mass_density_g_cm3, dtype=np.float64)
    if np.any(density <= 0.0):
        raise ValueError("mass_density_g_cm3 must be positive")
    helium_mass_fraction = 1.0 - hydrogen_mass_fraction
    shape = density.shape
    if len(shape) != 3:
        raise ValueError("mass_density_g_cm3 must be a three-dimensional array")
    dust = np.broadcast_to(np.asarray(dust_relative_abundance, dtype=np.float64), shape).copy()
    zeros = np.zeros(shape, dtype=np.float64)
    return StaticRTInput(
        grid=grid,
        hydrogen_number_density_cm3=density * hydrogen_mass_fraction / PROTON_MASS_G,
        helium_number_density_cm3=density * helium_mass_fraction / (4.0 * PROTON_MASS_G),
        temperature_k=temperature_k,
        dust_relative_abundance=dust,
        x_hii=zeros,
        x_heii=zeros,
        x_heiii=zeros,
        sources=sources,
    )


def write_static_rt_input(path: str | Path, snapshot: StaticRTInput) -> None:
    """Write a self-describing canonical P4 staging file."""

    import h5py

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["format"] = FORMAT_NAME
        handle.attrs["format_version"] = FORMAT_VERSION
        handle.attrs["cell_width_cm"] = snapshot.grid.cell_width_cm
        handle.create_dataset("grid/left_edge_cm", data=snapshot.grid.left_edge_cm)
        gas = handle.create_group("gas")
        gas.create_dataset("hydrogen_number_density_cm3", data=snapshot.hydrogen_number_density_cm3)
        gas.create_dataset("helium_number_density_cm3", data=snapshot.helium_number_density_cm3)
        gas.create_dataset("temperature_k", data=snapshot.temperature_k)
        gas.create_dataset("dust_relative_abundance", data=snapshot.dust_relative_abundance)
        ionization = handle.create_group("ionization")
        ionization.create_dataset("x_hii", data=snapshot.x_hii)
        ionization.create_dataset("x_heii", data=snapshot.x_heii)
        ionization.create_dataset("x_heiii", data=snapshot.x_heiii)
        if snapshot.sources is not None:
            sources = handle.create_group("sources")
            sources.create_dataset("cell_index", data=snapshot.sources.cell_index)
            sources.create_dataset("photon_luminosity_s", data=snapshot.sources.photon_luminosity_s)


def read_static_rt_input(path: str | Path) -> StaticRTInput:
    """Read a canonical P4 staging file and validate its physical contract."""

    import h5py

    with h5py.File(Path(path), "r") as handle:
        if handle.attrs.get("format", "") != FORMAT_NAME:
            raise ValueError("not an SNRT static RT input file")
        if int(handle.attrs.get("format_version", -1)) != FORMAT_VERSION:
            raise ValueError("unsupported SNRT static RT input format version")
        sources = None
        if "sources" in handle:
            sources = SourceCatalog(
                cell_index=np.asarray(handle["sources/cell_index"]),
                photon_luminosity_s=np.asarray(handle["sources/photon_luminosity_s"]),
            )
        return StaticRTInput(
            grid=GridSpec(
                cell_width_cm=float(handle.attrs["cell_width_cm"]),
                left_edge_cm=np.asarray(handle["grid/left_edge_cm"]),
            ),
            hydrogen_number_density_cm3=np.asarray(handle["gas/hydrogen_number_density_cm3"]),
            helium_number_density_cm3=np.asarray(handle["gas/helium_number_density_cm3"]),
            temperature_k=np.asarray(handle["gas/temperature_k"]),
            dust_relative_abundance=np.asarray(handle["gas/dust_relative_abundance"]),
            x_hii=np.asarray(handle["ionization/x_hii"]),
            x_heii=np.asarray(handle["ionization/x_heii"]),
            x_heiii=np.asarray(handle["ionization/x_heiii"]),
            sources=sources,
        )


@dataclass(frozen=True)
class RamsesFieldMap:
    """Explicit yt fields and thermal conversion assumptions.

    Supply either a temperature field, or a thermal-pressure field together
    with its declared mean molecular weight.  No field name, unit, or thermal
    composition is inferred by the adapter.
    """

    density: Any
    temperature: Any | None = None
    thermal_pressure: Any | None = None
    mean_molecular_weight: float | None = None
    dust_relative_abundance: Any | None = None
    equilibrium_temperature: bool = False

    def __post_init__(self) -> None:
        thermal_field_count = int(self.temperature is not None) + int(self.thermal_pressure is not None)
        if thermal_field_count != 1 and not (thermal_field_count == 0 and self.equilibrium_temperature):
            raise ValueError("supply one thermal field, or request equilibrium_temperature fallback")
        if self.mean_molecular_weight is not None and self.mean_molecular_weight <= 0.0:
            raise ValueError("mean_molecular_weight must be positive")


def _stage_yt_dataset(
    dataset: Any,
    output_path: str | Path,
    *,
    level: int,
    dimensions: tuple[int, int, int],
    fields: RamsesFieldMap,
    left_edge_code: tuple[float, float, float],
    right_edge_code: tuple[float, float, float] | None,
    sources: SourceCatalog | None,
    hydrogen_mass_fraction: float,
    mu_table: PrimordialMuTable | None,
    grackle_table: GrackleEquilibriumTable | None,
    grackle_metallicity_solar: float,
    thermal_atlas: ThermalAtlas | None,
    thermal_scale_factor: float | None,
) -> StaticRTInput:
    left_edge = dataset.arr(left_edge_code, "code_length")
    if right_edge_code is None:
        sampled_grid = dataset.covering_grid(level=level, left_edge=left_edge, dims=dimensions)
        cell_width_cm = float(sampled_grid.dds[0].to_value("cm"))
    else:
        right_edge = dataset.arr(right_edge_code, "code_length")
        if np.any(right_edge <= left_edge):
            raise ValueError("right_edge_code must be greater than left_edge_code in every dimension")
        sampled_grid = dataset.arbitrary_grid(left_edge, right_edge, dims=dimensions)
        cell_width_cm = float(((right_edge[0] - left_edge[0]) / dimensions[0]).to_value("cm"))
    density = np.asarray(sampled_grid[fields.density].to_value("g/cm**3"))
    pressure = None
    if fields.thermal_pressure is not None:
        pressure = np.asarray(sampled_grid[fields.thermal_pressure].to_value("dyn/cm**2"))
        valid_hydro = np.isfinite(density) & np.isfinite(pressure) & (density > 0.0) & (pressure > 0.0)
        if not np.any(valid_hydro):
            raise ValueError("selected RAMSES region contains no valid density/pressure cells")
        if thermal_atlas is not None:
            n_h_floor = 10.0 ** thermal_atlas.log_hydrogen_number_density_cm3[0]
            temperature_floor = 10.0 ** thermal_atlas.log_temperature_k[0]
        elif grackle_table is not None:
            n_h_floor = 10.0 ** grackle_table.log_hydrogen_number_density_cm3[0]
            temperature_floor = 10.0 ** grackle_table.log_temperature_k[0]
        else:
            n_h_floor = 1.0e-12
            temperature_floor = 1.0
        density_floor = n_h_floor * PROTON_MASS_G / hydrogen_mass_fraction
        density = np.where(valid_hydro, density, density_floor)
        pressure_floor = density * BOLTZMANN_ERG_K * temperature_floor / (1.3 * PROTON_MASS_G)
        pressure = np.where(valid_hydro, pressure, pressure_floor)
    hydrogen_number_density = density * hydrogen_mass_fraction / PROTON_MASS_G
    if fields.temperature is not None:
        temperature = np.asarray(sampled_grid[fields.temperature].to_value("K"))
    elif fields.equilibrium_temperature:
        if thermal_atlas is None or thermal_scale_factor is None:
            raise ValueError("equilibrium_temperature fallback requires thermal_atlas and thermal_scale_factor")
        temperature = thermal_atlas.equilibrium_temperature(
            thermal_scale_factor, hydrogen_number_density, grackle_metallicity_solar
        )
    else:
        if thermal_atlas is not None:
            if thermal_scale_factor is None:
                raise ValueError("thermal_atlas staging requires thermal_scale_factor")
            temperature = thermal_atlas.temperature_from_pressure(
                thermal_scale_factor,
                pressure,
                density,
                grackle_metallicity_solar,
                hydrogen_mass_fraction=hydrogen_mass_fraction,
            )
        elif grackle_table is not None:
            temperature = grackle_table.temperature_from_pressure(
                pressure, density, grackle_metallicity_solar, hydrogen_mass_fraction=hydrogen_mass_fraction
            )
        elif mu_table is None:
            if fields.mean_molecular_weight is None:
                raise ValueError("thermal_pressure requires mean_molecular_weight or mu_table")
            temperature = pressure * fields.mean_molecular_weight * PROTON_MASS_G / (density * BOLTZMANN_ERG_K)
        else:
            temperature = mu_table.temperature_from_pressure(pressure, density, hydrogen_number_density)
    if fields.dust_relative_abundance is None:
        dust = np.zeros_like(density)
    else:
        dust = np.asarray(sampled_grid[fields.dust_relative_abundance].to_value(""))
    if mu_table is None:
        x_hii = np.zeros(density.shape, dtype=np.float64)
        x_heii = np.zeros(density.shape, dtype=np.float64)
        x_heiii = np.zeros(density.shape, dtype=np.float64)
    else:
        _, x_hii, x_heii, x_heiii = mu_table.state(temperature, hydrogen_number_density)
    helium_number_density = density * (1.0 - hydrogen_mass_fraction) / (4.0 * PROTON_MASS_G)
    staged = StaticRTInput(
        grid=GridSpec(
            cell_width_cm=cell_width_cm,
            left_edge_cm=np.asarray(left_edge.to_value("cm")),
        ),
        hydrogen_number_density_cm3=hydrogen_number_density,
        helium_number_density_cm3=helium_number_density,
        temperature_k=temperature,
        dust_relative_abundance=dust,
        x_hii=x_hii,
        x_heii=x_heii,
        x_heiii=x_heiii,
        sources=sources,
    )
    write_static_rt_input(output_path, staged)
    return staged


def stage_ramses_with_yt(
    info_path: str | Path,
    output_path: str | Path,
    *,
    level: int,
    dimensions: tuple[int, int, int],
    fields: RamsesFieldMap,
    left_edge_code: tuple[float, float, float] = (0.0, 0.0, 0.0),
    right_edge_code: tuple[float, float, float] | None = None,
    sources: SourceCatalog | None = None,
    hydrogen_mass_fraction: float = HYDROGEN_MASS_FRACTION,
    mu_table: PrimordialMuTable | None = None,
    grackle_table: GrackleEquilibriumTable | None = None,
    grackle_metallicity_solar: float = 1.0e-6,
    thermal_atlas: ThermalAtlas | None = None,
    thermal_scale_factor: float | None = None,
) -> StaticRTInput:
    """Stage an explicitly mapped uniform RAMSES subvolume through yt.

    This adapter intentionally relies on yt's RAMSES frontend rather than
    decoding custom Fortran records.  If a legacy particle header prevents yt
    from constructing an index, the raised RuntimeError preserves that fact;
    it must be resolved by a hydro-only reader or an audited frontend patch,
    never by silently changing the particle interpretation.
    """

    try:
        import yt

        return _stage_yt_dataset(
            yt.load(str(info_path)),
            output_path,
            level=level,
            dimensions=dimensions,
            fields=fields,
            left_edge_code=left_edge_code,
            right_edge_code=right_edge_code,
            sources=sources,
            hydrogen_mass_fraction=hydrogen_mass_fraction,
            mu_table=mu_table,
            grackle_table=grackle_table,
            grackle_metallicity_solar=grackle_metallicity_solar,
            thermal_atlas=thermal_atlas,
            thermal_scale_factor=thermal_scale_factor,
        )
    except Exception as exc:
        raise RuntimeError(
            "RAMSES-to-SNRT staging failed before a canonical input was written. "
            "Inspect the RAMSES field map and particle/header compatibility; do not "
            "substitute guessed fields or particle records."
        ) from exc


def stage_ramses_hydro_only(
    info_path: str | Path,
    output_path: str | Path,
    *,
    level: int,
    dimensions: tuple[int, int, int],
    fields: RamsesFieldMap,
    scratch_directory: str | Path,
    left_edge_code: tuple[float, float, float] = (0.0, 0.0, 0.0),
    right_edge_code: tuple[float, float, float] | None = None,
    sources: SourceCatalog | None = None,
    hydrogen_mass_fraction: float = HYDROGEN_MASS_FRACTION,
    mu_table: PrimordialMuTable | None = None,
    grackle_table: GrackleEquilibriumTable | None = None,
    grackle_metallicity_solar: float = 1.0e-6,
    thermal_atlas: ThermalAtlas | None = None,
    thermal_scale_factor: float | None = None,
) -> StaticRTInput:
    """Stage RAMSES hydro while intentionally excluding incompatible particles.

    A temporary `output_XXXXX` view is created below ``scratch_directory``.
    It copies only text metadata and symlinks `amr`/`hydro` rank files, so the
    stock yt AMR reader never opens a legacy particle header.  The source
    snapshot is never modified.
    """

    info = Path(info_path).resolve()
    match = re.fullmatch(r"info_(\d{5})\.txt", info.name)
    if match is None:
        raise ValueError("info_path must be a RAMSES info_XXXXX.txt file")
    output_number = match.group(1)
    source_directory = info.parent
    amr_files = sorted(source_directory.glob(f"amr_{output_number}.out*"))
    hydro_files = sorted(source_directory.glob(f"hydro_{output_number}.out*"))
    if not amr_files or not hydro_files:
        raise FileNotFoundError("RAMSES hydro-only staging requires both amr_ and hydro_ rank files")
    scratch = Path(scratch_directory)
    scratch.mkdir(parents=True, exist_ok=True)
    try:
        from yt.frontends.ramses.data_structures import RAMSESDataset

        with TemporaryDirectory(prefix="snrt_hydro_", dir=scratch) as temporary_root:
            view = Path(temporary_root) / f"output_{output_number}"
            view.mkdir()
            shutil.copy2(info, view / info.name)
            for metadata_name in ("hydro_file_descriptor.txt", "namelist.txt"):
                metadata = source_directory / metadata_name
                if metadata.is_file():
                    shutil.copy2(metadata, view / metadata_name)
            for source in amr_files + hydro_files:
                os.symlink(source, view / source.name)
            dataset = RAMSESDataset(str(view / info.name))
            return _stage_yt_dataset(
                dataset,
                output_path,
                level=level,
                dimensions=dimensions,
                fields=fields,
                left_edge_code=left_edge_code,
                right_edge_code=right_edge_code,
                sources=sources,
                hydrogen_mass_fraction=hydrogen_mass_fraction,
                mu_table=mu_table,
                grackle_table=grackle_table,
                grackle_metallicity_solar=grackle_metallicity_solar,
                thermal_atlas=thermal_atlas,
                thermal_scale_factor=thermal_scale_factor,
            )
    except Exception as exc:
        raise RuntimeError(
            "hydro-only RAMSES staging failed before a canonical input was written. "
            "Preserve the source snapshot and inspect the explicit field map, selected "
            "AMR level, and hydro output compatibility."
        ) from exc
