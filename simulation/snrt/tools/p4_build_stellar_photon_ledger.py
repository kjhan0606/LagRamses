#!/usr/bin/env python3
"""Convert a tabulated stellar SED into a source photon ledger.

The SED table contract is deliberately narrow: each row contains

    age_myr, metallicity_solar, energy_ev,
    photon_rate_per_msun_per_ev_s

The last column is a photon-number spectrum per *initial stellar mass* in
photons s^-1 eV^-1 Msun^-1.  Energy/luminosity spectra are not accepted by
this tool, which prevents an unrecorded energy-to-photon conversion from
entering the source ledger.

Interpolation is linear in age and metallicity after transforming both
coordinates to log10.  The spectrum itself is interpolated linearly in
photon-number space, preserving exact zeroes.  Sources outside the table
range are rejected unless --clamp-table-range is supplied explicitly.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


TOOL_PATH = Path(__file__).resolve()
PROJECT_ROOT = TOOL_PATH.parents[3]
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_GROUP_EDGES = SNRT_ROOT / "config" / "p0_photon_group_edges_ev.txt"
REQUIRED_SED_COLUMNS = (
    "age_myr",
    "metallicity_solar",
    "energy_ev",
    "photon_rate_per_msun_per_ev_s",
)
REQUIRED_CATALOGUE_COLUMNS = (
    "source_id",
    "position_x_code",
    "position_y_code",
    "position_z_code",
    "age_myr",
    "birth_metallicity_mass_fraction",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_group_edges(path: Path) -> np.ndarray:
    values: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            fields = line.replace(",", " ").split()
            if len(fields) != 1:
                raise ValueError(
                    f"{path}:{line_number}: expected one edge value, got {raw_line.rstrip()!r}"
                )
            values.append(float(fields[0]))
    edges = np.asarray(values, dtype=np.float64)
    if edges.size < 2 or not np.all(np.isfinite(edges)) or np.any(edges <= 0.0):
        raise ValueError(f"{path}: group edges must be finite and strictly positive")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError(f"{path}: group edges must be strictly increasing")
    return edges


def _read_sed_table(path: Path) -> dict[str, np.ndarray | int]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing CSV header")
        missing = [name for name in REQUIRED_SED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path}: missing SED columns: {', '.join(missing)}")

        rows: list[tuple[float, float, float, float]] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                values = tuple(float(row[name]) for name in REQUIRED_SED_COLUMNS)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_number}: non-numeric SED row") from exc
            rows.append(values)

    if not rows:
        raise ValueError(f"{path}: SED table is empty")
    data = np.asarray(rows, dtype=np.float64)
    ages = np.unique(data[:, 0])
    metallicities = np.unique(data[:, 1])
    energies = np.unique(data[:, 2])
    if ages.size < 2 or metallicities.size < 2 or energies.size < 2:
        raise ValueError("SED table needs at least two age, metallicity, and energy samples")
    if (
        not np.all(np.isfinite(data))
        or np.any(ages <= 0.0)
        or np.any(metallicities <= 0.0)
        or np.any(energies <= 0.0)
        or np.any(data[:, 3] < 0.0)
    ):
        raise ValueError("SED table has invalid age, metallicity, energy, or photon-rate values")
    expected_rows = ages.size * metallicities.size * energies.size
    if data.shape[0] != expected_rows:
        raise ValueError(
            "SED table is not a complete rectangular grid: "
            f"rows={data.shape[0]} expected={expected_rows}"
        )

    spectrum = np.full((ages.size, metallicities.size, energies.size), np.nan, dtype=np.float64)
    age_index = {value: index for index, value in enumerate(ages)}
    metallicity_index = {value: index for index, value in enumerate(metallicities)}
    energy_index = {value: index for index, value in enumerate(energies)}
    for age, metallicity, energy, rate in data:
        index = (age_index[age], metallicity_index[metallicity], energy_index[energy])
        if np.isfinite(spectrum[index]):
            raise ValueError("SED table contains duplicate age/metallicity/energy rows")
        spectrum[index] = rate
    if not np.all(np.isfinite(spectrum)):
        raise ValueError("SED table is missing one or more age/metallicity/energy cells")
    if np.any(np.diff(energies) <= 0.0):
        raise ValueError("SED energy samples must be strictly increasing")
    return {
        "ages_myr": ages,
        "metallicities_solar": metallicities,
        "energies_ev": energies,
        "spectrum": spectrum,
        "row_count": int(data.shape[0]),
    }


def _read_catalogue(path: Path, solar_mass_fraction: float, metallicity_floor: float | None) -> dict[str, np.ndarray | int]:
    if solar_mass_fraction <= 0.0 or not np.isfinite(solar_mass_fraction):
        raise ValueError("solar metallicity mass fraction must be finite and positive")
    if metallicity_floor is not None and (
        metallicity_floor <= 0.0 or not np.isfinite(metallicity_floor)
    ):
        raise ValueError("metallicity floor must be finite and positive")

    source_ids: list[int] = []
    x_code: list[float] = []
    y_code: list[float] = []
    z_code: list[float] = []
    ages_myr: list[float] = []
    metallicity_mass_fraction: list[float] = []
    masses: dict[str, list[float]] = {
        "mass_msun": [],
        "initial_mass_msun": [],
    }
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing catalogue CSV header")
        missing = [name for name in REQUIRED_CATALOGUE_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path}: missing catalogue columns: {', '.join(missing)}")
        for name in masses:
            if name not in reader.fieldnames:
                raise ValueError(f"{path}: missing catalogue mass column {name!r}")
        for line_number, row in enumerate(reader, start=2):
            try:
                source_ids.append(int(row["source_id"]))
                x_code.append(float(row["position_x_code"]))
                y_code.append(float(row["position_y_code"]))
                z_code.append(float(row["position_z_code"]))
                ages_myr.append(float(row["age_myr"]))
                metallicity_mass_fraction.append(float(row["birth_metallicity_mass_fraction"]))
                for name in masses:
                    masses[name].append(float(row[name]))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_number}: non-numeric catalogue row") from exc

    if not source_ids:
        raise ValueError(f"{path}: catalogue is empty")
    ids = np.asarray(source_ids, dtype=np.int64)
    positions = np.column_stack((x_code, y_code, z_code)).astype(np.float64)
    ages = np.asarray(ages_myr, dtype=np.float64)
    metallicities = np.asarray(metallicity_mass_fraction, dtype=np.float64)
    mass_arrays = {name: np.asarray(values, dtype=np.float64) for name, values in masses.items()}
    if np.unique(ids).size != ids.size:
        raise ValueError(f"{path}: source_id values are not unique")
    if (
        not np.all(np.isfinite(positions))
        or np.any(positions < 0.0)
        or np.any(positions > 1.0)
        or not np.all(np.isfinite(ages))
        or np.any(ages <= 0.0)
        or not np.all(np.isfinite(metallicities))
        or np.any(metallicities < 0.0)
    ):
        raise ValueError(f"{path}: invalid source positions, ages, or metallicities")
    for name, values in mass_arrays.items():
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError(f"{path}: {name} must be finite and positive")

    metallicities_solar = metallicities / solar_mass_fraction
    floored = metallicities_solar <= 0.0
    if np.any(floored):
        if metallicity_floor is None:
            raise ValueError(
                f"{int(np.count_nonzero(floored))} sources have zero metallicity; "
                "supply --metallicity-floor-solar explicitly"
            )
        metallicities_solar[floored] = metallicity_floor
    return {
        "source_ids": ids,
        "positions": positions,
        "ages_myr": ages,
        "metallicities_solar": metallicities_solar,
        "mass_msun": mass_arrays["mass_msun"],
        "initial_mass_msun": mass_arrays["initial_mass_msun"],
        "row_count": int(ids.size),
        "metallicity_floor_count": int(np.count_nonzero(floored)),
    }


def _brackets(axis: np.ndarray, values: np.ndarray, name: str, clamp: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    log_axis = np.log10(axis)
    log_values = np.log10(values)
    outside = (log_values < log_axis[0]) | (log_values > log_axis[-1])
    outside_count = int(np.count_nonzero(outside))
    if outside_count and not clamp:
        raise ValueError(
            f"{outside_count} sources lie outside SED {name} range "
            f"[{axis[0]:.8g}, {axis[-1]:.8g}]; use --clamp-table-range explicitly to clamp"
        )
    clipped = np.clip(log_values, log_axis[0], log_axis[-1])
    lower = np.searchsorted(log_axis, clipped, side="right") - 1
    lower = np.clip(lower, 0, axis.size - 2)
    upper = lower + 1
    weight = (clipped - log_axis[lower]) / (log_axis[upper] - log_axis[lower])
    return lower, upper, weight, outside_count


def _interpolate_spectrum(
    table: dict[str, np.ndarray | int],
    ages_myr: np.ndarray,
    metallicities_solar: np.ndarray,
    clamp: bool,
) -> tuple[np.ndarray, dict[str, int]]:
    table_ages = np.asarray(table["ages_myr"])
    table_metallicities = np.asarray(table["metallicities_solar"])
    table_spectrum = np.asarray(table["spectrum"])
    age_lo, age_hi, age_weight, age_clamped = _brackets(table_ages, ages_myr, "age", clamp)
    metal_lo, metal_hi, metal_weight, metal_clamped = _brackets(
        table_metallicities, metallicities_solar, "metallicity", clamp
    )

    result = np.empty((ages_myr.size, table_spectrum.shape[2]), dtype=np.float64)
    for index in range(ages_myr.size):
        wa = age_weight[index]
        wz = metal_weight[index]
        result[index] = (
            (1.0 - wa) * (1.0 - wz) * table_spectrum[age_lo[index], metal_lo[index]]
            + wa * (1.0 - wz) * table_spectrum[age_hi[index], metal_lo[index]]
            + (1.0 - wa) * wz * table_spectrum[age_lo[index], metal_hi[index]]
            + wa * wz * table_spectrum[age_hi[index], metal_hi[index]]
        )
    return result, {"age": age_clamped, "metallicity": metal_clamped}


def _load_primordial_closure():
    sys.path.insert(0, str(SNRT_ROOT))
    from snrt_core.primordial import sed_weighted_group_closure

    return sed_weighted_group_closure


def build_ledger(
    catalogue_path: Path,
    sed_table_path: Path,
    group_edges_path: Path,
    output_csv: Path,
    output_metadata: Path,
    solar_mass_fraction: float,
    scale_factor: float,
    mass_field: str = "initial_mass_msun",
    metallicity_floor_solar: float | None = None,
    clamp_table_range: bool = False,
    escape_fraction: float = 1.0,
) -> dict[str, object]:
    if not (0.0 < scale_factor <= 1.0) or not np.isfinite(scale_factor):
        raise ValueError("source scale factor must be finite and in (0, 1]")
    if not (0.0 < escape_fraction <= 1.0) or not np.isfinite(escape_fraction):
        raise ValueError("escape fraction must be finite and in (0, 1]")
    if mass_field not in {"mass_msun", "initial_mass_msun"}:
        raise ValueError("mass field must be mass_msun or initial_mass_msun")
    for path in (catalogue_path, sed_table_path, group_edges_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output_csv.exists() or output_metadata.exists():
        raise FileExistsError("refusing to overwrite an existing stellar ledger or metadata file")

    edges = _read_group_edges(group_edges_path)
    table = _read_sed_table(sed_table_path)
    energies = np.asarray(table["energies_ev"])
    if energies[0] > edges[0] or energies[-1] < edges[-1]:
        raise ValueError(
            "SED energy range does not cover all group boundaries: "
            f"table=[{energies[0]:.8g}, {energies[-1]:.8g}] eV, "
            f"groups=[{edges[0]:.8g}, {edges[-1]:.8g}] eV"
        )
    catalogue = _read_catalogue(catalogue_path, solar_mass_fraction, metallicity_floor_solar)
    spectra, clamped = _interpolate_spectrum(
        table,
        np.asarray(catalogue["ages_myr"]),
        np.asarray(catalogue["metallicities_solar"]),
        clamp_table_range,
    )
    masses = np.asarray(catalogue[mass_field])
    integration_grid = np.unique(np.concatenate((energies, edges)))
    group_masks = [
        (integration_grid >= edges[index]) & (integration_grid <= edges[index + 1])
        for index in range(edges.size - 1)
    ]
    q_groups = np.empty((masses.size, edges.size - 1), dtype=np.float64)
    aggregate_spectrum = np.zeros(energies.size, dtype=np.float64)
    for index, spectrum in enumerate(spectra):
        aggregate_spectrum += masses[index] * spectrum
        grid_spectrum = np.interp(integration_grid, energies, spectrum)
        for group_index, mask in enumerate(group_masks):
            q_groups[index, group_index] = (
                masses[index]
                * escape_fraction
                * np.trapezoid(grid_spectrum[mask], integration_grid[mask])
            )
    if not np.all(np.isfinite(q_groups)) or np.any(q_groups < 0.0):
        raise ValueError("integrated stellar group photon rates are invalid")
    group_totals = np.sum(q_groups, axis=0)
    if not np.any(group_totals > 0.0):
        raise ValueError("stellar SED emits no photons in any requested group")

    sed_weighted_group_closure = _load_primordial_closure()
    # A stellar-only SED may have no photons in the hard X-ray groups.  Build
    # the physical closure independently for populated groups.  Empty groups
    # retain a monotonic representative energy and zero opacity; their values
    # cannot be reused if another source population later injects photons in
    # those groups.
    mean_energy_ev = np.sqrt(edges[:-1] * edges[1:])
    averaged_sigma = np.zeros((3, edges.size - 1), dtype=np.float64)
    excess_energy_ev = np.zeros_like(averaged_sigma)
    closure_status: list[str] = []
    for group_index, (lower, upper) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        if group_totals[group_index] <= 0.0:
            closure_status.append("empty_source_group_zero_photons")
            continue
        group_closure = sed_weighted_group_closure(
            np.asarray((lower, upper), dtype=np.float64), energies, aggregate_spectrum
        )
        mean_energy_ev[group_index] = float(group_closure.photon_weighted_energy_ev[0])
        averaged_sigma[0, group_index] = float(group_closure.cross_sections.hydrogen_i[0])
        averaged_sigma[1, group_index] = float(group_closure.cross_sections.helium_i[0])
        averaged_sigma[2, group_index] = float(group_closure.cross_sections.helium_ii[0])
        excess_energy_ev[:, group_index] = np.asarray(
            group_closure.photoelectron_excess_energy_ev[:, 0], dtype=np.float64
        )
        closure_status.append("sed_weighted_aggregate")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_metadata.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "source_id",
        "source_kind",
        "x_code",
        "y_code",
        "z_code",
        *[f"q_group_{index}_s" for index in range(edges.size - 1)],
        "age_myr",
        "metallicity_solar",
        "mass_msun",
        "initial_mass_msun",
        "source_scale_factor",
    ]
    with output_csv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index, source_id in enumerate(np.asarray(catalogue["source_ids"])):
            row: dict[str, object] = {
                "source_id": int(source_id),
                "source_kind": "star",
                "x_code": f"{catalogue['positions'][index, 0]:.17g}",
                "y_code": f"{catalogue['positions'][index, 1]:.17g}",
                "z_code": f"{catalogue['positions'][index, 2]:.17g}",
                "age_myr": f"{catalogue['ages_myr'][index]:.17g}",
                "metallicity_solar": f"{catalogue['metallicities_solar'][index]:.17g}",
                "mass_msun": f"{catalogue['mass_msun'][index]:.17g}",
                "initial_mass_msun": f"{catalogue['initial_mass_msun'][index]:.17g}",
                "source_scale_factor": f"{scale_factor:.17g}",
            }
            for group_index, value in enumerate(q_groups[index]):
                row[f"q_group_{group_index}_s"] = f"{value:.17g}"
            writer.writerow(row)

    cross_sections_cm2 = {
        "hydrogen_i": averaged_sigma[0].tolist(),
        "helium_i": averaged_sigma[1].tolist(),
        "helium_ii": averaged_sigma[2].tolist(),
    }
    excess_by_species = {
        "hydrogen_i": excess_energy_ev[0].tolist(),
        "helium_i": excess_energy_ev[1].tolist(),
        "helium_ii": excess_energy_ev[2].tolist(),
    }
    closure_metadata = {
        "method": "photon-number-weighted Verner cross sections; absorber-weighted photoelectron excess energy",
        "species_order": ["hydrogen_i", "helium_i", "helium_ii"],
        "cross_sections_cm2": cross_sections_cm2,
        "photoelectron_excess_energy_ev": excess_by_species,
    }
    groups = [
        {
            "index": int(index),
            "energy_interval_ev": [float(edges[index]), float(edges[index + 1])],
            "photon_weighted_mean_energy_ev": float(mean_energy_ev[index]),
            "total_photon_rate_s": float(group_totals[index]),
            "closure_status": closure_status[index],
        }
        for index in range(edges.size - 1)
    ]
    metadata: dict[str, object] = {
        "schema": "stellar_photon_source_ledger_v1",
        "status": "complete_stellar_photon_ledger",
        "source_kind": "star",
        "source_count": int(masses.size),
        "source_scale_factor": scale_factor,
        "input_catalogue": str(catalogue_path),
        "input_catalogue_sha256": _sha256(catalogue_path),
        "input_sed_table": str(sed_table_path),
        "input_sed_table_sha256": _sha256(sed_table_path),
        "group_edges_file": str(group_edges_path),
        "group_edges_sha256": _sha256(group_edges_path),
        "group_edges_ev": edges.tolist(),
        "sed_table_contract": {
            "columns": list(REQUIRED_SED_COLUMNS),
            "photon_rate_units": "photons s^-1 eV^-1 Msun^-1",
            "mass_normalization": "initial stellar mass",
            "age_interpolation_coordinate": "log10(age_myr)",
            "metallicity_interpolation_coordinate": "log10(metallicity_solar)",
            "spectrum_interpolation": "linear photon-number space",
            "table_range_policy": "clamp" if clamp_table_range else "reject",
            "age_grid_myr": np.asarray(table["ages_myr"]).tolist(),
            "metallicity_grid_solar": np.asarray(table["metallicities_solar"]).tolist(),
            "energy_grid_ev": np.asarray(table["energies_ev"]).tolist(),
            "row_count": int(table["row_count"]),
        },
        "normalization": {
            "mass_field_used": mass_field,
            "solar_mass_fraction": solar_mass_fraction,
            "metallicity_floor_solar": metallicity_floor_solar,
            "metallicity_floor_count": int(catalogue["metallicity_floor_count"]),
            "escape_fraction": escape_fraction,
            "escape_fraction_interpretation": "source-side multiplicative factor",
        },
        "interpolation_clamped_sources": clamped,
        "groups": groups,
        "group_photon_rate_total_s": group_totals.tolist(),
        "group_spectral_closure": {
            **closure_metadata,
            "group_status": closure_status,
            "empty_group_policy": "geometric-mean representative energy and zero opacity/excess; not reusable for later injected photons",
        },
        "closure_complete": all(status == "sed_weighted_aggregate" for status in closure_status),
        "output_csv_sha256": _sha256(output_csv),
        "algorithm": {
            "group_integration": "trapezoid on SED energy grid plus group boundaries",
            "closure_spectrum": "aggregate interpolated stellar spectrum weighted by selected source mass",
            "dust_attenuation": "not applied; this is an intrinsic/source-side ledger",
            "gas_transport": "delegated to RT solver",
        },
    }
    with output_metadata.open("x", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, required=True)
    parser.add_argument("--sed-table", type=Path, required=True)
    parser.add_argument("--group-edges", type=Path, default=DEFAULT_GROUP_EDGES)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-metadata", type=Path, required=True)
    parser.add_argument("--solar-metal-mass-fraction", type=float, required=True)
    parser.add_argument("--scale-factor", type=float, required=True)
    parser.add_argument(
        "--mass-field",
        choices=("mass_msun", "initial_mass_msun"),
        default="initial_mass_msun",
    )
    parser.add_argument("--metallicity-floor-solar", type=float)
    parser.add_argument("--clamp-table-range", action="store_true")
    parser.add_argument("--escape-fraction", type=float, default=1.0)
    args = parser.parse_args()
    try:
        metadata = build_ledger(
            catalogue_path=args.catalogue,
            sed_table_path=args.sed_table,
            group_edges_path=args.group_edges,
            output_csv=args.output_csv,
            output_metadata=args.output_metadata,
            solar_mass_fraction=args.solar_metal_mass_fraction,
            scale_factor=args.scale_factor,
            mass_field=args.mass_field,
            metallicity_floor_solar=args.metallicity_floor_solar,
            clamp_table_range=args.clamp_table_range,
            escape_fraction=args.escape_fraction,
        )
    except (FileExistsError, FileNotFoundError, OSError, ValueError) as exc:
        parser.error(str(exc))
    print(
        "STELLAR_PHOTON_LEDGER_OK "
        f"sources={metadata['source_count']} groups={len(metadata['group_edges_ev']) - 1}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
