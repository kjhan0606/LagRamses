#!/usr/bin/env python3
"""Convert a BPASS v2.2.1 HDF5 SSP grid into a P0 stellar photon ledger.

The Zenodo BPASS product stores ``L_nu`` in ``L_sun/Hz`` on a
(metallicity, age, wavelength) grid.  This adapter reduces every grid node to
photon-number group moments and only then interpolates those moments to the
native stellar catalogue.  It therefore avoids expanding the 100,000-sample
spectra into a multi-gigabyte CSV while retaining the group closure used by
SNRT.

Two assumptions must be explicit at the command line.  First, the selected
BPASS SSP spectrum is treated as normalized per initial stellar Msun.  The
HDF5 units attribute declares ``L_sun/Hz`` but does not spell out that mass
normalization.  Second, the product starts at 1 Angstrom and ends at 100,000
Angstroms, so the part of the lowest P0 group below 0.123984 eV may be
zero-padded only with ``--pad-low-energy-zero``.  Both choices are recorded as
candidate-ledger limitations.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import h5py
import numpy as np


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
if str(SNRT_ROOT) not in sys.path:
    sys.path.insert(0, str(SNRT_ROOT))

from snrt_core.primordial import (  # noqa: E402
    HE_I_FIT,
    HE_II_FIT,
    H_I_FIT,
    _verner_cross_section_numpy,
)
from tools import p4_build_stellar_photon_ledger as stellar_ledger  # noqa: E402


PLANCK_ERG_S = 6.62607015e-27
SOLAR_LUMINOSITY_ERG_S = 3.827e33
EV_PER_ANGSTROM = 12398.419843320026
DEFAULT_GROUP_EDGES = SNRT_ROOT / "config" / "p0_photon_group_edges_ev.txt"
REQUIRED_DATASETS = ("ages", "metallicities", "spectra", "wavelengths")


def _decode(value: object) -> object:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    return value


def _log_brackets(
    axis_log: np.ndarray,
    query_log: np.ndarray,
    name: str,
    clamp: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    outside = (query_log < axis_log[0]) | (query_log > axis_log[-1])
    outside_count = int(np.count_nonzero(outside))
    if outside_count and not clamp:
        raise ValueError(
            f"{outside_count} sources lie outside BPASS {name} range "
            f"[{10.0 ** axis_log[0]:.8g}, {10.0 ** axis_log[-1]:.8g}]; "
            "use --clamp-table-range explicitly to clamp"
        )
    clipped = np.clip(query_log, axis_log[0], axis_log[-1])
    lower = np.searchsorted(axis_log, clipped, side="right") - 1
    lower = np.clip(lower, 0, axis_log.size - 2)
    upper = lower + 1
    weight = (clipped - axis_log[lower]) / (axis_log[upper] - axis_log[lower])
    return lower, upper, weight, outside_count


def _interpolate_nodes(
    values: np.ndarray,
    metal_lower: np.ndarray,
    metal_upper: np.ndarray,
    metal_weight: np.ndarray,
    age_lower: np.ndarray,
    age_upper: np.ndarray,
    age_weight: np.ndarray,
) -> np.ndarray:
    """Bilinearly interpolate [metallicity, age, ...] node moments."""

    wm = metal_weight.reshape((-1,) + (1,) * (values.ndim - 2))
    wa = age_weight.reshape((-1,) + (1,) * (values.ndim - 2))
    lower_m = values[metal_lower, age_lower]
    upper_age = values[metal_lower, age_upper]
    upper_m = values[metal_upper, age_lower]
    upper_both = values[metal_upper, age_upper]
    return (
        (1.0 - wm) * (1.0 - wa) * lower_m
        + (1.0 - wm) * wa * upper_age
        + wm * (1.0 - wa) * upper_m
        + wm * wa * upper_both
    )


def _node_moments(
    source_path: Path,
    group_edges_ev: np.ndarray,
    pad_low_energy_zero: bool,
) -> dict[str, object]:
    with h5py.File(source_path, "r") as handle:
        missing = [name for name in REQUIRED_DATASETS if name not in handle]
        if missing:
            raise ValueError(f"BPASS HDF5 is missing datasets: {missing}")
        ages_gyr = np.asarray(handle["ages"], dtype=np.float64)
        metallicities_log_solar = np.asarray(handle["metallicities"], dtype=np.float64)
        wavelengths_angstrom = np.asarray(handle["wavelengths"], dtype=np.float64)
        spectra = handle["spectra"]
        if spectra.shape != (
            metallicities_log_solar.size,
            ages_gyr.size,
            wavelengths_angstrom.size,
        ):
            raise ValueError(
                "BPASS spectra shape does not match metallicity, age, and wavelength axes: "
                f"shape={spectra.shape} expected={(metallicities_log_solar.size, ages_gyr.size, wavelengths_angstrom.size)}"
            )
        if (
            ages_gyr.size < 2
            or metallicities_log_solar.size < 2
            or wavelengths_angstrom.size < 2
            or not np.isfinite(ages_gyr).all()
            or not np.isfinite(metallicities_log_solar).all()
            or not np.isfinite(wavelengths_angstrom).all()
            or np.any(ages_gyr <= 0.0)
            or np.any(wavelengths_angstrom <= 0.0)
            or np.any(np.diff(ages_gyr) <= 0.0)
            or np.any(np.diff(metallicities_log_solar) <= 0.0)
            or np.any(np.diff(wavelengths_angstrom) <= 0.0)
        ):
            raise ValueError("BPASS axes must be finite, positive, and strictly increasing")

        energy_ev = EV_PER_ANGSTROM / wavelengths_angstrom
        order = np.argsort(energy_ev)
        energy_ev = energy_ev[order]
        if group_edges_ev[0] < energy_ev[0] and not pad_low_energy_zero:
            raise ValueError(
                f"BPASS spectrum starts at {energy_ev[0]:.8g} eV below the requested "
                f"lowest group edge {group_edges_ev[0]:.8g} eV; supply --pad-low-energy-zero"
            )
        if group_edges_ev[-1] > energy_ev[-1]:
            raise ValueError(
                f"BPASS spectrum ends at {energy_ev[-1]:.8g} eV below the requested "
                f"highest group edge {group_edges_ev[-1]:.8g} eV"
            )
        integration_energy = np.unique(np.concatenate((energy_ev, group_edges_ev)))
        group_masks = [
            (integration_energy >= lower) & (integration_energy <= upper)
            for lower, upper in zip(group_edges_ev[:-1], group_edges_ev[1:], strict=True)
        ]
        fits = (H_I_FIT, HE_I_FIT, HE_II_FIT)
        thresholds = np.asarray([fit.threshold_ev for fit in fits], dtype=np.float64)
        shape = (metallicities_log_solar.size, ages_gyr.size, group_edges_ev.size - 1)
        photon_moment = np.zeros(shape, dtype=np.float64)
        energy_moment = np.zeros_like(photon_moment)
        sigma_moment = np.zeros((3, *shape), dtype=np.float64)
        excess_moment = np.zeros_like(sigma_moment)
        for metal_index in range(metallicities_log_solar.size):
            for age_index in range(ages_gyr.size):
                spectrum = np.asarray(spectra[metal_index, age_index, :], dtype=np.float64)[order]
                if not np.isfinite(spectrum).all() or np.any(spectrum < 0.0):
                    raise ValueError(
                        f"BPASS spectra[{metal_index},{age_index},:] contains invalid values"
                    )
                photon_per_ev = spectrum * SOLAR_LUMINOSITY_ERG_S / (PLANCK_ERG_S * energy_ev)
                integration_spectrum = np.interp(
                    integration_energy,
                    energy_ev,
                    photon_per_ev,
                    left=0.0 if pad_low_energy_zero else None,
                    right=0.0,
                )
                for group, mask in enumerate(group_masks):
                    energies = integration_energy[mask]
                    photons = integration_spectrum[mask]
                    photon_count = float(np.trapezoid(photons, energies))
                    if not np.isfinite(photon_count) or photon_count < 0.0:
                        raise ValueError(f"BPASS node {metal_index},{age_index} has invalid group photons")
                    photon_moment[metal_index, age_index, group] = photon_count
                    energy_moment[metal_index, age_index, group] = np.trapezoid(photons * energies, energies)
                    for species, fit in enumerate(fits):
                        sigma = _verner_cross_section_numpy(energies, fit)
                        weighted_sigma = float(np.trapezoid(photons * sigma, energies))
                        sigma_moment[species, metal_index, age_index, group] = weighted_sigma
                        excess_moment[species, metal_index, age_index, group] = np.trapezoid(
                            photons * sigma * np.maximum(energies - thresholds[species], 0.0),
                            energies,
                        )

        root_attributes = {str(key): _decode(value) for key, value in handle.attrs.items()}
    return {
        "ages_myr": ages_gyr * 1000.0,
        "metallicities_log_solar": metallicities_log_solar,
        "energy_range_ev": [float(energy_ev[0]), float(energy_ev[-1])],
        "photon_moment": photon_moment,
        "energy_moment": energy_moment,
        "sigma_moment": sigma_moment,
        "excess_moment": excess_moment,
        "root_attributes": root_attributes,
    }


def build_ledger(
    catalogue_path: Path,
    bpass_path: Path,
    group_edges_path: Path,
    output_csv: Path,
    output_metadata: Path,
    solar_mass_fraction: float,
    scale_factor: float,
    mass_field: str,
    metallicity_floor_solar: float | None,
    clamp_table_range: bool,
    escape_fraction: float,
    assume_per_initial_mass: bool,
    pad_low_energy_zero: bool,
) -> dict[str, object]:
    if not assume_per_initial_mass:
        raise ValueError("BPASS mass normalization must be acknowledged with --assume-per-initial-mass")
    if not np.isfinite(solar_mass_fraction) or solar_mass_fraction <= 0.0:
        raise ValueError("solar metallicity mass fraction must be finite and positive")
    if not np.isfinite(scale_factor) or not 0.0 < scale_factor <= 1.0:
        raise ValueError("source scale factor must be finite and in (0, 1]")
    if not (0.0 < escape_fraction <= 1.0) or not np.isfinite(escape_fraction):
        raise ValueError("escape fraction must be finite and in (0, 1]")
    for path in (catalogue_path, bpass_path, group_edges_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output_csv.exists() or output_metadata.exists():
        raise FileExistsError("refusing to overwrite an existing BPASS stellar ledger or metadata file")

    group_edges = stellar_ledger._read_group_edges(group_edges_path)
    table = _node_moments(bpass_path, group_edges, pad_low_energy_zero)
    catalogue = stellar_ledger._read_catalogue(
        catalogue_path,
        solar_mass_fraction,
        metallicity_floor_solar,
    )
    age_axis = np.asarray(table["ages_myr"], dtype=np.float64)
    metal_axis_log = np.asarray(table["metallicities_log_solar"], dtype=np.float64)
    query_age_log = np.log10(np.asarray(catalogue["ages_myr"], dtype=np.float64))
    query_metal_log = np.log10(np.asarray(catalogue["metallicities_solar"], dtype=np.float64))
    age_lower, age_upper, age_weight, age_clamped = _log_brackets(
        np.log10(age_axis), query_age_log, "age [Myr]", clamp_table_range
    )
    metal_lower, metal_upper, metal_weight, metal_clamped = _log_brackets(
        metal_axis_log, query_metal_log, "metallicity relative to solar", clamp_table_range
    )

    photon_nodes = np.asarray(table["photon_moment"], dtype=np.float64)
    energy_nodes = np.asarray(table["energy_moment"], dtype=np.float64)
    sigma_nodes = np.asarray(table["sigma_moment"], dtype=np.float64)
    excess_nodes = np.asarray(table["excess_moment"], dtype=np.float64)
    source_photons_per_mass = _interpolate_nodes(
        photon_nodes,
        metal_lower,
        metal_upper,
        metal_weight,
        age_lower,
        age_upper,
        age_weight,
    )
    source_energy_moment = _interpolate_nodes(
        energy_nodes,
        metal_lower,
        metal_upper,
        metal_weight,
        age_lower,
        age_upper,
        age_weight,
    )
    source_sigma_moment = np.stack(
        [
            _interpolate_nodes(
                sigma_nodes[species],
                metal_lower,
                metal_upper,
                metal_weight,
                age_lower,
                age_upper,
                age_weight,
            )
            for species in range(3)
        ],
        axis=1,
    )
    source_excess_moment = np.stack(
        [
            _interpolate_nodes(
                excess_nodes[species],
                metal_lower,
                metal_upper,
                metal_weight,
                age_lower,
                age_upper,
                age_weight,
            )
            for species in range(3)
        ],
        axis=1,
    )
    masses = np.asarray(catalogue[mass_field], dtype=np.float64)
    source_weights = masses * escape_fraction
    q_groups = source_photons_per_mass * source_weights[:, None]
    if not np.isfinite(q_groups).all() or np.any(q_groups < 0.0):
        raise ValueError("interpolated BPASS photon rates are invalid")
    group_totals = q_groups.sum(axis=0, dtype=np.float64)
    if not np.any(group_totals > 0.0):
        raise ValueError("BPASS catalogue emits no photons in the requested groups")

    aggregate_energy_moment = np.sum(source_energy_moment * source_weights[:, None], axis=0, dtype=np.float64)
    aggregate_sigma_moment = np.sum(source_sigma_moment * source_weights[:, None, None], axis=0, dtype=np.float64)
    aggregate_excess_moment = np.sum(source_excess_moment * source_weights[:, None, None], axis=0, dtype=np.float64)
    mean_energy = np.sqrt(group_edges[:-1] * group_edges[1:])
    mean_energy = np.divide(
        aggregate_energy_moment,
        group_totals,
        out=mean_energy.copy(),
        where=group_totals > 0.0,
    )
    averaged_sigma = np.divide(
        aggregate_sigma_moment,
        group_totals[None, :],
        out=np.zeros_like(aggregate_sigma_moment),
        where=group_totals[None, :] > 0.0,
    )
    excess_energy = np.divide(
        aggregate_excess_moment,
        aggregate_sigma_moment,
        out=np.zeros_like(aggregate_excess_moment),
        where=aggregate_sigma_moment > 0.0,
    )
    closure_status = [
        "bpass_node_moment_interpolation" if total > 0.0 else "empty_source_group_zero_photons"
        for total in group_totals
    ]

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_metadata.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source_id",
        "source_kind",
        "x_code",
        "y_code",
        "z_code",
        *[f"q_group_{index}_s" for index in range(group_totals.size)],
        "age_myr",
        "metallicity_solar",
        "mass_msun",
        "initial_mass_msun",
        "source_scale_factor",
    ]
    with output_csv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
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
            for group, value in enumerate(q_groups[index]):
                row[f"q_group_{group}_s"] = f"{value:.17g}"
            writer.writerow(row)

    root_attributes = table["root_attributes"]
    metadata: dict[str, object] = {
        "schema": "stellar_photon_source_ledger_v1",
        "status": "candidate_bpass_stellar_photon_ledger",
        "source_kind": "star",
        "source_count": int(masses.size),
        "source_scale_factor": scale_factor,
        "input_catalogue": str(catalogue_path.resolve()),
        "input_catalogue_sha256": stellar_ledger._sha256(catalogue_path),
        "input_bpass_hdf5": str(bpass_path.resolve()),
        "input_bpass_hdf5_sha256": stellar_ledger._sha256(bpass_path),
        "group_edges_file": str(group_edges_path.resolve()),
        "group_edges_sha256": stellar_ledger._sha256(group_edges_path),
        "group_edges_ev": group_edges.tolist(),
        "bpass_hdf5_contract": {
            "datasets": list(REQUIRED_DATASETS),
            "ages_units": "Gigayears (converted to Myr)",
            "metallicities_units": "log10 relative to Solar (dex)",
            "wavelength_units": "Angstroms (converted to eV)",
            "spectra_units": "L_sun/Hz as declared by HDF5",
            "spectra_mass_normalization": "assumed per initial stellar Msun by explicit CLI acknowledgement",
            "root_attributes": root_attributes,
            "age_grid_myr": age_axis.tolist(),
            "metallicity_log10_solar_grid": metal_axis_log.tolist(),
            "energy_range_ev": table["energy_range_ev"],
        },
        "normalization": {
            "mass_field_used": mass_field,
            "solar_mass_fraction": solar_mass_fraction,
            "metallicity_floor_solar": metallicity_floor_solar,
            "metallicity_floor_count": int(catalogue["metallicity_floor_count"]),
            "escape_fraction": escape_fraction,
            "escape_fraction_interpretation": "source-side multiplicative factor",
            "mass_normalization_acknowledged": assume_per_initial_mass,
        },
        "interpolation_clamped_sources": {
            "age": age_clamped,
            "metallicity": metal_clamped,
        },
        "low_energy_coverage": {
            "policy": "zero_below_table_support" if pad_low_energy_zero else "reject_outside_support",
            "table_minimum_energy_ev": float(np.asarray(table["energy_range_ev"])[0]),
            "requested_minimum_energy_ev": float(group_edges[0]),
        },
        "groups": [
            {
                "index": int(index),
                "energy_interval_ev": [float(group_edges[index]), float(group_edges[index + 1])],
                "photon_weighted_mean_energy_ev": float(mean_energy[index]),
                "total_photon_rate_s": float(group_totals[index]),
                "closure_status": closure_status[index],
            }
            for index in range(group_totals.size)
        ],
        "group_photon_rate_total_s": group_totals.tolist(),
        "group_spectral_closure": {
            "method": "BPASS L_nu-to-photon conversion, node group moments, and log-age/log-metallicity bilinear interpolation",
            "photon_conversion": "q_E = L_nu * L_sun / (h * E_eV)",
            "species_order": ["hydrogen_i", "helium_i", "helium_ii"],
            "cross_sections_cm2": {
                "hydrogen_i": averaged_sigma[0].tolist(),
                "helium_i": averaged_sigma[1].tolist(),
                "helium_ii": averaged_sigma[2].tolist(),
            },
            "photoelectron_excess_energy_ev": {
                "hydrogen_i": excess_energy[0].tolist(),
                "helium_i": excess_energy[1].tolist(),
                "helium_ii": excess_energy[2].tolist(),
            },
            "group_status": closure_status,
            "empty_group_policy": "geometric-mean representative energy and zero closure for groups without BPASS photons",
        },
        "closure_complete": all(status == "bpass_node_moment_interpolation" for status in closure_status),
        "output_csv_sha256": stellar_ledger._sha256(output_csv),
        "algorithm": {
            "group_integration": "trapezoid on the BPASS energy grid plus P0 group boundaries",
            "node_storage": "photon, energy, absorber, and excess-energy moments; no spectral extrapolation",
            "dust_attenuation": "not applied; this is an intrinsic/source-side ledger",
            "gas_transport": "delegated to RT solver",
        },
        "limits": [
            "The BPASS HDF5 spectra mass normalization is an explicit assumption because its units attribute is L_sun/Hz.",
            "The missing spectrum below the table minimum is zero-padded only under the explicit candidate policy.",
            "Nebular emission, source-specific stellar escape, and live stellar evolution are not included.",
        ],
    }
    with output_metadata.open("x", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, required=True)
    parser.add_argument("--bpass-hdf5", type=Path, required=True)
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
    parser.add_argument(
        "--assume-per-initial-mass",
        action="store_true",
        help="acknowledge the candidate assumption that BPASS spectra are per initial stellar Msun",
    )
    parser.add_argument(
        "--pad-low-energy-zero",
        action="store_true",
        help="explicitly set photon spectrum to zero below the BPASS wavelength support",
    )
    args = parser.parse_args()
    try:
        metadata = build_ledger(
            catalogue_path=args.catalogue,
            bpass_path=args.bpass_hdf5,
            group_edges_path=args.group_edges,
            output_csv=args.output_csv,
            output_metadata=args.output_metadata,
            solar_mass_fraction=args.solar_metal_mass_fraction,
            scale_factor=args.scale_factor,
            mass_field=args.mass_field,
            metallicity_floor_solar=args.metallicity_floor_solar,
            clamp_table_range=args.clamp_table_range,
            escape_fraction=args.escape_fraction,
            assume_per_initial_mass=args.assume_per_initial_mass,
            pad_low_energy_zero=args.pad_low_energy_zero,
        )
    except (FileExistsError, FileNotFoundError, OSError, ValueError) as exc:
        parser.error(str(exc))
    print(
        "BPASS_STELLAR_PHOTON_LEDGER_OK "
        f"sources={metadata['source_count']} groups={len(metadata['group_edges_ev']) - 1} "
        f"output={args.output_csv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
