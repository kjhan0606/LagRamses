#!/usr/bin/env python3
"""Build a Kirchhoff-equilibrium dust thermal/IR sidecar from Draine data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
if str(SNRT_ROOT) not in sys.path:
    sys.path.insert(0, str(SNRT_ROOT))

from snrt_core.dust import DUST_THERMAL_CODE_MANIFEST  # noqa: E402
from snrt_core.provenance import (  # noqa: E402
    PAYLOAD_HASH_SCHEME,
    build_code_manifest,
    canonical_payload_sha256,
)

from tools.build_draine_dust_opacity import (  # noqa: E402
    DEFAULT_GROUP_EDGES,
    DEFAULT_SOURCE_URL,
    EV_PER_MICRON,
    read_draine_table,
    read_group_edges,
)


EV_ERG = 1.602176634e-12
PLANCK_ERG_S = 6.62607015e-27
LIGHT_SPEED_CM_S = 2.99792458e10
BOLTZMANN_EV_K = 8.617333262145e-5
DEFAULT_TEMPERATURE_GRID_K = np.geomspace(5.0, 300.0, 64)


def _sha256(path: Path) -> str:
    from snrt_core.provenance import sha256_file

    return sha256_file(path)


def _planck_power_density(
    energy_ev: np.ndarray,
    absorption_per_h: np.ndarray,
    temperature_k: float,
) -> np.ndarray:
    """Return 4 pi C_abs B_E dE integrand in erg s^-1 H^-1 per eV."""

    energy_erg = energy_ev * EV_ERG
    x = energy_ev / (BOLTZMANN_EV_K * temperature_k)
    occupation = 1.0 / np.expm1(np.minimum(x, 700.0))
    spectral_radiance_per_erg = (
        2.0 * energy_erg**3 / (PLANCK_ERG_S**3 * LIGHT_SPEED_CM_S**2) * occupation
    )
    return 4.0 * np.pi * absorption_per_h * spectral_radiance_per_erg * EV_ERG


def _thermal_rows(
    edges: np.ndarray,
    table: dict[str, object],
    ir_group_indices: np.ndarray,
    temperature_grid_k: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    energy = np.asarray(table["energy_ev"], dtype=np.float64)
    absorption = np.asarray(table["absorption_per_h_cm2"], dtype=np.float64)
    integration_energy = np.unique(np.concatenate((energy, edges)))
    log_absorption = np.interp(
        np.log(integration_energy), np.log(energy), np.log(absorption)
    )
    absorption_grid = np.exp(log_absorption)
    total_power = np.empty(temperature_grid_k.size, dtype=np.float64)
    group_power = np.empty((temperature_grid_k.size, ir_group_indices.size), dtype=np.float64)
    group_photon_energy = np.empty_like(group_power)
    for row, temperature_k in enumerate(temperature_grid_k):
        total_energy = integration_energy[(integration_energy >= energy[0]) & (integration_energy <= energy[-1])]
        total_absorption = np.interp(total_energy, integration_energy, absorption_grid)
        total_integrand = _planck_power_density(total_energy, total_absorption, float(temperature_k))
        total_power[row] = np.trapezoid(total_integrand, total_energy)
        for output, group in enumerate(ir_group_indices):
            lower, upper = edges[group], edges[group + 1]
            selected = (integration_energy >= lower) & (integration_energy <= upper)
            group_energy = integration_energy[selected]
            group_absorption = absorption_grid[selected]
            group_integrand = _planck_power_density(
                group_energy, group_absorption, float(temperature_k)
            )
            power = float(np.trapezoid(group_integrand, group_energy))
            photon_power = float(np.trapezoid(group_integrand / group_energy, group_energy))
            if power <= 0.0 or photon_power <= 0.0:
                raise ValueError(f"IR group {group} has no finite thermal emission")
            group_power[row, output] = power
            group_photon_energy[row, output] = power / photon_power
    fractions = group_power / total_power[:, None]
    untracked = 1.0 - fractions.sum(axis=1)
    if np.any(untracked < -1.0e-10):
        raise ValueError("configured IR groups exceed the full Draine thermal power")
    return total_power, fractions, group_photon_energy, np.maximum(untracked, 0.0)


def build_thermal_metadata(
    source_path: Path,
    group_edges_path: Path = DEFAULT_GROUP_EDGES,
    *,
    ir_group_indices: list[int] | None = None,
    temperature_grid_k: np.ndarray = DEFAULT_TEMPERATURE_GRID_K,
    source_url: str = DEFAULT_SOURCE_URL,
) -> dict[str, object]:
    """Build a candidate single-temperature Kirchhoff thermal sidecar."""

    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    edges = read_group_edges(group_edges_path)
    table = read_draine_table(source_path)
    energy = np.asarray(table["energy_ev"], dtype=np.float64)
    temperatures = np.asarray(temperature_grid_k, dtype=np.float64)
    if temperatures.ndim != 1 or temperatures.size < 2 or not np.isfinite(temperatures).all():
        raise ValueError("temperature grid must be a finite one-dimensional array")
    if np.any(temperatures <= 0.0) or np.any(np.diff(temperatures) <= 0.0):
        raise ValueError("temperature grid must be strictly increasing and positive")
    if ir_group_indices is None:
        ir_group_indices_array = np.flatnonzero(edges[1:] <= 1.0).astype(np.int64)
    else:
        ir_group_indices_array = np.asarray(ir_group_indices, dtype=np.int64)
    if (
        ir_group_indices_array.ndim != 1
        or ir_group_indices_array.size == 0
        or np.any(ir_group_indices_array < 0)
        or np.any(ir_group_indices_array >= edges.size - 1)
        or np.unique(ir_group_indices_array).size != ir_group_indices_array.size
    ):
        raise ValueError("IR group indices are invalid or duplicated")
    if edges[ir_group_indices_array[0]] < energy[0] or np.any(
        edges[ir_group_indices_array + 1] > energy[-1]
    ):
        raise ValueError("Draine table does not cover all configured IR groups")

    total_power, fractions, photon_energy, untracked = _thermal_rows(
        edges, table, ir_group_indices_array, temperatures
    )
    if not np.all(np.diff(total_power) > 0.0):
        raise ValueError("Draine thermal power curve is not strictly increasing")
    metadata: dict[str, object] = {
        "schema": "snrt_dust_thermal_v1",
        "schema_version": 1,
        "status": "candidate_kirchhoff_equilibrium",
        "group_edges_ev": edges.tolist(),
        "group_edges_path": str(group_edges_path.resolve()),
        "group_edges_sha256": _sha256(group_edges_path.resolve()),
        "ir_group_indices": [int(value) for value in ir_group_indices_array],
        "temperature_k": temperatures.tolist(),
        "emitted_power_per_h_erg_s": total_power.tolist(),
        "ir_energy_fraction": fractions.tolist(),
        "ir_mean_photon_energy_ev": photon_energy.tolist(),
        "untracked_energy_fraction": untracked.tolist(),
        "fraction_tolerance": 1.0e-10,
        "reference_mixture": (
            "Draine/Weingartner-Draine carbonaceous-silicate Milky Way model, "
            "R_V=3.1, D03-renormalized, b_C=55.8 ppm"
        ),
        "thermal_source": (
            "Kirchhoff equilibrium derived from Draine absorption cross section "
            "and Planck B_E; no stochastic PAH heating"
        ),
        "single_temperature_assumption": (
            "one equilibrium temperature per cell and reference mixture; "
            "small-grain stochastic heating is deferred"
        ),
        "source_table": {
            "path": str(source_path.resolve()),
            "sha256": _sha256(source_path),
            "url": source_url,
            "row_count": int(table["row_count"]),
            "energy_range_ev": [float(energy[0]), float(energy[-1])],
            "dust_mass_per_h_g": float(table["dust_mass_per_h_g"]),
            "gas_to_dust_mass_ratio": float(table["gas_to_dust_mass_ratio"]),
            "derivation": "4*pi*C_abs(E)*B_E(T), integrated over source-table energy range",
        },
        "algorithm": {
            "opacity_interpolation": "log-log in energy",
            "thermal_quadrature": "trapezoid on the raw Draine grid plus group edges",
            "temperature_interpolation": (
                "linear-in-power interpolation on a log-temperature grid; "
                "32-step log-temperature bisection at runtime"
            ),
            "tracked_energy": "configured IR group energy fractions",
            "untracked_energy": "explicit complement outside configured groups",
            "photon_energy": "emission-power integral divided by photon-number integral per group",
            "cmb": (
                "runtime background power is the same table curve evaluated at T_CMB; "
                "net excess is split with the emitting-temperature fractions, "
                "so tracked IR is conservatively low near T_CMB"
            ),
        },
        "limits": [
            "This is a candidate thermal/emission closure, not a dust-mixture approval.",
            "The IR source is recorded one-pass and is not recursively re-transported.",
            "Stochastic heating, PAH features, dust-gas exchange, and source obscuration are deferred.",
        ],
    }
    metadata["builder"] = {
        "path": str(TOOL_PATH.resolve()),
        "sha256": _sha256(TOOL_PATH.resolve()),
    }
    metadata["closure_code_manifest"] = build_code_manifest(DUST_THERMAL_CODE_MANIFEST)
    metadata["payload_hash_scheme"] = PAYLOAD_HASH_SCHEME
    metadata["payload_sha256"] = canonical_payload_sha256(metadata)
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--group-edges", type=Path, default=DEFAULT_GROUP_EDGES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ir-groups", type=int, nargs="+", help="configured IR group indices")
    parser.add_argument("--temperature-min-k", type=float, default=5.0)
    parser.add_argument("--temperature-max-k", type=float, default=300.0)
    parser.add_argument("--temperature-count", type=int, default=64)
    parser.add_argument("--source-url", default=DEFAULT_SOURCE_URL)
    args = parser.parse_args()
    if args.temperature_min_k <= 0.0 or args.temperature_max_k <= args.temperature_min_k:
        raise ValueError("invalid thermal temperature bounds")
    if args.temperature_count < 2:
        raise ValueError("temperature-count must be at least two")
    metadata = build_thermal_metadata(
        args.source,
        args.group_edges,
        ir_group_indices=args.ir_groups,
        temperature_grid_k=np.geomspace(
            args.temperature_min_k, args.temperature_max_k, args.temperature_count
        ),
        source_url=args.source_url,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(
        "DUST_THERMAL_METADATA_OK "
        f"groups={len(metadata['group_edges_ev']) - 1} "
        f"ir_groups={len(metadata['ir_group_indices'])} temperatures={len(metadata['temperature_k'])} "
        f"source_sha256={metadata['source_table']['sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
