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
    parser.add_argument("--native-output", type=Path,
                        help="also write a reference-only native v3/v4 dust namelist")
    parser.add_argument("--native-source-ledger", type=Path,
                        help="AGN photon ledger supplying the actual group representative energies")
    parser.add_argument("--native-background-k", type=float, default=10.0)
    parser.add_argument("--native-scattering", choices=("none", "isotropic_elastic"), default="none",
                        help="optional primary elastic isotropization comparison; no recoil or IR scattering")
    parser.add_argument("--native-gas-exchange", choices=("none", "hydrogen_accommodation"), default="none")
    parser.add_argument("--collision-area-per-h", type=float,
                        help="explicit geometric grain cross section per reference H, cm2; not optical opacity")
    parser.add_argument("--accommodation", type=float, help="explicit thermal accommodation fraction (0,1]")
    parser.add_argument("--reference-heat-capacity-per-h", type=float,
                        help="explicit constant test capacity, erg/H/K; NOT inferred from opacity")
    parser.add_argument("--material-energy-table", type=Path,
                        help="v4 material JSON: temperature_k, internal_energy_erg_g, source_id, source_url, composition")
    args = parser.parse_args()
    if args.native_scattering != "none" and args.native_output is None:
        parser.error("native scattering selection requires --native-output")
    if (args.native_gas_exchange != "none" or args.collision_area_per_h is not None or
            args.accommodation is not None) and args.native_output is None:
        parser.error("gas-exchange selection requires --native-output")
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
    native = None
    if args.native_output is not None:
        if args.native_output.resolve() == args.output.resolve():
            parser.error("JSON and native output must use different paths")
        if args.native_source_ledger is None:
            parser.error("native output requires a source ledger")
        if (args.reference_heat_capacity_per_h is None) == (args.material_energy_table is None):
            parser.error("choose exactly one: reference heat capacity or material energy table")
        if args.native_output.exists():
            raise FileExistsError(args.native_output)
        native = build_native_reference_namelist(
            args.source, args.native_source_ledger, args.group_edges,
            np.asarray(metadata["temperature_k"]), args.native_background_k,
            args.reference_heat_capacity_per_h, material_path=args.material_energy_table,
            scattering=args.native_scattering, gas_exchange=args.native_gas_exchange,
            collision_area=args.collision_area_per_h, accommodation=args.accommodation)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    if native is not None:
        args.native_output.parent.mkdir(parents=True, exist_ok=True)
        with args.native_output.open("x", encoding="utf-8") as handle:
            handle.write(native)
        print(f"DUST_NATIVE_REFERENCE_OK path={args.native_output} production_approved=false")
    print(
        "DUST_THERMAL_METADATA_OK "
        f"groups={len(metadata['group_edges_ev']) - 1} "
        f"ir_groups={len(metadata['ir_group_indices'])} temperatures={len(metadata['temperature_k'])} "
        f"source_sha256={metadata['source_table']['sha256']}"
    )
    return 0


def build_native_reference_namelist(source: Path, ledger_path: Path, edges_path: Path,
                                   temperatures: np.ndarray, background: float,
                                   reference_capacity: float | None,
                                   *, material_path: Path | None = None, scattering: str = "none",
                                   gas_exchange: str = "none", collision_area: float | None = None,
                                   accommodation: float | None = None) -> str:
    """Export actual Draine optics with constant (v3) or tabulated U(T) (v4).

    Primary groups are monochromatic at the source ledger's representative
    energy: absorption and deposited energy use that same energy. This is
    an explicit grey approximation, not a source-spectrum weighted opacity.
    IR uses the full-domain spectral quadrature. Material identity and
    composition must be explicit; admitting a U(T) table does NOT establish
    its compatibility with WD01 optics or confer production approval.
    """
    import jax
    from snrt_core.dust_ir import prepare_spectral_table

    if (reference_capacity is None) == (material_path is None):
        raise ValueError("choose exactly one material representation")
    if scattering not in ("none", "isotropic_elastic"):
        raise ValueError("unknown primary scattering model")
    if scattering != "none" and material_path is None:
        raise ValueError("primary scattering comparison requires native v4")
    if gas_exchange == "none":
        if collision_area is not None or accommodation is not None:
            raise ValueError("collision parameters require an enabled gas-exchange model")
    elif gas_exchange == "hydrogen_accommodation":
        if material_path is None:
            raise ValueError("gas-exchange comparison requires native v4")
        if (collision_area is None or accommodation is None or not np.isfinite(collision_area) or
                collision_area <= 0 or not np.isfinite(accommodation) or not 0 < accommodation <= 1):
            raise ValueError("gas exchange requires explicit positive geometric area and accommodation in (0,1]")
    else:
        raise ValueError("unknown gas-exchange model")
    if reference_capacity is not None:
        if not np.isfinite(reference_capacity) or reference_capacity <= 0:
            raise ValueError("reference heat capacity must be finite and positive")
    edges = read_group_edges(edges_path)
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    if not np.array_equal(np.asarray(ledger["group_edges_ev"]), edges):
        raise ValueError("source ledger group edges differ from configured edges")
    if ledger["group_edges_sha256"] != _sha256(edges_path):
        raise ValueError("source ledger group-edge checksum does not match")
    groups = ledger["groups"]
    if [g["index"] for g in groups] != list(range(len(edges)-1)):
        raise ValueError("source ledger groups are missing, duplicated or reordered")
    means = np.asarray([g["photon_weighted_mean_energy_ev"] for g in groups])
    if (not np.isfinite(means).all() or np.any(means < edges[:-1]) or
            np.any(means > edges[1:])):
        raise ValueError("source ledger representative energy lies outside its group")
    raw = read_draine_table(source)
    energy = np.asarray(raw["energy_ev"])
    opacity = np.asarray(raw["absorption_per_h_cm2"])
    if edges[0] < energy[0] or edges[-1] > energy[-1]:
        raise ValueError("raw dust opacity does not cover the primary groups")
    primary = np.exp(np.interp(np.log(means), np.log(energy), np.log(opacity)))
    with jax.enable_x64(True):
        spectral = prepare_spectral_table(energy, opacity, temperatures, background)
    temperatures = np.unique(np.append(temperatures, background))
    material_u = None
    material_description = "Draine Kirchhoff IR; constant reference heat capacity is NOT physical input"
    if material_path is not None:
        material = json.loads(material_path.read_text(encoding="utf-8"))
        if material.get("schema") != "snrt_dust_material_energy_v1":
            raise ValueError("unknown material energy table schema")
        for key in ("source_id", "source_url", "composition"):
            if not isinstance(material.get(key), str) or not material[key].strip():
                raise ValueError(f"missing material identity: {key}")
        if material.get("energy_zero") != "U(0)=0; no zero-point energy":
            raise ValueError("material internal-energy zero must be explicit")
        mt = np.asarray(material["temperature_k"], dtype=float)
        mu = np.asarray(material["internal_energy_erg_g"], dtype=float)
        if (mt.ndim != 1 or mu.shape != mt.shape or len(mt) < 2 or
                not np.isfinite(mt).all() or not np.isfinite(mu).all() or
                np.any(mt <= 0) or np.any(mu <= 0) or
                np.any(np.diff(mt) <= 0) or np.any(np.diff(mu) <= 0)):
            raise ValueError("material T and U must be finite, positive and increasing")
        if temperatures[0] < mt[0] or temperatures[-1] > mt[-1]:
            raise ValueError("material table cannot extrapolate to the requested temperature/bath")
        # Use the union: discarding material knots would change the supplied
        # U(log T) function and its inverse. Rebuild the Planck table there.
        temperatures = np.unique(np.concatenate((
            temperatures, mt[(mt >= temperatures[0]) & (mt <= temperatures[-1])])))
        if len(temperatures) > 256:
            raise ValueError("union of material and emission knots exceeds native bounds")
        material_u = np.interp(np.log(temperatures), np.log(mt), mu) * raw["dust_mass_per_h_g"]
        with jax.enable_x64(True):
            spectral = prepare_spectral_table(energy, opacity, temperatures, background)
        material_description = "Draine Kirchhoff IR; tabulated U(log T); " + material["source_id"]
        if len(material_description) > 256 or any(c in material_description for c in "\n\r'\""):
            raise ValueError("material source_id cannot be encoded safely in native namelist")
    nodes = np.asarray(spectral.energy_ev)
    weights = np.asarray(spectral.weights_ev)
    sigma = np.asarray(spectral.absorption_per_h_cm2)
    if len(groups) > 32 or len(temperatures) > 256 or len(nodes) > 256:
        raise ValueError("native dust table dimensions exceed the compiled contract bounds")
    rows = [
        ("! Draine WD01/RV3.1 optics + tabulated material U(T)." if material_path is not None
         else "! Draine WD01/RV3.1 optics + explicit constant TEST heat capacity."),
        "! Primary opacity is sampled at each source representative energy (grey approximation).",
        "! IR uses full raw-domain quadrature; scattering and stochastic heating are absent.",
        "! This file does not confer production approval.",
        f"! Raw optical source: {DEFAULT_SOURCE_URL}",
        "&snrt_dust_contract",
        f" contract_version={4 if material_path is not None else 3}, ngroups_input={len(groups)}, ntemperature_input={len(temperatures)},",
        " opacity_status='reference_control', thermal_status='reference_thermal_control',",
        " source_id='draine_wd01_rv31_monochromatic_source_groups',",
        f" source_sha256='{_sha256(ledger_path)}',",
        f" source_table_sha256='{_sha256(source)}',",
        f" group_edges_sha256='{_sha256(edges_path)}', approval_id='',",
        f" thermal_source='{material_description}',",
        f" mass_per_h_input={raw['dust_mass_per_h_g']:.17e},",
        f" heat_capacity_per_h_erg_k_input={0.0 if material_path is not None else reference_capacity:.17e},",
        f" nir_input={len(nodes)}, ir_status='reference_ir_control', ir_background_input={background:.17e},",
    ]
    arrays = dict(edges_input=edges, absorption_input=primary, mean_energy_input=means,
                  temperature_input=temperatures, power_input=np.asarray(spectral.table.power),
                  ir_energy_input=nodes, ir_weight_input=weights, ir_absorption_input=sigma)
    if material_path is not None:
        rows.append(f" material_sha256='{_sha256(material_path)}',")
        arrays["internal_energy_per_h_erg_input"] = material_u
    if scattering == "isotropic_elastic":
        # C_sca = C_ext * albedo. Keep exact zero endpoints without taking
        # log(0); positive segments use the same log-log sampling as absorption.
        sca = np.asarray(raw["scattering_per_h_cm2"])
        hi = np.clip(np.searchsorted(energy, means, side="right"), 1, len(energy)-1)
        lo = hi-1
        f = np.log(means/energy[lo])/np.log(energy[hi]/energy[lo])
        sampled = (1-f)*sca[lo]+f*sca[hi]
        positive = (sca[lo] > 0) & (sca[hi] > 0)
        sampled[positive] = np.exp((1-f[positive])*np.log(sca[lo[positive]]) +
                                  f[positive]*np.log(sca[hi[positive]]))
        rows[2] = "! Primary isotropic elastic scattering ON; IR scattering and stochastic heating absent."
        rows.insert(3, "! Isotropic is a comparison approximation, NOT the measured Draine phase function.")
        rows.insert(4, "! First-order split; no radiation pressure/recoil; not asymptotic-preserving diffusion.")
        rows.append(" scattering_model='isotropic_elastic',")
        arrays["scattering_input"] = sampled
    if gas_exchange != "none":
        rows.insert(1, "! Hydrogen-equivalent geometric gas/dust accommodation; no electron/ion Coulomb model.")
        rows.insert(2, "! Joint implicit gas/dust/IR thermal solve; frozen speed/Cv per IR substep, split from chemistry.")
        rows.insert(3, "! Geometric collision area is an EXPLICIT comparison input, not derived from Draine optical opacity.")
        rows.append(" gas_exchange_model='hydrogen_accommodation',")
        rows.append(f" collision_area_per_h_input={collision_area:.17e}, accommodation_input={accommodation:.17e},")
    for key, values in arrays.items():
        for start in range(0, len(values), 3):
            rows.append((f" {key}=" if start == 0 else " ") +
                        ",".join(f"{v:.17e}" for v in values[start:start+3]) + ",")
    return "\n".join(rows + ["/", ""])


if __name__ == "__main__":
    raise SystemExit(main())
