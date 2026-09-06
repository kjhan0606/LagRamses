#!/usr/bin/env python3
"""Frozen-primary-heating IR energy transport study; opt-in, float64 CPU."""

import argparse
import json
from pathlib import Path
import sys

import h5py
import jax
import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.dust import read_dust_opacity_metadata, read_dust_thermal_metadata
from snrt_core.dust_ir import build_ir_step, prepare_excess_table
from snrt_core.provenance import sha256_file
from snrt_core.quadrature import level_symmetric_quadrature
from snrt_core.snapshot import read_static_rt_input
from snrt_core.transport import TransportConfig, angular_integral, initial_intensity
from tools.build_draine_dust_opacity import read_draine_table
from tools.build_draine_dust_thermal import _planck_power_density


def planck_band_opacity(source, edges, indices, reference_k):
    """Integral sigma_abs B_E / integral B_E in each complete IR band."""
    if not np.isfinite(reference_k) or reference_k <= 0:
        raise ValueError("reference temperature must be positive")
    data = read_draine_table(source)
    energy = data["energy_ev"]
    result = []
    for group in indices:
        lo, hi = edges[group:group + 2]
        if lo < energy[0] or hi > energy[-1]:
            raise ValueError("source table does not cover IR opacity group")
        grid = np.unique(np.concatenate(([lo, hi], energy[(energy > lo) & (energy < hi)])))
        sigma = np.exp(np.interp(np.log(grid), np.log(energy), np.log(data["absorption_per_h_cm2"])))
        planck = _planck_power_density(grid, np.ones_like(grid), reference_k)
        denominator = np.trapezoid(planck, grid)
        value = np.trapezoid(sigma * planck, grid) / denominator
        if not np.isfinite(value) or value <= 0:
            raise ValueError("nonfinite or zero IR Planck mean")
        result.append(value)
    return np.asarray(result)


def evolve(table, density, primary, absorption, *, width, duration, light_speed,
           order=4, courant=0.4, tolerance=1e-9, max_iterations=128):
    """Evolve a frozen source; return state and compact energy accounting."""
    controls = (width, duration, light_speed, courant, tolerance)
    if not np.isfinite(controls).all() or min(controls) <= 0 or courant > 1:
        raise ValueError("invalid IR evolution controls")
    if (np.asarray(primary).shape != np.asarray(density).shape
            or not np.isfinite(primary).all() or np.any(np.asarray(primary) < 0)):
        raise ValueError("invalid frozen primary heating")
    directions, weights = level_symmetric_quadrature(order, dtype=jnp.float64)
    base_dt = courant * width / (light_speed * np.abs(np.asarray(directions)).sum(axis=1).max())
    steps = max(1, int(np.ceil(duration / base_dt)))
    dt = duration / steps
    config = TransportConfig((width,) * 3, dt, light_speed)
    step = build_ir_step(config, directions, weights, table, density, absorption,
                         tolerance=tolerance, max_iterations=max_iterations)
    energy = initial_intensity(len(absorption), len(directions), np.asarray(density).shape, jnp.float64)
    primary = jnp.asarray(primary)
    escaped = outside = reprocessed = 0.0
    emitted_photons = jnp.zeros_like(jnp.asarray(absorption))
    max_iterations_used = 0
    max_balance = 0.0
    stationarity = 0.0
    for _ in range(steps):
        result = step(energy, primary)
        if not bool(result.valid):
            if bool(result.thermal_invalid):
                reason = "thermal input/range invalid"
            elif not np.isfinite(result.energy).all() or np.any(np.asarray(result.energy) < 0):
                reason = "invalid radiation state"
            else:
                reason = "nonconvergence"
            raise RuntimeError(
                f"IR step rejected ({reason}): iterations={int(result.iterations)} "
                f"local={float(result.local_relative):.6g} balance={float(result.balance_relative):.6g}")
        injected_step = float(primary.sum()) * dt * width**3
        old_total = float(angular_integral(energy, weights).sum()) * width**3
        new_total = float(angular_integral(result.energy, weights).sum()) * width**3
        stationarity = abs(new_total - old_total) / max(injected_step, old_total, np.finfo(float).tiny)
        energy = result.energy
        emitted_photons += result.emitted_photons
        escaped += float(result.escaped_energy)
        outside += float(result.outside_energy.sum()) * width**3
        reprocessed += float(result.absorbed_energy.sum()) * width**3
        max_iterations_used = max(max_iterations_used, int(result.iterations))
        max_balance = max(max_balance, float(result.balance_relative))
    injected = float(primary.sum()) * duration * width**3
    stored = float(angular_integral(energy, weights).sum()) * width**3
    balance = stored + escaped + outside - injected
    relative = abs(balance) / max(injected, np.finfo(float).tiny)
    if relative > tolerance or not np.isfinite(relative):
        raise RuntimeError(f"IR global energy closure failed: {relative:.6g}")
    tau_step = light_speed * dt * np.asarray(absorption)
    source_response = np.ones_like(tau_step)
    np.divide(-np.expm1(-tau_step), tau_step, out=source_response, where=tau_step > 0)
    return {
        "energy_density": np.asarray(angular_integral(energy, weights)),
        "grain_temperature_k": np.asarray(result.temperature),
        "emitted_photons_cm3": np.asarray(emitted_photons),
        "primary_energy_erg": injected, "stored_energy_erg": stored,
        "escaped_energy_erg": escaped, "outside_energy_erg": outside,
        "reprocessed_energy_erg": reprocessed, "balance_residual_erg": balance,
        "balance_relative": relative, "max_step_balance_relative": max_balance,
        "stationarity_relative": stationarity, "max_iterations": max_iterations_used,
        "max_in_step_self_absorption_fraction": float(np.max(1 - source_response)),
        "outside_fraction_of_emitted": outside / max(injected + reprocessed, np.finfo(float).tiny),
        "max_cell_tau": float(np.max(np.asarray(absorption) * width)),
        "dt_s": dt, "steps": steps,
    }


def run(args):
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(output)
    jax.config.update("jax_enable_x64", True)
    static = read_static_rt_input(args.input)
    opacity = read_dust_opacity_metadata(args.dust_opacity_metadata)
    if opacity.schema != "snrt_dust_opacity_v3" or opacity.source_table_sha256 is None:
        raise ValueError("IR study requires pinned v3 dust opacity")
    closure = read_dust_thermal_metadata(
        args.dust_thermal_metadata, expected_group_edges_ev=opacity.group_edges_ev,
        expected_group_edges_sha256=opacity.group_edges_sha256,
        expected_source_table_sha256=opacity.source_table_sha256,
        expected_dust_mass_per_h_g=opacity.dust_mass_per_h_g)
    with h5py.File(args.p5_heating, "r") as handle:
        if handle.attrs.get("format") != "snrt_p5_thermochemical_pilot" or not handle.attrs.get("validation_passed", False):
            raise ValueError("IR heating must come from a validated P5 snapshot")
        for key, path in (("static_input_sha256", args.input),
                          ("dust_opacity_metadata_sha256", args.dust_opacity_metadata),
                          ("dust_thermal_metadata_sha256", args.dust_thermal_metadata)):
            if handle.attrs.get(key) != sha256_file(path):
                raise ValueError(f"P5 heating snapshot binding mismatch: {key}")
        scale_factor = float(handle.attrs["scale_factor"])
        primary = np.asarray(handle["rates/dust_heating_erg_cm3_s"], dtype=float)
    table = prepare_excess_table(closure, 2.7255 / scale_factor)
    source = Path(json.loads(Path(args.dust_thermal_metadata).read_text())["source_table"]["path"])
    density = static.hydrogen_number_density_cm3 * static.dust_relative_abundance
    if primary.shape != density.shape:
        raise ValueError("P5 heating shape does not match static input")
    if not 0 < args.reduced_light_fraction <= 1:
        raise ValueError("invalid reduced light fraction")
    cross_section = planck_band_opacity(source, closure.group_edges_ev,
                                       closure.ir_group_indices, args.opacity_temperature_k)
    absorption = cross_section[:, None, None, None] * density
    result = evolve(table, density, primary, absorption, width=static.grid.cell_width_cm,
                    duration=args.duration_s, light_speed=2.99792458e10 * args.reduced_light_fraction,
                    order=args.sn_order, courant=args.courant, tolerance=args.tolerance,
                    max_iterations=args.max_iterations)
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "x") as handle:
        handle.attrs["format"] = "snrt_dust_ir_transport_study_v1"
        handle.attrs["status"] = "static_gray_candidate"
        handle.attrs["primary_heating_semantics"] = "frozen_final_P5_rate"
        handle.attrs["radiation_semantics"] = "energy_density_excess_above_CMB"
        handle.attrs["outside_escape_semantics"] = (
            "below_IR_lower_edge_plus_above_IR_upper_edge; cold dust dominated by below 0.01 eV; "
            "free complement escape underestimates reabsorption/trapping under the fixed-opacity model")
        handle.attrs["temperature_semantics"] = (
            "diagnostic; zero is inactive sentinel; tiny increments round to CMB temperature; "
            "use differential emission for source rates, not temperature for SED reconstruction")
        handle.attrs["scattering"] = "omitted_IR_absorption_only_study"
        handle.attrs["opacity_temperature_k"] = args.opacity_temperature_k
        handle.attrs["opacity_weighting"] = "integral_sigma_B_over_integral_B"
        handle.attrs["opacity_per_h_cm2"] = cross_section
        handle.attrs["ir_group_indices"] = closure.ir_group_indices
        handle.attrs["group_edges_ev"] = closure.group_edges_ev
        handle.attrs["cmb_temperature_k"] = 2.7255 / scale_factor
        handle.attrs["sn_order"] = args.sn_order
        handle.attrs["reduced_light_fraction"] = args.reduced_light_fraction
        handle.attrs["cell_width_cm"] = static.grid.cell_width_cm
        handle.attrs["tolerance"] = args.tolerance
        for name, path in (("static_input", args.input), ("p5_heating", args.p5_heating),
                           ("opacity", args.dust_opacity_metadata), ("thermal", args.dust_thermal_metadata),
                           ("source_table", source)):
            handle.attrs[name + "_sha256"] = sha256_file(path)
        code_files = (Path(__file__), ROOT / "snrt_core/dust_ir.py", ROOT / "snrt_core/transport.py",
                      ROOT / "snrt_core/dust.py", ROOT / "tools/build_draine_dust_opacity.py",
                      ROOT / "tools/build_draine_dust_thermal.py")
        handle.attrs["code_sha256"] = json.dumps({str(p.resolve()): sha256_file(p) for p in code_files})
        for key, value in result.items():
            if isinstance(value, np.ndarray):
                dataset = handle.create_dataset(key, data=value)
                dataset.attrs["units"] = {"energy_density": "erg cm^-3",
                                          "grain_temperature_k": "K",
                                          "emitted_photons_cm3": "photons cm^-3"}[key]
            else:
                handle.attrs[key] = value
        handle.attrs["validation_passed"] = True
    print(f"DUST_IR_STUDY_OK balance={result['balance_relative']:.3e} "
          f"reprocessed={result['reprocessed_energy_erg']:.6g} steps={result['steps']}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("input", "p5-heating", "dust-opacity-metadata", "dust-thermal-metadata", "output"):
        parser.add_argument("--" + arg, required=True)
    parser.add_argument("--duration-s", required=True, type=float)
    parser.add_argument("--opacity-temperature-k", type=float, default=20)
    parser.add_argument("--reduced-light-fraction", type=float, default=0.01)
    parser.add_argument("--sn-order", type=int, choices=(4, 6, 8), default=4)
    parser.add_argument("--courant", type=float, default=0.4)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--max-iterations", type=int, default=128)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
