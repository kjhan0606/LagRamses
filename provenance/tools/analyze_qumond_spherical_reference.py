#!/usr/bin/env python3
"""Compare lagRamses force snapshots with its uniform-grid QUMOND operators."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np


SOFTENING = 0.05
CENTER = 0.5
UNITS_LENGTH_CM = 3.0856775814913673e22
UNITS_TIME_S = 3.8707679655184575e15
UNITS_DENSITY_CGS = 1.0e-24
G_CGS = 6.67430e-8
MSUN_G = 1.98847e33
A0_CGS = 1.2e-8
FOURPI = 4.0 * math.pi


def read_uniform_snapshot(path: Path) -> tuple[np.ndarray, np.ndarray, float, int]:
    with h5py.File(path, "r") as handle:
        header = handle["header"].attrs
        levelmin = int(np.asarray(header["levelmin"]).item())
        levelmax = int(np.asarray(header["nlevelmax"]).item())
        boxlen = float(np.asarray(header["boxlen"]).item())
        if levelmin != levelmax:
            raise ValueError(f"expected a uniform grid, got levels {levelmin}:{levelmax}")
        level = levelmin
        nmesh = 2**level
        dx = boxlen / nmesh
        amr = handle[f"amr/level_{level}"]
        gravity = handle[f"gravity/level_{level}"]
        xg = np.column_stack([np.asarray(amr[f"xg_{axis}"]) for axis in (1, 2, 3)])
        ngrid = len(xg)
        if 8 * ngrid != nmesh**3:
            raise ValueError(f"grid coverage has {8 * ngrid} cells, expected {nmesh**3}")
        son = np.asarray(amr["son_flag"])
        if np.any(son != 0):
            raise ValueError("reference run unexpectedly contains refined cells")

        child = np.arange(8)
        child_xyz = np.column_stack((child % 2, (child // 2) % 2, child // 4))
        positions = xg[:, None, :] + (child_xyz[None, :, :] - 0.5) * dx
        positions %= boxlen
        ijk = np.rint(positions / dx - 0.5).astype(np.int64) % nmesh
        flat_indices = np.ravel_multi_index(ijk.reshape(-1, 3).T, (nmesh,) * 3)
        counts = np.bincount(flat_indices, minlength=nmesh**3)
        if not np.all(counts == 1):
            raise ValueError("AMR coordinates do not map one-to-one onto the uniform grid")

        forces = np.empty((3, nmesh, nmesh, nmesh), dtype=np.float64)
        coordinate_grid = np.empty((nmesh, nmesh, nmesh, 3), dtype=np.float64)
        coordinate_grid[ijk[..., 0], ijk[..., 1], ijk[..., 2]] = positions
        for axis in range(3):
            values = np.asarray(gravity[f"f_{axis + 1}"]).reshape(ngrid, 8)
            grid = np.empty((nmesh, nmesh, nmesh), dtype=np.float64)
            grid[ijk[..., 0], ijk[..., 1], ijk[..., 2]] = values
            forces[axis] = grid

    if not np.all(np.isfinite(forces)):
        raise ValueError(f"non-finite force values in {path}")
    return forces, coordinate_grid, boxlen, nmesh


def source_density(radius: np.ndarray) -> np.ndarray:
    density_cap = 1.0 / (SOFTENING * (1.0 + SOFTENING) ** 2)
    safe_radius = np.maximum(radius, np.finfo(np.float64).tiny)
    return np.minimum(1.0 / (safe_radius * (1.0 + safe_radius) ** 2), density_cap)


def periodic_poisson(density: np.ndarray, dx: float) -> np.ndarray:
    nmesh = density.shape[0]
    ix = np.arange(nmesh)
    iz = np.arange(nmesh // 2 + 1)
    eigenvalue = 2.0 * (
        np.cos(2.0 * np.pi * ix[:, None, None] / nmesh)
        + np.cos(2.0 * np.pi * ix[None, :, None] / nmesh)
        + np.cos(2.0 * np.pi * iz[None, None, :] / nmesh)
        - 3.0
    )
    green = np.zeros_like(eigenvalue)
    nonzero = eigenvalue != 0.0
    green[nonzero] = dx**2 / eigenvalue[nonzero]
    rhs_k = np.fft.rfftn(FOURPI * (density - np.mean(density)))
    phi_k = rhs_k * green
    return np.fft.irfftn(phi_k, s=density.shape).real


def fourth_order_force(phi: np.ndarray, dx: float) -> np.ndarray:
    force = np.empty((3, *phi.shape), dtype=np.float64)
    for axis in range(3):
        near_left = np.roll(phi, 1, axis=axis)
        near_right = np.roll(phi, -1, axis=axis)
        far_left = np.roll(phi, 2, axis=axis)
        far_right = np.roll(phi, -2, axis=axis)
        force[axis] = (2.0 / (3.0 * dx)) * (near_left - near_right) - (
            1.0 / (12.0 * dx)
        ) * (far_left - far_right)
    return force


def central_divergence(vector: np.ndarray, dx: float) -> np.ndarray:
    result = np.zeros(vector.shape[1:], dtype=np.float64)
    for axis in range(3):
        result += (np.roll(vector[axis], -1, axis=axis) - np.roll(vector[axis], 1, axis=axis)) / (
            2.0 * dx
        )
    return result


def discrete_reference(density: np.ndarray, dx: float, a0_code: float) -> tuple[np.ndarray, np.ndarray]:
    newtonian = fourth_order_force(periodic_poisson(density, dx), dx)
    gnorm = np.sqrt(np.sum(newtonian**2, axis=0))
    x = np.maximum(gnorm / a0_code, 1.0e-30)
    nu_minus_one = 0.5 * (np.sqrt(1.0 + 4.0 / x) - 1.0)
    phantom_density = -central_divergence(nu_minus_one[None, ...] * newtonian, dx) / FOURPI
    qumond = fourth_order_force(periodic_poisson(density + phantom_density, dx), dx)
    return newtonian, qumond


def enclosed_source_mass(radius: np.ndarray) -> np.ndarray:
    rho_core = source_density(np.asarray(SOFTENING))
    core_mass = FOURPI * rho_core * SOFTENING**3 / 3.0
    f_r = np.log1p(radius) + 1.0 / (1.0 + radius)
    f_core = math.log1p(SOFTENING) + 1.0 / (1.0 + SOFTENING)
    mass_outer = core_mass + FOURPI * (f_r - f_core)
    mass_core = FOURPI * rho_core * radius**3 / 3.0
    return np.where(radius <= SOFTENING, mass_core, mass_outer)


def error_statistics(measured: np.ndarray, reference: np.ndarray, radius: np.ndarray) -> dict[str, float]:
    mask = (radius > 0.06) & (radius < 0.35)
    denominator = np.sqrt(np.sum(reference**2, axis=0))
    difference = np.sqrt(np.sum((measured - reference) ** 2, axis=0))
    valid = mask & (denominator > np.max(denominator[mask]) * 1.0e-4)
    relative = difference[valid] / denominator[valid]
    return {
        "cells": int(relative.size),
        "median_relative_vector_force_error": float(np.median(relative)),
        "p90_relative_vector_force_error": float(np.quantile(relative, 0.90)),
        "max_relative_vector_force_error": float(np.max(relative)),
    }


def radial_acceleration(
    force: np.ndarray, radius: np.ndarray, coordinates: np.ndarray, boxlen: float
) -> np.ndarray:
    half_box = 0.5 * boxlen
    displacement = coordinates - half_box
    displacement = (displacement + half_box) % boxlen - half_box
    unit_radius = np.zeros_like(displacement)
    np.divide(displacement, radius[..., None], out=unit_radius, where=radius[..., None] > 0)
    return -np.sum(force * np.moveaxis(unit_radius, -1, 0), axis=0)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: analyze_qumond_spherical_reference.py RUN_DIRECTORY")
    root = Path(sys.argv[1]).resolve()
    cases = {}
    positions = None
    boxlen = None
    nmesh = None
    for model in ("newtonian", "qumond"):
        case = root / model
        snapshots = sorted(case.glob("output_[0-9][0-9][0-9][0-9][0-9]/data_*.h5"))
        if not snapshots:
            raise ValueError(f"{model} has no HDF5 output")
        final_force = None
        model_positions = None
        model_boxlen = None
        model_nmesh = None
        snapshot_records = []
        valid_force_snapshots = []
        for snapshot in snapshots:
            force, output_positions, output_boxlen, output_nmesh = read_uniform_snapshot(snapshot)
            if model_positions is None:
                model_positions, model_boxlen, model_nmesh = output_positions, output_boxlen, output_nmesh
            elif not np.array_equal(model_positions, output_positions):
                raise ValueError(f"{model} output meshes differ across snapshots")
            max_abs_force = float(np.max(np.abs(force)))
            has_force = max_abs_force > 0.0
            snapshot_record = {
                "path": str(snapshot),
                "bytes": snapshot.stat().st_size,
                "max_abs_force_code": max_abs_force,
                "has_nonzero_force": has_force,
            }
            snapshot_records.append(snapshot_record)
            if has_force:
                valid_force_snapshots.append((snapshot, force))
                final_force = force
        if model == "newtonian":
            positions, boxlen, nmesh = model_positions, model_boxlen, model_nmesh
        elif not np.array_equal(positions, model_positions):
            raise ValueError("Newtonian and QUMOND output meshes do not match")
        if model_positions is None or model_boxlen is None:
            raise ValueError(f"{model} output is incomplete")
        if final_force is None:
            raise ValueError(
                f"{model} has no nonzero force snapshot; check that the analytic source was deposited "
                "and that the selected output follows the Poisson solve"
            )
        model_displacement = model_positions - 0.5 * model_boxlen
        model_displacement = (model_displacement + 0.5 * model_boxlen) % model_boxlen - 0.5 * model_boxlen
        model_radius = np.sqrt(np.sum(model_displacement**2, axis=-1))
        cases[model] = {
            "force": final_force,
            "snapshot_records": snapshot_records,
            "snapshot_count": len(snapshots),
            "force_snapshot_count": len(valid_force_snapshots),
            "zero_force_snapshot_count": len(snapshots) - len(valid_force_snapshots),
            "force_snapshot_paths": [str(item[0]) for item in valid_force_snapshots],
        }

    assert positions is not None and boxlen is not None and nmesh is not None
    dx = boxlen / nmesh
    displacement = positions - CENTER
    displacement = (displacement + 0.5 * boxlen) % boxlen - 0.5 * boxlen
    radius = np.sqrt(np.sum(displacement**2, axis=-1))
    density = source_density(radius)
    acceleration_unit = UNITS_LENGTH_CM / UNITS_TIME_S**2
    a0_code = A0_CGS / acceleration_unit
    a0_code_from_g_units = G_CGS * UNITS_DENSITY_CGS * UNITS_TIME_S**2
    if not math.isclose(a0_code_from_g_units, 1.0, rel_tol=2.0e-5):
        raise ValueError(f"chosen density/time units give G_code={a0_code_from_g_units}")

    reference_newtonian, reference_qumond = discrete_reference(density, dx, a0_code)
    metrics = {
        "newtonian_discrete_reference": error_statistics(cases["newtonian"]["force"], reference_newtonian, radius),
        "qumond_discrete_reference": error_statistics(cases["qumond"]["force"], reference_qumond, radius),
    }

    gbar_isolated = enclosed_source_mass(np.maximum(radius, 1.0e-30)) / np.maximum(radius, 1.0e-30) ** 2
    x_iso = np.maximum(gbar_isolated / a0_code, 1.0e-30)
    nu_iso = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 / x_iso))
    gqumond_isolated = nu_iso * gbar_isolated
    radius_kpc = radius * boxlen * UNITS_LENGTH_CM / 3.0856775814913673e21
    gbar_cgs = gbar_isolated * acceleration_unit
    gqumond_cgs = radial_acceleration(cases["qumond"]["force"], radius, positions, boxlen) * acceleration_unit
    gnewton_cgs = radial_acceleration(cases["newtonian"]["force"], radius, positions, boxlen) * acceleration_unit
    gqumond_isolated_cgs = gqumond_isolated * acceleration_unit
    empirical = gbar_cgs / (-np.expm1(-np.sqrt(gbar_cgs / A0_CGS)))
    shell_rows = []
    edges = np.arange(0.06, 0.361, 0.02)
    for lower, upper in zip(edges[:-1], edges[1:]):
        selected = (radius > lower) & (radius <= upper)
        if np.count_nonzero(selected) < 10:
            continue
        median_gbar = float(np.median(gbar_cgs[selected]))
        median_gq = float(np.median(gqumond_cgs[selected]))
        median_empirical = float(np.median(empirical[selected]))
        shell_rows.append({
            "radius_kpc": float(np.median(radius_kpc[selected])),
            "cell_count": int(np.count_nonzero(selected)),
            "gbar_cgs": median_gbar,
            "g_newtonian_code_cgs": float(np.median(gnewton_cgs[selected])),
            "g_qumond_code_cgs": median_gq,
            "g_qumond_isolated_spherical_cgs": float(np.median(gqumond_isolated_cgs[selected])),
            "g_empirical_rar_cgs": median_empirical,
            "delta_log10_qumond_to_empirical_rar_dex": float(np.log10(median_gq / median_empirical)),
        })

    with (root / "rar_profile.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(shell_rows[0]))
        writer.writeheader()
        writer.writerows(shell_rows)

    summary = {
        "status": "ANALYZED",
        "scope": "compiled uniform-grid spherical source benchmark; not a disc simulation or a Lambda CDM versus QUMOND galaxy comparison",
        "source": {
            "analytic_profile": "rho=min[1/(r*(1+r)^2), 1/(epsilon*(1+epsilon)^2)]",
            "epsilon_code": SOFTENING,
            "box_length_code": boxlen,
            "box_length_kpc": boxlen * UNITS_LENGTH_CM / 3.0856775814913673e21,
            "grid": f"{nmesh}^3 uniform, level {LEVEL}",
            "gas_and_particle_dynamics": "disabled; fixed analytic density source",
            "number_of_snapshots": {model: item["snapshot_count"] for model, item in cases.items()},
        },
        "units": {
            "mass_unit_msun": UNITS_DENSITY_CGS * UNITS_LENGTH_CM**3 / MSUN_G,
            "acceleration_unit_cgs": acceleration_unit,
            "G_code": a0_code_from_g_units,
            "a0_cgs": A0_CGS,
            "a0_code": a0_code,
        },
        "execution": {
            model: {
                "snapshot_count": item["snapshot_count"],
                "snapshots": item["snapshot_records"],
                "force_snapshot_count": item["force_snapshot_count"],
                "zero_force_snapshot_count": item["zero_force_snapshot_count"],
                "force_snapshot_paths": item["force_snapshot_paths"],
            }
            for model, item in cases.items()
        },
        "discrete_operator_agreement": metrics,
        "radial_bins": shell_rows,
        "interpretation_limit": "The isolated spherical relation is a reference in the central region. Periodic-image and finite-grid errors require a box-size and resolution comparison before treating the RAR residual as a physical result.",
    }
    (root / "analysis.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
