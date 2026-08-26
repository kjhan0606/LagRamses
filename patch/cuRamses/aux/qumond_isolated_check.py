#!/usr/bin/env python3
"""Isolated QUMOND two-Poisson benchmark on a periodic uniform mesh.

The source is a narrow spherical Gaussian with its box mean removed.  The
test follows the production QUMOND ordering exactly: solve Newtonian
Poisson, form ``rho_ph=-div[(nu-1) g_N]``, and solve Poisson again.  In the
central region periodic-image corrections are negligible and spherical
QUMOND requires ``g_Q=nu(|g_N|/a0) g_N``.  This exercises the discrete
gradient, divergence, phantom-density sign, and second Poisson solve rather
than merely evaluating the interpolation function.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def solve_poisson(source: np.ndarray, k2: np.ndarray) -> np.ndarray:
    source_k = np.fft.fftn(source)
    potential_k = np.zeros_like(source_k)
    np.divide(-source_k, k2, out=potential_k, where=k2 > 0.0)
    return np.fft.ifftn(potential_k).real


def gradient(field: np.ndarray, dx: float) -> tuple[np.ndarray, ...]:
    return tuple(
        (np.roll(field, -1, axis) - np.roll(field, 1, axis)) / (2.0 * dx)
        for axis in range(3)
    )


def divergence(vector: tuple[np.ndarray, ...], dx: float) -> np.ndarray:
    return sum(
        (np.roll(vector[axis], -1, axis) - np.roll(vector[axis], 1, axis))
        / (2.0 * dx)
        for axis in range(3)
    )


def benchmark(nmesh: int = 64) -> dict[str, float | int | str | bool]:
    dx = 1.0 / nmesh
    phase = 2.0 * np.pi * np.fft.fftfreq(nmesh)
    eigenvalue = 2.0 - 2.0 * np.cos(phase)
    k2 = (
        eigenvalue[:, None, None]
        + eigenvalue[None, :, None]
        + eigenvalue[None, None, :]
    ) / dx**2

    x = (np.arange(nmesh) + 0.5) * dx - 0.5
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    radius = np.sqrt(xx**2 + yy**2 + zz**2)
    sigma = 0.035
    density = np.exp(-radius**2 / (2.0 * sigma**2)) / (2.0 * np.pi * sigma**2) ** 1.5
    density -= density.mean()

    phi_newton = solve_poisson(density, k2)
    force_newton = tuple(-component for component in gradient(phi_newton, dx))
    g_newton = np.sqrt(sum(component**2 for component in force_newton))

    calibration_shell = (radius > 0.08) & (radius < 0.30)
    a0 = float(np.median(g_newton[calibration_shell]))
    nu = 0.5 + 0.5 * np.sqrt(1.0 + 4.0 * a0 / np.maximum(g_newton, 1.0e-30))
    phantom_flux = tuple((nu - 1.0) * component for component in force_newton)
    phantom_density = -divergence(phantom_flux, dx)

    phi_qumond = solve_poisson(density + phantom_density, k2)
    force_qumond = tuple(-component for component in gradient(phi_qumond, dx))
    radial_unit = (xx, yy, zz)
    radial_newton = sum(
        force_newton[axis] * radial_unit[axis] for axis in range(3)
    ) / np.maximum(radius, 1.0e-30)
    radial_qumond = sum(
        force_qumond[axis] * radial_unit[axis] for axis in range(3)
    ) / np.maximum(radius, 1.0e-30)
    target = nu * radial_newton

    comparison_shell = (radius > 0.10) & (radius < 0.22)
    comparison_shell &= np.abs(target) > 0.03 * np.max(np.abs(target[comparison_shell]))
    relative_error = np.abs(radial_qumond[comparison_shell] - target[comparison_shell]) / np.abs(
        target[comparison_shell]
    )
    median_error = float(np.median(relative_error))
    p90_error = float(np.quantile(relative_error, 0.90))
    max_error = float(np.max(relative_error))
    force_ratio = float(np.median(radial_qumond[comparison_shell] / target[comparison_shell]))
    mass_closure = float(abs(phantom_density.mean()))

    passed = (
        median_error < 3.0e-3
        and p90_error < 1.0e-2
        and max_error < 2.0e-2
        and mass_closure < 1.0e-11
    )
    return {
        "model": "QUMOND",
        "test": "periodic spherical two-Poisson isolated solver",
        "nmesh": nmesh,
        "a0_code": a0,
        "median_relative_force_error": median_error,
        "p90_relative_force_error": p90_error,
        "max_relative_force_error": max_error,
        "median_measured_to_spherical_target": force_ratio,
        "phantom_density_mean_abs": mass_closure,
        "passed": passed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nmesh", type=int, default=64)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = benchmark(args.nmesh)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(payload, end="")
    if args.json is not None:
        args.json.write_text(payload)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
