#!/usr/bin/env python3
"""Check that the HPM pressure branch creates a finite, directional kick."""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np


def latest_particles(case: Path) -> tuple[np.ndarray, np.ndarray]:
    outputs = sorted(case.glob("output_*/data_*.h5"))
    if not outputs:
        raise RuntimeError(f"{case}: no HDF5 output")
    with h5py.File(outputs[-1], "r") as handle:
        particles = handle["particles"]
        velocity = np.column_stack(
            [np.asarray(particles[f"v_{axis}"], dtype=float) for axis in (1, 2, 3)]
        )
        return velocity, np.asarray(particles["mass"], dtype=float)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: analyze_adm_hpm_gradient.py ROOT")
    root = Path(sys.argv[1])
    baseline_velocity, mass = latest_particles(root / "baseline")
    hpm_velocity, hpm_mass = latest_particles(root / "hpm")
    if not np.array_equal(mass, hpm_mass):
        raise RuntimeError("HPM altered macro-particle masses")
    delta = hpm_velocity - baseline_velocity
    if not np.all(np.isfinite(delta)):
        raise RuntimeError("HPM produced non-finite velocity")
    if np.max(np.abs(delta[:, 0])) <= 1.0e-12:
        raise RuntimeError("HPM pressure branch produced no x-direction kick")
    if np.max(np.abs(delta[:, 1:])) > 1.0e-11:
        raise RuntimeError("one-dimensional pressure test produced transverse acceleration")
    momentum_delta = np.sum(mass[:, None] * delta, axis=0)
    scale = np.sum(mass) * np.max(np.abs(delta[:, 0]))
    if abs(momentum_delta[0]) > 5.0e-8 * scale:
        raise RuntimeError(f"HPM pressure kick is not momentum balanced: {momentum_delta[0]:.3e}")
    print("ADM_HPM_SMOKE_RESULT=pressure-gradient-kick-passed")
    print(f"max_delta_vx={np.max(np.abs(delta[:, 0])):.8e}")
    print(f"momentum_delta_x={momentum_delta[0]:.8e}")


if __name__ == "__main__":
    main()
