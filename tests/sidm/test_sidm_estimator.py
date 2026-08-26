#!/usr/bin/env python3
"""Regression and resolution diagnostics for the production SIDM estimator."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "patch" / "cuRamses" / "sidm_estimator.f90"
DRIVER = Path(__file__).with_suffix(".f90")


def compile_and_run() -> np.ndarray:
    compiler = shutil.which("gfortran")
    if compiler is None:
        raise SystemExit("gfortran is required for the SIDM estimator regression test")
    with tempfile.TemporaryDirectory(prefix="sidm-estimator-") as temporary:
        executable = Path(temporary) / "test_sidm_estimator"
        subprocess.run(
            [compiler, "-O2", "-std=f2008", "-Wall", "-Wextra", str(SOURCE),
             str(DRIVER), "-o", str(executable)],
            check=True,
        )
        result = subprocess.run(
            [str(executable)], check=True, text=True, capture_output=True
        )
    return np.fromstring(result.stdout, sep=" ").reshape(-1, 4)


def partition_ratio(total_particles: int, cells: int) -> float:
    if total_particles % cells != 0:
        raise ValueError("This diagnostic requires an equal cell occupancy")
    return (total_particles-cells)/(total_particles-1)


def main() -> None:
    rows = compile_and_run()
    probability = rows[:, 2]
    expected = rows[:, 3]
    if not np.allclose(probability, expected, rtol=2.0e-15, atol=0.0):
        raise AssertionError("Production pair probability differs from its exact estimator")
    print(f"pair-probability regression: {rows.shape[0]} occupancies passed")

    for occupancy in (2, 4, 8, 16, 32, 64, 128):
        cells = 64
        total_particles = occupancy*cells
        ratio = partition_ratio(total_particles, cells)
        print(
            f"equal-cell occupancy={occupancy:3d} "
            f"resolved/global-rate={ratio:.6f} "
            f"finite-cell deficit={1.0-ratio:.2%}"
        )

    for p_max in (0.1, 0.05, 0.02, 0.01):
        linearised_bias = 1.0-(1.0-math.exp(-p_max))/p_max
        print(f"Pmax={p_max:.3f} one-event linearisation bias={linearised_bias:.2%}")


if __name__ == "__main__":
    main()
