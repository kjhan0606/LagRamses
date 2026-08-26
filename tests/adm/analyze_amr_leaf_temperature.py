#!/usr/bin/env python3
"""Validate the ADM AMR leaf-density and new-run temperature smoke test."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py
import numpy as np


K_B_CGS = 1.380649e-16
GEV_TO_G = 1.78266192e-24
LEAF_COUNT = 16**3
INITIAL_TEMPERATURE = 1.0e3
DARK_PROTON_MASS_GEV = 40.0
SCALE_V = 3.08567758e21 / 5.0e9
INITIALISATION = re.compile(
    r"ADM new-run temperature initialized:\s*T_D=\s*"
    r"(?P<temperature>[+-]?\d\.\d+E[+-]\d+)\s*K, edp="
    r"\s*(?P<energy>[+-]?\d\.\d+E[+-]\d+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    return parser.parse_args()


def temperature_from_energy(energy: np.ndarray) -> np.ndarray:
    return (
        (2.0 / 3.0) * energy * SCALE_V**2
        * DARK_PROTON_MASS_GEV * GEV_TO_G / K_B_CGS
    )


def main() -> None:
    root = parse_args().root
    print("multiplicity particles mass_sum T_mean[K] T_span[K]")
    post_step_temperatures: list[float] = []
    for multiplicity in (1, 2, 4, 8):
        case = root / f"multiplicity_{multiplicity}"
        log = (case / "run.log").read_text(encoding="utf-8")
        match = INITIALISATION.search(log)
        if match is None:
            raise RuntimeError(f"{case}: ADM initialisation message missing")
        if "Entering dark_cooling_fine for level" not in log:
            raise RuntimeError(f"{case}: dark-cooling AMR path was not entered")

        initial_temperature = float(match.group("temperature"))
        initial_energy = float(match.group("energy"))
        expected_energy = (
            1.5 * K_B_CGS * INITIAL_TEMPERATURE
            / (DARK_PROTON_MASS_GEV * GEV_TO_G * SCALE_V**2)
        )
        if not np.isclose(initial_temperature, INITIAL_TEMPERATURE, rtol=2.0e-12):
            raise RuntimeError(
                f"{case}: initial temperature is {initial_temperature} K, not "
                f"{INITIAL_TEMPERATURE} K"
            )
        if not np.isclose(initial_energy, expected_energy, rtol=3.0e-5):
            raise RuntimeError(f"{case}: initial specific energy has the wrong unit conversion")

        # foutput=1 writes output_00001 at the initial state and output_00002
        # after mesh refinement.  With nstepmax=2, output_00003 follows a
        # fully refined leaf update and contains the cooling response.
        output = case / "output_00003" / "data_00003.h5"
        with h5py.File(output, "r") as handle:
            particles = handle["particles"]
            energy = np.asarray(particles["dark_energy_int"], dtype=float)
            mass = np.asarray(particles["mass"], dtype=float)

        expected_particles = LEAF_COUNT * multiplicity
        if energy.size != expected_particles:
            raise RuntimeError(
                f"{case}: expected {expected_particles} particles, got {energy.size}"
            )
        if not np.all(np.isfinite(energy)):
            raise RuntimeError(f"{case}: non-finite ADM internal energy")
        if not np.isclose(mass.sum(), LEAF_COUNT * 1.0e-6, rtol=2.0e-14):
            raise RuntimeError(f"{case}: total macro-particle mass changed")

        temperature = temperature_from_energy(energy)
        mean_temperature = float(np.average(temperature, weights=mass))
        span_temperature = float(temperature.max() - temperature.min())
        if mean_temperature < 1.0 or mean_temperature > INITIAL_TEMPERATURE:
            raise RuntimeError(
                f"{case}: post-step temperature {mean_temperature} K is outside [1, 1000] K"
            )
        if mean_temperature >= 0.98 * INITIAL_TEMPERATURE:
            raise RuntimeError(f"{case}: fully refined cooling response was not resolved")
        post_step_temperatures.append(mean_temperature)
        print(
            f"{multiplicity:12d} {energy.size:9d} {mass.sum():.8e} "
            f"{mean_temperature:.8f} {span_temperature:.3e}"
        )
    if not np.allclose(post_step_temperatures, post_step_temperatures[0], rtol=1.0e-4):
        raise RuntimeError(
            "post-step ADM temperature depends on the number of macro-particles per leaf"
        )
    print("ADM_AMR_SMOKE_RESULT=initial-temperature-and-leaf-occupancy-passed")


if __name__ == "__main__":
    main()
