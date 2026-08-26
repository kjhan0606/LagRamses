#!/usr/bin/env python3
"""Report the small-box response to physically labelled ADM initial temperatures."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import h5py
import numpy as np


K_B_CGS = 1.380649e-16
GEV_TO_G = 1.78266192e-24
DARK_PROTON_MASS_GEV = 40.0
SCALE_V = 3.08567758e21 / 5.0e9
LEAF_COUNT = 16**3
MULTIPLICITY = 4
INITIALISATION = re.compile(
    r"ADM new-run temperature initialized:\s*T_D=\s*"
    r"(?P<temperature>[+-]?\d\.\d+E[+-]\d+)\s*K, edp="
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
    manifest = json.loads((root / "thermal_history.json").read_text(encoding="utf-8"))
    results: list[dict[str, float | str]] = []
    print("label z_kd T_init[K] T_after_two_steps[K] fractional_drop")

    for candidate in manifest["candidates"]:
        if not candidate["above_code_floor"]:
            raise RuntimeError(f"{candidate['label']}: temperature is below the 1 K code floor")
        label = str(candidate["label"])
        expected_temperature = float(candidate["temperature_k"])
        case = root / label
        log = (case / "run.log").read_text(encoding="utf-8")
        match = INITIALISATION.search(log)
        if match is None:
            raise RuntimeError(f"{label}: ADM initialisation message missing")
        if "Entering dark_cooling_fine for level" not in log:
            raise RuntimeError(f"{label}: dark-cooling AMR path was not entered")
        initialized_temperature = float(match.group("temperature"))
        if not np.isclose(initialized_temperature, expected_temperature, rtol=5.0e-5):
            raise RuntimeError(
                f"{label}: initialized at {initialized_temperature} K, expected {expected_temperature} K"
            )

        outputs = sorted(case.glob("output_*/data_*.h5"))
        if not outputs:
            raise RuntimeError(f"{label}: no HDF5 particle output was written")
        output = outputs[-1]
        with h5py.File(output, "r") as handle:
            particles = handle["particles"]
            energy = np.asarray(particles["dark_energy_int"], dtype=float)
            mass = np.asarray(particles["mass"], dtype=float)
        if energy.size != LEAF_COUNT * MULTIPLICITY:
            raise RuntimeError(f"{label}: unexpected particle count {energy.size}")
        if not np.all(np.isfinite(energy)):
            raise RuntimeError(f"{label}: non-finite ADM internal energy")
        if not np.isclose(mass.sum(), LEAF_COUNT * 1.0e-6, rtol=2.0e-14):
            raise RuntimeError(f"{label}: total macro-particle mass changed")

        temperature = temperature_from_energy(energy)
        final_temperature = float(np.average(temperature, weights=mass))
        if final_temperature < 1.0 or final_temperature > expected_temperature * (1.0 + 1.0e-10):
            raise RuntimeError(f"{label}: cooling-only update left the physical temperature range")
        fractional_drop = 1.0 - final_temperature / expected_temperature
        results.append(
            {
                "label": label,
                "z_kd": float(candidate["z_kd"]),
                "initial_temperature_k": expected_temperature,
                "final_temperature_k": final_temperature,
                "fractional_drop": fractional_drop,
                "particle_output": output.parent.name,
            }
        )
        print(
            f"{label} {float(candidate['z_kd']):.8g} {expected_temperature:.8e} "
            f"{final_temperature:.8e} {fractional_drop:.8e}"
        )

    (root / "temperature_sensitivity.json").write_text(
        json.dumps({"results": results}, indent=2) + "\n", encoding="utf-8"
    )
    print("ADM_TEMPERATURE_SENSITIVITY_RESULT=all-candidates-completed")


if __name__ == "__main__":
    main()
