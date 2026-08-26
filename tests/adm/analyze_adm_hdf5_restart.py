#!/usr/bin/env python3
"""Compare continuous and HDF5-restarted ADM particle states."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


PARTICLE_FIELDS = (
    "x_1",
    "x_2",
    "x_3",
    "v_1",
    "v_2",
    "v_3",
    "mass",
    "levelp",
    "ptypep",
    "dark_energy_int",
    "dark_h2_frac",
)
EXPECTED_PARTICLES = 16**3 * 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    return parser.parse_args()


def read_particles(path: Path) -> dict[str, np.ndarray]:
    with h5py.File(path, "r") as handle:
        particles = handle["particles"]
        identity = np.asarray(particles["identity"])
        order = np.argsort(identity, kind="stable")
        state = {"identity": identity[order]}
        for field in PARTICLE_FIELDS:
            state[field] = np.asarray(particles[field])[order]
    return state


def main() -> None:
    root = parse_args().root
    continuous_log = (root / "continuous" / "run.log").read_text(encoding="utf-8")
    restart_log = (root / "restart" / "restart.log").read_text(encoding="utf-8")
    if "HDF5 PART backup files read completed" not in restart_log:
        raise RuntimeError("HDF5 particle restore was not reached")
    if "ADM new-run temperature initialized:" in restart_log:
        raise RuntimeError("restart incorrectly reinitialised ADM internal energy")
    if "ADM new-run temperature initialized:" not in continuous_log:
        raise RuntimeError("fresh ADM initialisation message missing")

    continuous = read_particles(
        root / "continuous" / "output_00004" / "data_00004.h5"
    )
    restarted = read_particles(root / "restart" / "output_00004" / "data_00004.h5")
    if continuous["identity"].size != EXPECTED_PARTICLES:
        raise RuntimeError("unexpected continuous particle count")
    if not np.array_equal(continuous["identity"], restarted["identity"]):
        raise RuntimeError("particle identities differ after HDF5 restart")

    max_relative_difference = 0.0
    for field in PARTICLE_FIELDS:
        reference = continuous[field]
        candidate = restarted[field]
        if field in ("levelp", "ptypep"):
            if not np.array_equal(reference, candidate):
                raise RuntimeError(f"{field} differs after HDF5 restart")
            continue
        # The ADM fields are copied through HDF5 exactly.  Position and
        # velocity differ only through the non-bitwise force reconstruction
        # after an otherwise identical HDF5 restart.
        if field in ("dark_energy_int", "dark_h2_frac", "mass"):
            if not np.array_equal(candidate, reference):
                raise RuntimeError(f"{field} is not bitwise preserved by HDF5 restart")
        else:
            np.testing.assert_allclose(candidate, reference, rtol=5.0e-7, atol=2.0e-10)
        scale = max(float(np.max(np.abs(reference))), 1.0e-300)
        max_relative_difference = max(
            max_relative_difference,
            float(np.max(np.abs(candidate - reference)) / scale),
        )

    energy = continuous["dark_energy_int"]
    print(f"particles={energy.size}")
    print(f"dark_energy_min={energy.min():.16e}")
    print(f"dark_energy_max={energy.max():.16e}")
    print(f"max_relative_difference={max_relative_difference:.3e}")
    print("ADM_HDF5_RESTART_RESULT=continuous-and-restarted-states-match")


if __name__ == "__main__":
    main()
