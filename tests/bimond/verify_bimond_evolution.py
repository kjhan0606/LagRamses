#!/usr/bin/env python3
"""Check BIMOND kick/drift identities from successive native HDF5 outputs."""

import json
import sys
from pathlib import Path

import h5py
import numpy as np

from verify_bimond_two_sector import assess, interpolate


def evolution(directory):
    initial = np.loadtxt(directory / "ic_part")
    previous_x = initial[:, :3] + 0.5
    previous_v = initial[:, 3:6]
    previous_dt = 0.0
    errors = []
    files = sorted(directory.glob("output_*/data_*.h5"))
    assert len(files) == 4
    for filename in files:
        with h5py.File(filename) as handle:
            header = handle["header"].attrs
            level = int(np.asarray(header["levelmin"]).item())
            n = 2**level
            dt = float(np.asarray(header["dtnew"])[level - 1])
            old_dt = float(np.asarray(header["dtold"])[level - 1])
            time = float(np.asarray(header["time"]).item())
            assert np.isclose(old_dt, previous_dt, rtol=1e-12, atol=1e-30)
            particle = handle["particles"]
            order = np.argsort(particle["identity"][:])
            identity = particle["identity"][:][order]
            assert np.array_equal(identity, np.arange(1, len(initial) + 1))
            position = np.column_stack([particle[f"x_{i}"][:][order] for i in (1, 2, 3)])
            velocity = np.column_stack([particle[f"v_{i}"][:][order] for i in (1, 2, 3)])
            mass = particle["mass"][:][order]
            types = particle["ptypep"][:][order]
            assert np.array_equal(types, initial[:, 7].astype(int))
            grid = handle[f"amr/level_{level}"]
            gravity = handle[f"gravity/level_{level}"]
            xg = np.column_stack([grid[f"xg_{i}"][:] for i in (1, 2, 3)])
            child = np.array([[c % 2, c // 2 % 2, c // 4] for c in range(8)])
            index = np.rint(xg[:, None, :] * n + child - 1).astype(int) % n
            acceleration = np.zeros_like(position)
            for tag, prefix in [(1, "f_"), (30, "bimond_f_twin_")]:
                fields = np.empty((3, n, n, n))
                for axis in range(3):
                    fields[axis, index[..., 0], index[..., 1], index[..., 2]] = gravity[f"{prefix}{axis + 1}"][:].reshape(-1, 8)
                mask = types == tag
                acceleration[mask] = interpolate(previous_x[mask], fields)
        expected_v = previous_v + 0.5 * (old_dt + dt) * acceleration
        expected_x = (previous_x + dt * expected_v) % 1
        velocity_error = np.sqrt(np.sum(mass[:, None] * (velocity - expected_v)**2) / np.sum(mass[:, None] * expected_v**2))
        position_error = np.max(np.abs((position - expected_x + 0.5) % 1 - 0.5))
        assert velocity_error < 1e-10
        assert position_error < 1e-12
        signed_p = np.sum(np.where(types == 1, 1, -1)[:, None] * mass[:, None] * velocity, axis=0)
        momentum_scale = np.sum(mass * np.linalg.norm(velocity, axis=1))
        momentum_error = np.linalg.norm(signed_p) / momentum_scale
        assert momentum_error < 1e-2
        if directory.name.startswith("equal"):
            assert np.array_equal(position[0::2], position[1::2])
            assert np.array_equal(velocity[0::2], velocity[1::2])
        errors.append({"time": time,
                       "dt": dt, "velocity_mass_weighted_rms": float(velocity_error),
                       "position_max": float(position_error), "signed_momentum_relative": float(momentum_error)})
        previous_x, previous_v, previous_dt = position, velocity, dt
    return errors


def main():
    root = Path(sys.argv[1])
    results = {"equal_evolution": evolution(root / "equal16"),
               "pair_evolution": evolution(root / "pair16")}
    if len(sys.argv) > 2:
        reference_root = Path(sys.argv[2])
        results.update(shift64=assess(root / "shift64"), mpi2_pair32=assess(root / "mpi2_pair32"))
        previous = json.loads((reference_root / "summary.json").read_text())
        results["shifted_signed_force_decreases"] = results["shift64"]["signed_force_relative"] < previous["shift32"]["signed_force_relative"]
        assert results["shifted_signed_force_decreases"]
        assert results["shift64"]["signed_force_relative"] < 1e-2
    results["passed"] = True
    (root / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
