#!/usr/bin/env python3
"""Prepare and assess native BIMOND static particle/field checks on Slurm."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from analyze_qumond_spherical_reference import (
    A0_CGS, FOURPI, UNITS_LENGTH_CM, UNITS_TIME_S,
    central_divergence, fourth_order_force, periodic_poisson,
)


def deposit(position, mass, n):
    scaled = position * n - 0.5
    left = np.floor(scaled).astype(int)
    right_weight = scaled - left
    density = np.zeros((n, n, n))
    for child in range(8):
        offset = np.array([child % 2, child // 2 % 2, child // 4])
        index = (left + offset) % n
        weight = np.prod(np.where(offset, right_weight, 1 - right_weight), axis=1)
        np.add.at(density, tuple(index.T), mass * weight * n**3)
    return density


def interpolate(position, field):
    n = field.shape[1]
    scaled = position * n - 0.5
    left = np.floor(scaled).astype(int)
    right_weight = scaled - left
    result = np.zeros_like(position)
    for child in range(8):
        offset = np.array([child % 2, child // 2 % 2, child // 4])
        index = (left + offset) % n
        weight = np.prod(np.where(offset, right_weight, 1 - right_weight), axis=1)
        result += field[:, index[:, 0], index[:, 1], index[:, 2]].T * weight[:, None]
    return result


def prepare(root, template, cases=None):
    if cases is None:
        cases = [("equal32", 5, False, 0), ("pair32", 5, True, 0),
                 ("pair64", 6, True, 0), ("shift32", 5, True, 0.27)]
    for name, level, displaced, shift in cases:
        directory = root / name
        directory.mkdir()
        n = 2**level
        grid = np.moveaxis(np.indices((n,) * 3), 0, -1).reshape(-1, 3)
        position = (grid + 0.5) / n
        masses = []
        for center_x in [0.4, 0.6 if displaced else 0.4]:
            delta = position - [center_x, 0.5, 0.5]
            profile = np.exp(np.sum(np.cos(2 * np.pi * delta) - 1, axis=1) / (2 * np.pi * 0.08)**2)
            masses.append(1e-3 * profile / np.sum(profile))
        position = (position + shift / n) % 1
        records = np.empty((2 * len(position), 8))
        for sector, ptype in enumerate([1, 30]):
            records[sector::2, :3] = position - 0.5
            records[sector::2, 3:6] = 0
            records[sector::2, 6] = masses[sector]
            records[sector::2, 7] = ptype
        np.savetxt(directory / "ic_part", records, fmt=["%.17e"] * 7 + ["%d"])
        text = template.read_text()
        replacements = {"pic=.false.": "pic=.true.", "hydro=.true.": "hydro=.false.",
                        "gravity_type=-1": "gravity_type=0", "nparttot=1": "nparttot=1000000",
                        "levelmin=6": f"levelmin={level}", "levelmax=6": f"levelmax={level}",
                        "&AMR_PARAMS": "&AMR_PARAMS\n ngridmax_auto=.false.",
                        "&INIT_PARAMS": "&INIT_PARAMS\n initfile(1)='.'", 
                        "&MOND_PARAMS": "&MOND_PARAMS\n use_bimond=.true.\n bimond_twin_matter=.true.\n bimond_output_fields=.true."}
        # read_params shifts the input list from slot 1 to levelmin.
        for old, new in replacements.items():
            if text.count(old) != 1:
                raise ValueError(old)
            text = text.replace(old, new)
        (directory / "effective.nml").write_text(text)
        (directory / "jobcontrol.txt").touch()
    (root / "output_policy.json").write_text(json.dumps({
        "purpose": "static two-sector particle-mesh force verification", "cases": [c[0] for c in cases],
        "noutput_per_case": 1, "aout": [1.1], "tout": ["1e100"], "foutput": 1,
        "fbackup": 1000000000, "nstepmax": 1, "expected_dumps": 4,
        "estimated_maximum_bytes_per_dump": 150000000, "reserved_total_bytes": 1073741824,
        "physics": "collisionless static sources; hydro, feedback, metals, dust and RT inactive",
    }, indent=2) + "\n")


def assess(directory):
    with h5py.File(directory / "output_00001/data_00001.h5") as handle:
        level = int(np.asarray(handle["header"].attrs["levelmin"]).item())
        n = 2**level
        grid = handle[f"amr/level_{level}"]
        force_group = handle[f"gravity/level_{level}"]
        xg = np.column_stack([grid[f"xg_{i}"][:] for i in (1, 2, 3)])
        children = np.array([[c % 2, c // 2 % 2, c // 4] for c in range(8)])
        index = np.rint((xg[:, None, :] + (children - 0.5) / n) * n - 0.5).astype(int) % n
        assert len(xg) * 8 == n**3
        assert np.unique(np.ravel_multi_index(index.reshape(-1, 3).T, (n,) * 3)).size == n**3

        def field(name):
            values = np.empty((n,) * 3)
            values[index[..., 0], index[..., 1], index[..., 2]] = force_group[name][:].reshape(-1, 8)
            return values

        measured = [np.stack([field(f"f_{i}") for i in (1, 2, 3)]),
                    np.stack([field(f"bimond_f_twin_{i}") for i in (1, 2, 3)])]
        densities = [field("bimond_rho_ordinary"), field("bimond_rho_twin")]
        phantom = field("bimond_rho_phantom")
        particle = handle["particles"]
        position = np.column_stack([particle[f"x_{i}"][:] for i in (1, 2, 3)])
        velocity = np.column_stack([particle[f"v_{i}"][:] for i in (1, 2, 3)])
        mass = particle["mass"][:]
        types = particle["ptypep"][:]
    masks = [types == 1, types == 30]
    assert set(types) == {1, 30}
    assert all(np.count_nonzero(mask) == n**3 for mask in masks)
    source_errors = []
    for density, mask in zip(densities, masks):
        independent = deposit(position[mask], mass[mask], n)
        error = np.linalg.norm(density - independent) / np.linalg.norm(independent)
        source_errors.append(float(error))
        assert error < 1e-12
        assert abs(np.sum(density) / n**3 / np.sum(mass[mask]) - 1) < 1e-12
    difference_force = fourth_order_force(periodic_poisson(densities[0] - densities[1], 1 / n), 1 / n)
    magnitude = np.linalg.norm(difference_force, axis=0)
    y = magnitude / (A0_CGS / (UNITS_LENGTH_CM / UNITS_TIME_S**2))
    nu = np.zeros_like(y)
    active = y >= 1e-30
    nu[active] = 0.5 * (np.sqrt(1 + 4 / y[active]) - 1)
    reference_phantom = -central_divergence(nu * difference_force, 1 / n) / FOURPI
    reference = [fourth_order_force(periodic_poisson(d + reference_phantom, 1 / n), 1 / n) for d in densities]
    common_norm = max(np.linalg.norm(f) for f in reference)
    errors = [float(np.linalg.norm(f - ref) / common_norm) for f, ref in zip(measured, reference)]
    assert max(errors) < 1e-9
    phantom_error = np.linalg.norm(phantom - reference_phantom) / max(np.linalg.norm(densities[0]), np.linalg.norm(reference_phantom))
    assert phantom_error < 1e-9
    # static=True writes interpolated acceleration into vp, testing species routing.
    kick_errors = []
    net_force = []
    for field_values, mask in zip(measured, masks):
        acceleration = interpolate(position[mask], field_values)
        kick_errors.append(float(np.linalg.norm(velocity[mask] - acceleration) / np.linalg.norm(acceleration)))
        net_force.append(np.sum(mass[mask, None] * acceleration, axis=0))
    assert max(kick_errors) < 1e-10
    signed_force = np.linalg.norm(net_force[0] - net_force[1])
    net_scale = sum(np.linalg.norm(f) for f in net_force)
    result = {"grid": n, "density_rms_errors": source_errors,
              "force_rms_errors_common_norm": errors, "phantom_error": float(phantom_error),
              "static_kick_rms_errors": kick_errors, "net_sector_forces": [f.tolist() for f in net_force],
              "signed_force_relative": float(signed_force / net_scale) if net_scale else 0.0}
    if directory.name.startswith("equal"):
        assert np.array_equal(densities[0], densities[1])
        assert np.count_nonzero(phantom) == 0
        result["equal_sector_cancellation_exact"] = True
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["prepare", "assess"])
    parser.add_argument("root", type=Path)
    parser.add_argument("--template", type=Path)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.root, args.template)
    else:
        results = {name: assess(args.root / name) for name in ["equal32", "pair32", "pair64", "shift32"]}
        results["signed_force_resolution_decreases"] = results["pair64"]["signed_force_relative"] < results["pair32"]["signed_force_relative"]
        results["passed_field_and_static_kick_checks"] = True
        (args.root / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
