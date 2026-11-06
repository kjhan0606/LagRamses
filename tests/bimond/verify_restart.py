#!/usr/bin/env python3
"""Prepare and compare native uniform-grid BIMOND checkpoint continuations."""

import argparse
import json
import shutil
from pathlib import Path

import h5py
import numpy as np

from verify_bimond_two_sector import prepare


def prepare_full(root, template):
    prepare(root, template, cases=[("full", 4, True, 0.27)])
    path = root / "full/effective.nml"
    text = path.read_text().replace("static=.true.", "static=.false.").replace("nstepmax=1", "nstepmax=4")
    text = text.replace("&MOND_PARAMS", "&MOND_PARAMS\n courant_factor=0.05")
    path.write_text(text)
    (root / "output_policy.json").write_text(json.dumps({
        "purpose": "four-step particle restart equivalence and model-identity rejection",
        "noutput": 1, "aout": [1.1], "tout": ["1e100"], "foutput": 1, "fbackup": 1000000000,
        "nstepmax": 4, "restart_step": 2, "expected_new_full_dumps": 4,
        "maximum_new_resume_dumps": 3, "invalid_restart_dumps": 0,
        "estimated_maximum_bytes_per_dump": 3000000, "reserved_total_bytes": 268435456,
        "physics": "collisionless ordinary and twin gravity; hydro and subgrid physics inactive",
    }, indent=2) + "\n")


def prepare_resumes(root):
    for name in ["resume", "wrong_a0", "missing_identity", "wrapped_tag"]:
        directory = root / name
        directory.mkdir()
        shutil.copytree(root / "full/output_00002", directory / "output_00002")
        text = (root / "full/effective.nml").read_text().replace("nrestart=0", "nrestart=2")
        text = text.replace("&OUTPUT_PARAMS", "&OUTPUT_PARAMS\n informat='hdf5'")
        if name == "wrong_a0":
            text = text.replace("a0_mond=1.2d-8", "a0_mond=1.3d-8")
        path = directory / "output_00002/data_00002.h5"
        if name == "missing_identity":
            with h5py.File(path, "r+") as handle:
                del handle["header"].attrs["bimond_enabled"]
        if name == "wrapped_tag":
            with h5py.File(path, "r+") as handle:
                types = handle["particles/ptypep"]
                index = np.flatnonzero(types[:] == 30)[0]
                types[index] = 286  # int8 conversion would silently wrap to 30.
        (directory / "effective.nml").write_text(text)
        (directory / "jobcontrol.txt").touch()


def read_final(directory):
    candidates = []
    for path in directory.glob("output_*/data_*.h5"):
        with h5py.File(path) as handle:
            step = int(np.asarray(handle["header"].attrs["nstep_coarse"]).item())
        candidates.append((step, str(path)))
    step, path = max(candidates)
    assert step == 4
    with h5py.File(path) as handle:
        header = handle["header"].attrs
        level = int(np.asarray(header["levelmin"]).item())
        n = 2**level
        particles = handle["particles"]
        order = np.argsort(particles["identity"][:])
        inventory = {name: particles[name][:][order] for name in [
            "identity", "ptypep", "mass", "x_1", "x_2", "x_3", "v_1", "v_2", "v_3"]}
        grid = handle[f"amr/level_{level}"]
        gravity = handle[f"gravity/level_{level}"]
        xg = np.column_stack([grid[f"xg_{i}"][:] for i in (1, 2, 3)])
        child = np.array([[i % 2, i // 2 % 2, i // 4] for i in range(8)])
        cell = np.rint(xg[:, None, :] * n + child - 1).astype(int) % n
        keys = np.ravel_multi_index(cell.reshape(-1, 3).T, (n, n, n))
        assert len(np.unique(keys)) == n**3
        cell_order = np.argsort(keys)
        fields = {name: gravity[name][:].reshape(-1)[cell_order] for name in gravity.keys()}
        time = float(np.asarray(header["time"]).item())
        schema = int(np.asarray(header["bimond_nr_schema"]).item())
    return inventory, fields, time, schema


def assess(root):
    full, f_full, t_full, version = read_final(root / "full")
    resumed, f_resumed, t_resumed, resumed_version = read_final(root / "resume")
    assert version == resumed_version == 1
    assert t_full == t_resumed
    for name in ["identity", "ptypep", "mass"]:
        assert np.array_equal(full[name], resumed[name]), name
    ordinary = full["ptypep"] == 1
    twin = full["ptypep"] == 30
    assert np.count_nonzero(ordinary) == np.count_nonzero(twin) == 16**3
    state_errors = {name: float(np.max(np.abs(full[name] - resumed[name])))
                    for name in full if name.startswith(("x_", "v_"))}
    assert max(state_errors.values()) < 1e-12
    assert set(f_full) == set(f_resumed)
    field_errors = {}
    for name, expected in f_full.items():
        measured = f_resumed[name]
        denominator = max(np.linalg.norm(expected), np.linalg.norm(f_full["bimond_rho_ordinary"]) * 1e-14)
        field_errors[name] = float(np.linalg.norm(measured - expected) / denominator)
    assert max(field_errors.values()) < 1e-6
    result = {"passed": True, "final_step": 4, "time": t_full,
              "ordinary_count": int(ordinary.sum()), "twin_count": int(twin.sum()),
              "state_max_absolute_errors": state_errors, "field_normalized_rms_errors": field_errors,
              "scope": "same-rank uniform periodic CPU restart; both recomputed sector forces and particle state"}
    (root / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["prepare", "resume", "assess"])
    parser.add_argument("root", type=Path)
    parser.add_argument("--template", type=Path)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare_full(args.root, args.template)
    elif args.mode == "resume":
        prepare_resumes(args.root)
    else:
        assess(args.root)
