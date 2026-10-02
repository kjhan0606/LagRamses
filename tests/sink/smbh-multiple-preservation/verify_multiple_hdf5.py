#!/usr/bin/env python3
"""Check five-seed SMBH preservation and full HDF5 sink restart state."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from verify_multiple_smoke import EXPECTED_IDS, INITIAL_MASS_BY_ID, check_ledger


def read_sinks(root: Path, output: int) -> dict[str, np.ndarray]:
    path = root / f"output_{output:05d}" / f"data_{output:05d}.h5"
    if not (path.parent / "COMPLETE").is_file():
        raise ValueError(f"incomplete HDF5 checkpoint: {path}")
    with h5py.File(path, "r") as stream:
        group = stream["sinks"]
        if int(group.attrs["nsink"]) != 4 or int(group.attrs["nindsink"]) != 5:
            raise ValueError("HDF5 live SMBH counts are wrong")
        if int(group.attrs["sink_stat_format"]) != 20261003:
            raise ValueError("HDF5 global sink-stat marker is missing")
        values = {name: group[name][...] for name in group.keys()}
    ids = values["idsink"]
    if set(map(int, ids)) != EXPECTED_IDS:
        raise ValueError(f"HDF5 live SMBH IDs changed: {ids}")
    for name, field in values.items():
        if field.shape != (4,) or not np.all(np.isfinite(field)):
            raise ValueError(f"invalid HDF5 SMBH field {name}")
    statistics = [name for name in values if name.startswith("sink_stat_")]
    if len(statistics) != 7 or not any(np.any(values[name] != 0) for name in statistics):
        raise ValueError("global HDF5 sink statistics are absent")
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("continuous_directory", type=Path)
    parser.add_argument("restart_directory", type=Path, nargs="?")
    args = parser.parse_args()
    first = read_sinks(args.continuous_directory, 1)
    second = read_sinks(args.continuous_directory, 2)
    masses = dict(zip(map(int, first["idsink"]), first["msink"], strict=True))
    for sink_id, expected_mass in INITIAL_MASS_BY_ID.items():
        if not np.isclose(masses[sink_id], expected_mass, rtol=1e-12, atol=1e-12):
            raise ValueError(f"HDF5 binary compaction changed mass of sink {sink_id}")
    if not np.isclose(sum(masses.values()), 0.150, rtol=1e-12, atol=1e-12):
        raise ValueError("HDF5 binary compaction does not conserve seed mass")
    if args.restart_directory is None:
        binary_events, multiple_events = check_ledger(
            args.continuous_directory / "smbh_capture_ledger_v1.jsonl"
        )
        print(json.dumps({"status": "hdf5_multiple_smoke_passed",
                          "binary_events": binary_events,
                          "multiple_events": multiple_events,
                          "dataset_count": len(second)}, sort_keys=True))
        return

    restarted = read_sinks(args.restart_directory, 2)
    if set(second) != set(restarted):
        raise ValueError("HDF5 sink dataset sets differ after restart")
    largest = 0.0
    for name in second:
        expected, actual = second[name], restarted[name]
        largest = max(largest, float(np.max(np.abs(expected - actual))))
        if np.issubdtype(expected.dtype, np.integer):
            matches = np.array_equal(expected, actual)
        else:
            matches = np.allclose(expected, actual, rtol=1e-8, atol=1e-10)
        if not matches:
            raise ValueError(f"HDF5 restart state differs in {name}: "
                             f"continuous={expected}, restarted={actual}")
    print(json.dumps({"status": "hdf5_multiple_restart_passed",
                      "compared_datasets": len(second),
                      "largest_absolute_difference": largest}, sort_keys=True))


if __name__ == "__main__":
    main()
