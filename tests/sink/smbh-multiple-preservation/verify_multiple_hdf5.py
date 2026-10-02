#!/usr/bin/env python3
"""Check five-seed SMBH preservation and full HDF5 sink restart state."""

from __future__ import annotations

import argparse
from fractions import Fraction
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
        domain = stream["domain"]
        if int(np.asarray(domain.attrs["bound_key_split_format"]).item()) != 1:
            raise ValueError("HDF5 Hilbert bounds are not losslessly split")
        high = domain["bound_key_hi"][...]
        low = domain["bound_key_lo"][...]
        legacy = domain["bound_key"][...]
        if high.shape != low.shape or high.shape != legacy.shape or \
                high.size < 2 or not np.array_equal(high, legacy) or \
                not np.all(np.isfinite(high)) or not np.all(np.isfinite(low)):
            raise ValueError("invalid HDF5 split Hilbert bounds")
        exact_bounds = [Fraction(float(hi)) + Fraction(float(lo))
                        for hi, lo in zip(high, low, strict=True)]
        if any(left >= right for left, right in
               zip(exact_bounds, exact_bounds[1:])):
            raise ValueError("non-monotone HDF5 Hilbert bounds")
        group = stream["sinks"]
        if int(np.asarray(group.attrs["nsink"]).item()) != 4 or \
                int(np.asarray(group.attrs["nindsink"]).item()) != 5:
            raise ValueError("HDF5 live SMBH counts are wrong")
        if int(np.asarray(group.attrs["sink_stat_format"]).item()) != 20261003:
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


def compare_full_checkpoint(continuous: Path, restarted: Path, output: int) -> tuple[int, int]:
    name = f"output_{output:05d}/data_{output:05d}.h5"
    with h5py.File(continuous / name, "r") as reference, \
            h5py.File(restarted / name, "r") as actual:
        reference_names = ["/"]
        actual_names = ["/"]
        reference.visit(reference_names.append)
        actual.visit(actual_names.append)
        if reference_names != actual_names:
            raise ValueError("HDF5 restart object topology differs")
        dataset_count = 0
        attribute_count = 0
        for path in reference_names:
            left, right = reference[path], actual[path]
            if type(left) is not type(right) or set(left.attrs) != set(right.attrs):
                raise ValueError(f"HDF5 restart object or attributes differ at {path}")
            for key in left.attrs:
                attribute_count += 1
                if not np.array_equal(left.attrs[key], right.attrs[key]):
                    raise ValueError(f"HDF5 restart attribute differs at {path}:{key}")
            if isinstance(left, h5py.Dataset):
                dataset_count += 1
                if left.dtype != right.dtype or not np.array_equal(left[...], right[...]):
                    raise ValueError(f"HDF5 restart dataset differs at {path}")
    return dataset_count, attribute_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("continuous_directory", type=Path)
    parser.add_argument("restart_directory", type=Path, nargs="?")
    parser.add_argument("--target-output", type=int, default=2)
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
        if args.target_output > 2:
            read_sinks(args.continuous_directory, args.target_output)
        binary_events, multiple_events = check_ledger(
            args.continuous_directory / "smbh_capture_ledger_v1.jsonl"
        )
        print(json.dumps({"status": "hdf5_multiple_smoke_passed",
                          "binary_events": binary_events,
                          "multiple_events": multiple_events,
                          "dataset_count": len(second)}, sort_keys=True))
        return

    target = second if args.target_output == 2 else read_sinks(
        args.continuous_directory, args.target_output
    )
    restarted = read_sinks(args.restart_directory, args.target_output)
    if set(target) != set(restarted):
        raise ValueError("HDF5 sink dataset sets differ after restart")
    for name in target:
        expected, actual = target[name], restarted[name]
        if not np.array_equal(expected, actual):
            raise ValueError(f"HDF5 restart state differs in {name}: "
                             f"continuous={expected}, restarted={actual}")
    dataset_count, attribute_count = compare_full_checkpoint(
        args.continuous_directory, args.restart_directory, args.target_output
    )
    print(json.dumps({"status": "hdf5_multiple_restart_passed",
                      "compared_sink_datasets": len(target),
                      "compared_checkpoint_datasets": dataset_count,
                      "compared_checkpoint_attributes": attribute_count}, sort_keys=True))


if __name__ == "__main__":
    main()
