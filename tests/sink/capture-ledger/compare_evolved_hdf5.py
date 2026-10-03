#!/usr/bin/env python3
"""Quantify HDF5 state replay; do not turn a structural pass into a physics pass."""

import json
import sys

import h5py
import numpy as np


def datasets(handle):
    result = {}

    def collect(name, obj):
        if isinstance(obj, h5py.Dataset):
            result[name] = obj

    handle.visititems(collect)
    return result


def main(reference_path, replay_path):
    with h5py.File(reference_path) as reference, h5py.File(replay_path) as replay:
        reference_sets = datasets(reference)
        replay_sets = datasets(replay)
        if reference_sets.keys() != replay_sets.keys():
            raise ValueError(
                f"dataset mismatch: reference-only={sorted(reference_sets.keys() - replay_sets.keys())}; "
                f"replay-only={sorted(replay_sets.keys() - reference_sets.keys())}"
            )

        groups = {}
        for name in sorted(reference_sets):
            original = reference_sets[name][...]
            resumed = replay_sets[name][...]
            if original.shape != resumed.shape or original.dtype != resumed.dtype:
                raise ValueError(f"shape/dtype mismatch in {name}")
            group = name.split("/", 1)[0]
            if np.issubdtype(original.dtype, np.number):
                finite_original = np.isfinite(original)
                finite_resumed = np.isfinite(resumed)
                if not np.array_equal(finite_original, finite_resumed):
                    raise ValueError(f"nonfinite mask mismatch in {name}")
                nonfinite = int(original.size - np.count_nonzero(finite_original))
                if nonfinite and group in ("hydro", "gravity", "particles"):
                    raise ValueError(f"nonfinite active state in {name}")
                if not np.array_equal(np.isnan(original), np.isnan(resumed)):
                    raise ValueError(f"NaN mask mismatch in {name}")
                difference = np.abs(original[finite_original].astype(np.float64) - resumed[finite_resumed].astype(np.float64))
                maximum = float(np.max(difference)) if difference.size else 0.0
                scale = float(np.max(np.abs(original[finite_original]))) if np.any(finite_original) else 0.0
                relative = maximum / max(scale, np.finfo(np.float64).tiny)
                equal = bool(np.array_equal(original[finite_original], resumed[finite_resumed]))
            else:
                nonfinite = 0
                equal = bool(np.array_equal(original, resumed))
                maximum = 0.0 if equal else None
                relative = maximum
            summary = groups.setdefault(
                group,
                {"datasets": 0, "unequal_datasets": 0, "max_abs": 0.0, "worst_abs_dataset": None,
                 "max_relative_to_reference_peak": 0.0, "worst_relative_dataset": None,
                 "nonfinite_values": 0, "nonfinite_datasets": []},
            )
            summary["datasets"] += 1
            if not equal:
                summary["unequal_datasets"] += 1
            if nonfinite:
                summary["nonfinite_values"] += nonfinite
                summary["nonfinite_datasets"].append(name)
            if maximum is not None and maximum > summary["max_abs"]:
                summary["max_abs"] = maximum
                summary["worst_abs_dataset"] = name
            if relative is not None and relative > summary["max_relative_to_reference_peak"]:
                summary["max_relative_to_reference_peak"] = relative
                summary["worst_relative_dataset"] = name

        header = {}
        for name in ("nstep_coarse", "const", "mass_tot_0"):
            if name not in reference["/header"].attrs or name not in replay["/header"].attrs:
                raise ValueError(f"missing checkpoint header attribute {name}")
            left = np.asarray(reference["/header"].attrs[name])
            right = np.asarray(replay["/header"].attrs[name])
            header[name] = {
                "reference": left.tolist(),
                "replay": right.tolist(),
                "exact": bool(np.array_equal(left, right)),
            }

    report = {
        "reference": reference_path,
        "replay": replay_path,
        "dataset_count": len(reference_sets),
        "all_datasets_exact": all(item["unequal_datasets"] == 0 for item in groups.values()),
        "active_state_exact": all(
            groups[group]["unequal_datasets"] == 0
            for group in ("hydro", "gravity", "sinks")
        ),
        "groups": groups,
        "header": header,
        "interpretation": "Exact active-state replay is a strict regression check, not a physical calibration claim.",
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare_evolved_hdf5.py FRESH.h5 RESTARTED.h5")
    main(sys.argv[1], sys.argv[2])
