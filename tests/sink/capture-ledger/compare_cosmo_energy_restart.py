#!/usr/bin/env python3
"""Bounded two-rank cosmological DMO restart comparison with a roundoff gate."""

import json
import sys

import h5py
import numpy as np


RELATIVE_PEAK_TOLERANCE = 1.0e-12
HEADER_FIELDS = (
    "nstep_coarse", "const", "mass_tot_0", "aexp_old",
    "epot_tot_old", "epot_tot_int",
)
ENERGY_FIELDS = frozenset(("epot_tot_old", "epot_tot_int"))


def datasets(handle):
    result = {}

    def collect(name, obj):
        if isinstance(obj, h5py.Dataset):
            result[name] = obj

    handle.visititems(collect)
    return result


def compare_array(left, right, allow_roundoff):
    if left.shape != right.shape or left.dtype != right.dtype:
        return {"pass": False, "reason": "shape_or_dtype_mismatch"}
    if np.issubdtype(left.dtype, np.number):
        if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
            return {"pass": False, "reason": "nonfinite_value"}
    if np.array_equal(left, right):
        return {"pass": True, "max_relative_to_reference_peak": 0.0}
    if not allow_roundoff or not np.issubdtype(left.dtype, np.floating):
        return {"pass": False, "reason": "nonidentical_discrete_or_exact_field"}
    if left.size == 0:
        return {"pass": False, "reason": "empty_mismatch"}
    maximum = float(np.max(np.abs(left - right)))
    peak = float(np.max(np.abs(left)))
    relative = maximum / max(peak, np.finfo(np.float64).tiny)
    return {
        "pass": relative <= RELATIVE_PEAK_TOLERANCE,
        "max_abs": maximum,
        "reference_peak": peak,
        "max_relative_to_reference_peak": relative,
    }


def compare(reference_path, replay_path):
    with h5py.File(reference_path) as reference, h5py.File(replay_path) as replay:
        original = datasets(reference)
        resumed = datasets(replay)
        if original.keys() != resumed.keys():
            raise ValueError("HDF5 dataset-key mismatch")
        ids_left = reference["/particles/identity"][...]
        ids_right = replay["/particles/identity"][...]
        if np.unique(ids_left).size != ids_left.size or np.unique(ids_right).size != ids_right.size:
            raise ValueError("cosmological test requires unique particle identities")
        order_left = np.argsort(ids_left)
        order_right = np.argsort(ids_right)
        if not np.array_equal(ids_left[order_left], ids_right[order_right]):
            raise ValueError("particle identity set differs")

        fields = {}
        for name in sorted(original):
            left = original[name][...]
            right = resumed[name][...]
            if name.startswith("particles/") and name != "particles/npart_per_cpu":
                left = left[order_left]
                right = right[order_right]
            allow_roundoff = name.startswith("gravity/") or name.startswith("particles/")
            fields[name] = compare_array(left, right, allow_roundoff)

        headers = {}
        for name in HEADER_FIELDS:
            left = np.asarray(reference["/header"].attrs[name])
            right = np.asarray(replay["/header"].attrs[name])
            headers[name] = compare_array(left, right, name in ENERGY_FIELDS)
            headers[name]["reference"] = left.tolist()
            headers[name]["replay"] = right.tolist()
        integral = float(np.asarray(reference["/header"].attrs["epot_tot_int"]).item())
        expansion = float(np.asarray(reference["/header"].attrs["aexp_old"]).item())

    return {
        "reference": reference_path,
        "replay": replay_path,
        "relative_peak_tolerance": RELATIVE_PEAK_TOLERANCE,
        "nonzero_cosmological_integral": integral != 0.0,
        "expansion_advanced": expansion != 0.1,
        "fields": fields,
        "headers": headers,
        "replay_within_roundoff": all(item["pass"] for item in fields.values())
        and all(item["pass"] for item in headers.values())
        and integral != 0.0 and expansion != 0.1,
    }


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare_cosmo_energy_restart.py FRESH.h5 RESTARTED.h5")
    print(json.dumps(compare(sys.argv[1], sys.argv[2]), indent=2, sort_keys=True))
