#!/usr/bin/env python3
"""Gate 256^3/512^3 density--velocity estimator response convergence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile

import numpy as np


APPROVED_ESTIMATOR_SHA256 = (
    "f82279067ea69c26f0d9a85c8ff5b14b70bce16bae0a6ccde1626480db106f27"
)
REQUIRED_CONTROL_MODEL = "lcdm_phase_matched"
PRIMARY_FIELDS = ("p_delta_theta", "p_theta_theta")
GROWTH_FIELDS = ("f_cross", "f_auto")
VARIANTS = {
    "empty_fill": {
        "p_delta_theta": "p_delta_theta_zero_fill",
        "p_theta_theta": "p_theta_theta_zero_fill",
    },
    "velocity_window": {
        "p_delta_theta": "p_delta_theta_velocity_not_deconvolved",
        "p_theta_theta": "p_theta_theta_velocity_not_deconvolved",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_product(path: Path) -> tuple[dict[str, np.ndarray], dict]:
    sidecar = path.with_suffix(".json")
    manifest_path = path.with_suffix(".manifest.json")
    complete_path = path.with_suffix(".COMPLETE")
    for product in (path, sidecar, manifest_path, complete_path):
        if not product.is_file():
            raise FileNotFoundError(f"incomplete product set: {product}")
    manifest = json.loads(manifest_path.read_text())
    complete = json.loads(complete_path.read_text())
    if complete.get("manifest_sha256") != sha256(manifest_path):
        raise ValueError(f"COMPLETE marker does not match manifest: {path}")
    expected = {
        path.name: sha256(path), sidecar.name: sha256(sidecar),
    }
    if manifest.get("files") != expected:
        raise ValueError(f"artifact hashes do not match manifest: {path}")
    with np.load(path, allow_pickle=False) as archive:
        arrays = {key: archive[key].copy() for key in archive.files}
    required = {"k_h_mpc", *PRIMARY_FIELDS, *GROWTH_FIELDS}
    required.update(value for mapping in VARIANTS.values() for value in mapping.values())
    missing = required - arrays.keys()
    if missing:
        raise ValueError(f"missing arrays in {path}: {sorted(missing)}")
    k = arrays["k_h_mpc"]
    if k.ndim != 1 or not k.size or not np.all(np.isfinite(k)) or np.any(np.diff(k) <= 0):
        raise ValueError(f"invalid k grid in {path}")
    for field in required - {"k_h_mpc"}:
        if arrays[field].shape != k.shape:
            raise ValueError(f"shape mismatch for {field} in {path}")
    metadata = json.loads(sidecar.read_text())
    return arrays, metadata


def exact_common_shells(products: list[dict[str, np.ndarray]]) -> np.ndarray:
    reference = products[0]["k_h_mpc"]
    for product in products[1:]:
        if not np.array_equal(product["k_h_mpc"], reference):
            raise ValueError("products do not use exact common k shells")
    return reference


def response(model: np.ndarray, control: np.ndarray, label: str) -> np.ndarray:
    valid = np.isfinite(model) & np.isfinite(control) & (control != 0.0)
    if not np.all(valid):
        bad = int(valid.size - np.count_nonzero(valid))
        raise ValueError(f"{label} has {bad} non-finite or zero-control shells")
    return model / control


def maximum_absolute(values: np.ndarray, mask: np.ndarray, label: str) -> float:
    selected = values[mask]
    if not selected.size or not np.all(np.isfinite(selected)):
        raise ValueError(f"no finite trusted shells for {label}")
    return float(np.max(np.abs(selected)))


def occupancy_direct_pass(metadata: dict) -> bool:
    occupancy = metadata.get("occupancy", {})
    fractions = (
        occupancy.get("empty_fraction_unshifted"),
        occupancy.get("empty_fraction_shifted"),
    )
    passes = (
        occupancy.get("fill_passes_unshifted"),
        occupancy.get("fill_passes_shifted"),
    )
    return all(
        isinstance(value, (int, float)) and math.isfinite(value) and value <= 1.0e-3
        for value in fractions
    ) and all(isinstance(value, int) and value <= 2 for value in passes)


def product_identity(metadata: dict, label: str) -> dict:
    snapshot = metadata.get("snapshot_preflight", {})
    identity = {
        "model": snapshot.get("model"), "nmesh": metadata.get("nmesh"),
        "aexp": metadata.get("aexp"), "boxlen_mpc_h": metadata.get("boxlen_mpc_h"),
        "source_output": metadata.get("source_output"),
        "script_sha256": metadata.get("script_sha256"),
    }
    if not isinstance(identity["model"], str) or not identity["model"]:
        raise ValueError(f"missing model identity in {label}")
    if not isinstance(identity["nmesh"], int):
        raise ValueError(f"missing integer nmesh in {label}")
    for key in ("aexp", "boxlen_mpc_h"):
        value = identity[key]
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"missing finite {key} in {label}")
    if not isinstance(identity["source_output"], str) or not identity["source_output"]:
        raise ValueError(f"missing source_output in {label}")
    if not isinstance(identity["script_sha256"], str) or len(identity["script_sha256"]) != 64:
        raise ValueError(f"missing estimator SHA256 in {label}")
    identity["output_name"] = Path(identity["source_output"]).name
    return identity


def validate_identities(
    metadata: dict[str, dict], requested_model: str, control_model: str,
    estimator_sha256: str, output_number: int, expected_aexp: float,
    boxlen_mpc_h: float,
) -> dict:
    identities = {
        name: product_identity(value, name) for name, value in metadata.items()
    }
    for resolution in (256, 512):
        if identities[f"model_{resolution}"]["nmesh"] != resolution:
            raise ValueError(f"model_{resolution} metadata has wrong nmesh")
        if identities[f"control_{resolution}"]["nmesh"] != resolution:
            raise ValueError(f"control_{resolution} metadata has wrong nmesh")
    model_names = {identities[name]["model"] for name in ("model_256", "model_512")}
    control_names = {
        identities[name]["model"] for name in ("control_256", "control_512")
    }
    if model_names != {requested_model}:
        raise ValueError(f"model inputs disagree with --model: {sorted(model_names)}")
    if control_names != {control_model}:
        raise ValueError(f"invalid control identity: {sorted(control_names)}")
    expected_output = f"output_{output_number:05d}"
    for name, identity in identities.items():
        if not math.isclose(
            identity["aexp"], expected_aexp, rel_tol=0.0, abs_tol=1e-12,
        ):
            raise ValueError(f"redshift/aexp mismatch in {name}")
        if not math.isclose(
            identity["boxlen_mpc_h"], boxlen_mpc_h,
            rel_tol=0.0, abs_tol=1e-12,
        ):
            raise ValueError(f"box length mismatch in {name}")
        if identity["output_name"] != expected_output:
            raise ValueError(f"snapshot output-number mismatch in {name}")
        if identity["script_sha256"] != estimator_sha256:
            raise ValueError(f"estimator SHA mismatch in {name}")
    return identities


def evaluate(args: argparse.Namespace) -> dict:
    paths = {
        "model_256": args.model_256.resolve(), "control_256": args.control_256.resolve(),
        "model_512": args.model_512.resolve(), "control_512": args.control_512.resolve(),
    }
    loaded = {name: load_product(path) for name, path in paths.items()}
    arrays = {name: value[0] for name, value in loaded.items()}
    metadata = {name: value[1] for name, value in loaded.items()}
    identities = validate_identities(
        metadata, args.model, args.control_model, args.estimator_sha256,
        args.output_number, args.expected_aexp, args.boxlen_mpc_h,
    )
    k = exact_common_shells(list(arrays.values()))
    trusted = k <= args.kmax
    if not np.any(trusted):
        raise ValueError(f"no common shells at k <= {args.kmax}")

    convergence = {}
    for field in (*PRIMARY_FIELDS, *GROWTH_FIELDS):
        response_256 = response(
            arrays["model_256"][field], arrays["control_256"][field], f"256 {field}"
        )
        response_512 = response(
            arrays["model_512"][field], arrays["control_512"][field], f"512 {field}"
        )
        maximum = maximum_absolute(
            response_512 / response_256 - 1.0, trusted, f"mesh {field}"
        )
        convergence[field] = {
            "max_abs_fractional_difference": maximum,
            "pass": maximum < args.mesh_tolerance,
        }

    sensitivity = {}
    for variant_name, mapping in VARIANTS.items():
        sensitivity[variant_name] = {}
        for resolution in (256, 512):
            model = arrays[f"model_{resolution}"]
            control = arrays[f"control_{resolution}"]
            sensitivity[variant_name][str(resolution)] = {}
            for primary, variant in mapping.items():
                primary_response = response(
                    model[primary], control[primary], f"{resolution} primary {primary}"
                )
                variant_response = response(
                    model[variant], control[variant], f"{resolution} {variant}"
                )
                maximum = maximum_absolute(
                    variant_response / primary_response - 1.0,
                    trusted, f"{resolution} {variant_name} {primary}",
                )
                sensitivity[variant_name][str(resolution)][primary] = {
                    "max_abs_fractional_difference": maximum,
                    "pass": maximum < args.systematic_tolerance,
                }

    occupancy = {
        name: {
            "direct_pass": occupancy_direct_pass(meta),
            "values": meta.get("occupancy", {}),
        }
        for name, meta in metadata.items()
    }
    convergence_pass = all(item["pass"] for item in convergence.values())
    velocity_window_pass = all(
        item["pass"]
        for resolution in sensitivity["velocity_window"].values()
        for item in resolution.values()
    )
    empty_sensitivity_pass = all(
        item["pass"]
        for resolution in sensitivity["empty_fill"].values()
        for item in resolution.values()
    )
    empty_direct_pass = all(item["direct_pass"] for item in occupancy.values())
    empty_fill_pass = empty_sensitivity_pass
    overall = convergence_pass and velocity_window_pass and empty_fill_pass
    return {
        "status": "PASS" if overall else "HOLD",
        "model": args.model,
        "control_model": args.control_model,
        "expected_identity": {
            "estimator_sha256": args.estimator_sha256,
            "output_number": args.output_number,
            "aexp": args.expected_aexp,
            "boxlen_mpc_h": args.boxlen_mpc_h,
        },
        "comparator": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "thresholds": {
            "kmax_h_mpc": args.kmax,
            "mesh_fraction": args.mesh_tolerance,
            "systematic_fraction": args.systematic_tolerance,
        },
        "trusted_shells": int(np.count_nonzero(trusted)),
        "k_range_h_mpc": [float(k[trusted][0]), float(k[trusted][-1])],
        "convergence": convergence,
        "sensitivity": sensitivity,
        "occupancy": occupancy,
        "identities": identities,
        "gates": {
            "mesh_convergence": convergence_pass,
            "velocity_window": velocity_window_pass,
            "empty_fill": empty_fill_pass,
            "empty_fill_direct_occupancy": empty_direct_pass,
            "empty_fill_response_sensitivity_required": empty_sensitivity_pass,
        },
        "inputs": {
            name: {
                "path": str(path), "sha256": sha256(path),
                "sidecar_sha256": sha256(path.with_suffix(".json")),
            }
            for name, path in paths.items()
        },
    }


def write_atomic(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True)
    parser.add_argument("--control-model", required=True)
    parser.add_argument("--estimator-sha256", required=True)
    parser.add_argument("--output-number", type=int, required=True)
    parser.add_argument("--expected-aexp", type=float, required=True)
    parser.add_argument("--boxlen-mpc-h", type=float, required=True)
    parser.add_argument("--model-256", type=Path, required=True)
    parser.add_argument("--control-256", type=Path, required=True)
    parser.add_argument("--model-512", type=Path, required=True)
    parser.add_argument("--control-512", type=Path, required=True)
    parser.add_argument("--kmax", type=float, default=0.1)
    parser.add_argument("--mesh-tolerance", type=float, default=0.01)
    parser.add_argument("--systematic-tolerance", type=float, default=0.005)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{64}", args.estimator_sha256):
        parser.error("--estimator-sha256 must be a lowercase 64-digit SHA256")
    if args.estimator_sha256 != APPROVED_ESTIMATOR_SHA256:
        parser.error("--estimator-sha256 is not the audited estimator SHA256")
    if args.control_model != REQUIRED_CONTROL_MODEL:
        parser.error(f"--control-model must be {REQUIRED_CONTROL_MODEL}")
    if args.output_number < 0:
        parser.error("--output-number must be non-negative")
    for name in (
        "expected_aexp", "boxlen_mpc_h", "kmax", "mesh_tolerance",
        "systematic_tolerance",
    ):
        value = getattr(args, name)
        if not math.isfinite(value) or value <= 0.0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    result = evaluate(args)
    write_atomic(args.output.resolve(), result)
    print(f"{result['status']} {args.output}")


if __name__ == "__main__":
    main()
