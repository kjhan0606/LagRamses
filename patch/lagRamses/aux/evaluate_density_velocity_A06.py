#!/usr/bin/env python3
"""Evaluate the Paper-II A06 density--velocity consistency gate."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import tempfile

import numpy as np


APPROVED_ESTIMATOR_SHA256 = (
    "6638ac68b095d43be0edc6a71b8441d9078558d775b97fdb6a9095984fe18b18"
)
APPROVED_COMPARATOR_SHA256 = (
    "416de952cb36648e5438810597791f1d9201e8944cc6ff199a4a250e3f763a55"
)
CONTROL_MODEL = "lcdm_phase_matched"
MODELS = ("cpl_cluster_m09_p02", "f5", "f6")
TARGET_AEXP = 0.333333333000002
BOXLEN_MPC_H = 512.0
KMAX_H_MPC = 0.1
MIN_R_DELTA_THETA = 0.995
MAX_F_CROSS_AUTO_FRACTION = 0.02
MAX_MESH_RESPONSE_FRACTION = 0.01
APPROVED_PK_DEPENDENCY_SHA256 = (
    "b78c72666be27de1db5e93474829a110b5f142cef68912ca73c6e57fcef6af53"
)
PROVENANCE_CONTRACT_PATH = Path(__file__).with_name("A06_PROVENANCE_CONTRACT.json")
PROVENANCE_CONTRACT_SHA256 = (
    "2afcb52d86486fcaaa8cb987b6dedd44b4ef5062b52748529ee74734a57d6402"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


COMPARATOR_PATH = Path(__file__).with_name(
    "compare_density_velocity_convergence.py"
).resolve()
if sha256(COMPARATOR_PATH) != APPROVED_COMPARATOR_SHA256:
    raise RuntimeError(f"unaudited comparator SHA256: {COMPARATOR_PATH}")
SPEC = importlib.util.spec_from_file_location("a06_convergence", COMPARATOR_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import comparator: {COMPARATOR_PATH}")
convergence = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(convergence)

if sha256(PROVENANCE_CONTRACT_PATH) != PROVENANCE_CONTRACT_SHA256:
    raise RuntimeError("unaudited A06 provenance contract")
PROVENANCE_CONTRACT = json.loads(PROVENANCE_CONTRACT_PATH.read_text())


def validate_external_provenance() -> dict:
    for label in ("phase_audit", "extension_manifest"):
        record = PROVENANCE_CONTRACT[label]
        path = Path(record["path"])
        if sha256(path) != record["sha256"]:
            raise ValueError(f"A06 {label} SHA mismatch")
    phase = json.loads(Path(PROVENANCE_CONTRACT["phase_audit"]["path"]).read_text())
    if phase.get("complete") is not True or phase.get("failures") != []:
        raise ValueError("IC phase audit is not complete/PASS")
    minimum = PROVENANCE_CONTRACT["phase_audit"]["minimum_correlation"]
    for model in (CONTROL_MODEL, *MODELS):
        value = phase.get("correlations", {}).get(model, {}).get(
            "versus_lcdm_phase_matched"
        )
        if not isinstance(value, (int, float)) or value < minimum:
            raise ValueError(f"phase-matched gate failed for {model}: {value}")
    campaign_manifest = Path(PROVENANCE_CONTRACT["campaign"]) / "campaign.json"
    if sha256(campaign_manifest) != PROVENANCE_CONTRACT["campaign_manifest_sha256"]:
        raise ValueError("campaign manifest SHA mismatch")
    for model, contract in PROVENANCE_CONTRACT["models"].items():
        completion = Path(contract["completion_path"])
        if sha256(completion) != contract["completion_sha256"]:
            raise ValueError(f"completion evidence changed for {model}")
        run_namelist = Path(contract["source_output"]).parent / "run.nml"
        if sha256(run_namelist) != contract["run_namelist_sha256"]:
            raise ValueError(f"run.nml changed for {model}")
        if "original_run_log" in contract and sha256(
            Path(contract["original_run_log"])
        ) != contract["original_run_log_sha256"]:
            raise ValueError(f"original run log changed for {model}")
    return {
        "contract_path": str(PROVENANCE_CONTRACT_PATH.resolve()),
        "contract_sha256": PROVENANCE_CONTRACT_SHA256,
        "phase_audit_path": str(
            Path(PROVENANCE_CONTRACT["phase_audit"]["path"]).resolve()
        ),
        "phase_audit_sha256": PROVENANCE_CONTRACT["phase_audit"]["sha256"],
        "extension_manifest_sha256": PROVENANCE_CONTRACT["extension_manifest"][
            "sha256"
        ],
    }


def require_close(value: object, expected: float, label: str) -> None:
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"missing finite {label}")
    if not math.isclose(float(value), expected, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"wrong {label}: {value}, expected {expected}")


def load_validated_product(path: Path, expected_model: str, nmesh: int):
    arrays, metadata = convergence.load_product(path.resolve())
    identity = convergence.product_identity(metadata, str(path))
    if identity["model"] != expected_model:
        raise ValueError(f"wrong model in {path}: {identity['model']}")
    if identity["nmesh"] != nmesh:
        raise ValueError(f"wrong nmesh in {path}: {identity['nmesh']}")
    if identity["script_sha256"] != APPROVED_ESTIMATOR_SHA256:
        raise ValueError(f"unaudited estimator in {path}")
    if identity["output_name"] != "output_00002":
        raise ValueError(f"wrong output in {path}: {identity['output_name']}")
    require_close(identity["aexp"], TARGET_AEXP, f"aexp in {path}")
    require_close(identity["boxlen_mpc_h"], BOXLEN_MPC_H, f"box length in {path}")
    contract = PROVENANCE_CONTRACT["models"][expected_model]
    if identity["source_output"] != contract["source_output"]:
        raise ValueError(f"wrong exact source output in {path}")
    preflight = metadata.get("snapshot_preflight", {})
    if preflight.get("source_commit") != contract["source_commit"]:
        raise ValueError(f"source commit mismatch in {path}")
    completion = preflight.get("completion_contract", {})
    if completion.get("path") != contract["completion_path"] or completion.get(
        "sha256"
    ) != contract["completion_sha256"] or completion.get(
        "binary_sha256"
    ) != contract["binary_sha256"]:
        raise ValueError(f"completion provenance mismatch in {path}")
    provenance = preflight.get("provenance_files", {})
    if provenance.get("run.nml", {}).get("sha256") != contract[
        "run_namelist_sha256"
    ]:
        raise ValueError(f"run.nml SHA mismatch in {path}")
    if provenance.get("campaign.json", {}).get("sha256") != PROVENANCE_CONTRACT[
        "campaign_manifest_sha256"
    ]:
        raise ValueError(f"campaign manifest mismatch in {path}")
    parameters = preflight.get("model_parameters", {})
    for key, expected in contract["parameters"].items():
        if parameters.get(key) != expected:
            raise ValueError(
                f"model parameter mismatch in {path}: {key}={parameters.get(key)}"
            )
    dependency = preflight.get("python_dependencies", {}).get("measure_dmo_pk", {})
    if dependency.get("sha256") != APPROVED_PK_DEPENDENCY_SHA256:
        raise ValueError(f"measure_dmo_pk dependency mismatch in {path}")
    if "original_run_log" in contract:
        original = Path(contract["original_run_log"])
        if sha256(original) != contract["original_run_log_sha256"]:
            raise ValueError("F5 original run log changed")
    required = {"k_h_mpc", "nmodes", "r_delta_theta", "f_cross", "f_auto"}
    missing = required - arrays.keys()
    if missing:
        raise ValueError(f"missing arrays in {path}: {sorted(missing)}")
    k = arrays["k_h_mpc"]
    for field in required - {"k_h_mpc"}:
        if arrays[field].shape != k.shape:
            raise ValueError(f"shape mismatch for {field} in {path}")
    if (
        arrays["nmodes"].ndim != 1
        or not np.issubdtype(arrays["nmodes"].dtype, np.integer)
        or np.any(arrays["nmodes"] <= 0)
    ):
        raise ValueError(f"nmodes must be positive integers in {path}")
    return arrays, metadata, identity


def product_provenance(path: Path) -> dict:
    files = {
        "npz": path,
        "sidecar": path.with_suffix(".json"),
        "manifest": path.with_suffix(".manifest.json"),
        "complete": path.with_suffix(".COMPLETE"),
    }
    return {
        name: {"path": str(item.resolve()), "sha256": sha256(item.resolve())}
        for name, item in files.items()
    }


def evaluate_model(analysis_root: Path, model: str) -> dict:
    paths = {
        f"model_{nmesh}": analysis_root / model /
        f"density_velocity_00002_n{nmesh}.npz"
        for nmesh in (256, 512)
    }
    paths.update({
        f"control_{nmesh}": analysis_root / CONTROL_MODEL /
        f"density_velocity_00002_n{nmesh}.npz"
        for nmesh in (256, 512)
    })
    loaded = {}
    for role, path in paths.items():
        expected_model = model if role.startswith("model") else CONTROL_MODEL
        nmesh = int(role.rsplit("_", 1)[1])
        loaded[role] = load_validated_product(path, expected_model, nmesh)

    reference_k = loaded["model_256"][0]["k_h_mpc"]
    reference_modes = loaded["model_256"][0]["nmodes"]
    for role, (arrays, _, _) in loaded.items():
        if not np.array_equal(arrays["k_h_mpc"], reference_k):
            raise ValueError(f"non-identical k shells for {model}: {role}")
        if not np.array_equal(arrays["nmodes"], reference_modes):
            raise ValueError(f"non-identical mode counts for {model}: {role}")
    trusted = (
        (reference_k > 0.0) & (reference_k <= KMAX_H_MPC) &
        (reference_modes > 0)
    )
    if not np.any(trusted):
        raise ValueError(f"no trusted shells for {model}")

    field_gates = {}
    for role, (arrays, _, _) in loaded.items():
        selected = {
            field: arrays[field][trusted]
            for field in ("r_delta_theta", "f_cross", "f_auto")
        }
        if not all(np.all(np.isfinite(values)) for values in selected.values()):
            raise ValueError(f"non-finite A06 values for {model}: {role}")
        if np.any(selected["f_auto"] == 0.0):
            raise ValueError(f"zero f_auto for {model}: {role}")
        min_r = float(np.min(selected["r_delta_theta"]))
        max_f_difference = float(np.max(np.abs(
            selected["f_cross"] / selected["f_auto"] - 1.0
        )))
        field_gates[role] = {
            "min_r_delta_theta": min_r,
            "max_abs_f_cross_over_f_auto_minus_one": max_f_difference,
            "r_delta_theta_pass": min_r > MIN_R_DELTA_THETA,
            "f_cross_auto_pass": max_f_difference < MAX_F_CROSS_AUTO_FRACTION,
        }

    required_mesh = ("p_delta_theta", "p_theta_theta", "f_cross", "f_auto")
    mesh_values = {}
    for field in required_mesh:
        model_256 = loaded["model_256"][0][field][trusted]
        control_256 = loaded["control_256"][0][field][trusted]
        model_512 = loaded["model_512"][0][field][trusted]
        control_512 = loaded["control_512"][0][field][trusted]
        if np.any(control_256 == 0.0) or np.any(control_512 == 0.0):
            raise ValueError(f"zero control response denominator for {model}: {field}")
        response_256 = model_256 / control_256
        response_512 = model_512 / control_512
        if np.any(response_256 == 0.0):
            raise ValueError(f"zero N256 response for {model}: {field}")
        mesh_values[field] = float(
            np.max(np.abs(response_512 / response_256 - 1.0))
        )

    comparator_path = analysis_root / model / "convergence_00002_n256_n512.json"
    comparator_bytes = comparator_path.read_bytes()
    comparator = json.loads(comparator_bytes)
    if comparator.get("status") != "PASS" or comparator.get("model") != model:
        raise ValueError(f"comparator is not PASS for {model}")
    if comparator.get("control_model") != CONTROL_MODEL:
        raise ValueError(f"wrong comparator control for {model}")
    if comparator.get("comparator", {}).get("sha256") != APPROVED_COMPARATOR_SHA256:
        raise ValueError(f"unaudited comparator artifact for {model}")
    expected_thresholds = {
        "kmax_h_mpc": KMAX_H_MPC,
        "mesh_fraction": MAX_MESH_RESPONSE_FRACTION,
        "systematic_fraction": 0.005,
    }
    if comparator.get("thresholds") != expected_thresholds:
        raise ValueError(f"comparator thresholds changed for {model}")
    if comparator.get("trusted_shells") != int(np.count_nonzero(trusted)) or (
        comparator.get("k_range_h_mpc")
        != [float(reference_k[trusted][0]), float(reference_k[trusted][-1])]
    ):
        raise ValueError(f"comparator trusted k range mismatch for {model}")
    mesh = comparator.get("convergence", {})
    for field in required_mesh:
        entry = mesh.get(field, {})
        value = entry.get("max_abs_fractional_difference")
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"missing mesh statistic {field} for {model}")
        if entry.get("pass") is not True or value >= MAX_MESH_RESPONSE_FRACTION:
            raise ValueError(f"mesh response gate failed for {model}: {field}={value}")
        if not math.isclose(
            float(value), mesh_values[field], rel_tol=1.0e-12, abs_tol=1.0e-15
        ):
            raise ValueError(f"comparator/direct mesh mismatch for {model}: {field}")
    for role, path in paths.items():
        record = comparator.get("inputs", {}).get(role, {})
        if record.get("sha256") != sha256(path) or record.get(
            "sidecar_sha256"
        ) != sha256(path.with_suffix(".json")):
            raise ValueError(f"comparator/product hash mismatch for {model}: {role}")

    gate_values = [
        gate[name]
        for gate in field_gates.values()
        for name in ("r_delta_theta_pass", "f_cross_auto_pass")
    ]
    overall = all(gate_values) and max(mesh_values.values()) < MAX_MESH_RESPONSE_FRACTION
    return {
        "status": "PASS" if overall else "HOLD",
        "trusted_shells": int(np.count_nonzero(trusted)),
        "k_range_h_mpc": [
            float(reference_k[trusted][0]), float(reference_k[trusted][-1])
        ],
        "field_consistency": field_gates,
        "mesh_response_convergence": mesh_values,
        "comparator": {
            "path": str(comparator_path.resolve()),
            "sha256": hashlib.sha256(comparator_bytes).hexdigest(),
        },
        "inputs": {role: product_provenance(path) for role, path in paths.items()},
        "identities": {role: loaded[role][2] for role in loaded},
    }


def evaluate(analysis_root: Path) -> dict:
    external_provenance = validate_external_provenance()
    results = {model: evaluate_model(analysis_root, model) for model in MODELS}
    overall = all(item["status"] == "PASS" for item in results.values())
    return {
        "status": "PASS" if overall else "HOLD",
        "scope": "Paper-II A06 z=2 density--velocity estimator consistency",
        "science_interpretation": (
            "Estimator consistency and mesh-response convergence only; this is not "
            "an external modified-gravity growth benchmark."
        ),
        "models": results,
        "thresholds": {
            "kmax_h_mpc": KMAX_H_MPC,
            "min_r_delta_theta_strict": MIN_R_DELTA_THETA,
            "max_abs_f_cross_over_f_auto_minus_one_strict":
                MAX_F_CROSS_AUTO_FRACTION,
            "max_mesh_response_fraction_strict": MAX_MESH_RESPONSE_FRACTION,
        },
        "provenance": {
            "script_path": str(Path(__file__).resolve()),
            "script_sha256": sha256(Path(__file__).resolve()),
            "comparator_path": str(COMPARATOR_PATH),
            "comparator_sha256": APPROVED_COMPARATOR_SHA256,
            "estimator_sha256": APPROVED_ESTIMATOR_SHA256,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "external": external_provenance,
        },
    }


def write_atomic(path: Path, result: dict) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    directory_fd = os.open(path.parent, os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--complete-output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(args.analysis_root.resolve())
    write_atomic(args.output, result)
    if result["status"] == "PASS":
        write_atomic(
            args.complete_output,
            {
                "status": "A06_COMPLETE",
                "evaluation": str(args.output.resolve()),
                "evaluation_sha256": sha256(args.output.resolve()),
                "evaluator_sha256": sha256(Path(__file__).resolve()),
                "provenance_contract_sha256": PROVENANCE_CONTRACT_SHA256,
            },
        )
    print(f"{result['status']} {args.output.resolve()}")


if __name__ == "__main__":
    main()
