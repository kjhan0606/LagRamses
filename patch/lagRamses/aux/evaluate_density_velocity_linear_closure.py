#!/usr/bin/env python3
"""Evaluate the Paper-II z=2 N256 density--velocity linear-closure gate."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import re
import tempfile

import numpy as np
import scipy
from scipy.integrate import solve_ivp

APPROVED_ESTIMATOR_SHA256 = (
    "6638ac68b095d43be0edc6a71b8441d9078558d775b97fdb6a9095984fe18b18"
)
APPROVED_PK_DEPENDENCY_SHA256 = (
    "b78c72666be27de1db5e93474829a110b5f142cef68912ca73c6e57fcef6af53"
)
APPROVED_COMPARATOR_SHA256 = (
    "416de952cb36648e5438810597791f1d9201e8944cc6ff199a4a250e3f763a55"
)
LCDM_MODEL = "lcdm_phase_matched"
CPL_MODEL = "cpl_m09_p02"
TARGET_A = 0.333333333000002
KMAX_H_MPC = 0.08
CAMPAIGN = Path(
    "/gpfs/kjhan/Hydro/DE_nonstd/DMO_production_L512_N1024_20260729"
)
LCDM_NAMELIST = CAMPAIGN / LCDM_MODEL / "run.nml"
CPL_NAMELIST = CAMPAIGN / CPL_MODEL / "run.nml"
APPROVED_NAMELIST_SHA256 = {
    "lcdm": "aff3e5317f0a1b229626207335a4b69b1d05bacaccb917a8424f5108f6edd4f7",
    "cpl": "85d64021e6981b0530fa3eae1ccd9724668636d5749294068a28bc01b7e4c630",
}


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
COMPARATOR_SPEC = importlib.util.spec_from_file_location(
    "paper2_density_velocity_convergence", COMPARATOR_PATH
)
if COMPARATOR_SPEC is None or COMPARATOR_SPEC.loader is None:
    raise RuntimeError(f"cannot load source-bound comparator: {COMPARATOR_PATH}")
convergence = importlib.util.module_from_spec(COMPARATOR_SPEC)
COMPARATOR_SPEC.loader.exec_module(convergence)


def namelist_parameters(path: Path) -> dict[str, float]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing or empty namelist: {path}")
    values: dict[str, float] = {}
    for raw in path.read_text(errors="strict").splitlines():
        line = raw.split("!", 1)[0]
        match = re.match(r"\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([^,\s/]+)", line)
        if not match:
            continue
        key, token = match.group(1).lower(), match.group(2).replace("D", "E").replace("d", "e")
        try:
            values[key] = float(token)
        except ValueError:
            continue
    for key in ("omega_m", "omega_l"):
        if key not in values or not math.isfinite(values[key]):
            raise ValueError(f"missing finite {key} in {path}")
    values.setdefault("w0", -1.0)
    values.setdefault("wa", 0.0)
    return values


def growth_rate(
    parameters: dict[str, float], target_a: float, initial_a: float,
    rtol: float, atol: float,
) -> float:
    omega_m, omega_de = parameters["omega_m"], parameters["omega_l"]
    omega_k = 1.0 - omega_m - omega_de
    w0, wa = parameters["w0"], parameters["wa"]
    if not 0.0 < initial_a < target_a <= 1.0:
        raise ValueError("growth integration requires 0 < initial_a < target_a <= 1")

    def derivative(log_a: float, state: np.ndarray) -> tuple[float, float]:
        a = math.exp(log_a)
        matter = omega_m * a**-3
        curvature = omega_k * a**-2
        dark_energy = omega_de * a**(-3.0 * (1.0 + w0 + wa)) * math.exp(
            3.0 * wa * (a - 1.0)
        )
        e2 = matter + curvature + dark_energy
        if not math.isfinite(e2) or e2 <= 0.0:
            raise ValueError(f"non-positive E(a)^2 at a={a}")
        dark_slope = -3.0 * (1.0 + w0 + wa) + 3.0 * wa * a
        dlnh = 0.5 * (-3.0 * matter - 2.0 * curvature + dark_slope * dark_energy) / e2
        growth, growth_prime = state
        return growth_prime, -(2.0 + dlnh) * growth_prime + 1.5 * matter / e2 * growth

    solution = solve_ivp(
        derivative, (math.log(initial_a), math.log(target_a)),
        (initial_a, initial_a), method="DOP853", rtol=rtol, atol=atol,
    )
    if not solution.success or solution.y.shape[1] == 0:
        raise RuntimeError(f"growth ODE failed: {solution.message}")
    growth, growth_prime = solution.y[:, -1]
    result = float(growth_prime / growth)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"invalid growth rate: {result}")
    return result


def expected_growth(path: Path, expected_sha256: str, target_a: float) -> dict:
    observed_sha256 = sha256(path)
    if observed_sha256 != expected_sha256:
        raise ValueError(f"unaudited run.nml SHA256: {path}: {observed_sha256}")
    parameters = namelist_parameters(path)
    baseline = growth_rate(parameters, target_a, 1.0e-4, 1.0e-11, 1.0e-13)
    early = growth_rate(parameters, target_a, 1.0e-5, 1.0e-11, 1.0e-13)
    tight = growth_rate(parameters, target_a, 1.0e-4, 1.0e-12, 1.0e-14)
    early_difference = abs(early / baseline - 1.0)
    tight_difference = abs(tight / baseline - 1.0)
    if early_difference >= 1.0e-5 or tight_difference >= 1.0e-5:
        raise ValueError(
            f"growth ODE convergence failed: early={early_difference}, tight={tight_difference}"
        )
    return {
        "f_linear": baseline, "initial_a": 1.0e-4, "target_a": target_a,
        "method": "DOP853", "rtol": 1.0e-11, "atol": 1.0e-13,
        "early_initial_a": 1.0e-5, "early_relative_difference": early_difference,
        "tight_rtol": 1.0e-12, "tight_atol": 1.0e-14,
        "tight_relative_difference": tight_difference,
        "parameters": parameters, "namelist_path": str(path.resolve()),
        "namelist_sha256": observed_sha256,
    }


def validate_product(metadata: dict, model: str) -> dict:
    identity = convergence.product_identity(metadata, model)
    if identity["model"] != model or identity["nmesh"] != 256:
        raise ValueError(f"wrong model or mesh for {model}: {identity}")
    if identity["script_sha256"] != APPROVED_ESTIMATOR_SHA256:
        raise ValueError(f"unaudited estimator SHA for {model}")
    if identity["output_name"] != "output_00002":
        raise ValueError(f"wrong output number for {model}")
    if not math.isclose(identity["aexp"], TARGET_A, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"wrong aexp for {model}: {identity['aexp']}")
    if not math.isclose(identity["boxlen_mpc_h"], 512.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"wrong box length for {model}: {identity['boxlen_mpc_h']}")
    dependency = metadata.get("snapshot_preflight", {}).get(
        "python_dependencies", {}
    ).get("measure_dmo_pk", {})
    if dependency.get("sha256") != APPROVED_PK_DEPENDENCY_SHA256:
        raise ValueError(f"unaudited measure_dmo_pk dependency for {model}")
    measured_kmax = metadata.get("kmax_h_mpc")
    if not isinstance(measured_kmax, (int, float)) or not math.isfinite(measured_kmax):
        raise ValueError(f"missing finite measured kmax for {model}")
    if measured_kmax < KMAX_H_MPC:
        raise ValueError(f"artifact kmax is below A05 range for {model}")
    return identity


def validate_arrays(arrays: dict[str, np.ndarray], model: str) -> None:
    required = ("k_h_mpc", "nmodes", "f_cross", "f_auto", "r_delta_theta")
    missing = set(required) - arrays.keys()
    if missing:
        raise ValueError(f"missing A05 arrays for {model}: {sorted(missing)}")
    k = arrays["k_h_mpc"]
    if k.ndim != 1 or not k.size:
        raise ValueError(f"invalid k array for {model}")
    for field in required[1:]:
        if arrays[field].ndim != 1 or arrays[field].shape != k.shape:
            raise ValueError(f"invalid {field} shape for {model}")
    nmodes = arrays["nmodes"]
    if not np.issubdtype(nmodes.dtype, np.integer) or np.any(nmodes <= 0):
        raise ValueError(f"nmodes must be positive integers for {model}")
    if not np.all(np.isfinite(k)) or np.any(np.diff(k) <= 0.0):
        raise ValueError(f"invalid k values for {model}")
    if k[-1] < KMAX_H_MPC:
        raise ValueError(f"artifact shells do not reach A05 kmax for {model}")


def artifact_provenance(path: Path) -> dict:
    artifacts = {
        "npz": path, "sidecar": path.with_suffix(".json"),
        "manifest": path.with_suffix(".manifest.json"),
        "complete": path.with_suffix(".COMPLETE"),
    }
    return {
        name: {"path": str(item), "sha256": sha256(item)}
        for name, item in artifacts.items()
    }


def weighted_statistics(
    observable: np.ndarray, expected: float, weights: np.ndarray,
    k: np.ndarray,
) -> dict:
    residual = observable / expected - 1.0
    rms = float(np.sqrt(np.sum(weights * residual**2) / np.sum(weights)))
    mean = float(np.sum(weights * residual) / np.sum(weights))
    index = int(np.argmax(np.abs(residual)))
    diagnostics = (rms, mean, float(abs(residual[index])), float(k[index]))
    if not all(math.isfinite(value) for value in diagnostics):
        raise ValueError("non-finite closure statistic")
    return {
        "weighted_rms": rms, "weighted_signed_mean": mean,
        "max_abs_shell_residual": float(abs(residual[index])),
        "max_abs_shell_k_h_mpc": float(k[index]),
        "max_abs_shell_nmodes": int(weights[index]), "pass": rms < 0.02,
    }


def evaluate(args: argparse.Namespace) -> dict:
    products = {
        "lcdm": convergence.load_product(args.lcdm.resolve()),
        "cpl": convergence.load_product(args.cpl.resolve()),
    }
    metadata = {key: value[1] for key, value in products.items()}
    identities = {
        "lcdm": validate_product(metadata["lcdm"], LCDM_MODEL),
        "cpl": validate_product(metadata["cpl"], CPL_MODEL),
    }
    lcdm, cpl = products["lcdm"][0], products["cpl"][0]
    validate_arrays(lcdm, LCDM_MODEL)
    validate_arrays(cpl, CPL_MODEL)
    if not np.array_equal(lcdm["k_h_mpc"], cpl["k_h_mpc"]):
        raise ValueError("LCDM and CPL k shells differ")
    if not np.array_equal(lcdm["nmodes"], cpl["nmodes"]):
        raise ValueError("LCDM and CPL mode counts differ")
    k, nmodes = lcdm["k_h_mpc"], lcdm["nmodes"]
    trusted = (k > 0.0) & (k <= KMAX_H_MPC) & (nmodes > 0)
    if not np.any(trusted):
        raise ValueError("no trusted low-k shells")
    for model, arrays in (("lcdm", lcdm), ("cpl", cpl)):
        for field in ("f_cross", "f_auto", "r_delta_theta"):
            if not np.all(np.isfinite(arrays[field][trusted])):
                raise ValueError(f"non-finite {model} {field} in trusted shells")
    growth = {
        "lcdm": expected_growth(
            LCDM_NAMELIST, APPROVED_NAMELIST_SHA256["lcdm"], TARGET_A
        ),
        "cpl": expected_growth(
            CPL_NAMELIST, APPROVED_NAMELIST_SHA256["cpl"], TARGET_A
        ),
    }
    weights, selected_k = nmodes[trusted].astype(float), k[trusted]
    closure: dict[str, dict] = {}
    differential: dict[str, dict] = {}
    for model, arrays in (("lcdm", lcdm), ("cpl", cpl)):
        closure[model] = {}
        for field in ("f_cross", "f_auto"):
            closure[model][field] = weighted_statistics(
                arrays[field][trusted], growth[model]["f_linear"], weights, selected_k
            )
        closure[model]["min_r_delta_theta"] = float(
            np.min(arrays["r_delta_theta"][trusted])
        )
    expected_ratio = growth["cpl"]["f_linear"] / growth["lcdm"]["f_linear"]
    for field in ("f_cross", "f_auto"):
        if np.any(lcdm[field][trusted] == 0.0):
            raise ValueError(f"zero LCDM denominator in {field}")
        residual = (
            cpl[field][trusted] / lcdm[field][trusted] / expected_ratio - 1.0
        )
        rms = float(np.sqrt(np.sum(weights * residual**2) / np.sum(weights)))
        mean = float(np.sum(weights * residual) / np.sum(weights))
        index = int(np.argmax(np.abs(residual)))
        diagnostics = (rms, mean, float(abs(residual[index])), float(selected_k[index]))
        if not all(math.isfinite(value) for value in diagnostics):
            raise ValueError(f"non-finite differential statistic for {field}")
        differential[field] = {
            "weighted_rms": rms, "weighted_signed_mean": mean,
            "max_abs_shell_residual": float(abs(residual[index])),
            "max_abs_shell_k_h_mpc": float(selected_k[index]),
            "max_abs_shell_nmodes": int(weights[index]), "pass": rms < 0.01,
        }
    gates = [
        closure[model][field]["pass"]
        for model in ("lcdm", "cpl") for field in ("f_cross", "f_auto")
    ] + [differential[field]["pass"] for field in ("f_cross", "f_auto")]
    return {
        "status": "PASS" if all(gates) else "HOLD",
        "scope": "A05 N256 linear-closure diagnostic only",
        "science_eligible": False,
        "science_hold_reason": "A04 N256/N512 mesh convergence is pending",
        "trusted_shells": int(np.count_nonzero(trusted)),
        "k_range_h_mpc": [float(selected_k[0]), float(selected_k[-1])],
        "thresholds": {"closure_weighted_rms": 0.02, "differential_weighted_rms": 0.01},
        "expected_growth": growth, "closure": closure, "differential": differential,
        "inputs": {
            key: {
                "artifacts": artifact_provenance(path.resolve()),
                "validated_identity": identities[key],
            }
            for key, path in (("lcdm", args.lcdm), ("cpl", args.cpl))
        },
        "provenance": {
            "script_path": str(Path(__file__).resolve()),
            "script_sha256": sha256(Path(__file__).resolve()),
            "comparator_path": str(COMPARATOR_PATH),
            "comparator_sha256": sha256(COMPARATOR_PATH),
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__,
            "equation": "D'' + (2+dlnH/dlna)D' - 3 Omega_m(a) D / 2 = 0",
        },
    }


def write_atomic(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lcdm", type=Path, required=True)
    parser.add_argument("--cpl", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(args)
    write_atomic(args.output.resolve(), result)
    print(f"{result['status']} {args.output}")


if __name__ == "__main__":
    main()
