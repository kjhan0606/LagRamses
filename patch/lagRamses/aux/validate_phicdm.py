#!/usr/bin/env python3
"""Validate the production Ratra--Peebles phiCDM benchmark.

The gate combines independent checks that are easy to accidentally conflate:

1. parse the background values printed by the compiled lagRamses executable;
2. compare them with the parameter-matched lagCAMB scalar-field background;
3. check lagCAMB shooting, transfer, A_s and velocity diagnostics;
4. require exact requested RAMSES output expansion factors through z=0; and
5. compare the common-phase initial P_phiCDM/P_LCDM ratio with linear lagCAMB.

CAMB retains radiation while the N-body background intentionally starts from
the usual matter+DE RAMSES expansion.  At z_start=49 this leaves an explicitly
budgeted, sub-0.2% difference in rho_phi; later differences are much smaller.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

from dmo_benchmark_setup import configure_camb_dark_energy, load_local_camb


BG_PATTERN = re.compile(
    r"PHICDM_BG a=\s*(?P<a>[0-9.E+-]+)"
    r"\s+fde=\s*(?P<fde>[0-9.E+-]+)"
    r"\s+w=\s*(?P<w>[0-9.E+-]+)"
    r"\s+E2=\s*(?P<e2>[0-9.E+-]+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--model", default="phicdm_a01")
    parser.add_argument("--camb-dir", type=Path)
    parser.add_argument("--background-fde-rtol", type=float, default=2.0e-3)
    parser.add_argument("--background-w-atol", type=float, default=5.0e-4)
    parser.add_argument("--initial-pk-ratio-rtol", type=float, default=1.0e-3)
    parser.add_argument("--initial-pk-kmax", type=float, default=0.3)
    parser.add_argument(
        "--linear-ratio-accuracy-rtol", type=float, default=5.0e-4
    )
    parser.add_argument("--report", type=Path)
    return parser.parse_args()


def read_background(log_path: Path) -> dict[float, dict[str, float]]:
    samples: dict[float, dict[str, float]] = {}
    for match in BG_PATTERN.finditer(log_path.read_text()):
        row = {key: float(match.group(key)) for key in ("a", "fde", "w", "e2")}
        samples[row["a"]] = row
    if not samples:
        raise RuntimeError(f"no PHICDM_BG records in {log_path}")
    return samples


def camb_results(
    metadata: dict, model: str, camb, accuracy_boost: float = 1.0
):
    h = metadata["H0"] / 100.0
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=metadata["H0"],
        ombh2=metadata["omega_b"] * h**2,
        omch2=(metadata["omega_m"] - metadata["omega_b"]) * h**2,
        mnu=0.0,
    )
    pars.InitPower.set_params(ns=metadata["n_s"], As=metadata["A_s"])
    configure_camb_dark_energy(pars, model, camb)
    if model == "phicdm_a01" and accuracy_boost > 1.0:
        pars.DarkEnergy.npoints = 10_000
    pars.set_accuracy(
        AccuracyBoost=accuracy_boost, lAccuracyBoost=accuracy_boost
    )
    pars.set_matter_power(
        redshifts=[metadata["zstart"], *metadata["output_redshifts"]],
        kmax=2.0,
    )
    return camb.get_results(pars)


def linear_pk(results, redshift: float, k: np.ndarray) -> np.ndarray:
    kh, redshifts, spectra = results.get_linear_matter_power_spectrum(
        hubble_units=True, k_hunit=True
    )
    iz = int(np.argmin(np.abs(redshifts - redshift)))
    if not math.isclose(redshifts[iz], redshift, rel_tol=0.0, abs_tol=1.0e-8):
        raise RuntimeError(f"lagCAMB did not return z={redshift:g}")
    return np.exp(np.interp(np.log(k), np.log(kh), np.log(spectra[iz])))


def exact_output_a(model_dir: Path) -> list[float]:
    values = []
    for info in sorted(model_dir.glob("output_*/info_*.txt")):
        for line in info.read_text().splitlines():
            if line.lstrip().startswith("aexp"):
                values.append(float(line.split("=", 1)[1]))
                break
    return values


def main() -> int:
    args = parse_args()
    campaign = args.campaign.resolve()
    metadata = json.loads((campaign / "campaign.json").read_text())
    model_blocks = metadata["models"][args.model]["blocks"]
    alpha_match = re.search(
        r"(?im)^\s*quint_alpha\s*=\s*([0-9.EeDd+-]+)", model_blocks
    )
    tracker_match = re.search(
        r"(?im)^\s*quint_ic_mode\s*=\s*1\s*$", model_blocks
    )
    if not alpha_match or not tracker_match:
        raise ValueError(
            f"{args.model} is not a configured Ratra-Peebles tracker campaign"
        )
    alpha = float(alpha_match.group(1).replace("d", "e").replace("D", "E"))
    camb_dir = args.camb_dir or Path(metadata["camb_dir"])
    camb = load_local_camb(camb_dir)

    model_dir = campaign / args.model
    log_path = model_dir / "manual.out"
    if not log_path.is_file():
        logs = sorted(model_dir.glob("run-*.out"))
        if not logs:
            raise FileNotFoundError(f"no RAMSES log below {model_dir}")
        log_path = logs[-1]
    background = read_background(log_path)

    model_results = camb_results(metadata, args.model, camb)
    lcdm_results = camb_results(metadata, "lcdm", camb)
    model_results_high = camb_results(metadata, args.model, camb, 2.0)
    lcdm_results_high = camb_results(metadata, "lcdm", camb, 2.0)
    a_values = np.array(sorted(a for a in background if a >= 1.0 / 50.0))
    camb_fde, camb_w = model_results.get_dark_energy_rho_w(a_values)
    ramses_fde = np.array([background[a]["fde"] for a in a_values])
    ramses_w = np.array([background[a]["w"] for a in a_values])
    fde_residual = np.abs(ramses_fde / camb_fde - 1.0)
    w_residual = np.abs(ramses_w - camb_w)

    tracker_target = -2.0 / (alpha + 2.0)
    earliest = background[min(background)]
    tracker_w_residual = abs(earliest["w"] - tracker_target)
    closure_fde = abs(background[1.0]["fde"] - 1.0)
    closure_e2 = abs(background[1.0]["e2"] - 1.0)

    omega_de = model_results.get_Omega("de", 0.0)
    shooting_residual = abs(
        model_results.Params.DarkEnergy.omega_solved / omega_de - 1.0
    )

    transfer = metadata["transfer_diagnostics"][args.model]
    transfer_checks = {
        "force_pnorm_relative_scatter": transfer[
            "force_pnorm_relative_scatter"
        ],
        "vfact_scale_maximum_deviation": transfer[
            "vfact_scale_maximum_deviation"
        ],
        "sigma8_z0": transfer["sigma8_z0"],
    }

    requested_a = [
        1.0 / (1.0 + metadata["zstart"]),
        *(1.0 / (1.0 + z) for z in metadata["output_redshifts"]),
    ]
    actual_a = exact_output_a(model_dir)
    outputs_complete = len(actual_a) == len(requested_a) and all(
        math.isclose(got, want, rel_tol=0.0, abs_tol=1.0e-8)
        for got, want in zip(actual_a, requested_a)
    )

    first_output = model_dir / "output_00001" / "pk_cic_00001.dat"
    lcdm_first = campaign / "lcdm" / "output_00001" / "pk_cic_00001.dat"
    if not first_output.is_file() or not lcdm_first.is_file():
        raise FileNotFoundError("initial CIC spectra are required")
    model_pk = np.loadtxt(first_output)
    lcdm_pk = np.loadtxt(lcdm_first)
    if not np.array_equal(model_pk[:, 0], lcdm_pk[:, 0]):
        raise RuntimeError("model and LCDM CIC k grids differ")
    k = model_pk[:, 0]
    mask = (k <= args.initial_pk_kmax) & (lcdm_pk[:, 1] > 0.0)
    measured_ratio = model_pk[mask, 1] / lcdm_pk[mask, 1]
    theory_ratio = linear_pk(
        model_results, metadata["zstart"], k[mask]
    ) / linear_pk(lcdm_results, metadata["zstart"], k[mask])
    pk_ratio_residual = np.abs(measured_ratio / theory_ratio - 1.0)

    accuracy_k = np.geomspace(1.0e-3, 1.0, 500)
    ratio_accuracy_residuals = []
    for redshift in (metadata["zstart"], 0.0):
        standard_ratio = linear_pk(
            model_results, redshift, accuracy_k
        ) / linear_pk(lcdm_results, redshift, accuracy_k)
        high_ratio = linear_pk(
            model_results_high, redshift, accuracy_k
        ) / linear_pk(lcdm_results_high, redshift, accuracy_k)
        ratio_accuracy_residuals.append(
            float(np.max(np.abs(standard_ratio / high_ratio - 1.0)))
        )
    linear_ratio_accuracy_residual = max(ratio_accuracy_residuals)

    checks = {
        "lagramses_present_fde_closure": bool(closure_fde < 1.0e-10),
        "lagramses_present_E2_closure": bool(closure_e2 < 1.0e-10),
        "lagramses_tracker_limit": bool(tracker_w_residual < 1.0e-7),
        "lagcamb_shooting_closure": bool(shooting_residual < 1.0e-8),
        "background_fde_match": bool(
            float(fde_residual.max()) < args.background_fde_rtol
        ),
        "background_w_match": bool(
            float(w_residual.max()) < args.background_w_atol
        ),
        "transfer_amplitude_constant": bool(
            transfer_checks["force_pnorm_relative_scatter"] < 1.0e-6
        ),
        "velocity_correction_scale_independent": bool(
            transfer_checks["vfact_scale_maximum_deviation"] < 1.0e-4
        ),
        "outputs_complete_through_z0": bool(outputs_complete),
        "initial_pk_ratio_matches_lagcamb": bool(
            float(pk_ratio_residual.max()) < args.initial_pk_ratio_rtol
        ),
        "lagcamb_linear_ratio_accuracy_converged": bool(
            linear_ratio_accuracy_residual
            < args.linear_ratio_accuracy_rtol
        ),
    }
    report = {
        "complete": all(checks.values()),
        "checks": checks,
        "thresholds": {
            "background_fde_rtol": args.background_fde_rtol,
            "background_w_atol": args.background_w_atol,
            "initial_pk_ratio_rtol": args.initial_pk_ratio_rtol,
            "initial_pk_kmax_h_mpc": args.initial_pk_kmax,
            "linear_ratio_accuracy_rtol": args.linear_ratio_accuracy_rtol,
        },
        "metrics": {
            "maximum_background_fde_relative_residual": float(
                fde_residual.max()
            ),
            "maximum_background_w_absolute_residual": float(w_residual.max()),
            "tracker_w_absolute_residual": tracker_w_residual,
            "lagcamb_shooting_relative_residual": shooting_residual,
            "maximum_initial_pk_ratio_relative_residual": float(
                pk_ratio_residual.max()
            ),
            "maximum_lagcamb_linear_ratio_accuracy_residual": (
                linear_ratio_accuracy_residual
            ),
            **transfer_checks,
        },
        "campaign": str(campaign),
        "model": args.model,
        "alpha": alpha,
        "lagramses_log": str(log_path),
        "camb_module": str(Path(camb.__file__).resolve()),
    }
    report_path = args.report or campaign / f"{args.model}_validation.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
