#!/usr/bin/env python3
"""Fail closed unless all five DE extension gate runs are production-ready."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np


MODELS = (
    "cpl_cluster_m09_p02",
    "hs10_m01",
    "nu_lcdm",
    "f6_nu",
    "ede03",
)
REQUIRED_MARKERS = {
    "cpl_cluster_m09_p02": ("DE perturbation (table)", "CPL Dark Energy"),
    "hs10_m01": ("Horndeski mu(a,k) gravity enabled",),
    "nu_lcdm": ("Neutrino linear response",),
    "f6_nu": ("Neutrino linear response", "f(R) Hu-Sawicki gravity enabled"),
    "ede03": ("Early Dark Energy",),
}
BAD_MARKERS = re.compile(
    r"\b(?:FATAL|NaN|segmentation|SIGSEGV|MPI_ABORT|out of memory|killed)\b",
    re.IGNORECASE,
)
MCONS = re.compile(r"mcons=\s*([+-]?[0-9.]+(?:[Ee][+-]?[0-9]+)?)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_aexp(info: Path) -> float:
    for line in info.read_text().splitlines():
        fields = line.split()
        if fields and fields[0].lower().startswith("aexp"):
            return float(fields[-1])
    raise RuntimeError(f"aexp not found in {info}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    args = parser.parse_args()
    campaign = args.campaign.resolve()
    control = campaign / "extension5_20260805"
    validation = control / "validation"
    failures: list[str] = []
    audit: dict[str, object] = {"complete": False, "models": {}}

    reference_pnorm: float | None = None
    for model_name in MODELS:
        diag_path = campaign / "transfers" / f"transfer_{model_name}_z49.json"
        transfer_path = diag_path.with_suffix(".dat")
        if not diag_path.is_file() or not transfer_path.is_file():
            failures.append(f"{model_name}: transfer or diagnostics missing")
            continue
        diag = json.loads(diag_path.read_text())
        pnorm = float(diag["force_pnorm"])
        if reference_pnorm is None:
            reference_pnorm = pnorm
        if abs(pnorm / reference_pnorm - 1.0) > 5.0e-5:
            failures.append(f"{model_name}: common-As pnorm mismatch")
        if float(diag["force_pnorm_relative_scatter"]) > 1.0e-6:
            failures.append(f"{model_name}: non-constant CAMB/MUSIC normalization")
        transfer = np.loadtxt(transfer_path)
        if transfer.ndim != 2 or transfer.shape[1] < 13 or not np.isfinite(transfer).all():
            failures.append(f"{model_name}: invalid transfer table")

        ic_dir = validation / f"ics_{model_name}" / "level_006"
        components = (
            "ic_deltab",
            "ic_poscx",
            "ic_poscy",
            "ic_poscz",
            "ic_velcx",
            "ic_velcy",
            "ic_velcz",
        )
        for component in components:
            path = ic_dir / component
            if not path.is_file() or path.stat().st_size < 1_000_000:
                failures.append(f"{model_name}: incomplete IC component {component}")

        music_log_path = validation / f"music_{model_name}.conf_log.txt"
        music_log = music_log_path.read_text() if music_log_path.is_file() else ""
        if model_name == MODELS[0]:
            if "RNG-slab:" not in music_log:
                failures.append(f"{model_name}: distributed numeric-seed RNG path not exercised")
        else:
            phase_markers = (
                "white-noise filename is an explicit phase anchor",
                "Reading compact MUSIC white noise file",
            )
            for marker in phase_markers:
                if marker not in music_log:
                    failures.append(f"{model_name}: phase-anchor reader marker missing: {marker}")

        model_dir = validation / model_name
        stdout = model_dir / "run-gate.out"
        stderr = model_dir / "run-gate.err"
        text = stdout.read_text() if stdout.is_file() else ""
        if stderr.is_file() and stderr.stat().st_size:
            failures.append(f"{model_name}: non-empty stderr")
        if BAD_MARKERS.search(text):
            failures.append(f"{model_name}: fatal marker in log")
        for marker in REQUIRED_MARKERS[model_name]:
            if marker not in text:
                failures.append(f"{model_name}: activation marker missing: {marker}")
        values = [float(match.group(1)) for match in MCONS.finditer(text)]
        if not values or max(abs(value) for value in values) > 1.0e-12:
            failures.append(f"{model_name}: mass conservation gate failed")
        final_info = model_dir / "output_00005" / "info_00005.txt"
        if not final_info.is_file() or not math.isclose(read_aexp(final_info), 1.0, abs_tol=1.0e-12):
            failures.append(f"{model_name}: exact z=0 output missing")
        audit["models"][model_name] = {
            "transfer_sha256": sha256(transfer_path),
            "force_pnorm": pnorm,
            "omega_nu": float(diag.get("omega_nu", 0.0)),
            "z0": final_info.is_file(),
            "max_abs_mcons": max((abs(value) for value in values), default=None),
        }

    de_table = campaign / "tables" / "de_cpl_cluster_m09_p02.dat"
    nu_table = campaign / "tables" / "neutrino_mnu006.dat"
    for path in (de_table, nu_table):
        if not path.is_file() or path.stat().st_size < 10_000:
            failures.append(f"response table missing/incomplete: {path}")
    if (campaign / "wnoise_0010.bin").stat().st_size != 8_589_934_604:
        failures.append("production white-noise size changed")

    audit["failures"] = failures
    audit["complete"] = not failures
    output = control / "GATE_AUDIT.json"
    output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    (control / "GATE_PASSED").write_text("all five low-resolution gates passed\n")
    print("DE extension-5 gate: PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
