#!/usr/bin/env python3
"""Validate CAMB vtotal initial velocities for the coupled-DE DMO benchmark."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import struct

import numpy as np

from dmo_benchmark_setup import OMEGA_B, OMEGA_M


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--model", default="cde10")
    parser.add_argument("--velocity-ratio-rtol", type=float, default=1.0e-3)
    parser.add_argument(
        "--initial-density-ratio-rtol", type=float, default=1.0e-3
    )
    parser.add_argument("--kmax", type=float, default=1.0)
    parser.add_argument("--initial-pk-kmax", type=float, default=0.3)
    parser.add_argument("--report", type=Path)
    return parser.parse_args()


def read_grafic(path: Path) -> tuple[dict[str, float | int], np.ndarray]:
    """Read one single-level GRAFIC2 field as an array ordered (z,y,x)."""
    with path.open("rb") as stream:
        marker = struct.unpack("=i", stream.read(4))[0]
        header_raw = stream.read(marker)
        if struct.unpack("=i", stream.read(4))[0] != marker:
            raise RuntimeError(f"broken GRAFIC header record: {path}")
        if marker != struct.calcsize("=3i8f"):
            raise RuntimeError(f"unexpected GRAFIC header size {marker}: {path}")
        values = struct.unpack("=3i8f", header_raw)
        nx, ny, nz = values[:3]
        planes = []
        for _ in range(nz):
            plane_bytes = struct.unpack("=i", stream.read(4))[0]
            expected = nx * ny * np.dtype(np.float32).itemsize
            if plane_bytes != expected:
                raise RuntimeError(f"unexpected GRAFIC plane size: {path}")
            plane = np.frombuffer(
                stream.read(plane_bytes), dtype=np.float32
            ).copy()
            if struct.unpack("=i", stream.read(4))[0] != plane_bytes:
                raise RuntimeError(f"broken GRAFIC plane record: {path}")
            planes.append(plane)
        if stream.read(1):
            raise RuntimeError(f"trailing bytes in GRAFIC field: {path}")
    header = {
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "dx": values[3],
        "astart": values[7],
        "omega_m": values[8],
        "omega_l": values[9],
        "H0": values[10],
    }
    return header, np.asarray(planes).reshape(nz, ny, nx)


def velocity_divergence(
    level_dir: Path, boxlen: float
) -> tuple[np.ndarray, np.ndarray]:
    fields = []
    reference_header = None
    for component in "xyz":
        header, field = read_grafic(level_dir / f"ic_velc{component}")
        if reference_header is None:
            reference_header = header
        elif header != reference_header:
            raise RuntimeError("GRAFIC velocity headers differ")
        fields.append(field)
    nz, ny, nx = fields[0].shape
    kx = 2.0 * np.pi * np.fft.rfftfreq(nx, d=boxlen / nx)[None, None, :]
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=boxlen / ny)[None, :, None]
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=boxlen / nz)[:, None, None]
    divergence = 1j * (
        kx * np.fft.rfftn(fields[0])
        + ky * np.fft.rfftn(fields[1])
        + kz * np.fft.rfftn(fields[2])
    )
    return np.sqrt(kx * kx + ky * ky + kz * kz), divergence


def mass_weighted_velocity(table: np.ndarray) -> np.ndarray:
    omega_c = OMEGA_M - OMEGA_B
    return (omega_c * table[:, 10] + OMEGA_B * table[:, 11]) / OMEGA_M


def interpolate_positive(
    source_k: np.ndarray, source_value: np.ndarray, target_k: np.ndarray
) -> np.ndarray:
    return np.exp(
        np.interp(np.log(target_k), np.log(source_k), np.log(source_value))
    )


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
    model = args.model
    diagnostics = metadata["transfer_diagnostics"][model]
    if diagnostics.get("dmo_velocity_mode") != "camb_vtotal":
        raise RuntimeError(f"{model} does not use the CAMB vtotal path")

    level = int(metadata["levelmin"])
    boxlen = float(metadata["boxlen_mpc_h"])
    lcdm_level = campaign / "ics_lcdm" / f"level_{level:03d}"
    model_level = campaign / f"ics_{model}" / f"level_{level:03d}"
    k, theta_lcdm = velocity_divergence(lcdm_level, boxlen)
    model_k, theta_model = velocity_divergence(model_level, boxlen)
    if not np.array_equal(k, model_k):
        raise RuntimeError("model and LCDM velocity Fourier grids differ")

    zstart = float(metadata["zstart"])
    lcdm_transfer = np.loadtxt(
        campaign / "transfers" / f"transfer_lcdm_z{zstart:g}.dat"
    )
    model_transfer = np.loadtxt(
        campaign / "transfers" / f"transfer_{model}_z{zstart:g}.dat"
    )
    lcdm_v = mass_weighted_velocity(lcdm_transfer)
    model_v = mass_weighted_velocity(model_transfer)
    lcdm_v_on_model = interpolate_positive(
        lcdm_transfer[:, 0], lcdm_v, model_transfer[:, 0]
    )
    transfer_ratio = model_v / lcdm_v_on_model

    fundamental = 2.0 * np.pi / boxlen
    shell = np.rint(k / fundamental).astype(int)
    shell_rows = []
    for index in range(1, int(args.kmax / fundamental) + 1):
        mask = shell == index
        if np.count_nonzero(mask) < 10:
            continue
        weights = np.abs(theta_lcdm[mask]) ** 2
        if not np.all(np.isfinite(weights)) or weights.sum() <= 0.0:
            continue
        measured = float(
            np.sum(
                (theta_model[mask] * np.conj(theta_lcdm[mask])).real
            )
            / weights.sum()
        )
        theory_modes = interpolate_positive(
            model_transfer[:, 0], transfer_ratio, k[mask]
        )
        theory = float(np.sum(theory_modes * weights) / weights.sum())
        coherence = float(
            abs(np.sum(theta_model[mask] * np.conj(theta_lcdm[mask])))
            / math.sqrt(
                float(np.sum(np.abs(theta_model[mask]) ** 2) * weights.sum())
            )
        )
        shell_rows.append(
            {
                "shell": index,
                "k_weighted_h_mpc": float(
                    np.sum(k[mask] * weights) / weights.sum()
                ),
                "modes": int(np.count_nonzero(mask)),
                "measured_ratio": measured,
                "lagcamb_ratio": theory,
                "relative_residual": measured / theory - 1.0,
                "coherence": coherence,
            }
        )
    if not shell_rows:
        raise RuntimeError("no velocity shells were measured")
    velocity_residual = max(abs(row["relative_residual"]) for row in shell_rows)
    minimum_coherence = min(row["coherence"] for row in shell_rows)

    _, density_lcdm = read_grafic(lcdm_level / "ic_deltab")
    _, density_model = read_grafic(model_level / "ic_deltab")
    delta_lcdm = np.fft.rfftn(density_lcdm)
    delta_model = np.fft.rfftn(density_model)
    model_total = model_transfer[:, 6]
    lcdm_total_on_model = interpolate_positive(
        lcdm_transfer[:, 0], lcdm_transfer[:, 6], model_transfer[:, 0]
    )
    density_transfer_ratio = model_total / lcdm_total_on_model
    density_shell_rows = []
    for index in range(1, int(args.kmax / fundamental) + 1):
        mask = shell == index
        if np.count_nonzero(mask) < 10:
            continue
        weights = np.abs(delta_lcdm[mask]) ** 2
        if weights.sum() <= 0.0:
            continue
        measured = float(
            np.sum((delta_model[mask] * np.conj(delta_lcdm[mask])).real)
            / weights.sum()
        )
        theory_modes = interpolate_positive(
            model_transfer[:, 0], density_transfer_ratio, k[mask]
        )
        theory = float(np.sum(theory_modes * weights) / weights.sum())
        coherence = float(
            abs(np.sum(delta_model[mask] * np.conj(delta_lcdm[mask])))
            / math.sqrt(
                float(np.sum(np.abs(delta_model[mask]) ** 2) * weights.sum())
            )
        )
        density_shell_rows.append(
            {
                "shell": index,
                "modes": int(np.count_nonzero(mask)),
                "measured_ratio": measured,
                "lagcamb_ratio": theory,
                "relative_residual": measured / theory - 1.0,
                "coherence": coherence,
            }
        )
    density_residual = max(
        abs(row["relative_residual"]) for row in density_shell_rows
    )
    minimum_density_coherence = min(
        row["coherence"] for row in density_shell_rows
    )

    # The particle CIC field already contains the intended 2LPT displacement
    # and is therefore not expected to equal linear CAMB at sub-per-mille
    # precision. Retain that comparison as a diagnostic, while gating the
    # actual linear source field stored in ic_deltab above.
    model_pk = np.loadtxt(
        campaign / model / "output_00001" / "pk_cic_00001.dat"
    )
    lcdm_pk = np.loadtxt(
        campaign / "lcdm" / "output_00001" / "pk_cic_00001.dat"
    )
    if not np.array_equal(model_pk[:, 0], lcdm_pk[:, 0]):
        raise RuntimeError("model and LCDM CIC k grids differ")
    pk_mask = (
        (model_pk[:, 0] <= args.initial_pk_kmax)
        & (model_pk[:, 1] > 0.0)
        & (lcdm_pk[:, 1] > 0.0)
    )
    measured_pk_ratio = model_pk[pk_mask, 1] / lcdm_pk[pk_mask, 1]
    transfer_pk_ratio = (model_total / lcdm_total_on_model) ** 2
    theory_pk_ratio = interpolate_positive(
        model_transfer[:, 0], transfer_pk_ratio, model_pk[pk_mask, 0]
    )
    initial_pk_residual = float(
        np.max(np.abs(measured_pk_ratio / theory_pk_ratio - 1.0))
    )

    requested_a = [
        1.0 / (1.0 + zstart),
        *(1.0 / (1.0 + z) for z in metadata["output_redshifts"]),
    ]
    actual_a = exact_output_a(campaign / model)
    outputs_complete = len(actual_a) == len(requested_a) and all(
        math.isclose(got, want, rel_tol=0.0, abs_tol=1.0e-8)
        for got, want in zip(actual_a, requested_a)
    )
    log_path = campaign / model / "manual.out"
    log_text = log_path.read_text()
    clean_log = (
        "Run completed" in log_text
        and not re.search(r"(?i)(\\bnan\\b|\\babort\\b|\\bfatal\\b)", log_text)
    )

    checks = {
        "camb_vtotal_enabled": bool(
            diagnostics.get("dmo_velocity_transfer_enabled")
            and math.isclose(diagnostics["vfact_scale"], 1.0)
        ),
        "scale_dependence_detected": bool(
            diagnostics["velocity_growth_ratio_maximum_deviation"] > 1.0e-3
        ),
        "velocity_ratio_matches_lagcamb": bool(
            velocity_residual < args.velocity_ratio_rtol
        ),
        "velocity_phase_coherent": bool(minimum_coherence > 0.999),
        "initial_density_transfer_matches_lagcamb": bool(
            density_residual < args.initial_density_ratio_rtol
        ),
        "initial_density_phase_coherent": bool(
            minimum_density_coherence > 0.999
        ),
        "outputs_complete_through_z0": bool(outputs_complete),
        "lagramses_log_clean": bool(clean_log),
    }
    report = {
        "complete": all(checks.values()),
        "checks": checks,
        "metrics": {
            "maximum_velocity_ratio_relative_residual": velocity_residual,
            "minimum_velocity_shell_coherence": minimum_coherence,
            "maximum_initial_density_ratio_relative_residual": (
                density_residual
            ),
            "minimum_initial_density_shell_coherence": (
                minimum_density_coherence
            ),
            "initial_particle_2lpt_pk_ratio_residual_diagnostic": (
                initial_pk_residual
            ),
            "input_velocity_growth_scale_dependence": diagnostics[
                "velocity_growth_ratio_maximum_deviation"
            ],
        },
        "thresholds": {
            "velocity_ratio_rtol": args.velocity_ratio_rtol,
            "initial_density_ratio_rtol": args.initial_density_ratio_rtol,
            "velocity_kmax_h_mpc": args.kmax,
            "initial_pk_kmax_h_mpc": args.initial_pk_kmax,
        },
        "velocity_shells": shell_rows,
        "density_shells": density_shell_rows,
        "campaign": str(campaign),
        "model": model,
        "lagramses_log": str(log_path),
    }
    report_path = args.report or campaign / f"{model}_velocity_validation.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
