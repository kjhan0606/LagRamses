#!/usr/bin/env python3
"""Compare review-only shell-boundary estimators for staged HESMA profiles.

HESMA supplies density values at shell-centre velocities, but not the exact
shell edges needed to reproduce a mass integral from the spherical average.
This tool makes that numerical ambiguity explicit by comparing two documented
edge policies. Neither result is promoted to an event mass or energy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any
from zipfile import BadZipFile, ZipFile


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_ROOT = PROJECT_ROOT / "assets" / "review_only" / "fp2_snia" / "hesma_yysd4_xap92"
DEFAULT_MANIFEST = PROJECT_ROOT / "manifests" / "fp2_snia_hesma_yysd4_review_v1.json"
DEFAULT_AUDIT = SNRT_ROOT / "data" / "fp2_snia_hesma_source_audit.json"
SOLAR_MASS_G = 1.98847e33
KM_S_TO_CM_S = 1.0e5
ESTIMATOR_NAMES = (
    "inner_zero_outer_half_bin",
    "half_bin_extrapolated_both_ends",
)


class ProfileEstimatorError(ValueError):
    """The review-only profile estimator input is malformed."""


def _sha256(path: Path) -> str:
    if not path.is_file():
        raise ProfileEstimatorError(f"file is missing: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProfileEstimatorError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ProfileEstimatorError(f"JSON object required: {path}")
    return value


def _read_profile(archive: Path, member: str) -> list[tuple[float, float]]:
    try:
        with ZipFile(archive) as zfile:
            text = zfile.read(member).decode("utf-8")
    except (BadZipFile, KeyError, OSError, UnicodeDecodeError) as exc:
        raise ProfileEstimatorError(f"cannot read density profile {member}") from exc
    rows: list[tuple[float, float]] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.strip().startswith("Velocity [km/s]"):
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ProfileEstimatorError(f"density profile {member} line {line_number} is not 2-column")
        try:
            velocity_km_s, density = (float(value) for value in fields)
        except ValueError as exc:
            raise ProfileEstimatorError(f"density profile {member} line {line_number} is non-numeric") from exc
        if not math.isfinite(velocity_km_s) or not math.isfinite(density) or velocity_km_s < 0.0 or density < 0.0:
            raise ProfileEstimatorError(f"density profile {member} line {line_number} is invalid")
        rows.append((velocity_km_s, density))
    if len(rows) < 2 or any(rows[index][0] >= rows[index + 1][0] for index in range(len(rows) - 1)):
        raise ProfileEstimatorError(f"density profile {member} has an invalid velocity grid")
    return rows


def _edges(velocities: list[float], estimator: str) -> list[float]:
    if estimator not in ESTIMATOR_NAMES:
        raise ProfileEstimatorError(f"unknown profile estimator: {estimator}")
    spacings = [velocities[index + 1] - velocities[index] for index in range(len(velocities) - 1)]
    edges = [(velocities[index] + velocities[index + 1]) / 2.0 for index in range(len(velocities) - 1)]
    if estimator == "inner_zero_outer_half_bin":
        return [0.0, *edges, velocities[-1] + spacings[-1] / 2.0]
    return [max(0.0, velocities[0] - spacings[0] / 2.0), *edges, velocities[-1] + spacings[-1] / 2.0]


def _integrate(rows: list[tuple[float, float]], epoch_s: float, estimator: str) -> dict[str, Any]:
    velocities_cm_s = [velocity * KM_S_TO_CM_S for velocity, _ in rows]
    edges_cm_s = _edges(velocities_cm_s, estimator)
    mass_g = 0.0
    energy_erg = 0.0
    for index, (_, density) in enumerate(rows):
        shell_volume = 4.0 * math.pi * epoch_s**3 * (
            edges_cm_s[index + 1] ** 3 - edges_cm_s[index] ** 3
        ) / 3.0
        shell_mass = density * shell_volume
        mass_g += shell_mass
        energy_erg += 0.5 * shell_mass * velocities_cm_s[index] ** 2
    return {
        "estimator": estimator,
        "mass_estimate_msun": mass_g / SOLAR_MASS_G,
        "kinetic_energy_estimate_erg": energy_erg,
        "edge_policy": (
            "inner edge fixed at zero; interior midpoints; outer edge half-bin extrapolated"
            if estimator == "inner_zero_outer_half_bin"
            else "inner and outer edges half-bin extrapolated; interior midpoints; inner edge clamped at zero"
        ),
        "event_scalar_approval_status": "review_only_not_authoritative",
    }


def compare_estimators(
    *,
    root: Path = DEFAULT_ROOT,
    manifest_path: Path = DEFAULT_MANIFEST,
    audit_path: Path = DEFAULT_AUDIT,
) -> dict[str, Any]:
    root = Path(root).resolve()
    manifest_path = Path(manifest_path).resolve()
    audit_path = Path(audit_path).resolve()
    manifest = _read_json(manifest_path)
    audit = _read_json(audit_path)
    if audit.get("status") != "review_only_source_format_passed":
        raise ProfileEstimatorError("HESMA source audit is not clean")
    if audit.get("package_sha256") != manifest.get("package_sha256"):
        raise ProfileEstimatorError("HESMA audit and manifest package hashes disagree")

    rows: list[dict[str, Any]] = []
    density_archive = root / "ddt_2013_density.zip"
    for model_id in audit["model_ids"]:
        model_report = audit["model_reports"][model_id]
        member = model_report["density"]["path"]
        profile_rows = _read_profile(density_archive, member)
        epoch_s = model_report["density"]["profile_epoch_s"]
        integrated_mass = model_report["abundances"]["stable_element_mass_msun"]
        estimates = {
            estimator: _integrate(
                profile_rows,
                epoch_s,
                estimator,
            )
            for estimator in ESTIMATOR_NAMES
        }
        for estimate in estimates.values():
            estimate["relative_difference_from_integrated_stable_mass"] = abs(
                estimate["mass_estimate_msun"] - integrated_mass
            ) / max(integrated_mass, 1.0e-30)
        rows.append(
            {
                "model_id": model_id,
                "density_member": member,
                "row_count": len(profile_rows),
                "profile_epoch_s": epoch_s,
                "integrated_stable_element_mass_msun": integrated_mass,
                "estimators": estimates,
                "source_audit_physical_warnings": [
                    warning for warning in audit.get("physical_warnings", [])
                    if warning.get("model") == model_id
                ],
            }
        )

    return {
        "schema": "snrt-fp2-snia-hesma-profile-estimator-comparison",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "review_only_diagnostic",
        "source": {
            "record_id": manifest["source_record_id"],
            "package_sha256": manifest["package_sha256"],
            "manifest_sha256": _sha256(manifest_path),
            "source_audit_sha256": _sha256(audit_path),
        },
        "interpretation": {
            "purpose": "quantify shell-edge sensitivity of spherical-average profile diagnostics",
            "selection": "no estimator selected for production",
            "mass_and_energy": "diagnostic estimates only; not authoritative event fields",
            "momentum": "signed vector remains unavailable from spherical-average profiles",
            "source_epoch": "use exact source-reported per-model epoch from the HESMA record",
        },
        "model_count": len(rows),
        "estimators": list(ESTIMATOR_NAMES),
        "models": rows,
        "admission": {
            "canonical_conversion_allowed": False,
            "runtime_activation_allowed": False,
            "selected_estimator": None,
        },
        "tool_sha256": _sha256(TOOL_PATH),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    try:
        report = compare_estimators(root=args.root, manifest_path=args.manifest, audit_path=args.audit)
        payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(payload, encoding="utf-8")
    except (ProfileEstimatorError, OSError) as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, indent=2), file=sys.stderr)
        return 2
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
