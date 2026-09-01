#!/usr/bin/env python3
"""Audit the Stockinger et al. (2020) low-mass CCSN release."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable

import h5py
import numpy as np


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_stockinger2020_candidate_contract_v1.json"
_TRACKED_ELEMENTS = {"H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"}


class StockingerAuditError(ValueError):
    """The staged Stockinger release violates its review contract."""


def _hash(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise StockingerAuditError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StockingerAuditError(f"cannot read contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-stockinger2020-candidate-contract"
        or contract.get("schema_version") != 1
    ):
        raise StockingerAuditError("unsupported Stockinger candidate contract")
    policy = contract.get("audit_policy", {})
    approval = contract.get("approval", {})
    if policy.get("canonical_rows_emitted") != 0:
        raise StockingerAuditError("review contract unexpectedly emits canonical rows")
    if policy.get("vsh_dataset_use_allowed") is not False:
        raise StockingerAuditError("known-bad vsh metadata must be quarantined")
    if policy.get("terminal_yields_as_stable_decay_products_allowed") is not False:
        raise StockingerAuditError("radioactive event yields cannot masquerade as stable yields")
    if approval.get("canonical_conversion_allowed") is not False:
        raise StockingerAuditError("review contract unexpectedly permits conversion")
    return contract


def _finite(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise StockingerAuditError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise StockingerAuditError(f"{field} is not finite: {token!r}")
    return value


def _parse_yield(path: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    total: tuple[float, float] | None = None
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise StockingerAuditError(f"cannot read {path}: {exc}") from exc
    for line_number, raw in enumerate(lines, start=1):
        fields = raw.split()
        if len(fields) != 3 or raw.lstrip().startswith("#"):
            continue
        values = (_finite(fields[1], "grid mass"), _finite(fields[2], "outflow mass"))
        if any(value < 0.0 for value in values):
            raise StockingerAuditError(f"{path.name}:{line_number}: negative yield")
        if fields[0] == "total":
            total = values
        else:
            rows.append(
                {
                    "species": fields[0],
                    "grid_mass_msun": values[0],
                    "positive_radial_velocity_mass_msun": values[1],
                }
            )
    if total is None or not rows:
        raise StockingerAuditError(f"malformed yield table: {path}")
    names = [row["species"] for row in rows]
    if len(names) != len(set(names)):
        raise StockingerAuditError(f"duplicate species in {path}")
    sums = (
        sum(row["grid_mass_msun"] for row in rows),
        sum(row["positive_radial_velocity_mass_msun"] for row in rows),
    )
    return {
        "file": path.name,
        "species": names,
        "species_count": len(rows),
        "reported_total_grid_mass_msun": total[0],
        "reported_total_positive_radial_velocity_mass_msun": total[1],
        "species_sum_grid_mass_msun": sums[0],
        "species_sum_positive_radial_velocity_mass_msun": sums[1],
        "grid_mass_closure_residual_msun": sums[0] - total[0],
        "positive_radial_velocity_mass_closure_residual_msun": sums[1] - total[1],
    }


def _element_for_species(species: str) -> str | None:
    if species == "p":
        return "H"
    if species == "n" or species == "Tr":
        return None
    match = re.match(r"([A-Z][a-z]?)", species)
    return match.group(1) if match else None


def _decode_attribute(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _audit_energy(path: Path, contract: dict[str, Any]) -> dict[str, Any]:
    models: dict[str, Any] = {}
    try:
        with h5py.File(path, "r") as handle:
            time = np.asarray(handle["time"][:], dtype=float)
            time_unit = _decode_attribute(handle["time"].attrs.get("description"))
            if time.ndim != 1 or len(time) != 300 or not np.all(np.isfinite(time)):
                raise StockingerAuditError("energy time axis is malformed")
            if np.any(np.diff(time) <= 0.0) or time_unit != "Postbounce time in s":
                raise StockingerAuditError("energy time semantics drifted")
            for model_name, model_contract in contract["models"].items():
                group_name = model_contract["hdf5_group"]
                if group_name not in handle:
                    raise StockingerAuditError(f"missing HDF5 model group: {group_name}")
                group = handle[group_name]
                energy = np.asarray(group["energy"][:], dtype=float)
                if energy.shape != time.shape or np.any(np.isinf(energy)):
                    raise StockingerAuditError(f"{group_name}: invalid energy array")
                valid = np.flatnonzero(np.isfinite(energy))
                if len(valid) == 0 or np.any(energy[valid] < 0.0):
                    raise StockingerAuditError(f"{group_name}: no valid diagnostic energy")
                energy_unit = _decode_attribute(group["energy"].attrs.get("description"))
                if energy_unit != contract["audit_policy"]["energy_unit_attribute"]:
                    raise StockingerAuditError(f"{group_name}: energy unit drifted")
                last = int(valid[-1])
                maximum = int(valid[np.argmax(energy[valid])])
                vsh_attribute = _decode_attribute(group["vsh"].attrs.get("description"))
                if vsh_attribute != contract["known_release_findings"]["vsh_attribute_value"]:
                    raise StockingerAuditError(f"{group_name}: expected vsh metadata finding drifted")
                models[model_name] = {
                    "hdf5_group": group_name,
                    "finite_sample_count": int(len(valid)),
                    "last_finite_time_postbounce_s": float(time[last]),
                    "last_finite_diagnostic_energy_1e50_erg": float(energy[last]),
                    "last_finite_diagnostic_energy_erg": float(energy[last] * 1.0e50),
                    "maximum_time_postbounce_s": float(time[maximum]),
                    "maximum_diagnostic_energy_1e50_erg": float(energy[maximum]),
                    "vsh_attribute": vsh_attribute,
                    "vsh_used": False,
                }
    except (OSError, KeyError) as exc:
        raise StockingerAuditError(f"cannot audit energy HDF5 {path}: {exc}") from exc
    return {
        "time_sample_count": 300,
        "time_minimum_postbounce_s": 0.0,
        "time_maximum_postbounce_s": 3.14,
        "models": models,
        "vsh_metadata_semantics_pass": False,
        "vsh_quarantined": True,
    }


def audit_stockinger2020_candidate(
    *, root: Path = DEFAULT_ROOT, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    root = Path(root).resolve()
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    source = contract["source"]
    base = root / source["release_root_relative_path"]

    fingerprints: dict[str, dict[str, Any]] = {}
    for name, expected in source["files"].items():
        size, sha256 = _hash(base / name)
        if size != expected["bytes"] or sha256 != expected["sha256"]:
            raise StockingerAuditError(f"staged Stockinger source fingerprint drifted: {name}")
        fingerprints[name] = {"bytes": size, "sha256": sha256}

    snapshots: dict[str, Any] = {}
    final_reports: dict[str, Any] = {}
    max_closure = 0.0
    for model_name, model_contract in contract["models"].items():
        model_snapshots = []
        for days in model_contract["yield_snapshot_days"]:
            matches = sorted(base.glob(f"{model_name}_yields_{days:07.3f}d.txt"))
            if len(matches) != 1:
                raise StockingerAuditError(
                    f"{model_name}: expected one yield snapshot for {days} d, found {len(matches)}"
                )
            report = _parse_yield(matches[0])
            if report["species_count"] != model_contract["expected_species_count"]:
                raise StockingerAuditError(f"{model_name}: species count drifted")
            for field in (
                "grid_mass_closure_residual_msun",
                "positive_radial_velocity_mass_closure_residual_msun",
            ):
                max_closure = max(max_closure, abs(report[field]))
            model_snapshots.append({"time_days": days, **report})
        if [item["species"] for item in model_snapshots] != [
            model_snapshots[0]["species"]
        ] * len(model_snapshots):
            raise StockingerAuditError(f"{model_name}: species sequence drifted over snapshots")
        snapshots[model_name] = model_snapshots
        final = next(
            report for report in model_snapshots
            if report["file"] == model_contract["final_yield_file"]
        )
        elements = {
            element for species in final["species"]
            if (element := _element_for_species(species)) in _TRACKED_ELEMENTS
        }
        final_reports[model_name] = {
            "zams_mass_msun": model_contract["zams_mass_msun"],
            "metallicity": model_contract["metallicity"],
            "collapse_type": model_contract["collapse_type"],
            "energy_model": model_contract["energy_model"],
            "final_yield_file": final["file"],
            "final_time_days": next(
                item["time_days"] for item in model_snapshots if item["file"] == final["file"]
            ),
            "final_total_ejecta_msun": final["reported_total_grid_mass_msun"],
            "tracked_elements_present": sorted(elements),
            "tracked_elements_absent": sorted(_TRACKED_ELEMENTS - elements),
            "radioactive_or_unresolved_tracer_species": [
                species for species in final["species"]
                if species in {"Ti44", "Cr48", "Fe52", "Ni56", "Co56", "Tr"}
            ],
        }
    tolerance = contract["audit_policy"]["yield_total_closure_tolerance_msun"]
    if max_closure > tolerance:
        raise StockingerAuditError("yield species sum does not close to reported total")
    if any(report["tracked_elements_absent"] != ["N"] for report in final_reports.values()):
        raise StockingerAuditError("tracked-element coverage drifted")

    energy = _audit_energy(base / "shock_energy_kick.h5", contract)
    return {
        "schema": "snrt-g2-stockinger2020-candidate-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "candidate_acquired_energy_yields_audited_license_unresolved",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "source_identity": {
            "candidate_id": source["candidate_id"],
            "article_doi": source["article_doi"],
            "source_page": source["source_page"],
            "research_access_verified": True,
            "redistribution_license_verified": False,
            "file_count": len(fingerprints),
            "files": fingerprints,
        },
        "model_grid": {
            "model_count": len(final_reports),
            "zams_mass_msun": sorted(report["zams_mass_msun"] for report in final_reports.values()),
            "models": final_reports,
            "cross_model_interpolation_allowed": False,
        },
        "yield_snapshots": snapshots,
        "yield_mass_closure": {
            "maximum_absolute_species_sum_residual_msun": max_closure,
            "tolerance_msun": tolerance,
            "pass": True,
        },
        "diagnostic_explosion_energy": energy,
        "quality_findings": {
            "nitrogen_missing_from_all_event_yield_vectors": True,
            "radioactive_decay_projection_complete": False,
            "vsh_metadata_semantics_pass": False,
            "vsh_dataset_used": False,
            "precollapse_wind_in_event_yield_tables": False,
        },
        "semantic_firewalls": {
            "precollapse_wind_and_terminal_event_double_count_forbidden": True,
            "event_snapshots_as_stellar_age_history_allowed": False,
            "cross_model_interpolation_allowed": False,
            "terminal_yields_as_stable_decay_products_allowed": False,
        },
        "blockers": [
            "no_explicit_redistribution_license_identified",
            "three_models_do_not_define_a_continuous_8_to_11_msun_grid",
            "solar_and_zero_metallicity_models_cannot_be_interpolated_as_one_sequence",
            "nitrogen_missing_from_reduced_chemistry_vector",
            "radioactive_decay_and_tracer_projection_not_approved",
            "no_canonical_terminal_momentum_field",
            "source_selection_and_cross_source_seam_not_approved",
        ],
        "contract_path": str(contract_path),
        "contract_sha256": _hash(contract_path)[1],
        "audit_code_sha256": _hash(TOOL_PATH)[1],
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_stockinger2020_candidate(root=args.root, contract_path=args.contract)
    except StockingerAuditError as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, indent=2), file=sys.stderr)
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
