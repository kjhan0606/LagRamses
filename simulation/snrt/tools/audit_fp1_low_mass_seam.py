#!/usr/bin/env python3
"""Audit the F-P1 0.8--1.0 Msun lifetime seam without promoting yields.

The Huscher et al. candidate has single-star, lifetime-integrated ejecta at
the two seam endpoints.  Endpoint presence is not a lifetime model: the
production resolver needs an age convention and an age-resolved release
history (or an explicitly approved terminal lumping rule).  This audit makes
that distinction machine-checkable and intentionally emits no canonical row.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_FATE_MAP = SNRT_ROOT / "config" / "fp1_population_fate_map_v1.json"
DEFAULT_HUSCHER_CONTRACT = SNRT_ROOT / "config" / "g2_huscher2025_candidate_contract_v1.json"
DEFAULT_HUSCHER_AUDIT = SNRT_ROOT / "data" / "g2_huscher2025_candidate_audit.json"
DEFAULT_JSON_OUT = SNRT_ROOT / "data" / "fp1_low_mass_seam_review.json"
SEAM_ID = "low_mass_lifetime_seam"
SEAM = [0.8, 1.0]


class LowMassSeamAuditError(ValueError):
    """Review inputs are inconsistent or overclaim the lifetime seam."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise LowMassSeamAuditError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise LowMassSeamAuditError(f"{label} must be a JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise LowMassSeamAuditError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def audit_low_mass_seam(
    *,
    fate_map_path: Path = DEFAULT_FATE_MAP,
    huscher_contract_path: Path = DEFAULT_HUSCHER_CONTRACT,
    huscher_audit_path: Path = DEFAULT_HUSCHER_AUDIT,
) -> dict[str, Any]:
    fate_map_path = Path(fate_map_path).resolve()
    huscher_contract_path = Path(huscher_contract_path).resolve()
    huscher_audit_path = Path(huscher_audit_path).resolve()
    fate_map = _read_json(fate_map_path, "F-P1 fate map")
    contract = _read_json(huscher_contract_path, "Huscher candidate contract")
    candidate_audit = _read_json(huscher_audit_path, "Huscher candidate audit")

    intervals = fate_map.get("intervals")
    if not isinstance(intervals, list):
        raise LowMassSeamAuditError("F-P1 fate map intervals are missing")
    seam = next((item for item in intervals if isinstance(item, dict) and item.get("id") == SEAM_ID), None)
    if seam is None or seam.get("mass_msun") != SEAM:
        raise LowMassSeamAuditError("F-P1 low-mass lifetime seam is missing or changed")
    if seam.get("fate_class") != "unresolved":
        raise LowMassSeamAuditError("low-mass lifetime seam must remain unresolved until approved")

    source = contract.get("source")
    grid = contract.get("single_star_grid")
    projection = contract.get("canonical_projection_policy")
    approval = contract.get("approval")
    if not all(isinstance(item, dict) for item in (source, grid, projection, approval)):
        raise LowMassSeamAuditError("Huscher candidate contract is incomplete")
    masses = grid.get("mass_msun")
    if not isinstance(masses, list) or not all(isinstance(value, (int, float)) for value in masses):
        raise LowMassSeamAuditError("Huscher single-star mass grid is malformed")
    endpoint_coverage = {str(mass): mass in masses for mass in SEAM}
    if not all(endpoint_coverage.values()):
        raise LowMassSeamAuditError("Huscher candidate does not cover both seam endpoints")

    release_semantics = grid.get("yield_column_semantics")
    if release_semantics != "gross isotope mass ejected over the stellar lifetime in Msun per initial star":
        raise LowMassSeamAuditError("Huscher lifetime-integrated yield semantics changed")
    population_tables = contract.get("population_tables")
    if not isinstance(population_tables, dict):
        raise LowMassSeamAuditError("Huscher population-table contract is missing")
    runtime_use = population_tables.get("runtime_use")
    if runtime_use != "independent population-level validation only; never convolve a second IMF or mix with per-star rows":
        raise LowMassSeamAuditError("Huscher population-table firewall changed")

    audit_blockers = candidate_audit.get("blockers")
    required_blocker = "single_star_files_are_lifetime_integrated_and_have_no_per_star_age_resolved_release_history"
    if not isinstance(audit_blockers, list) or required_blocker not in audit_blockers:
        raise LowMassSeamAuditError("Huscher audit no longer records the missing age-resolved history")
    if projection.get("canonical_rows_emitted") != 0 or approval.get("canonical_conversion_allowed") is not False:
        raise LowMassSeamAuditError("Huscher candidate unexpectedly permits canonical conversion")

    return {
        "schema": "snrt-fp1-low-mass-lifetime-seam-review",
        "schema_version": 1,
        "gate": "F-P1",
        "status": "review_only_candidate_covers_endpoints_lifetime_unresolved",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "runtime_activation_allowed": False,
        "seam": {
            "id": SEAM_ID,
            "mass_msun": SEAM,
            "fate_map_status": seam.get("evidence_status"),
            "endpoint_coverage": endpoint_coverage,
        },
        "candidate": {
            "candidate_id": source.get("candidate_id"),
            "article_doi": source.get("article_doi"),
            "data_doi": source.get("data_doi"),
            "license": source.get("license"),
            "contract": {"path": str(huscher_contract_path), "sha256": _sha256(huscher_contract_path)},
            "audit": {"path": str(huscher_audit_path), "sha256": _sha256(huscher_audit_path)},
            "mass_grid_msun": masses,
            "yield_semantics": release_semantics,
        },
        "resolved": False,
        "blockers": [
            "no_approved_lifetime_source_or_age_resolved_release_history",
            "lifetime_integrated_ejecta_are_not_an_age_resolved_source",
            "terminal_fate_and_remnant_ownership_are_not_approved",
            "Huscher_population_tables_are_IMF_weighted_and_second_convolution_is_forbidden",
        ],
        "required_next_inputs": [
            "immutable_lifetime_source_and_age_convention",
            "per_star_or_explicitly_approved_terminal_lumping_release_history",
            "non_overlapping_AGB_wind_and_terminal_remnant_ownership",
            "named_physics_approval_id_and_source_package_fingerprint",
        ],
        "interpretation": (
            "The candidate reaches both 0.8 and 1.0 Msun endpoints, but its files contain "
            "lifetime-integrated ejecta rather than an approved age-resolved lifetime/fate map. "
            "The F-P1 seam therefore remains unresolved and contributes no canonical feedback."
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fate-map", type=Path, default=DEFAULT_FATE_MAP)
    parser.add_argument("--huscher-contract", type=Path, default=DEFAULT_HUSCHER_CONTRACT)
    parser.add_argument("--huscher-audit", type=Path, default=DEFAULT_HUSCHER_AUDIT)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    args = parser.parse_args(argv)
    try:
        report = audit_low_mass_seam(
            fate_map_path=args.fate_map,
            huscher_contract_path=args.huscher_contract,
            huscher_audit_path=args.huscher_audit,
        )
    except LowMassSeamAuditError as exc:
        print(f"F-P1 low-mass seam audit ERROR: {exc}", file=sys.stderr)
        return 2
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("F-P1 low-mass seam: review_only_lifetime_unresolved")
    print("endpoint_coverage=0.8,1.0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
