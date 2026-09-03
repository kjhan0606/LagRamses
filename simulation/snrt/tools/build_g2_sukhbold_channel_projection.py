#!/usr/bin/env python3
"""Build a lossless review-only component projection of Sukhbold 2016."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

from audit_g2_sukhbold2016_candidate import (
    DEFAULT_ROOT,
    SukhboldAuditError,
    audit_sukhbold2016_candidate,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_sukhbold_channel_projection_contract_v1.json"
_EXPECTED_FIREWALL_KEYS = (
    "component_channel_assignment_approved",
    "wind_and_terminal_may_be_summed_before_channel_ownership_approval",
    "stable_and_radioactive_segments_may_be_summed",
    "selected_radioactive_inventory_is_decay_complete",
    "stable_component_sum_is_exact_total_returned_mass",
    "source_age_history_may_be_invented",
    "source_lifetime_may_be_invented",
    "exact_launch_momentum_may_be_derived",
    "cross_source_interpolation_allowed",
    "out_of_domain_extrapolation_allowed",
)
_EXPECTED_COMPONENT_CONTRACT = {
    "wind": ("wind", 1, "integrated presupernova wind only"),
    "terminal_ccsn": ("ejecta", 3, "terminal supernova ejecta only"),
}


class SukhboldProjectionError(ValueError):
    """The Sukhbold component projection violates its fail-closed contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise SukhboldProjectionError(f"cannot read {path}: {exc}") from exc
    return digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SukhboldProjectionError(f"cannot read projection contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-sukhbold-channel-projection-contract"
        or contract.get("schema_version") != 1
    ):
        raise SukhboldProjectionError("unsupported Sukhbold projection contract")
    firewalls = contract.get("semantic_firewalls", {})
    if set(firewalls) != set(_EXPECTED_FIREWALL_KEYS) or any(
        firewalls.get(key) is not False for key in _EXPECTED_FIREWALL_KEYS
    ):
        raise SukhboldProjectionError("projection semantic firewalls are not fail-closed")
    approval = contract.get("approval", {})
    if approval.get("canonical_rows_emitted") != 0:
        raise SukhboldProjectionError("projection unexpectedly emits canonical rows")
    if approval.get("canonical_conversion_allowed") is not False:
        raise SukhboldProjectionError("projection unexpectedly permits conversion")
    if approval.get("runtime_deposition_allowed") is not False:
        raise SukhboldProjectionError("projection unexpectedly permits deposition")
    required = approval.get("required_before_approval")
    if not isinstance(required, list) or not required or any(
        not isinstance(value, str) or not value for value in required
    ):
        raise SukhboldProjectionError("projection approval prerequisites are missing")
    components = contract.get("component_projection")
    if set(components or {}) != set(_EXPECTED_COMPONENT_CONTRACT):
        raise SukhboldProjectionError("projection component map is incomplete or has unknown entries")
    return contract


def _component_record(
    *,
    mass: float,
    source: dict[str, Any],
    component: str,
    component_contract: dict[str, Any],
    tracked_elements: list[str],
    engine: str = "Z9.6",
    source_branch: str = "Z9.6",
) -> dict[str, Any]:
    expected_contract = _EXPECTED_COMPONENT_CONTRACT.get(component)
    if expected_contract is None:
        raise SukhboldProjectionError(f"unsupported source component: {component}")
    expected_source_column, expected_channel, expected_ownership = expected_contract
    actual_contract = (
        component_contract.get("source_column"),
        component_contract.get("proposed_runtime_channel"),
        component_contract.get("ownership"),
    )
    if actual_contract != expected_contract:
        raise SukhboldProjectionError(
            f"{component}: contract identity {actual_contract!r} does not match "
            f"{expected_contract!r}"
        )
    is_wind = component == "wind"
    prefix = expected_source_column
    tracked = source[f"stable_{prefix}_by_tracked_element_msun"]
    if list(tracked) != tracked_elements:
        raise SukhboldProjectionError(f"s{mass}: tracked-element order drifted")
    stable_mass = float(source[f"stable_{prefix}_sum_msun"])
    untracked = float(source[f"untracked_stable_{prefix}_msun"])
    if not math.isfinite(stable_mass) or not math.isfinite(untracked):
        raise SukhboldProjectionError(f"s{mass}: non-finite stable source component")
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tracked.values()):
        raise SukhboldProjectionError(f"s{mass}: invalid tracked stable source component")
    radioactive = {
        isotope: float(values[f"{prefix}_msun"])
        for isotope, values in source["selected_radioactive_inventory"].items()
    }
    if len(radioactive) != 20 or any(value < 0.0 for value in radioactive.values()):
        raise SukhboldProjectionError(f"s{mass}: radioactive sidecar drifted")
    return {
        "zams_mass_msun": mass,
        "metallicity": "source_labelled_solar",
        "engine": engine,
        "source_branch": source_branch,
        "source_component": expected_source_column,
        "proposed_runtime_channel": expected_channel,
        "component_ownership": expected_ownership,
        "source_lifetime_yr": None,
        "release_age_yr": None,
        "stable_component_mass_msun": stable_mass,
        "stable_mass_by_tracked_element_msun": tracked,
        "untracked_stable_component_mass_msun": untracked,
        "source_stable_segment_mass_budget_residual_msun": source.get(
            "stable_segment_mass_budget_residual_msun"
        ),
        "source_stable_segment_mass_budget_exact_closure_claimed": source.get(
            "exact_mass_closure_claimed", False
        ),
        "source_cross_segment_duplicate_isotopes": source.get(
            "cross_segment_duplicate_isotopes", []
        ),
        "selected_radioactive_inventory_msun": radioactive,
        "selected_radioactive_inventory_sum_msun": math.fsum(radioactive.values()),
        "decay_complete_returned_mass_msun": None,
        "final_kinetic_energy_erg": None if is_wind else source["final_kinetic_energy_erg"],
        "baryonic_mass_cut_after_fallback_msun": None if is_wind else source["baryonic_mass_cut_after_fallback_msun"],
        "fallback_mass_msun": None if is_wind else source["fallback_mass_msun"],
        "canonical_scalar_launch_momentum_g_cm_s": None,
        "canonical_source_frame_vector_momentum_g_cm_s": None,
        "canonical_row_emitted": False,
        "runtime_channel_assignment_approved": False,
    }


def build_sukhbold_channel_projection(
    *,
    root: Path = DEFAULT_ROOT,
    contract_path: Path = DEFAULT_CONTRACT,
    include_records: bool = True,
) -> dict[str, Any]:
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    source = audit_sukhbold2016_candidate(root=Path(root))
    if source["canonical_rows_emitted"] != 0 or source["production_ready"] is not False:
        raise SukhboldProjectionError("source audit unexpectedly permits production use")
    if source["source_identity"]["candidate_id"] != contract["source_candidate_id"]:
        raise SukhboldProjectionError("source candidate identity drifted")
    tracked_elements = list(contract["tracked_elements"])
    components = contract["component_projection"]
    records: list[dict[str, Any]] = []
    grid = source["z96_grid"]
    for mass in grid["zams_mass_msun"]:
        model = grid["models"][str(mass)]
        for component in ("wind", "terminal_ccsn"):
            records.append(
                _component_record(
                    mass=mass,
                    source=model,
                    component=component,
                    component_contract=components[component],
                    tracked_elements=tracked_elements,
                )
            )
    high_mass_records: list[dict[str, Any]] = []
    for engine, engine_report in source["high_mass_engine_evidence"]["engines"].items():
        for mass_key, model in engine_report.get("high_mass_yields", {}).items():
            mass = float(mass_key)
            for component in ("wind", "terminal_ccsn"):
                high_mass_records.append(
                    _component_record(
                        mass=mass,
                        source=model,
                        component=component,
                        component_contract=components[component],
                        tracked_elements=tracked_elements,
                        engine=engine,
                        source_branch=engine,
                    )
                )
        for mass_key, model in engine_report.get("high_mass_implosion_winds", {}).items():
            high_mass_records.append(
                _component_record(
                    mass=float(mass_key),
                    source=model,
                    component="wind",
                    component_contract=components["wind"],
                    tracked_elements=tracked_elements,
                    engine=engine,
                    source_branch="implosions_W18",
                )
            )
    report: dict[str, Any] = {
        "schema": "snrt-g2-sukhbold-channel-projection-review",
        "schema_version": 1,
        "gate": "G2",
        "status": "review_only_blocked_decay_age_boundary_and_approval",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "runtime_deposition_allowed": False,
        "canonical_rows_emitted": 0,
        "source_candidate_id": source["source_identity"]["candidate_id"],
        "source_identity": source["source_identity"],
        "model_count": grid["model_count"],
        "record_count": len(records),
        "high_mass_record_count": len(high_mass_records),
        "record_count_by_source_component": {
            "wind": sum(record["source_component"] == "wind" for record in records),
            "ejecta": sum(record["source_component"] == "ejecta" for record in records),
        },
        "high_mass_record_count_by_source_component": {
            "wind": sum(record["source_component"] == "wind" for record in high_mass_records),
            "ejecta": sum(record["source_component"] == "ejecta" for record in high_mass_records),
        },
        "high_mass_record_count_by_engine_and_source_component": {
            engine: {
                "wind": sum(
                    record["engine"] == engine and record["source_component"] == "wind"
                    for record in high_mass_records
                ),
                "ejecta": sum(
                    record["engine"] == engine and record["source_component"] == "ejecta"
                    for record in high_mass_records
                ),
            }
            for engine in source["high_mass_engine_evidence"]["engines"]
        },
        "high_mass_missing_source_masses_by_engine": {
            engine: engine_report.get("missing_high_mass_yield_masses", [])
            for engine, engine_report in source["high_mass_engine_evidence"]["engines"].items()
        },
        "tracked_elements": tracked_elements,
        "source_nulls_preserved": {
            "lifetime_and_release_age": True,
            "decay_complete_returned_mass": True,
            "canonical_launch_momentum": True,
        },
        "semantic_firewalls": contract["semantic_firewalls"],
        "blockers": contract["approval"]["required_before_approval"],
        "contract_path": str(contract_path),
        "contract_sha256": _sha256(contract_path),
        "source_audit_code_sha256": source["audit_code_sha256"],
        "builder_code_sha256": _sha256(TOOL_PATH),
    }
    if include_records:
        report["records"] = records
        report["high_mass_records"] = high_mass_records
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--summary-only", action="store_true")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = build_sukhbold_channel_projection(
            root=args.root,
            contract_path=args.contract,
            include_records=not args.summary_only,
        )
    except (SukhboldProjectionError, SukhboldAuditError) as exc:
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
