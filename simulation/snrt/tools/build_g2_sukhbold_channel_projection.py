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
    if not firewalls or any(value is not False for value in firewalls.values()):
        raise SukhboldProjectionError("projection semantic firewalls are not fail-closed")
    approval = contract.get("approval", {})
    if approval.get("canonical_rows_emitted") != 0:
        raise SukhboldProjectionError("projection unexpectedly emits canonical rows")
    if approval.get("canonical_conversion_allowed") is not False:
        raise SukhboldProjectionError("projection unexpectedly permits conversion")
    if approval.get("runtime_deposition_allowed") is not False:
        raise SukhboldProjectionError("projection unexpectedly permits deposition")
    return contract


def _component_record(
    *,
    mass: float,
    source: dict[str, Any],
    component: str,
    component_contract: dict[str, Any],
    tracked_elements: list[str],
) -> dict[str, Any]:
    is_wind = component == "wind"
    prefix = "wind" if is_wind else "ejecta"
    tracked = source[f"stable_{prefix}_by_tracked_element_msun"]
    if list(tracked) != tracked_elements:
        raise SukhboldProjectionError(f"s{mass}: tracked-element order drifted")
    stable_mass = float(source[f"stable_{prefix}_sum_msun"])
    untracked = float(source[f"untracked_stable_{prefix}_msun"])
    if not math.isclose(
        math.fsum(float(value) for value in tracked.values()) + untracked,
        stable_mass,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise SukhboldProjectionError(f"s{mass}: reduced stable vector does not close")
    radioactive = {
        isotope: float(values[f"{prefix}_msun"])
        for isotope, values in source["selected_radioactive_inventory"].items()
    }
    if len(radioactive) != 20 or any(value < 0.0 for value in radioactive.values()):
        raise SukhboldProjectionError(f"s{mass}: radioactive sidecar drifted")
    return {
        "zams_mass_msun": mass,
        "metallicity": "source_labelled_solar",
        "engine": "Z9.6",
        "source_component": component_contract["source_column"],
        "proposed_runtime_channel": component_contract["proposed_runtime_channel"],
        "component_ownership": component_contract["ownership"],
        "source_lifetime_yr": None,
        "release_age_yr": None,
        "stable_component_mass_msun": stable_mass,
        "stable_mass_by_tracked_element_msun": tracked,
        "untracked_stable_component_mass_msun": untracked,
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
        "model_count": grid["model_count"],
        "record_count": len(records),
        "record_count_by_source_component": {
            "wind": sum(record["source_component"] == "wind" for record in records),
            "ejecta": sum(record["source_component"] == "ejecta" for record in records),
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
