#!/usr/bin/env python3
"""Audit the mass-only wind history recoverable from Limongi table5."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

from adapt_g2_candidate_sources import (
    DEFAULT_ROOT,
    LIMONGI_ID,
    SourceAdapterError,
    adapt_candidate,
)
from fp1_limongi_phase_history import build_phase_histories


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_limongi_phase_mass_history_contract_v1.json"
PHASE_HISTORY_TOOL_PATH = TOOL_PATH.with_name("fp1_limongi_phase_history.py")


class PhaseHistoryAuditError(ValueError):
    """Limongi phase history violates its fail-closed review contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise PhaseHistoryAuditError(f"cannot read {path}: {exc}") from exc
    return digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PhaseHistoryAuditError(f"cannot read phase-history contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-limongi-phase-mass-history-contract"
        or contract.get("schema_version") != 1
    ):
        raise PhaseHistoryAuditError("unsupported Limongi phase-history contract")
    limitations = contract.get("limitations", {})
    if limitations.get("canonical_conversion_allowed") is not False:
        raise PhaseHistoryAuditError("phase-history review unexpectedly permits conversion")
    if limitations.get("runtime_deposition_allowed") is not False:
        raise PhaseHistoryAuditError("phase-history review unexpectedly permits deposition")
    return contract


def _yield_coordinate(record: dict[str, Any]) -> tuple[int, int, float]:
    source = record["source_model_coordinate"]
    return (
        int(source["rotation_velocity_km_s"]),
        int(source["metallicity_feh"]),
        float(source["initial_mass_msun"]),
    )


def _extreme(residuals: list[dict[str, Any]], fn: Any) -> dict[str, Any]:
    return fn(residuals, key=lambda value: value["residual_msun"])


def audit_limongi_phase_mass_history(
    *, root: Path = DEFAULT_ROOT, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    source = adapt_candidate(LIMONGI_ID, root=Path(root), include_records=True)
    properties = source["source_components"]["evolutionary_properties"]
    if properties["all_duplicate_rows_physically_identical"] is not True:
        raise PhaseHistoryAuditError("non-identical Limongi phase duplicates cannot collapse")
    try:
        canonical_histories, phase_diagnostics = build_phase_histories(
            properties["records"],
            contract["phase_order"],
        )
    except ValueError as exc:
        raise PhaseHistoryAuditError(str(exc)) from exc
    histories: dict[tuple[int, int, float], dict[str, Any]] = {}
    for coordinate, canonical in canonical_histories.items():
        histories[coordinate] = {
            "terminal_age_yr": canonical["terminal_age_yr"],
            "phase_node_count_excluding_age_zero": canonical["unique_phase_count"],
            "terminal_cumulative_wind_mass_msun": canonical[
                "terminal_cumulative_wind_mass_msun"
            ],
            "nodes": [
                {
                    "phase": "age_zero",
                    "age_yr": 0.0,
                    "cumulative_wind_mass_msun": 0.0,
                }
            ]
            + [
                {
                    "phase": node["phase"],
                    "age_yr": node["cumulative_age_yr"],
                    "cumulative_wind_mass_msun": node[
                        "cumulative_wind_mass_msun"
                    ],
                }
                for node in canonical["nodes"]
            ],
        }
    collapsed_extra_rows = phase_diagnostics["collapsed_duplicate_row_count"]
    monotonic_mass_violation_count = phase_diagnostics[
        "nonincreasing_total_mass_violation_count"
    ]
    negative_cumulative_mass_count = phase_diagnostics[
        "negative_cumulative_mass_count"
    ]
    maximum_age_yr = phase_diagnostics["maximum_terminal_age_yr"]
    minimum_age_yr = phase_diagnostics["minimum_terminal_age_yr"]
    phase_node_counts = [
        history["phase_node_count_excluding_age_zero"]
        for history in histories.values()
    ]

    components = source["source_components"]
    wind_sums = {
        _yield_coordinate(record): math.fsum(
            record["source_reported_isotopic_yields"].values()
        )
        for record in components["wind_yields"]["records"]
    }
    recommended_sums = {
        _yield_coordinate(record): math.fsum(
            record["source_reported_isotopic_yields"].values()
        )
        for record in components["recommended_yields"]["records"]
    }
    residuals: list[dict[str, Any]] = []
    for coordinate, history in sorted(histories.items()):
        source_component = "table9_wind" if coordinate[2] <= 25.0 else "table8_setR_wind_only"
        integrated = wind_sums.get(coordinate) if coordinate[2] <= 25.0 else recommended_sums.get(coordinate)
        if integrated is None:
            raise PhaseHistoryAuditError(f"missing integrated wind yield at {coordinate}")
        phase_mass = history["terminal_cumulative_wind_mass_msun"]
        residuals.append(
            {
                "coordinate": {
                    "rotation_velocity_km_s": coordinate[0],
                    "metallicity_feh": coordinate[1],
                    "initial_mass_msun": coordinate[2],
                },
                "integrated_yield_component": source_component,
                "phase_endpoint_wind_mass_msun": phase_mass,
                "integrated_isotopic_wind_mass_msun": integrated,
                "residual_msun": phase_mass - integrated,
                "absolute_residual_msun": abs(phase_mass - integrated),
            }
        )
    maximum_absolute = max(residuals, key=lambda value: value["absolute_residual_msun"])
    quantization_half_width = 0.5 * float(
        contract["limitations"]["phase_endpoint_total_mass_precision_msun"]
    )
    over_quantization = [
        value for value in residuals if value["absolute_residual_msun"] > quantization_half_width
    ]
    return {
        "schema": "snrt-g2-limongi-phase-mass-history-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "mass_history_recoverable_composition_and_closure_blocked",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "runtime_deposition_allowed": False,
        "contract_path": str(contract_path),
        "contract_sha256": _sha256(contract_path),
        "audit_code_sha256": _sha256(TOOL_PATH),
        "phase_history_shared_code_sha256": _sha256(PHASE_HISTORY_TOOL_PATH),
        "source_adapter_code_sha256": source["adapter_code_sha256"],
        "duplicate_resolution": {
            "duplicate_coordinate_count": properties[
                "duplicate_model_phase_coordinate_count"
            ],
            "collapsed_extra_row_count": collapsed_extra_rows,
            "all_collapsed_rows_physically_identical": True,
        },
        "mass_history": {
            "model_count": len(histories),
            "phase_row_count_after_exact_collapse": sum(phase_node_counts),
            "minimum_phase_node_count_per_model": min(phase_node_counts),
            "maximum_phase_node_count_per_model": max(phase_node_counts),
            "age_zero_anchor_count": len(histories),
            "minimum_terminal_age_yr": minimum_age_yr,
            "maximum_terminal_age_yr": maximum_age_yr,
            "monotonic_mass_violation_count": monotonic_mass_violation_count,
            "negative_cumulative_mass_count": negative_cumulative_mass_count,
            "observed_source_order_matches_contract_rank": phase_diagnostics[
                "observed_source_order_matches_contract_rank"
            ],
            "phase_order_violation_count": phase_diagnostics[
                "phase_order_violation_count"
            ],
            "time_resolved_mass_available": True,
            "time_resolved_isotopic_composition_available": False,
        },
        "phase_order_provenance": {
            "source": "g2_limongi_phase_mass_history_contract_v1",
            "source_attested_for_intermediate_burning_order": False,
            "interpretation": (
                "MS/H/He/PSN labels are source-attested; the C/Ne/O/Si "
                "ordering is a project contract assumption."
            ),
        },
        "terminal_integrated_wind_closure": {
            "model_count": len(residuals),
            "residual_definition": "initial_minus_table5_terminal_mass minus integrated isotopic wind mass",
            "minimum_residual_msun": _extreme(residuals, min),
            "maximum_residual_msun": _extreme(residuals, max),
            "maximum_absolute_residual": maximum_absolute,
            "table5_total_mass_quantization_half_width_msun": quantization_half_width,
            "model_count_exceeding_quantization_half_width": len(over_quantization),
            "all_models_close_within_printed_quantization": not over_quantization,
        },
        "interpretation": (
            "Table5 supplies a monotonic phase-endpoint wind-mass history after "
            "exact duplicate collapse, but not a time-resolved isotopic composition. "
            "The endpoint mass and integrated isotope sums also disagree beyond "
            "the nominal printed-mass half-bin for some models, so the history "
            "cannot normalize canonical ejecta without a documented reconciliation."
        ),
        "blockers": contract["approval"]["required_before_approval"],
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
        report = audit_limongi_phase_mass_history(
            root=args.root, contract_path=args.contract
        )
    except (PhaseHistoryAuditError, SourceAdapterError) as exc:
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
