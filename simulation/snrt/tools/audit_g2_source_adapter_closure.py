#!/usr/bin/env python3
"""Audit source-internal closure without claiming canonical G2 completeness."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

from adapt_g2_candidate_sources import (
    DEFAULT_CONTRACT,
    DEFAULT_ROOT,
    LIMONGI_ID,
    NUGRID_ID,
    SourceAdapterError,
    adapt_candidate,
)


TOOL_PATH = Path(__file__).resolve()
DIAGNOSTIC_ROUNDING_TOLERANCE_MSUN = 1.0e-3


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _limongi_coordinate(record: dict[str, Any], field: str) -> tuple[int, int, float]:
    coordinate = record[field]
    return (
        coordinate["rotation_velocity_km_s"],
        coordinate["metallicity_feh"],
        coordinate["initial_mass_msun"],
    )


def _audit_limongi(report: dict[str, Any]) -> dict[str, Any]:
    components = report["source_components"]
    recommended = components["recommended_yields"]["records"]
    wind = components["wind_yields"]["records"]
    evolutionary = components["evolutionary_properties"]["records"]
    presupernova = components["presupernova_properties"]["records"]

    recommended_coordinates = {
        _limongi_coordinate(record, "source_model_coordinate") for record in recommended
    }
    evolutionary_coordinates = {
        _limongi_coordinate(record, "source_coordinate") for record in evolutionary
    }
    presupernova_coordinates = {
        _limongi_coordinate(record, "source_coordinate") for record in presupernova
    }

    def budget(component_records: list[dict[str, Any]]) -> dict[str, Any]:
        ratios = [
            record["source_reported_yield_sum"]
            / record["source_model_coordinate"]["initial_mass_msun"]
            for record in component_records
        ]
        over_budget = [ratio for ratio in ratios if ratio > 1.0]
        return {
            "model_count": len(component_records),
            "source_yield_sum_exceeds_initial_mass_count": len(over_budget),
            "maximum_source_yield_sum_to_initial_mass_ratio": max(ratios),
            "canonical_mass_closure_claim": "not_allowed; source units and physical returned/remnant definitions are unresolved",
        }

    return {
        "status": "source_internal_diagnostics_pass_canonical_closure_blocked",
        "recommended_yield_budget": budget(recommended),
        "wind_yield_budget": budget(wind),
        "recommended_minus_wind_diagnostic": report[
            "uninterpreted_component_relation_diagnostics"
        ],
        "coordinate_alignment": {
            "recommended_yield_model_count": len(recommended_coordinates),
            "evolutionary_property_model_count": len(evolutionary_coordinates),
            "presupernova_property_model_count": len(presupernova_coordinates),
            "recommended_missing_evolutionary_property_count": len(
                recommended_coordinates - evolutionary_coordinates
            ),
            "recommended_missing_presupernova_property_count": len(
                recommended_coordinates - presupernova_coordinates
            ),
            "recommended_missing_presupernova_coordinates": [
                {
                    "rotation_velocity_km_s": coordinate[0],
                    "metallicity_feh": coordinate[1],
                    "initial_mass_msun": coordinate[2],
                }
                for coordinate in sorted(recommended_coordinates - presupernova_coordinates)
            ],
        },
        "duplicate_evolutionary_model_phase_coordinate_count": components[
            "evolutionary_properties"
        ]["duplicate_model_phase_coordinate_count"],
        "canonical_closure_available": False,
        "canonical_closure_blockers": [
            "source_mass_unit_is_literature_supported_but_project_approval_is_missing",
            "isotope_decay_and_element_mapping_unresolved",
            "returned_mass_and_remnant_mass_semantics_unresolved",
            "age_resolved_cumulative_history_absent",
            "energy_and_momentum_absent",
        ],
    }


def _nugrid_coordinate(record: dict[str, Any]) -> tuple[float, float]:
    coordinate = record["source_coordinate"]
    return coordinate["initial_mass_msun"], coordinate["metallicity_mass_fraction"]


def _duplicate_diagnostics(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[float, float], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[_nugrid_coordinate(record)].append(record)
    diagnostics: list[dict[str, Any]] = []
    for coordinate, members in sorted(grouped.items()):
        if len(members) < 2:
            continue
        elements = members[0]["source_reported_element_yields_msun"].keys()
        maximum_difference = max(
            abs(
                members[first]["source_reported_element_yields_msun"][element]
                - members[second]["source_reported_element_yields_msun"][element]
            )
            for first in range(len(members))
            for second in range(first + 1, len(members))
            for element in elements
        )
        diagnostics.append(
            {
                "initial_mass_msun": coordinate[0],
                "metallicity_mass_fraction": coordinate[1],
                "multiplicity": len(members),
                "lifetimes_identical": len({member["lifetime_yr"] for member in members}) == 1,
                "final_masses_identical": len({member["final_mass_msun"] for member in members}) == 1,
                "maximum_element_yield_difference_msun": maximum_difference,
                "numerically_identical_source_records": maximum_difference == 0.0,
            }
        )
    return diagnostics


def _audit_nugrid(report: dict[str, Any]) -> dict[str, Any]:
    components: dict[str, Any] = report["source_components"]
    component_closure: dict[str, Any] = {}
    for name, component in components.items():
        records = component["records"]
        residuals = [
            record["diagnostic_initial_minus_final_minus_yields_msun"] for record in records
        ]
        abundance_residuals = [
            1.0 - math.fsum(record["source_initial_mass_fractions"].values())
            for record in records
        ]
        component_closure[name] = {
            "block_count": len(records),
            "nonnegative_source_yield_value_count": sum(
                value >= 0.0
                for record in records
                for value in record["source_reported_element_yields_msun"].values()
            ),
            "negative_source_yield_value_count": sum(
                value < 0.0
                for record in records
                for value in record["source_reported_element_yields_msun"].values()
            ),
            "maximum_absolute_initial_minus_final_minus_yields_msun": max(
                abs(value) for value in residuals
            ),
            "blocks_within_diagnostic_rounding_tolerance": sum(
                abs(value) <= DIAGNOSTIC_ROUNDING_TOLERANCE_MSUN for value in residuals
            ),
            "diagnostic_rounding_tolerance_msun": DIAGNOSTIC_ROUNDING_TOLERANCE_MSUN,
            "rounding_tolerance_is_not_a_canonical_approval_threshold": True,
            "maximum_absolute_initial_abundance_sum_residual": max(
                abs(value) for value in abundance_residuals
            ),
            "duplicate_coordinate_diagnostics": _duplicate_diagnostics(records),
        }
    return {
        "status": "source_internal_diagnostics_pass_canonical_closure_blocked",
        "component_closure": component_closure,
        "component_relation_diagnostics": report[
            "uninterpreted_component_relation_diagnostics"
        ],
        "canonical_closure_available": False,
        "canonical_closure_blockers": [
            "duplicate_mass_metallicity_coordinate_requires_source_selection_policy",
            "integrated_snapshots_are_not_cumulative_age_histories",
            "source_component_semantics_resolved_but_project_channel_ownership_unapproved",
            "runtime_mass_grid_incomplete",
            "energy_and_momentum_absent",
        ],
    }


def audit_source_adapter_closure(
    *,
    root: Path = DEFAULT_ROOT,
    manifest_path: Path | None = None,
    contract_path: Path = DEFAULT_CONTRACT,
) -> dict[str, Any]:
    limongi = adapt_candidate(
        LIMONGI_ID,
        root=root,
        manifest_path=manifest_path,
        contract_path=contract_path,
        include_records=True,
    )
    nugrid = adapt_candidate(
        NUGRID_ID,
        root=root,
        manifest_path=manifest_path,
        contract_path=contract_path,
        include_records=True,
    )
    return {
        "schema": "snrt-g2-source-adapter-closure-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "review_only_blocked",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "candidates": {
            LIMONGI_ID: _audit_limongi(limongi),
            NUGRID_ID: _audit_nugrid(nugrid),
        },
        "gate_blockers": [
            "source_internal_consistency_does_not_supply_missing_age_energy_momentum_or_channel_semantics",
            "candidate_sources_are_not_approved",
            "no_complete_runtime_grid_for_channels_1_to_3",
        ],
        "adapter_code_sha256": limongi["adapter_code_sha256"],
        "adapter_contract_sha256": limongi["adapter_contract_sha256"],
        "source_semantics_evidence_sha256": limongi["source_semantics_evidence_sha256"],
        "closure_audit_code_sha256": _sha256(TOOL_PATH),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, help="default: ROOT/acquisition_manifest_v1.json")
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_source_adapter_closure(
            root=args.root,
            manifest_path=args.manifest,
            contract_path=args.contract,
        )
    except SourceAdapterError as exc:
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
