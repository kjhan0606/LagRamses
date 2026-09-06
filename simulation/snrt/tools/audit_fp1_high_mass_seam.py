#!/usr/bin/env python3
"""Audit the F-P1 40--120 Msun seam as source-node review evidence.

The staged Sukhbold W18/N20 records demonstrate engine-dependent,
non-monotonic outcomes. They are useful for choosing a future fate resolver,
but they do not by themselves supply the project's complete age, decay,
wind-ownership, energy, momentum, or redistribution contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_FATE_MAP = SNRT_ROOT / "config" / "fp1_population_fate_map_v1.json"
DEFAULT_CANDIDATE_AUDIT = SNRT_ROOT / "data" / "g2_sukhbold2016_candidate_audit.json"
DEFAULT_PROJECTION = SNRT_ROOT / "data" / "g2_sukhbold_channel_projection_review.json"
DEFAULT_JSON_OUT = SNRT_ROOT / "data" / "fp1_high_mass_seam_review.json"
SEAM_ID = "massive_terminal_fate_seam"
SEAM = [40.0, 120.0]
EXPECTED_MASSES = [40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 100.0, 120.0]
EXPECTED_ENGINES = ["N20", "W18"]
ROUNDED_SOURCE_REVIEW_RELATIVE_TOLERANCE = 7.0e-3


class HighMassSeamAuditError(ValueError):
    """Review inputs are inconsistent or overclaim the high-mass seam."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HighMassSeamAuditError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise HighMassSeamAuditError(f"{label} must be a JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise HighMassSeamAuditError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def audit_high_mass_seam(
    *,
    fate_map_path: Path = DEFAULT_FATE_MAP,
    candidate_audit_path: Path = DEFAULT_CANDIDATE_AUDIT,
    projection_path: Path = DEFAULT_PROJECTION,
) -> dict[str, Any]:
    fate_map_path = Path(fate_map_path).resolve()
    candidate_audit_path = Path(candidate_audit_path).resolve()
    projection_path = Path(projection_path).resolve()
    fate_map = _read_json(fate_map_path, "F-P1 fate map")
    candidate_audit = _read_json(candidate_audit_path, "Sukhbold candidate audit")
    projection = _read_json(projection_path, "Sukhbold projection review")

    intervals = fate_map.get("intervals")
    if not isinstance(intervals, list):
        raise HighMassSeamAuditError("F-P1 fate map intervals are missing")
    seam = next((item for item in intervals if isinstance(item, dict) and item.get("id") == SEAM_ID), None)
    if seam is None or seam.get("mass_msun") != SEAM:
        raise HighMassSeamAuditError("F-P1 high-mass seam is missing or changed")
    if seam.get("fate_class") != "unresolved":
        raise HighMassSeamAuditError("high-mass seam must remain unresolved until source approval")

    evidence = candidate_audit.get("high_mass_engine_evidence")
    engines = evidence.get("engines") if isinstance(evidence, dict) else None
    if not isinstance(engines, dict) or sorted(engines) != EXPECTED_ENGINES:
        raise HighMassSeamAuditError("Sukhbold W18/N20 engine coverage changed")

    engine_summary: dict[str, Any] = {}
    for engine_name in EXPECTED_ENGINES:
        engine = engines[engine_name]
        if not isinstance(engine, dict):
            raise HighMassSeamAuditError(f"Sukhbold engine record is malformed: {engine_name}")
        results = engine.get("high_mass_results")
        if not isinstance(results, dict):
            raise HighMassSeamAuditError(f"high-mass results are missing: {engine_name}")
        masses = sorted(float(key) for key in results)
        if masses != EXPECTED_MASSES:
            raise HighMassSeamAuditError(f"high-mass mass nodes changed for {engine_name}")
        outcomes = [results[str(mass)]["outcome"] for mass in EXPECTED_MASSES]
        if "no_positive_explosion_energy" not in outcomes or "explosion_energy_positive" not in outcomes:
            raise HighMassSeamAuditError(f"engine evidence is not outcome-diverse: {engine_name}")
        engine_summary[engine_name] = {
            "mass_nodes_msun": EXPECTED_MASSES,
            "outcomes": dict(zip((str(mass) for mass in EXPECTED_MASSES), outcomes)),
            "positive_energy_mass_nodes_msun": [
                mass for mass in EXPECTED_MASSES if results[str(mass)]["outcome"] == "explosion_energy_positive"
            ],
            "nonpositive_energy_mass_nodes_msun": [
                mass for mass in EXPECTED_MASSES if results[str(mass)]["outcome"] == "no_positive_explosion_energy"
            ],
        }

    mass_budget = evidence.get("mass_budget_review")
    if not isinstance(mass_budget, dict):
        raise HighMassSeamAuditError("high-mass mass-budget review is missing")
    if mass_budget.get("records_evaluated") != 6:
        raise HighMassSeamAuditError("high-mass mass-budget record count changed")
    maximum_relative_residual = mass_budget.get("maximum_absolute_relative_residual")
    if not isinstance(maximum_relative_residual, (int, float)) or isinstance(
        maximum_relative_residual, bool
    ) or not math.isfinite(float(maximum_relative_residual)):
        raise HighMassSeamAuditError("high-mass mass-budget residual is malformed")
    within_review_tolerance = (
        float(maximum_relative_residual) <= ROUNDED_SOURCE_REVIEW_RELATIVE_TOLERANCE
    )
    if not within_review_tolerance:
        raise HighMassSeamAuditError("high-mass rounded-source review tolerance exceeded")
    if mass_budget.get("review_bound_applied") is not False or mass_budget.get(
        "exact_mass_closure_claimed"
    ) is not False:
        raise HighMassSeamAuditError("candidate audit overclaims high-mass closure")

    def wind_record(engine_name: str, mass: float) -> dict[str, Any] | None:
        engine = engines[engine_name]
        key = str(mass)
        yields = engine.get("high_mass_yields", {})
        implosions = engine.get("high_mass_implosion_winds", {})
        record = yields.get(key) if isinstance(yields, dict) else None
        if record is None and isinstance(implosions, dict):
            record = implosions.get(key)
        return record if isinstance(record, dict) else None

    common_wind_masses = [
        mass
        for mass in EXPECTED_MASSES
        if wind_record("N20", mass) is not None and wind_record("W18", mass) is not None
    ]
    if common_wind_masses != [60.0, 80.0, 100.0, 120.0]:
        raise HighMassSeamAuditError("cross-engine common wind coverage changed")
    wind_comparisons: list[dict[str, Any]] = []
    radioactive_epoch_warnings: list[dict[str, Any]] = []
    k40_duplicate_record_count = 0
    for mass in common_wind_masses:
        n20 = wind_record("N20", mass)
        w18 = wind_record("W18", mass)
        assert n20 is not None and w18 is not None
        n20_wind = float(n20["stable_wind_sum_msun"])
        w18_wind = float(w18["stable_wind_sum_msun"])
        n20_elements = n20.get("stable_wind_by_tracked_element_msun")
        w18_elements = w18.get("stable_wind_by_tracked_element_msun")
        if not isinstance(n20_elements, dict) or not isinstance(w18_elements, dict):
            raise HighMassSeamAuditError("cross-engine stable-wind element data are missing")
        if set(n20_elements) != set(w18_elements):
            raise HighMassSeamAuditError("cross-engine stable-wind element sets disagree")
        element_differences = {
            element: float(n20_elements[element]) - float(w18_elements[element])
            for element in sorted(n20_elements)
            if float(n20_elements[element]) != float(w18_elements[element])
        }
        wind_comparisons.append(
            {
                "zams_mass_msun": mass,
                "n20_stable_wind_msun": n20_wind,
                "w18_stable_wind_msun": w18_wind,
                "signed_n20_minus_w18_msun": n20_wind - w18_wind,
                "bit_identical": n20_wind == w18_wind and not element_differences,
                "element_differences_msun": element_differences,
            }
        )
        for record in (n20, w18):
            duplicates = record.get("cross_segment_duplicate_isotopes") or record.get(
                "duplicate_isotopes"
            )
            if isinstance(duplicates, list) and "k40" in duplicates:
                k40_duplicate_record_count += 1
        n20_radio = n20.get("selected_radioactive_inventory")
        w18_radio = w18.get("selected_radioactive_inventory")
        if not isinstance(n20_radio, dict) or not isinstance(w18_radio, dict):
            raise HighMassSeamAuditError("cross-engine radioactive wind data are missing")
        for isotope in sorted(set(n20_radio) & set(w18_radio)):
            n20_value = float(n20_radio[isotope]["wind_msun"])
            w18_value = float(w18_radio[isotope]["wind_msun"])
            if n20_value == w18_value:
                continue
            if n20_value > 0.0 and w18_value > 0.0:
                log10_spread = abs(math.log10(n20_value / w18_value))
                zero_vs_positive = False
            else:
                log10_spread = None
                zero_vs_positive = (n20_value == 0.0) != (w18_value == 0.0)
            if zero_vs_positive or (log10_spread is not None and log10_spread >= 6.0):
                radioactive_epoch_warnings.append(
                    {
                        "zams_mass_msun": mass,
                        "isotope": isotope,
                        "n20_wind_msun": n20_value,
                        "w18_wind_msun": w18_value,
                        "absolute_log10_spread_dex": log10_spread,
                        "zero_vs_positive": zero_vs_positive,
                    }
                )
    failed_records = [
        (engine_name, float(mass), result)
        for engine_name in EXPECTED_ENGINES
        for mass, result in engines[engine_name]["high_mass_results"].items()
        if result["outcome"] == "no_positive_explosion_energy"
    ]
    failed_outcome_count = len(failed_records)
    remnant_keys = {
        "baryonic_remnant_mass_msun",
        "gravitational_remnant_mass_msun",
        "remnant_mass_msun",
        "mass_cut_msun",
    }

    def has_source_remnant(record: dict[str, Any]) -> bool:
        return any(
            key in record
            and isinstance(record[key], (int, float))
            and not isinstance(record[key], bool)
            and math.isfinite(float(record[key]))
            and float(record[key]) >= 0.0
            for key in remnant_keys
        )

    failed_nodes_with_source_remnant_count = sum(
        has_source_remnant(result)
        or (
            (wind := wind_record(engine_name, mass)) is not None
            and has_source_remnant(wind)
        )
        for engine_name, mass, result in failed_records
    )
    terminal_yield_record_count = sum(
        len(engines[engine_name].get("high_mass_yields", {}))
        for engine_name in EXPECTED_ENGINES
    )
    failed_wind_record_count = sum(
        wind_record(engine_name, mass) is not None
        for engine_name, mass, _ in failed_records
    )
    common_wind_record_count = 2 * len(common_wind_masses)
    all_k40_duplicates = k40_duplicate_record_count == common_wind_record_count
    source_erratum_required = any(
        not record["bit_identical"] for record in wind_comparisons
    )
    radioactive_wind_values_admissible = (
        not radioactive_epoch_warnings and k40_duplicate_record_count == 0
    )

    projection_records = projection.get("records")
    if not isinstance(projection_records, list):
        raise HighMassSeamAuditError("Sukhbold projection records are missing")
    if projection.get("canonical_rows_emitted") != 0 or projection.get("production_ready") is not False:
        raise HighMassSeamAuditError("Sukhbold projection unexpectedly permits production")
    if any(record.get("canonical_row_emitted") for record in projection_records if isinstance(record, dict)):
        raise HighMassSeamAuditError("Sukhbold review projection emitted a canonical row")

    return {
        "schema": "snrt-fp1-high-mass-seam-review",
        "schema_version": 1,
        "gate": "F-P1",
        "status": "review_only_engine_comparison_fate_unresolved",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "runtime_activation_allowed": False,
        "seam": {
            "id": SEAM_ID,
            "mass_msun": SEAM,
            "fate_map_status": seam.get("evidence_status"),
            "mass_nodes_msun": EXPECTED_MASSES,
        },
        "candidate": {
            "candidate_id": candidate_audit.get("source_identity", {}).get("candidate_id"),
            "article_doi": candidate_audit.get("source_identity", {}).get("article_doi"),
            "data_doi": candidate_audit.get("source_identity", {}).get("data_doi"),
            "source_identity": candidate_audit.get("source_identity"),
            "candidate_audit": {"path": str(candidate_audit_path), "sha256": _sha256(candidate_audit_path)},
            "projection_review": {"path": str(projection_path), "sha256": _sha256(projection_path)},
            "engines": engine_summary,
        },
        "resolved": False,
        "review_closure": {
            "scope": "available rounded Sukhbold high-mass terminal-yield records only",
            "record_count": mass_budget["records_evaluated"],
            "maximum_absolute_residual_msun": mass_budget.get(
                "maximum_absolute_residual_msun"
            ),
            "maximum_absolute_relative_residual": float(maximum_relative_residual),
            "review_relative_tolerance": ROUNDED_SOURCE_REVIEW_RELATIVE_TOLERANCE,
            "within_review_tolerance": within_review_tolerance,
            "exact_mass_closure_claimed": False,
            "production_acceptance_claimed": False,
        },
        "source_node_completeness": {
            "outcome_record_count": len(EXPECTED_ENGINES) * len(EXPECTED_MASSES),
            "terminal_yield_record_count": terminal_yield_record_count,
            "failed_outcome_count": failed_outcome_count,
            "failed_wind_record_count": failed_wind_record_count,
            "failed_nodes_with_source_remnant_count": failed_nodes_with_source_remnant_count,
            "complete": (
                failed_nodes_with_source_remnant_count == failed_outcome_count
                and failed_wind_record_count == failed_outcome_count
                and terminal_yield_record_count + failed_outcome_count
                == len(EXPECTED_ENGINES) * len(EXPECTED_MASSES)
            ),
        },
        "cross_engine_wind_review": {
            "common_mass_nodes_msun": common_wind_masses,
            "stable_wind_comparisons": wind_comparisons,
            "all_common_stable_winds_bit_identical": all(
                record["bit_identical"] for record in wind_comparisons
            ),
            "source_erratum_or_explanation_required": source_erratum_required,
            "k40_cross_segment_duplicate_present_in_all_common_records": all_k40_duplicates,
            "k40_cross_segment_duplicate_record_count": k40_duplicate_record_count,
            "common_wind_record_count": common_wind_record_count,
            "radioactive_reference_epoch_warning_count": len(radioactive_epoch_warnings),
            "radioactive_reference_epoch_warnings": radioactive_epoch_warnings,
            "radioactive_wind_values_admissible": radioactive_wind_values_admissible,
        },
        "blockers": [
            "no_approved_source_node_fate_map_for_40_to_120_msun",
            "engine_branch_outcomes_are_comparison_evidence_not_a_project_fate_law",
            "age_resolved_wind_history_and_terminal_lumping_policy_missing",
            "complete_decay_projection_and_element_closure_missing",
            "canonical_event_energy_and_momentum_deposition_contract_missing",
            "wind_terminal_remnant_ownership_and_cross_source_seams_unapproved",
            "third_party_redistribution_permission_not_verified",
            *(
                ["failed_outcome_remnant_and_n20_wind_coverage_missing"]
                if failed_nodes_with_source_remnant_count < failed_outcome_count
                or failed_wind_record_count < failed_outcome_count
                else []
            ),
            *(
                ["cross_engine_stable_wind_discrepancy_requires_source_erratum_or_explanation"]
                if source_erratum_required else []
            ),
            *(
                ["radioactive_wind_reference_epoch_inconsistent_between_engine_branches"]
                if radioactive_epoch_warnings else []
            ),
            *(
                ["cross_segment_k40_requires_explicit_decay_projection_and_duplicate_resolution"]
                if k40_duplicate_record_count else []
            ),
        ],
        "required_next_inputs": [
            "licensed_source_node_fate_map_with_source_hull_and_approval",
            "age_resolved_wind_or_approved_terminal_lumping_contract",
            "decay_projection_and_reduced_chemistry_closure",
            "energy_momentum_and_deposition_semantics",
            "non_overlapping_wind_terminal_remnant_ownership",
        ],
        "interpretation": (
            "W18/N20 contain both positive- and non-positive-energy outcomes across the same "
            "40--120 Msun mass range. This proves mass-only direct-collapse fallback is unsafe; "
            "it does not select a production fate map or authorize any canonical feedback row."
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fate-map", type=Path, default=DEFAULT_FATE_MAP)
    parser.add_argument("--candidate-audit", type=Path, default=DEFAULT_CANDIDATE_AUDIT)
    parser.add_argument("--projection", type=Path, default=DEFAULT_PROJECTION)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    args = parser.parse_args(argv)
    try:
        report = audit_high_mass_seam(
            fate_map_path=args.fate_map,
            candidate_audit_path=args.candidate_audit,
            projection_path=args.projection,
        )
    except HighMassSeamAuditError as exc:
        print(f"F-P1 high-mass seam audit ERROR: {exc}", file=sys.stderr)
        return 2
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("F-P1 high-mass seam: review_only_engine_comparison_fate_unresolved")
    print("engines=N20,W18 mass_nodes=9")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
