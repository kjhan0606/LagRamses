#!/usr/bin/env python3
"""Prepare an explicit, non-ranking HESMA SNIa source-selection packet."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

from fp2_provenance import project_relative


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_MODEL_COMPARISON = SNRT_ROOT / "data" / "fp2_snia_hesma_model_comparison.json"
DEFAULT_PROFILE_COMPARISON = SNRT_ROOT / "data" / "fp2_snia_hesma_profile_estimator_comparison.json"


class SelectionPacketError(ValueError):
    """The source-selection packet inputs are not review-complete."""


def _sha256(path: Path) -> str:
    if not path.is_file():
        raise SelectionPacketError(f"file is missing: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SelectionPacketError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SelectionPacketError(f"JSON object required: {path}")
    return value


def build_packet(
    *,
    model_comparison_path: Path = DEFAULT_MODEL_COMPARISON,
    profile_comparison_path: Path = DEFAULT_PROFILE_COMPARISON,
) -> dict[str, Any]:
    model_comparison_path = Path(model_comparison_path).resolve()
    profile_comparison_path = Path(profile_comparison_path).resolve()
    model_comparison = _read_json(model_comparison_path)
    profile_comparison = _read_json(profile_comparison_path)
    if model_comparison.get("status") != "review_only_no_model_selected":
        raise SelectionPacketError("HESMA model comparison is not unselected review data")
    if profile_comparison.get("status") != "review_only_diagnostic":
        raise SelectionPacketError("HESMA profile comparison is not review-only diagnostic data")
    if model_comparison.get("model_count") != 15 or profile_comparison.get("model_count") != 15:
        raise SelectionPacketError("HESMA comparison model counts are not 15")
    if model_comparison.get("selection", {}).get("selected_model_id") is not None:
        raise SelectionPacketError("model comparison already contains a selected model")
    if profile_comparison.get("admission", {}).get("selected_estimator") is not None:
        raise SelectionPacketError("profile comparison already contains a selected estimator")

    profile_by_model = {row["model_id"]: row for row in profile_comparison["models"]}
    screening_rows: list[dict[str, Any]] = []
    for row in model_comparison["models"]:
        model_id = row["model_id"]
        profile_row = profile_by_model[model_id]
        estimator_differences = {
            name: value["relative_difference_from_integrated_stable_mass"]
            for name, value in profile_row["estimators"].items()
        }
        within_review_threshold = all(value <= 0.05 for value in estimator_differences.values())
        has_source_warning = row["physical_warning_count"] > 0
        screening_rows.append(
            {
                "model_id": model_id,
                "review_screen": (
                    "profile_consistent_review_candidate"
                    if within_review_threshold and not has_source_warning
                    else "physical_warning_requires_resolution"
                ),
                "screening_threshold_relative_mass_difference": 0.05,
                "screening_is_not_production_approval": True,
                "profile_mass_relative_difference_by_estimator": estimator_differences,
                "profile_kinetic_energy_by_estimator_erg": {
                    name: value["kinetic_energy_estimate_erg"]
                    for name, value in profile_row["estimators"].items()
                },
                "integrated_stable_element_mass_msun": row["all_source_stable_element_mass_msun"],
                "physical_warnings": row["physical_warnings"],
                "quarantine_required": any(
                    warning.get("requires_quarantine") is True
                    for warning in row["physical_warnings"]
                ),
                "source_review_classification": row.get("profile_review_classification"),
            }
        )

    return {
        "schema": "snrt-fp2-snia-hesma-source-selection-packet",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "review_only_selection_pending",
        "source": {
            "candidate_id": "hesma_model_archive_snia_profiles",
            "record_id": model_comparison["source"]["record_id"],
            "package_sha256": model_comparison["source"]["package_sha256"],
            "model_comparison_path": project_relative(model_comparison_path),
            "model_comparison_sha256": _sha256(model_comparison_path),
            "profile_comparison_path": project_relative(profile_comparison_path),
            "profile_comparison_sha256": _sha256(profile_comparison_path),
        },
        "screening_policy": {
            "purpose": "separate source-format/profile consistency review from physical source approval",
            "threshold": "existing 5% profile-mass warning threshold",
            "automatic_selection": False,
            "automatic_population_weighting": False,
            "automatic_event_energy": False,
        },
        "models": screening_rows,
        "selection": {
            "selected_model_id": None,
            "selected_population_mixture": None,
            "selected_profile_estimator": None,
            "approval_id": None,
            "status": "requires_explicit_physics_decision",
        },
        "physical_event_contract": {
            "decay_convention": None,
            "decay_horizon_yr": None,
            "isotope_to_project_element_policy": None,
            "returned_mass_msun_per_event": None,
            "terminal_remnant_msun_per_event": None,
            "energy_erg_per_event": None,
            "momentum_g_cm_s_per_event": None,
            "population_weight": None,
        },
        "admission": {
            "canonical_conversion_allowed": False,
            "runtime_activation_allowed": False,
            "canonical_rows_emitted": 0,
        },
        "blockers": [
            "explicitly select a model or population mixture",
            "resolve any selected-model profile warning",
            "approve the complete physical event contract before conversion",
        ],
        "tool_sha256": _sha256(TOOL_PATH),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-comparison", type=Path, default=DEFAULT_MODEL_COMPARISON)
    parser.add_argument("--profile-comparison", type=Path, default=DEFAULT_PROFILE_COMPARISON)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    try:
        report = build_packet(
            model_comparison_path=args.model_comparison,
            profile_comparison_path=args.profile_comparison,
        )
        payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(payload, encoding="utf-8")
    except (SelectionPacketError, OSError) as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, indent=2), file=sys.stderr)
        return 2
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
