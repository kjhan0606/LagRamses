#!/usr/bin/env python3
"""Checks for the review-only LC18 failed-wind cross-source audit."""

from __future__ import annotations

import argparse
import contextlib
import copy
import hashlib
import io
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from adapt_g2_candidate_sources import LIMONGI_ID, adapt_candidate  # noqa: E402
import audit_fp1_lc18_failed_wind_crosscheck as crosscheck_module  # noqa: E402
from audit_fp1_lc18_failed_wind_crosscheck import (  # noqa: E402
    Lc18FailedWindCrosscheckError,
    _build_phase_histories,
    audit_lc18_failed_wind_crosscheck,
)
from audit_g2_limongi_phase_mass_history import (  # noqa: E402
    audit_limongi_phase_mass_history,
)


SOURCE_FILES = (
    ROOT.parents[1]
    / "external"
    / "g2_candidates"
    / "boccioli_roberti2026_ccsn"
    / "LC18.zip",
    ROOT.parents[1]
    / "external"
    / "g2_candidates"
    / "limongi_chieffi_2018_cds"
    / "table5.dat",
    ROOT.parents[1]
    / "external"
    / "g2_candidates"
    / "limongi_chieffi_2018_cds"
    / "table7.dat",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _config_data_files() -> tuple[Path, ...]:
    paths = []
    for directory in (ROOT / "config", ROOT / "data"):
        paths.extend(path for path in directory.rglob("*") if path.is_file())
    return tuple(sorted(paths))


def _phase_fixture() -> tuple[list[dict], list[str]]:
    source = adapt_candidate(
        LIMONGI_ID,
        root=ROOT.parents[1] / "external" / "g2_candidates",
        include_records=True,
    )
    phase_contract = json.loads(
        (ROOT / "config" / "g2_limongi_phase_mass_history_contract_v1.json").read_text(
            encoding="utf-8"
        )
    )
    return (
        copy.deepcopy(source["source_components"]["evolutionary_properties"]["records"]),
        list(phase_contract["phase_order"]),
    )


def _first_model_records(records: list[dict]) -> list[dict]:
    coordinate = (
        records[0]["source_coordinate"]["rotation_velocity_km_s"],
        records[0]["source_coordinate"]["metallicity_feh"],
        records[0]["source_coordinate"]["initial_mass_msun"],
    )
    selected = [
        record
        for record in records
        if (
            record["source_coordinate"]["rotation_velocity_km_s"],
            record["source_coordinate"]["metallicity_feh"],
            record["source_coordinate"]["initial_mass_msun"],
        )
        == coordinate
        and int(record["source_coordinate"]["phase_occurrence"]) == 1
    ]
    assert selected
    return selected


def _expect_phase_error(
    records: list[dict], phase_order: list[str], diagnostic_key: str
) -> None:
    try:
        _build_phase_histories(records, phase_order)
    except Lc18FailedWindCrosscheckError as exc:
        assert "phase-history invariants violated" in str(exc), str(exc)
        assert exc.diagnostics[diagnostic_key] > 0, exc.diagnostics
    else:
        raise AssertionError(f"expected phase-history error for {diagnostic_key}")


def _phase_invariant_tests() -> None:
    records, phase_order = _phase_fixture()

    negative_duration = copy.deepcopy(records)
    _first_model_records(negative_duration)[0]["phase_duration_yr"] = -1.0
    _expect_phase_error(
        negative_duration,
        phase_order,
        "strictly_increasing_cumulative_age_violation_count",
    )

    increasing_mass = copy.deepcopy(records)
    model_records = _first_model_records(increasing_mass)
    ordered = sorted(
        model_records,
        key=lambda record: phase_order.index(record["source_coordinate"]["phase"]),
    )
    assert len(ordered) >= 2
    increasing_mass_value = float(ordered[0]["total_mass_msun"]) + 1.0
    ordered[1]["total_mass_msun"] = increasing_mass_value
    _expect_phase_error(
        increasing_mass,
        phase_order,
        "nonincreasing_total_mass_violation_count",
    )

    missing_terminal = copy.deepcopy(records)
    model_records = _first_model_records(missing_terminal)
    psn = next(
        record
        for record in model_records
        if record["source_coordinate"]["phase"] == "PSN"
    )
    missing_terminal.remove(psn)
    _expect_phase_error(
        missing_terminal,
        phase_order,
        "missing_psn_terminal_phase_count",
    )


def _differential_phase_audit_test() -> None:
    root = ROOT.parents[1] / "external" / "g2_candidates"
    legacy = audit_limongi_phase_mass_history(root=root)
    records, phase_order = _phase_fixture()
    _, current = _build_phase_histories(records, phase_order)
    assert current["model_count"] == legacy["mass_history"]["model_count"] == 108
    assert (
        current["unique_phase_row_count"]
        == legacy["mass_history"]["phase_row_count_after_exact_collapse"]
        == 845
    )
    assert (
        current["minimum_unique_phase_count_per_model"]
        == legacy["mass_history"]["minimum_phase_node_count_per_model"]
        == 3
    )
    assert (
        current["maximum_unique_phase_count_per_model"]
        == legacy["mass_history"]["maximum_phase_node_count_per_model"]
        == 8
    )
    assert (
        current["collapsed_duplicate_row_count"]
        == legacy["duplicate_resolution"]["collapsed_extra_row_count"]
        == 19
    )
    assert current["strictly_increasing_cumulative_age_violation_count"] == 0
    assert legacy["mass_history"]["monotonic_mass_violation_count"] == 0


def _cli_diagnostics_test() -> None:
    original = crosscheck_module.audit_lc18_failed_wind_crosscheck

    def fail_with_diagnostics(**_: object) -> dict:
        raise Lc18FailedWindCrosscheckError(
            "phase-history invariants violated",
            diagnostics={
                "strictly_increasing_cumulative_age_violation_count": 1,
                "nonincreasing_total_mass_violation_count": 0,
                "missing_psn_terminal_phase_count": 0,
            },
        )

    crosscheck_module.audit_lc18_failed_wind_crosscheck = fail_with_diagnostics
    stderr = io.StringIO()
    try:
        with contextlib.redirect_stderr(stderr):
            status = crosscheck_module.main([])
    finally:
        crosscheck_module.audit_lc18_failed_wind_crosscheck = original
    assert status == 2
    error_payload = json.loads(stderr.getvalue())
    assert error_payload["status"] == "error"
    assert error_payload["diagnostics"][
        "strictly_increasing_cumulative_age_violation_count"
    ] == 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    before = {str(path): _sha256(path) for path in SOURCE_FILES}
    config_data_before = {
        str(path): _sha256(path) for path in _config_data_files()
    }
    report = audit_lc18_failed_wind_crosscheck(
        root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    after = {str(path): _sha256(path) for path in SOURCE_FILES}
    config_data_after = {
        str(path): _sha256(path) for path in _config_data_files()
    }
    assert before == after
    assert config_data_before == config_data_after

    assert report["status"] == (
        "failed_wind_anomaly_independently_crosschecked_unresolved"
    )
    assert report["production_ready"] is False
    assert report["publication_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    assert report["runtime_deposition_allowed"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["physical_nodes_emitted"] == 0
    assert report["inquiry_sent"] is False

    joined = report["join"]
    assert joined["summary_row_count"] == 108
    assert joined["cds_model_count"] == 108
    assert joined["joined_model_count"] == 108
    assert joined["one_to_one"] is True
    assert joined["unmatched_release_coordinates"] == []
    assert joined["unmatched_cds_coordinates"] == []

    successful = report["successful_release_control"]
    assert successful["model_count"] == 52
    assert successful["summary_wind_positive_count"] == 52
    assert successful["release_wind_table_nonzero_count"] == 52
    assert successful["summary_minus_cds_above_nominal_half_bin_count"] == 52
    assert math.isclose(
        successful["maximum_absolute_summary_minus_release_wind_table_msun"],
        0.007183005956193367,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert math.isclose(
        successful["maximum_absolute_summary_minus_cds_terminal_wind_msun"],
        0.5842,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )

    failed = report["failed_wind_anomaly"]
    assert failed["model_count"] == 56
    assert failed["summary_wind_positive_count"] == 56
    assert failed["release_wind_table_exact_zero_count"] == 56
    assert failed["cds_terminal_wind_positive_count"] == 53
    assert failed["cds_terminal_wind_zero_count"] == 3
    assert failed["unresolved_count"] == 56
    assert report["quality_findings"]["lc18_readme_consistency_pass"] is False
    assert report["quality_findings"]["failed_wind_anomaly_resolved"] is False
    assert report["quality_findings"][
        "cross_source_difference_silently_reconciled"
    ] is False

    comparison = report["cross_source_wind_comparison"]
    assert comparison["model_count"] == 108
    assert comparison["above_nominal_half_bin_count"] == 108
    assert comparison["agreement_required_for_this_review"] is False
    assert math.isclose(
        comparison["maximum_absolute_difference_msun"],
        1.5476,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )

    history = report["phase_history"]
    assert history["model_count"] == 108
    assert history["unique_phase_row_count"] == 845
    assert history["minimum_unique_phase_count_per_model"] == 3
    assert history["maximum_unique_phase_count_per_model"] == 8
    assert history["collapsed_duplicate_row_count"] == 19
    assert history["strictly_increasing_cumulative_age_violation_count"] == 0
    assert history["nonincreasing_total_mass_violation_count"] == 0
    assert history["missing_psn_terminal_phase_count"] == 0

    structure = report["presupernova_structure"]
    assert structure["available_model_count"] == 96
    assert structure["explicit_null_model_count"] == 12
    assert len(structure["missing_coordinates"]) == 12
    assert structure["binding_energy_is_not_injected_explosion_energy"] is True

    rows = report["rows"]
    assert len(rows) == 108
    assert sum(row["exploded"] for row in rows) == 52
    assert sum(
        row["resolution"] == "unresolved_failed_wind_anomaly" for row in rows
    ) == 56
    assert sum(row["cds_presupernova_structure"] is None for row in rows) == 12
    assert report["admission"]["hard_blockers_unchanged"] is True
    assert report["admission"]["hard_blockers"] == [
        "failed_model_wind_summary_table_anomaly_requires_author_or_corrected_release",
        "age_resolved_wind_missing",
        "per_node_injected_energy_mapping_missing",
        "canonical_momentum_and_deposition_missing",
    ]
    assert report["admission"]["production_qualified"] is False
    assert report["source_identity"][
        "limongi_chieffi_cds_redistribution_status"
    ] == "no_explicit_catalogue_license_identified"
    assert report["source_identity"]["limongi_chieffi_cds_rights"][
        "authoritative_for_verdict"
    ] is False
    assert report["source_identity"]["review_use_only"] is True
    assert report["phase_history"]["phase_order_provenance"][
        "source_attested_for_intermediate_burning_order"
    ] is False

    residuals = report["cross_source_wind_comparison"]["residual_statistics"]
    assert (
        residuals["all_models"]["relative_difference_denominator"]
        == "summary_wind_mass_msun"
    )
    failed_residuals = residuals["failed_wind_anomaly"]
    assert failed_residuals["model_count"] == 56
    assert failed_residuals["positive_signed_difference_count"] == 56
    assert failed_residuals["negative_signed_difference_count"] == 0
    assert failed_residuals["zero_signed_difference_count"] == 0
    assert failed_residuals["relative_null_count_zero_denominator"] == 0

    _phase_invariant_tests()
    _differential_phase_audit_test()
    _cli_diagnostics_test()

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print("FP1_LC18_FAILED_WIND_CROSSCHECK_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
