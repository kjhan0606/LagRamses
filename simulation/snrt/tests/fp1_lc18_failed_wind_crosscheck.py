#!/usr/bin/env python3
"""Checks for the review-only LC18 failed-wind cross-source audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_fp1_lc18_failed_wind_crosscheck import (  # noqa: E402
    audit_lc18_failed_wind_crosscheck,
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    before = {str(path): _sha256(path) for path in SOURCE_FILES}
    report = audit_lc18_failed_wind_crosscheck(
        root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    after = {str(path): _sha256(path) for path in SOURCE_FILES}
    assert before == after

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
    assert len(report["admission"]["hard_blockers"]) == 4
    assert report["admission"]["production_qualified"] is False
    assert report["source_identity"][
        "limongi_chieffi_cds_redistribution_status"
    ] == "no_explicit_catalogue_license_identified"
    assert report["source_identity"]["review_use_only"] is True

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print("FP1_LC18_FAILED_WIND_CROSSCHECK_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
