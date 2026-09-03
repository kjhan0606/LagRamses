#!/usr/bin/env python3
"""Tests for the review-only HESMA profile estimator comparison."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_hesma_profile_estimators import compare_estimators  # noqa: E402


def main() -> int:
    project_root = ROOT.parents[1]
    report = compare_estimators(
        root=project_root / "assets" / "review_only" / "fp2_snia" / "hesma_yysd4_xap92",
        manifest_path=project_root / "manifests" / "fp2_snia_hesma_yysd4_review_v1.json",
    )
    assert report["status"] == "review_only_diagnostic"
    assert report["model_count"] == 15
    assert report["admission"]["canonical_conversion_allowed"] is False
    assert report["admission"]["runtime_activation_allowed"] is False
    assert report["admission"]["selected_estimator"] is None
    assert set(report["estimators"]) == {
        "inner_zero_outer_half_bin",
        "half_bin_extrapolated_both_ends",
    }
    n100 = next(row for row in report["models"] if row["model_id"] == "n100")
    assert n100["row_count"] == 93
    assert n100["estimators"]["inner_zero_outer_half_bin"]["mass_estimate_msun"] > 1.4
    assert n100["estimators"]["half_bin_extrapolated_both_ends"]["mass_estimate_msun"] > 1.4
    n300c = next(row for row in report["models"] if row["model_id"] == "n300c")
    assert n300c["source_audit_physical_warnings"]
    assert n300c["estimators"]["inner_zero_outer_half_bin"][
        "relative_difference_from_integrated_stable_mass"
    ] > 1.0
    print("FP2_SNIa_HESMA_PROFILE_ESTIMATOR_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
