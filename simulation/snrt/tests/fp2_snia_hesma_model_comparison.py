#!/usr/bin/env python3
"""Tests for the review-only HESMA model comparison matrix."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from build_hesma_snia_model_comparison import build_comparison  # noqa: E402


def main() -> int:
    project_root = ROOT.parents[1]
    report = build_comparison(
        root=project_root / "assets" / "review_only" / "fp2_snia" / "hesma_yysd4_xap92",
        manifest_path=project_root / "manifests" / "fp2_snia_hesma_yysd4_review_v1.json",
    )
    assert report["status"] == "review_only_no_model_selected"
    assert report["model_count"] == 15
    assert report["selection"]["selected_model_id"] is None
    assert report["selection"]["population_mixture"] is None
    assert report["admission"]["canonical_conversion_allowed"] is False
    assert report["admission"]["runtime_activation_allowed"] is False
    assert report["admission"]["canonical_rows_emitted"] == 0
    n100 = next(row for row in report["models"] if row["model_id"] == "n100")
    assert n100["profile_mass_closure_status"] == "within_5_percent"
    assert n100["physical_warning_count"] == 0
    n300c = next(row for row in report["models"] if row["model_id"] == "n300c")
    assert n300c["profile_mass_closure_status"] == "warning_profile_mass_mismatch"
    assert n300c["physical_warning_count"] == 1
    assert all(row["selection_status"] == "not_selected" for row in report["models"])
    print("FP2_SNIa_HESMA_MODEL_COMPARISON_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
