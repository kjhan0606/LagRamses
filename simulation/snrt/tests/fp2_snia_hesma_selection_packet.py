#!/usr/bin/env python3
"""Tests for the non-ranking HESMA source-selection packet."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from build_hesma_snia_selection_packet import build_packet  # noqa: E402


def main() -> int:
    report = build_packet(
        model_comparison_path=ROOT / "data" / "fp2_snia_hesma_model_comparison.json",
        profile_comparison_path=ROOT / "data" / "fp2_snia_hesma_profile_estimator_comparison.json",
    )
    assert report["status"] == "review_only_selection_pending"
    assert len(report["models"]) == 15
    assert report["selection"]["selected_model_id"] is None
    assert report["selection"]["selected_population_mixture"] is None
    assert report["physical_event_contract"]["energy_erg_per_event"] is None
    assert report["admission"]["canonical_conversion_allowed"] is False
    assert report["admission"]["runtime_activation_allowed"] is False
    candidate = next(row for row in report["models"] if row["model_id"] == "n100")
    assert candidate["review_screen"] == "profile_consistent_review_candidate"
    warning = next(row for row in report["models"] if row["model_id"] == "n300c")
    assert warning["review_screen"] == "physical_warning_requires_resolution"
    assert warning["quarantine_required"] is True
    assert warning["source_review_classification"] == "source_data_anomaly_requires_quarantine"
    assert sum(row["review_screen"] == "profile_consistent_review_candidate" for row in report["models"]) == 13
    print("FP2_SNIa_HESMA_SELECTION_PACKET_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
