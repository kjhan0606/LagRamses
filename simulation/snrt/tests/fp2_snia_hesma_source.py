#!/usr/bin/env python3
"""Admission checks for the review-only HESMA F-P2 source package."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_fp2_snia_hesma import audit_source  # noqa: E402


def main() -> int:
    root = ROOT.parents[1] / "assets" / "review_only" / "fp2_snia" / "hesma_yysd4_xap92"
    manifest = ROOT.parents[1] / "manifests" / "fp2_snia_hesma_yysd4_review_v1.json"
    report = audit_source(root, manifest)
    assert report["status"] == "review_only_source_format_passed", report
    assert report["source_integrity_passed"] is True
    assert report["record_id"] == "yysd4-xap92"
    assert report["record_access"]["record"] == "public"
    assert report["record_access"]["files"] == "public"
    assert report["model_count"] == 15
    assert report["canonical_conversion_allowed"] is False
    assert report["runtime_activation_allowed"] is False
    assert report["data_semantics"]["event_energy"].startswith("not supplied")
    assert report["data_semantics"]["event_momentum"].startswith("signed vector is not determined")
    n100 = report["model_reports"]["n100"]
    assert n100["abundances"]["missing_project_elements"] == []
    assert n100["isotopes"]["isotope_count"] == 384
    assert n100["isotopes"]["project_element_presence"]["H"] == "present_via_free_proton_column"
    assert n100["density"]["row_count"] == n100["isotopes"]["row_count"]
    assert n100["profile_mass_vs_integrated_abundance"]["status"] == "within_5_percent"
    assert any(warning["model"] == "n300c" for warning in report["physical_warnings"])
    assert report["physical_review_status"] == "review_only_with_physical_warnings"
    assert report["model_reports"]["n300c"]["profile_mass_vs_integrated_abundance"]["review_classification"] == "source_data_anomaly_requires_quarantine"
    assert any(
        warning["model"] == "n300c" and warning["requires_quarantine"] is True
        for warning in report["physical_warnings"]
    )

    payload = json.loads(manifest.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="snrt-fp2-hesma-") as directory:
        temporary = Path(directory)
        bad = copy.deepcopy(payload)
        bad["files"][0]["bytes"] += 1
        bad_manifest = temporary / "bad-manifest.json"
        bad_manifest.write_text(json.dumps(bad), encoding="utf-8")
        bad_report = audit_source(root, bad_manifest)
        assert bad_report["status"] == "blocked_source_format_integrity"
        assert "fingerprint_mismatch:ddt_2013_abundances.zip" in bad_report["audit_failures"]
    print("FP2_SNIa_HESMA_SOURCE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
