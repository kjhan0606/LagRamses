#!/usr/bin/env python3
"""Tests for the review-only HESMA SNIa source adapter."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from adapt_hesma_snia_source import HesmaAdapterError, adapt_source  # noqa: E402


def main() -> int:
    project_root = ROOT.parents[1]
    source_root = project_root / "assets" / "review_only" / "fp2_snia" / "hesma_yysd4_xap92"
    manifest = project_root / "manifests" / "fp2_snia_hesma_yysd4_review_v1.json"

    report = adapt_source("n100", root=source_root, manifest_path=manifest)
    assert report["status"] == "review_only_source_normalized_physics_blocked"
    assert report["source"]["selected_model_id"] == "n100"
    assert report["source"]["selection_policy"].startswith("explicit_model_argument_required")
    assert report["composition"]["element_order"] == [
        "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"
    ]
    assert report["composition"]["project_stable_element_mass_msun"] > 1.0
    assert report["composition"]["non_project_stable_element_mass_msun"] > 0.0
    assert report["profiles"]["density_isotope_row_count_match"] is True
    assert report["profiles"]["isotopes"]["isotope_count"] == 384
    assert report["profiles"]["physical_warnings"] == []
    event_contract = report["event_contract"]
    assert event_contract["returned_mass_msun_per_event"] is None
    assert event_contract["energy_erg_per_event"] is None
    assert event_contract["momentum_g_cm_s_per_event"] is None
    assert report["admission"]["canonical_conversion_allowed"] is False
    assert report["admission"]["runtime_activation_allowed"] is False
    assert report["admission"]["converter_input_emitted"] is False

    try:
        adapt_source("n300c", root=source_root, manifest_path=manifest)
    except HesmaAdapterError as exc:
        assert "quarantined" in str(exc)
    else:
        raise AssertionError("quarantined HESMA model was normalized")

    try:
        adapt_source("n1600c", root=source_root, manifest_path=manifest)
    except HesmaAdapterError as exc:
        assert "unresolved physical warnings" in str(exc)
    else:
        raise AssertionError("warning-bearing HESMA model was normalized")

    for bad_model in ("", "archive_default", "n999"):
        try:
            adapt_source(bad_model, root=source_root, manifest_path=manifest)
        except HesmaAdapterError:
            pass
        else:
            raise AssertionError(f"model selection unexpectedly accepted: {bad_model!r}")

    # The adapter returns a value and does not overwrite files as a side effect.
    with tempfile.TemporaryDirectory(prefix="snrt-fp2-hesma-adapter-") as directory:
        output = Path(directory) / "sentinel.json"
        output.write_text(json.dumps(report), encoding="utf-8")
        assert json.loads(output.read_text(encoding="utf-8"))["source"]["selected_model_id"] == "n100"

    print("FP2_SNIa_HESMA_ADAPTER_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
