#!/usr/bin/env python3
"""Tests for deterministic promotion of the approved HESMA n100 baseline."""

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

from promote_hesma_snia_source import PromotionError, promote  # noqa: E402


APPROVAL = "FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ"
COMMIT = "c6c8042b03406b9d69bc50434fe5d6af7f542be6"


def main() -> int:
    input_path = ROOT / "data" / "fp2_snia_hesma_n100_review_normalized.json"
    with tempfile.TemporaryDirectory(prefix="snrt-fp2-snia-promote-") as directory:
        output_path = Path(directory) / "approved.json"
        payload = promote(input_path, output_path, approval_id=APPROVAL, source_commit_binding=COMMIT)
        assert payload["status"] == "approved_physical_baseline_runtime_gated"
        assert payload["source"]["model_id"] == "n100"
        assert payload["event"]["returned_mass_msun_per_event"] == 1.4004633930489443
        assert payload["event"]["terminal_remnant_msun_per_event"] == 0.0
        assert payload["event"]["momentum_g_cm_s_per_event"] == [0.0, 0.0, 0.0]
        assert payload["event"]["untracked_ejecta_msun_per_event"] > 0.0
        assert payload["approval"]["runtime_activation_allowed"] is False
        assert payload["conversion"]["conversion_code_sha256"]

        bad_input = Path(directory) / "bad.json"
        bad = json.loads(input_path.read_text(encoding="utf-8"))
        bad["profiles"]["physical_warnings"] = [{"model": "n100", "severity": "test"}]
        bad_input.write_text(json.dumps(bad), encoding="utf-8")
        try:
            promote(bad_input, Path(directory) / "bad-output.json", approval_id=APPROVAL, source_commit_binding=COMMIT)
        except PromotionError as exc:
            assert "warnings" in str(exc)
        else:
            raise AssertionError("warning-bearing HESMA source was promoted")

        second_output = Path(directory) / "second.json"
        promote(input_path, second_output, approval_id=APPROVAL, source_commit_binding=COMMIT)
        assert output_path.read_bytes() == second_output.read_bytes()

    print("FP2_SNIa_HESMA_PROMOTION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
