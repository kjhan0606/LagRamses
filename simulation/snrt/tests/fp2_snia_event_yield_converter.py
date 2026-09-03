#!/usr/bin/env python3
"""Tests for the review-only SNIa event-yield converter."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from convert_snia_event_yields import ConversionError, convert  # noqa: E402


ELEMENTS = ["H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _document(source_path: Path) -> dict:
    return {
        "source": {
            "source_id": "unit_review_source",
            "citation": "unit-only converter fixture; not a physical source",
            "source_version": "unit-v1",
            "source_path": str(source_path),
            "source_sha256": _sha256(source_path),
            "license_status": "approved",
            "provenance_status": "approved",
            "approval_id": "UNIT_ONLY_NOT_PROJECT_APPROVAL",
            "decay_convention": "fully_decayed_for_unit_test",
            "decay_horizon_yr": 1.0e9,
            "metallicity_definition": "mass_fraction",
            "population_model": "unit-only",
            "model_selection": "unit-only",
            "element_order": ELEMENTS,
        },
        "rows": [
            {
                "model_id": "model_b",
                "metallicity_mass_fraction": 0.02,
                "returned_mass_msun_per_event": 1.4,
                "remnant_mass_msun_per_event": 0.0,
                "energy_erg_per_event": 1.0e51,
                "momentum_g_cm_s_per_event": [1.0e42, -2.0e42, 3.0e42],
                "ejecta_msun_per_event": [0.1] * 11,
                "net_yield_msun_per_event": [-0.01] + [0.02] * 10,
            },
            {
                "model_id": "model_a",
                "metallicity_mass_fraction": 0.001,
                "returned_mass_msun_per_event": 1.2,
                "remnant_mass_msun_per_event": 0.0,
                "energy_erg_per_event": 8.0e50,
                "momentum_g_cm_s_per_event": [0.0, 0.0, 0.0],
                "ejecta_msun_per_event": [0.05] * 11,
                "net_yield_msun_per_event": [0.01] * 11,
            },
        ],
    }


def _write_document(document: dict, path: Path) -> None:
    path.write_text(json.dumps(document, indent=2), encoding="utf-8")


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="snrt-fp2-event-yield-") as directory:
        temporary = Path(directory)
        source_path = temporary / "source.dat"
        source_path.write_bytes(b"unit source bytes\n")
        input_path = temporary / "input.json"
        _write_document(_document(source_path), input_path)
        output_path = temporary / "event_asset.json"
        sidecar_path = temporary / "event_asset.sidecar.json"

        first = convert(input_path, output_path, sidecar_path)
        assert first["row_count"] == 2
        assert first["asset_sha256"] == _sha256(output_path)
        asset = json.loads(output_path.read_text(encoding="utf-8"))
        assert [row["model_id"] for row in asset["rows"]] == ["model_a", "model_b"]
        assert asset["event_semantics"]["quantity_basis"] == "per_event"
        assert asset["event_semantics"]["terminal_remnant_msun_per_event"] == 0.0
        assert asset["rows"][0]["untracked_ejecta_msun_per_event"] > 0.0

        second_output = temporary / "event_asset_2.json"
        second_sidecar = temporary / "event_asset_2.sidecar.json"
        convert(input_path, second_output, second_sidecar)
        assert output_path.read_bytes() == second_output.read_bytes()
        assert sidecar_path.read_bytes() == second_sidecar.read_bytes()

        with tempfile.TemporaryDirectory(prefix="snrt-fp2-event-yield-bad-") as bad_directory:
            bad = Path(bad_directory)
            candidate = copy.deepcopy(_document(source_path))
            candidate["source"]["license_status"] = "candidate"
            candidate_input = bad / "candidate.json"
            _write_document(candidate, candidate_input)
            try:
                convert(candidate_input, bad / "out.json", bad / "out.sidecar.json")
            except ConversionError as exc:
                assert "license_status" in str(exc)
            else:
                raise AssertionError("unapproved license was accepted")

            mismatched = copy.deepcopy(_document(source_path))
            mismatched["source"]["source_sha256"] = "0" * 64
            mismatch_input = bad / "mismatch.json"
            _write_document(mismatched, mismatch_input)
            try:
                convert(mismatch_input, bad / "mismatch.out.json", bad / "mismatch.sidecar.json")
            except ConversionError as exc:
                assert "source_sha256" in str(exc)
            else:
                raise AssertionError("mismatched source checksum was accepted")

            over_return = copy.deepcopy(_document(source_path))
            over_return["rows"][0]["ejecta_msun_per_event"] = [0.2] * 11
            over_return_input = bad / "over_return.json"
            _write_document(over_return, over_return_input)
            try:
                convert(over_return_input, bad / "over.out.json", bad / "over.sidecar.json")
            except ConversionError as exc:
                assert "exceeding returned" in str(exc)
            else:
                raise AssertionError("over-return event row was accepted")

            remnant = copy.deepcopy(_document(source_path))
            remnant["rows"][0]["remnant_mass_msun_per_event"] = 0.1
            remnant_input = bad / "remnant.json"
            _write_document(remnant, remnant_input)
            try:
                convert(remnant_input, bad / "remnant.out.json", bad / "remnant.sidecar.json")
            except ConversionError as exc:
                assert "zero terminal remnant" in str(exc)
            else:
                raise AssertionError("SNIa remnant was accepted")

        try:
            convert(input_path, output_path, sidecar_path)
        except ConversionError as exc:
            assert "overwrite" in str(exc)
        else:
            raise AssertionError("existing output was overwritten")
    print("FP2_SNIa_EVENT_YIELD_CONVERTER_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
