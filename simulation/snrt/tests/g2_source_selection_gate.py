#!/usr/bin/env python3
"""Tests for the G2 review-only source-selection gate."""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import tempfile
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
import sys

if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_source_selection_gate import audit_selection  # noqa: E402


def _write(payload: dict, directory: Path, name: str) -> Path:
    path = directory / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-parked-agb", action="store_true",
                        help="Explicitly rerun parked KL16/CK22 source-reader checks; not physical approval")
    args = parser.parse_args()
    matrix_path = ROOT / "config" / "g2_source_selection_matrix_v1.json"
    fingerprint_path = ROOT / "data" / "g2_source_package_fingerprint_audit.json"
    report = audit_selection(matrix_path, fingerprint_path)
    assert report["status"] == "review_only_validation_branch_recorded", report
    assert report["production_ready"] is False
    assert report["runtime_activation_allowed"] is False
    assert report["review_validation_branch"]["candidate_id"] == "sukhbold2016_ccsn"
    assert report["review_validation_branch"]["approval_id"] is None
    assert report["production_source_id"] is None
    assert report["production_approval_id"] is None
    assert report["fingerprint_input_integrity_passed"] is True
    assert report["fingerprint_candidate_count"] == 11
    assert report["fingerprint_file_count"] == 65
    assert "review_validation_branch_is_not_a_production_source" in report["blockers"]

    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    if args.include_parked_agb:
        from read_karakas_lugaro2016 import read_karakas_lugaro2016, _blocks
        from read_cinquegrana_karakas2022 import read_cinquegrana_karakas2022, _read_table
        from adapt_g2_candidate_sources import SourceAdapterError

        ck22 = read_cinquegrana_karakas2022(matrix_path)
        assert len(ck22["records"]) == 103
        assert ck22["canonical_rows_emitted"] == 0 and not ck22["runtime_activation_allowed"]
        assert not any(ck22[k] for k in ("renormalization_applied", "negative_residue_clipping_applied", "decay_applied"))
        assert ck22["lifetime_yr"] is None and ck22["injected_energy_erg"] is None
        assert len(ck22["evolution_source_rows"]) == 121
        assert sum(len(m["evolution_mass_Z_matches"]) == 1 for m in ck22["records"]) == 102
        assert all(m["release_time_yr"] is None for m in ck22["records"])
        duplicate = next(m for m in ck22["records"] if m["initial_mass_msun"] == 5 and m["metallicity_mass_fraction"] == .08)
        assert [e["stellar_duration_myr"] for e in duplicate["evolution_mass_Z_matches"]] == ["68.0", "67.99"]
        assert duplicate["evolution_duration_candidate_yr"] is None
        slow = next(m for m in ck22["records"] if m["initial_mass_msun"] == 1 and m["metallicity_mass_fraction"] == .04)
        assert slow["evolution_duration_candidate_yr"] == 15.2e9  # No cosmological age clipping.
        assert slow["evolution_mass_Z_matches"][0]["core_mass_first_tp_msun"] == "0.573"
        assert slow["remnant_mass_msun"] is None  # Mc(1) is NOT a final remnant.
        no_tp = next(m for m in ck22["records"] if m["initial_mass_msun"] == 1 and m["metallicity_mass_fraction"] == .05)
        assert no_tp["evolution_mass_Z_matches"][0]["thermal_pulses"] == "0"
        assert no_tp["evolution_mass_Z_matches"][0]["core_mass_first_tp_msun"] == ""
        assert no_tp["evolution_endpoint_caveat"] is None  # Zero TP alone is not failure.
        incomplete = [m for m in ck22["records"] if m["evolution_endpoint_caveat"] == "KCJ22_section_4.2.1_explicit_early_termination"]
        assert {(m["initial_mass_msun"], m["metallicity_mass_fraction"]) for m in incomplete} == {(8., .09), (7., .1), (8., .1)}
        negative_gross = []
        for model in ck22["records"]:
            species = model["species_by_source_label"]
            assert model["missing_tracked_elements"] == ["Ca"]
            assert model["returned_mass_msun"] is None and model["remnant_mass_msun"] is None
            assert species["p"]["atomic_number"] == species["d"]["atomic_number"] == 1
            assert species["n"]["atomic_number"] is None and species["n14"]["atomic_number"] == 7
            assert species["g"]["atomic_number"] is None and species["p31"]["atomic_number"] == 15
            assert species["al-6"]["mass_number_label"] == species["al*6"]["mass_number_label"] == 26
            # Neither of these source-column denominators is silently chosen as M_ej.
            delta = model["H_gross_over_wind_fraction_msun"] - model["H_initial_wind_mass_over_X0_msun"]
            assert 9.9e-5 < delta < 4.01e-4
            negative_gross.extend(r["gross_mass_msun"] for r in species.values() if r["gross_mass_msun"] < 0)
            assert all(abs(r["net_identity_residual_msun"]) < 9.01e-8 for r in species.values())
        assert len(negative_gross) == 8 and min(negative_gross) == -6.9439942e-29
        ck_source = next(c for c in matrix["candidates"] if c["candidate_id"] == "cinquegrana_karakas2022_agb")
        sample_text = (Path(ck_source["source_asset_path"]) / "z04/yields_m1z04.dat").read_text()
        sample_ck = _read_table(sample_text)["species_by_source_label"]
        assert sample_ck["p"]["net_mass_msun"] == -7.6089774e-3
        assert sample_ck["p"]["gross_mass_msun"] == 2.5361665e-1
        read_bytes = Path.read_bytes
        def changed_source_bytes(path):
            data = read_bytes(path)
            return data + b"\n" if path.name == "yields_m1z04.dat" else data
        with patch.object(Path, "read_bytes", changed_source_bytes):
            try:
                read_cinquegrana_karakas2022(matrix_path)
            except SourceAdapterError as error:
                assert "fingerprint mismatch" in str(error)
            else:
                raise AssertionError("CK22 changed worktree bytes were accepted")
        # Atomic mass A is not atomic number; do not collapse separate Al26 states.
        for bad_text in (sample_text.replace("al*6", "al-6"),
                         sample_text.replace("2.5361665E-01", "nan"),
                         sample_text.replace("mass(i)_lost", "unknown")):
            try:
                _read_table(bad_text)
            except SourceAdapterError:
                pass
            else:
                raise AssertionError("invalid CK22 table was accepted")
        kl16 = read_karakas_lugaro2016(matrix_path)
        rows = kl16["records"]
        assert len(rows) == 62
        assert len(kl16["excluded_commented_records"]) == 2
        assert {(r["coordinate"]["initial_mass_msun"], r["coordinate"]["metallicity_mass_fraction"])
                for r in kl16["excluded_commented_records"]} == {(8.0, 0.014), (8.0, 0.03)}
        assert sum(r["net_yield_diagnostic_msun"] is not None for r in rows) == 46
        assert sum(value < 0 for r in rows if r["net_yield_diagnostic_msun"] is not None
                   for value in r["net_yield_diagnostic_msun"].values()) == 1201
        hydrogen = next(r for r in rows if r["coordinate"]["initial_mass_msun"] == 1.0
                        and r["coordinate"]["metallicity_mass_fraction"] == 0.014)
        assert math.isclose(hydrogen["net_yield_diagnostic_msun"][1], -0.008273595,
                            rel_tol=0, abs_tol=1e-12)
        assert kl16["canonical_rows_emitted"] == 0
        assert kl16["runtime_activation_allowed"] is False
        assert kl16["renormalization_applied"] is False
        assert kl16["lifetime_yr"] is None and kl16["injected_energy_erg"] is None
        assert all(len(r["elements_by_atomic_number"]) == 78 for r in rows)
        assert all(r["elements_by_atomic_number"][1]["source_symbol"] == "p"
                   and r["elements_by_atomic_number"][15]["source_symbol"] == "p" for r in rows)
        sample = next(r for r in rows if r["coordinate"]["initial_mass_msun"] == 3.5
                      and r["coordinate"]["metallicity_mass_fraction"] == 0.03)
        # Independently visible Table 7 values in Karakas & Lugaro (2016), p.14.
        for z, gross in {1: 1.76149, 2: 0.921855, 6: 0.0294504,
                         7: 0.0107646, 8: 0.0319316, 9: 9.20429e-6}.items():
            assert sample["elements_by_atomic_number"][z]["gross_mass_msun"] == gross
        disagreements = [r for r in rows if r["auxiliary_minus_header_final_mass_msun"] != 0.0]
        assert len(disagreements) == 1
        assert disagreements[0]["coordinate"]["initial_mass_msun"] == 4.0
        assert disagreements[0]["coordinate"]["metallicity_mass_fraction"] == 0.03
        assert disagreements[0]["final_mass_msun"] == 0.774
        assert disagreements[0]["auxiliary_final_mass_msun"] == 0.744
        assert math.isclose(max(r["gross_sum_minus_labelled_expelled_msun"] for r in rows),
                            0.03361670011913276, abs_tol=1e-12)
        for initial in kl16["initial_composition_records"]:
            assert len(initial["elements_by_atomic_number"]) == 81
        missing = next(r for r in rows if r["coordinate"]["initial_mass_msun"] == 7.5)
        assert missing["initial_composition_matching_lines"] == []
        assert missing["net_yield_diagnostic_msun"] is None
    fingerprints = json.loads(fingerprint_path.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="snrt-g2-selection-") as directory:
        temporary = Path(directory)
        if args.include_parked_agb:
            bad_evolution = copy.deepcopy(matrix)
            next(c for c in bad_evolution["candidates"] if c["candidate_id"] == "cinquegrana_karakas2022_agb")["evolution_reference"]["sha256"] = "0" * 64
            try:
                read_cinquegrana_karakas2022(_write(bad_evolution, temporary, "bad-evolution.json"))
            except SourceAdapterError as error:
                assert "evolution table fingerprint mismatch" in str(error)
            else:
                raise AssertionError("CK22 evolution fingerprint was not enforced")
            bad_ck = copy.deepcopy(matrix)
            next(c for c in bad_ck["candidates"] if c["candidate_id"] == "cinquegrana_karakas2022_agb")["source_file_sha256"]["readme_yields.dat"] = "0" * 64
            try:
                read_cinquegrana_karakas2022(_write(bad_ck, temporary, "bad-ck22.json"))
            except SourceAdapterError:
                pass
            else:
                raise AssertionError("CK22 input fingerprint was not enforced")
            tampered = copy.deepcopy(matrix)
            next(c for c in tampered["candidates"] if c["candidate_id"] == "karakas_lugaro2016_agb")["source_file_sha256"]["yield_z007.txt"] = "0" * 64
            try:
                read_karakas_lugaro2016(_write(tampered, temporary, "bad-kl16.json"))
            except SourceAdapterError:
                pass
            else:
                raise AssertionError("KL16 input fingerprint was not enforced")
            partial = "# Initial mass = 1.0, Z = 0.007, Y = 0.260, M_mix = 0\n# Final mass = 0.6, Mass expelled = 0.4\np 1 12 0 0 0.7 0.28\n# he 2 11 0 0 0.3 0.12\n"
            try:
                _blocks(partial, initial=False)
            except SourceAdapterError:
                pass
            else:
                raise AssertionError("partly commented KL16 node was activated")
        missing_branch = copy.deepcopy(matrix)
        missing_branch["review_selection"]["validation_branch"] = "not_staged"
        missing_matrix_path = _write(missing_branch, temporary, "matrix.json")
        failed = audit_selection(missing_matrix_path, fingerprint_path)
        assert failed["status"] == "review_selection_blocked_input_integrity"
        assert any(item["reason"] == "validation_branch_not_in_matrix" for item in failed["audit_failures"])

        overclaim = copy.deepcopy(matrix)
        overclaim["review_selection"]["production_source_id"] = "sukhbold2016_ccsn"
        overclaim["review_selection"]["production_approval_id"] = "UNAUTHORIZED"
        overclaim_path = _write(overclaim, temporary, "overclaim.json")
        failed = audit_selection(overclaim_path, fingerprint_path)
        assert failed["status"] == "review_selection_blocked_input_integrity"
        assert any(
            item["reason"] == "production_source_must_remain_unselected_until_physics_approval"
            for item in failed["audit_failures"]
        )

        bad_fingerprints = copy.deepcopy(fingerprints)
        bad_fingerprints["status"] = "candidate_fingerprint_blocked_input_integrity"
        bad_fingerprint_path = _write(bad_fingerprints, temporary, "bad-fingerprints.json")
        failed = audit_selection(matrix_path, bad_fingerprint_path)
        assert failed["status"] == "review_selection_blocked_input_integrity"
        assert any(item["reason"] == "fingerprint_audit_not_clean" for item in failed["audit_failures"])

    print("PARKED_AGB_SOURCE_CHECKS_" + ("RAN" if args.include_parked_agb else "NOT_REQUESTED"))
    print("G2_SOURCE_SELECTION_GATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
