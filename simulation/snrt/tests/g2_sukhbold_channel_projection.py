#!/usr/bin/env python3
"""Checks for the fail-closed Sukhbold component projection."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from build_g2_sukhbold_channel_projection import (  # noqa: E402
    SukhboldProjectionError,
    build_sukhbold_channel_projection,
)


def _assert_component_label_mutation_is_rejected() -> None:
    contract_path = ROOT / "config" / "g2_sukhbold_channel_projection_contract_v1.json"
    payload = json.loads(contract_path.read_text(encoding="utf-8"))
    payload["component_projection"]["wind"]["source_column"] = "ejecta"
    with tempfile.TemporaryDirectory(prefix="snrt-g2-projection-") as directory:
        mutated = Path(directory) / contract_path.name
        mutated.write_text(json.dumps(payload), encoding="utf-8")
        try:
            build_sukhbold_channel_projection(
                root=ROOT.parents[1] / "external" / "g2_candidates",
                contract_path=mutated,
            )
        except SukhboldProjectionError as exc:
            assert "contract identity" in str(exc)
        else:
            raise AssertionError("component/source-column mutation was accepted")


def _assert_contract_safety_mutations_are_rejected() -> None:
    contract_path = ROOT / "config" / "g2_sukhbold_channel_projection_contract_v1.json"
    mutations = []
    payload = json.loads(contract_path.read_text(encoding="utf-8"))
    missing_firewall = json.loads(json.dumps(payload))
    missing_firewall["semantic_firewalls"].pop("source_age_history_may_be_invented")
    mutations.append(missing_firewall)
    empty_prerequisites = json.loads(json.dumps(payload))
    empty_prerequisites["approval"]["required_before_approval"] = []
    mutations.append(empty_prerequisites)
    missing_component = json.loads(json.dumps(payload))
    missing_component["component_projection"].pop("wind")
    mutations.append(missing_component)
    bad_channel = json.loads(json.dumps(payload))
    bad_channel["component_projection"]["wind"]["proposed_runtime_channel"] = 3
    mutations.append(bad_channel)
    bad_ownership = json.loads(json.dumps(payload))
    bad_ownership["component_projection"]["wind"]["ownership"] = "terminal supernova ejecta only"
    mutations.append(bad_ownership)
    with tempfile.TemporaryDirectory(prefix="snrt-g2-projection-safety-") as directory:
        for index, mutated_payload in enumerate(mutations):
            mutated = Path(directory) / f"contract-{index}.json"
            mutated.write_text(json.dumps(mutated_payload), encoding="utf-8")
            try:
                build_sukhbold_channel_projection(
                    root=ROOT.parents[1] / "external" / "g2_candidates",
                    contract_path=mutated,
                )
            except SukhboldProjectionError:
                continue
            raise AssertionError(f"projection safety mutation {index} was accepted")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = build_sukhbold_channel_projection(
        root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    assert report["status"] == "review_only_blocked_decay_age_boundary_and_approval"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["model_count"] == 13
    assert report["record_count"] == 26
    assert report["record_count_by_source_component"] == {"wind": 13, "ejecta": 13}
    assert report["tracked_elements"] == [
        "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"
    ]
    assert all(report["source_nulls_preserved"].values())
    records = report["records"]
    assert all(record["canonical_row_emitted"] is False for record in records)
    assert all(record["release_age_yr"] is None for record in records)
    assert all(record["decay_complete_returned_mass_msun"] is None for record in records)
    assert all(record["canonical_scalar_launch_momentum_g_cm_s"] is None for record in records)
    for record in records:
        assert math.isclose(
            sum(record["stable_mass_by_tracked_element_msun"].values())
            + record["untracked_stable_component_mass_msun"],
            record["stable_component_mass_msun"],
            abs_tol=1.0e-12,
        )
        assert len(record["selected_radioactive_inventory_msun"]) == 20
        if record["source_component"] == "wind":
            assert record["final_kinetic_energy_erg"] is None
            assert record["fallback_mass_msun"] is None
        else:
            assert record["final_kinetic_energy_erg"] > 0.0
            assert record["fallback_mass_msun"] >= 0.0
    assert report["source_identity"]["file_count"] == 4
    assert all(len(value["sha256"]) == 64 for value in report["source_identity"]["files"].values())
    assert report["high_mass_record_count"] == 19
    assert report["high_mass_record_count_by_source_component"] == {"wind": 13, "ejecta": 6}
    assert all(record["canonical_row_emitted"] is False for record in report["high_mass_records"])
    assert all(record["runtime_channel_assignment_approved"] is False for record in report["high_mass_records"])
    assert {record["engine"] for record in report["high_mass_records"]} == {"W18", "N20"}
    assert report["high_mass_record_count_by_engine_and_source_component"] == {
        "W18": {"wind": 9, "ejecta": 2},
        "N20": {"wind": 4, "ejecta": 4},
    }
    assert report["high_mass_missing_source_masses_by_engine"] == {
        "W18": [40.0, 45.0, 50.0, 55.0, 70.0, 80.0, 100.0],
        "N20": [40.0, 45.0, 50.0, 55.0, 70.0],
    }
    for record in records + report["high_mass_records"]:
        assert record["source_cross_segment_duplicate_isotopes"] == ["k40"]
        assert math.isclose(
            sum(record["stable_mass_by_tracked_element_msun"].values())
            + record["untracked_stable_component_mass_msun"],
            record["stable_component_mass_msun"],
            abs_tol=1.0e-12,
        )
    assert all(record["source_branch"] in {"Z9.6", "W18", "N20", "implosions_W18"} for record in report["high_mass_records"])
    _assert_component_label_mutation_is_rejected()
    _assert_contract_safety_mutations_are_rejected()
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_SUKHBOLD_CHANNEL_PROJECTION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
