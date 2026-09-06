#!/usr/bin/env python3
"""Reject stale or non-passing canonical AGN nine-group provenance."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.primordial import H_I_FIT, HE_I_FIT, HE_II_FIT


ARTIFACT = ROOT / "data" / "agn_nine_group_validation.json"
LEDGER = ROOT / "data" / "p4_pilot_agn_photon_ledger.csv"
METADATA = ROOT / "data" / "p4_pilot_agn_photon_ledger.json"
GROUP_EDGES = ROOT / "config" / "p0_photon_group_edges_ev.txt"
GENERATOR = ROOT / "tools" / "p4_build_agn_photon_ledger.py"
VALIDATOR = ROOT / "tools" / "validate_agn_nine_group_ledger.py"
SOURCE_REBIND_TOOL = ROOT / "tools" / "p4_attach_pilot_sources.py"
P4_RUNNER = ROOT / "tools" / "p4_run_transport_pilot.py"
STATIC_INPUT = ROOT / "data" / "p4_coeval_static_rt_input_agn9.h5"
STATIC_METADATA = ROOT / "data" / "p4_coeval_static_rt_input_agn9.json"
TRANSPORT_CONTROL = ROOT / "data" / "p4_validation" / "p4_agn9_stage4_0p001myr.h5"
EXTERNAL_ASSET_MANIFEST = ROOT / "data" / "agn_nine_group_external_assets.json"
ATTESTATION_SCOPE = "simulation/snrt"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_finite(value: object) -> None:
    if isinstance(value, dict):
        for child in value.values():
            assert_finite(child)
    elif isinstance(value, list):
        for child in value:
            assert_finite(child)
    elif isinstance(value, float):
        assert math.isfinite(value)


def main() -> int:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    metadata = json.loads(METADATA.read_text(encoding="utf-8"))
    assert_finite(payload)
    assert payload["schema"] == "snrt_agn_nine_group_validation_v1"
    assert payload["passed"] is True
    assert all(payload["criteria"].values())

    configured_edges = np.loadtxt(GROUP_EDGES, comments="#")
    assert np.array_equal(
        payload["configuration"]["group_edges_ev"], configured_edges
    )
    assert np.array_equal(metadata["group_edges_ev"], configured_edges)
    assert payload["configuration"]["number_of_groups"] == 9
    assert (
        payload["configuration"]["interval_convention"]
        == "left_closed_right_open_except_final_closed"
    )
    assert metadata["group_interval_convention"] == payload["configuration"][
        "interval_convention"
    ]

    hard = payload["hard_xray_diagnostics"]
    assert hard["energy_interval_ev"] == [2000.0, 10000.0]
    assert hard["photon_rate_s"] > 0.0
    assert 0.1 < hard["photon_number_ratio_to_0p5_2kev"] < 0.4
    assert 0.5 < hard["energy_power_ratio_to_0p5_2kev"] < 2.0
    assert 0.0 < hard["fraction_of_total_photon_rate"] < 0.01
    assert 0.0 < hard["fraction_of_supported_sed_energy_power"] < 0.1
    assert 0.0 < hard["fraction_of_candidate_bolometric_luminosity"] < 0.02
    assert all(value > 0.0 for value in hard["cross_sections_cm2"])
    assert all(value > 0.0 for value in hard["photoelectron_excess_energy_ev"])

    transport = payload["transport_control"]
    assert transport["number_of_groups"] == 9
    assert transport["duration_myr"] == 0.001
    assert transport["validation_passed"] is True
    assert transport["photon_ledger_relative_error"] < 1.0e-12
    assert transport["hydrogen_ledger_l1_relative_error"] < 1.0e-3
    assert transport["maximum_fixed_point_residual"] < 1.0e-4
    assert transport["photoelectron_energy_ledger_l1_relative_error"] < 1.0e-12
    assert transport["electron_root_bracket_failure_count"] == 0

    intervals = np.column_stack((configured_edges[:-1], configured_edges[1:]))
    hard_matches = np.flatnonzero(np.all(intervals == (2000.0, 10000.0), axis=1))
    assert hard_matches.size == 1
    assert int(hard_matches[0]) == len(intervals) - 1
    closure = metadata["group_spectral_closure"]["cross_sections_cm2"]
    species_closure = (
        closure["hydrogen_i"],
        closure["helium_i"],
        closure["helium_ii"],
    )
    for values, threshold in zip(
        species_closure,
        (H_I_FIT.threshold_ev, HE_I_FIT.threshold_ev, HE_II_FIT.threshold_ev),
        strict=True,
    ):
        wholly_below = configured_edges[1:] <= threshold
        assert np.all(np.asarray(values)[wholly_below] == 0.0)

    expected_status = []
    for low, high in intervals:
        if high <= payload["configuration"]["sed_support_minimum_ev"]:
            expected_status.append("agn_sed_below_support_zero_photons")
        elif low < payload["configuration"]["sed_support_minimum_ev"]:
            expected_status.append("agn_sed_partially_supported_10ev_to_upper")
        else:
            expected_status.append("agn_sed_fully_supported")
    assert [group["closure_status"] for group in metadata["groups"]] == expected_status
    assert payload["ledger"]["source_count"] == 10
    assert payload["ledger"]["luminous_source_count"] == 5

    provenance = payload["provenance"]
    working_tree_status = subprocess.check_output(
        ("git", "status", "--short", "--untracked-files=all", "--", ATTESTATION_SCOPE),
        cwd=REPOSITORY_ROOT,
        text=True,
    )
    current_head = subprocess.check_output(
        ("git", "rev-parse", "HEAD"), cwd=REPOSITORY_ROOT, text=True
    ).strip()
    assert provenance["git_head"] == current_head
    assert provenance["working_tree_attestation_scope"] == ATTESTATION_SCOPE
    assert provenance["working_tree_clean"] is (working_tree_status == "")
    assert provenance["working_tree_status_sha256"] == hashlib.sha256(
        working_tree_status.encode("utf-8")
    ).hexdigest()
    for key, path in (
        ("validator_sha256", VALIDATOR),
        ("generator_sha256", GENERATOR),
        ("group_edges_sha256", GROUP_EDGES),
        ("ledger_sha256", LEDGER),
        ("metadata_sha256", METADATA),
        ("primordial_sha256", ROOT / "snrt_core" / "primordial.py"),
        ("source_ledger_sha256", ROOT / "snrt_core" / "source_ledger.py"),
        ("source_rebind_tool_sha256", SOURCE_REBIND_TOOL),
        ("p4_runner_sha256", P4_RUNNER),
        ("static_input_sha256", STATIC_INPUT),
        ("static_metadata_sha256", STATIC_METADATA),
        ("transport_control_sha256", TRANSPORT_CONTROL),
        ("external_asset_manifest_sha256", EXTERNAL_ASSET_MANIFEST),
    ):
        assert provenance[key] == sha256(path)
    candidates = Path(metadata["candidates"])
    assert provenance["candidates_sha256"] == sha256(candidates)

    with TemporaryDirectory(prefix="agn-pilot-artifact-test-") as directory:
        output = Path(directory) / "validation.json"
        result = subprocess.run(
            [
                sys.executable,
                str(VALIDATOR),
                "--output",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        fresh = json.loads(output.read_text(encoding="utf-8"))
        assert fresh["passed"] is True
        assert all(fresh["criteria"].values())
        assert fresh["configuration"]["source_mode"] == "pilot"
        assert fresh["criteria"] == payload["criteria"]

    print(
        "AGN_NINE_GROUP_ARTIFACT_OK "
        f"hard_q={hard['photon_rate_s']:.6g} "
        f"hard_to_soft_q={hard['photon_number_ratio_to_0p5_2kev']:.6g} "
        f"hard_supported_sed_fraction={hard['fraction_of_supported_sed_energy_power']:.6g} "
        f"hard_bolometric_fraction={hard['fraction_of_candidate_bolometric_luminosity']:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
