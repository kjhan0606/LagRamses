#!/usr/bin/env python3
"""Regression test for the fail-closed P0.1 source parity gate."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from validate_stellar_source_parity import (  # noqa: E402
    audit,
    binary_linkage_from_text,
    git_head,
    git_is_ancestor,
    production_build_log_contract,
    production_smoke_contract,
)


def main() -> int:
    payload = audit()
    assert payload["criteria"]["production_makefile_selects_declared_patch"]
    assert payload["criteria"]["production_objects_are_listed"]
    assert payload["criteria"]["production_build_inputs_resolve"]
    assert payload["criteria"]["production_stellar_objects_match_contract"]
    assert payload["criteria"]["production_sources_resolve_under_patch"]
    assert payload["criteria"]["production_source_manifest_exists_and_matches"]
    assert payload["criteria"]["runner_source_root_is_canonical_production"]
    assert payload["criteria"]["runner_production_sources_resolve"]
    assert payload["criteria"]["runner_production_objects_resolve"]
    assert payload["criteria"]["runner_production_sources_match_objects"]
    assert payload["criteria"]["runner_fixture_root_is_declared"]
    assert payload["criteria"]["runner_fixture_sources_resolve"]
    assert payload["criteria"]["runner_has_required_parity_gate"]
    assert payload["criteria"]["source_of_truth_is_canonical_production_tree"]
    assert payload["criteria"]["compile_parameters_match_contract"]
    assert payload["criteria"]["required_compile_flags_present"]
    assert payload["criteria"]["embedded_yields_disabled_by_default"]
    assert payload["criteria"]["production_linked_harness_exists"]
    assert payload["criteria"]["production_linked_harness_targets_bin_makefile"]
    assert payload["criteria"]["production_linked_harness_runs_smoke"]
    assert payload["criteria"]["production_linked_harness_object_coverage"]
    assert payload["criteria"]["production_source_manifest_order_matches_make_objects"]
    assert payload["criteria"]["production_makefile_consumes_source_manifest"]
    assert payload["criteria"]["production_make_database_resolves"]
    assert payload["criteria"]["production_shared_contract_profile_bounded"]
    assert payload["criteria"]["native_only_mirror_sources_declared"]
    assert payload["criteria"]["production_link_recipe_is_complete"]
    assert payload["criteria"]["runner_diagnostic_escape_is_fail_closed"]
    if payload["status"] == "blocked":
        assert not payload["criteria"]["production_linked_build_evidence"]
        assert "production_linked_build_evidence" in payload["blocking_reasons"]

    config = json.loads(
        (ROOT / "config/stellar_source_identity_v1.json").read_text(
            encoding="utf-8"
        )
    )
    production = config["production"]
    harness = config["production_linked_harness"]
    valid_lines = [f"P0_BUILD_COMMAND {harness['build_command']}"]
    for object_name in production["required_objects"]:
        source_name = object_name[:-2] + ".f90"
        valid_lines.append(
            "mpiifx -O3 -fpp -DNVECTOR=500 -DNDIM=3 -DNPRE=8 "
            "-DNENER=0 -DNVAR=18 -DSOLVER=hydro "
            "-DPHASE0_STELLAR_ENRICHMENT -DLONGINT -DQUADHILBERT "
            "-DOUTPUT_PARTICLE_POTENTIAL -DUSE_FFTW "
            f"-c ../patch/lagRamses/{source_name} -o {object_name}"
        )
    valid_lines.append("mpiifx stellar_objects ramses.o -o ramses_final3d -lfftw3")
    valid_log = "\n".join(valid_lines) + "\n"
    valid_contract = production_build_log_contract(
        valid_log,
        production["required_objects"],
        production["compile_parameters"],
        production["required_compile_flags"],
        production["embedded_yield_policy"]["macro"],
        harness["build_command"],
        str(Path("..") / production["patch_root"]),
        production["compile_policy"]["required_optimization_flag"],
        tuple(production["compile_policy"]["forbidden_compile_flags"]),
        harness["link_output"],
    )
    assert valid_contract["status"] == "pass"
    assert valid_contract["forced_rebuild"] is True
    assert valid_contract["compile_parameters"] == {
        name: str(value) for name, value in production["compile_parameters"].items()
    }
    binary_sha256 = "a" * 64
    assert (
        production_build_log_contract(
            valid_log + f"P0_BINARY_SHA256={binary_sha256}\n",
            production["required_objects"],
            production["compile_parameters"],
            production["required_compile_flags"],
            production["embedded_yield_policy"]["macro"],
            harness["build_command"],
            str(Path("..") / production["patch_root"]),
            production["compile_policy"]["required_optimization_flag"],
            tuple(production["compile_policy"]["forbidden_compile_flags"]),
            harness["link_output"],
            binary_sha256,
        )["status"]
        == "pass"
    )
    assert (
        production_build_log_contract(
            valid_log + f"P0_BINARY_SHA256={'b' * 64}\n",
            production["required_objects"],
            production["compile_parameters"],
            production["required_compile_flags"],
            production["embedded_yield_policy"]["macro"],
            harness["build_command"],
            str(Path("..") / production["patch_root"]),
            production["compile_policy"]["required_optimization_flag"],
            tuple(production["compile_policy"]["forbidden_compile_flags"]),
            harness["link_output"],
            binary_sha256,
        )["status"]
        == "blocked"
    )
    assert (
        production_build_log_contract(
            valid_log,
            production["required_objects"],
            production["compile_parameters"],
            production["required_compile_flags"],
            production["embedded_yield_policy"]["macro"],
            harness["build_command"],
            str(Path("..") / production["patch_root"]),
            production["compile_policy"]["required_optimization_flag"],
            tuple(production["compile_policy"]["forbidden_compile_flags"]),
            harness["link_output"],
            binary_sha256,
        )["status"]
        == "blocked"
    )

    negative_cases = {
        "missing_object": "stellar_ramses_runtime.o",
        "missing_build_command_marker": "P0_BUILD_COMMAND",
        "missing_phase0_define": "-DPHASE0_STELLAR_ENRICHMENT",
        "embedded_yields_enabled": "-DSTELLAR_EMBEDDED_YIELDS",
        "missing_required_flag": "-DOUTPUT_PARTICLE_POTENTIAL",
        "wrong_compile_parameter": "-DNDIM=3",
        "missing_link_output": "-o ramses_final3d",
    }
    for case, needle in negative_cases.items():
        negative_log = valid_log
        if case == "missing_object":
            negative_log = "\n".join(
                line
                for line in valid_log.splitlines()
                if needle not in line
            ) + "\n"
        elif case == "missing_build_command_marker":
            negative_log = valid_log.replace(needle + " " + harness["build_command"], "")
        elif case == "missing_phase0_define":
            negative_log = valid_log.replace(needle, "", 1)
        elif case == "embedded_yields_enabled":
            negative_log = valid_log.replace(
                "-DPHASE0_STELLAR_ENRICHMENT",
                "-DPHASE0_STELLAR_ENRICHMENT -DSTELLAR_EMBEDDED_YIELDS",
                1,
            )
        elif case == "missing_required_flag":
            negative_log = valid_log.replace(needle, "", 1)
        elif case == "wrong_compile_parameter":
            negative_log = valid_log.replace(needle, "-DNDIM=4", 1)
        elif case == "missing_link_output":
            negative_log = valid_log.replace(
                " -o ramses_final3d -lfftw3", " -o wrong_binary -lfftw3"
            )
        contract = production_build_log_contract(
            negative_log,
            production["required_objects"],
            production["compile_parameters"],
            production["required_compile_flags"],
            production["embedded_yield_policy"]["macro"],
            harness["build_command"],
            str(Path("..") / production["patch_root"]),
            production["compile_policy"]["required_optimization_flag"],
            tuple(production["compile_policy"]["forbidden_compile_flags"]),
            harness["link_output"],
        )
        assert contract["status"] == "blocked", case

    valid_smoke = (
        f"P0_SMOKE_COMMAND {harness['smoke_command']}\n"
        f"P0_BINARY_SHA256={binary_sha256}\n"
        " patch dir    = ../patch/lagRamses\n"
        f" last commit  = {git_head()}\n"
        "You should type: ramses3d input.nml [nrestart]\n"
        "P0_SMOKE_EXIT_CODE=0\n"
    )
    assert (
        production_smoke_contract(
            valid_smoke,
            harness["smoke_command"],
            harness["smoke_expected_output"],
            harness["smoke_expected_exit_code"],
            "../patch/lagRamses",
            git_head(),
            binary_sha256,
        )["status"]
        == "pass"
    )
    assert (
        production_smoke_contract(
            valid_smoke.replace("P0_SMOKE_EXIT_CODE=0", "P0_SMOKE_EXIT_CODE=1"),
            harness["smoke_command"],
            harness["smoke_expected_output"],
            harness["smoke_expected_exit_code"],
            "../patch/lagRamses",
            git_head(),
            binary_sha256,
        )["status"]
        == "blocked"
    )
    assert (
        production_smoke_contract(
            valid_smoke.replace(f"P0_BINARY_SHA256={binary_sha256}\n", ""),
            harness["smoke_command"],
            harness["smoke_expected_output"],
            harness["smoke_expected_exit_code"],
            "../patch/lagRamses",
            git_head(),
            binary_sha256,
        )["status"]
        == "blocked"
    )
    assert (
        production_smoke_contract(
            valid_smoke.replace("patch dir    = ../patch/lagRamses", "patch dir    = ../patch/wrong"),
            harness["smoke_command"],
            harness["smoke_expected_output"],
            harness["smoke_expected_exit_code"],
            "../patch/lagRamses",
            git_head(),
            binary_sha256,
        )["status"]
        == "blocked"
    )
    assert git_is_ancestor(git_head(), git_head())
    assert not git_is_ancestor("0" * 40, git_head())
    thermal_pattern = harness["linkage_symbol_patterns"][-1]
    assert re.search(thermal_pattern, "000000 T thermal_feedback_", re.IGNORECASE)
    assert not re.search(thermal_pattern, "000000 T unrelated_feedback_", re.IGNORECASE)
    assert binary_linkage_from_text(
        "000000 T thermal_feedback_\n", 0, [thermal_pattern]
    )["status"] == "pass"
    assert binary_linkage_from_text(
        "000000 D thermal_feedback_$format_pack\n", 0, [thermal_pattern]
    )["status"] == "blocked"

    print(
        "STELLAR_SOURCE_PARITY_GATE_OK "
        f"status={payload['status']} "
        f"differing_shared={len(payload['shared_contract']['differing_sources'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
