#!/usr/bin/env python3
"""Audit the AGN coarse-state ledger boundary and source transaction wiring.

This is a read-only, bounded audit.  It checks the source-owned algebra in a
JSONL ledger and the ordering/ownership markers in the production Fortran
patch.  It does not claim physical AGN SED closure, hydro closure, or a
durable cross-restart journal.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SNRT_ROOT = Path(__file__).resolve().parents[1]
if str(SNRT_ROOT) not in sys.path:
    sys.path.insert(0, str(SNRT_ROOT))

from snrt_core.sink_diagnostic import _canonicalize_agn_records


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_section(text: str, name: str, end_name: str) -> str:
    start = text.lower().find(f"subroutine {name.lower()}")
    end = text.lower().find(f"end subroutine {end_name.lower()}", start)
    if start < 0 or end < 0:
        return ""
    return text[start : end + len(f"end subroutine {end_name}")]


def _static_audit(
    source_path: Path,
    driver_path: Path,
    source_module_path: Path,
    makefile_path: Path,
    helper_path: Path,
    source_smoke_path: Path,
) -> dict[str, bool]:
    source_text = source_path.read_text(encoding="utf-8", errors="replace")
    driver_text = driver_path.read_text(encoding="utf-8", errors="replace")
    module_text = source_module_path.read_text(encoding="utf-8", errors="replace")
    makefile_text = makefile_path.read_text(encoding="utf-8", errors="replace")
    helper_text = helper_path.read_text(encoding="utf-8", errors="replace")
    source_smoke_text = source_smoke_path.read_text(encoding="utf-8", errors="replace")
    feedback = _source_section(source_text, "AGN_feedback", "AGN_feedback")
    writer = _source_section(source_text, "dump_agn_coarse_state", "dump_agn_coarse_state")
    driver = _source_section(driver_text, "snrt_ramses_advance_level", "snrt_ramses_advance_level")
    source_loop_match = re.search(
        r"do\s+isink\s*=\s*1\s*,\s*nsink(.*?)end\s+do",
        driver,
        flags=re.IGNORECASE | re.DOTALL,
    )
    source_loop = source_loop_match.group(1) if source_loop_match else ""

    dump_call_positions = [
        match.start()
        for match in re.finditer(r"\bcall\s+dump_agn_coarse_state\b", feedback, re.IGNORECASE)
    ]
    blast_position = re.search(r"\bcall\s+AGN_blast\b", feedback, re.IGNORECASE)
    mass_reset_position = re.search(r"dMsmbh\s*\(\s*isort\s*\)\s*=\s*0d0", feedback, re.IGNORECASE)
    coarse_reset_position = re.search(
        r"dMBH_coarse\s*=\s*0d0\s*;\s*dMEd_coarse\s*=\s*0d0", feedback, re.IGNORECASE
    )
    dump_before_blast = bool(
        dump_call_positions
        and blast_position
        and dump_call_positions[0] < blast_position.start()
    )
    dump_before_mass_reset = bool(
        dump_call_positions
        and mass_reset_position
        and dump_call_positions[0] < mass_reset_position.start()
    )
    dump_before_coarse_reset = bool(
        dump_call_positions
        and coarse_reset_position
        and dump_call_positions[0] < coarse_reset_position.start()
    )

    transaction_position = driver.find("call snrt_agn_deposit_transaction")
    transaction_suffix = driver[transaction_position:] if transaction_position >= 0 else ""
    successful_commit_block = re.search(
        r"if\s*\(\s*source_ok\s*\)\s*then(.*?)end\s*if",
        transaction_suffix,
        flags=re.IGNORECASE | re.DOTALL,
    )
    successful_commit_body = successful_commit_block.group(1) if successful_commit_block else ""
    accounting_commit_match = re.search(
        r"accounted_inflow\s*\(\s*isink\s*\)",
        successful_commit_body,
        flags=re.IGNORECASE,
    )
    accounting_commit_position = (
        transaction_position + successful_commit_block.start(1) + accounting_commit_match.start()
        if accounting_commit_match and successful_commit_block and transaction_position >= 0
        else -1
    )
    budget_position = driver.find("call snrt_agn_photon_budget")
    budget_region = driver[budget_position : budget_position + 500] if budget_position >= 0 else ""
    helper_transform = "epsilon_efficiency" in helper_text or "effective_efficiency" in helper_text
    static = {
        "one_coarse_writer_call_in_agn_feedback": len(dump_call_positions) == 1,
        "coarse_writer_precedes_agn_blast": dump_before_blast,
        "coarse_writer_precedes_dMsmbh_reset": dump_before_mass_reset,
        "coarse_writer_precedes_coarse_rate_reset": dump_before_coarse_reset,
        "rank1_writer_guard": "myid /= 1" in writer,
        "same_step_duplicate_guard_in_memory": "nstep_coarse == nstep_coarse_old" in writer,
        "same_step_duplicate_guard_in_writer": "nstep_coarse == last_dump_step" in writer,
        "pre_reset_marker_emitted": '"ledger_phase":"pre_feedback_pre_reset"' in writer,
        "instantaneous_marker_emitted": '"source_interval_kind":"instantaneous_pre_reset_state"' in writer,
        "raw_and_effective_efficiency_emitted": (
            '"raw_radiative_efficiency":' in writer
            and '"radiative_efficiency":' in writer
            and '"effective_radiative_efficiency":' in writer
        ),
        "efficiency_status_and_mode_emitted": (
            '"efficiency_status":' in writer
            and '"efficiency_status_name":' in writer
            and '"efficiency_mode":' in writer
            and '"efficiency_contract_ok":' in writer
        ),
        "shared_efficiency_helper_called_by_writer": (
            "call snrt_agn_resolve_efficiency" in writer
        ),
        "shared_efficiency_helper_called_by_driver": (
            "call snrt_agn_resolve_efficiency" in driver
        ),
        "helper_has_single_mode_transform": (
            "module snrt_agn_efficiency" in helper_text
            and "pure subroutine snrt_agn_resolve_efficiency" in helper_text
            and helper_transform
            and helper_text.count("effective_efficiency = resolved_base_efficiency *") == 1
        ),
        "writer_has_no_independent_mad_formula": (
            "epsilon_eff=epsilon_r*" not in writer
            and "epsilon_eff = epsilon_r *" not in writer
        ),
        "driver_has_no_independent_mad_formula": (
            "epsilon_eff=epsilon_r*" not in driver
            and "epsilon_eff = epsilon_r *" not in driver
        ),
        "driver_photon_budget_uses_effective_efficiency": (
            budget_position >= 0 and "epsilon_eff" in budget_region
        ),
        "driver_photon_budget_uses_supplied_inflow": (
            "call snrt_agn_photon_budget(delta_inflow" in driver
            and "accounted_inflow" in driver
        ),
        "driver_supplied_mass_is_bondi_eddington_minimum": (
            "supplied_mass=min(max(dMBH_coarse(isink),0.0d0)," in driver
            and "max(dMEd_coarse(isink),0.0d0)" in driver
        ),
        "driver_retained_mass_is_one_sided_check_only": (
            "retained_bound=(1.0d0-epsilon_eff)*delta_inflow" in driver
            and "dMsmbh(isink)" in driver
        ),
        "driver_rearms_retained_cursor_after_successful_commit": bool(
            successful_commit_block
            and re.search(
                r"retained_seen\s*\(\s*isink\s*\)\s*=\s*dMsmbh\s*\(\s*isink\s*\)",
                successful_commit_body,
                flags=re.IGNORECASE,
            )
        ),
        "driver_has_no_hidden_efficiency_clamp": (
            # This is a lexical negative guard, not a complete semantic proof.
            # The budget-region check also proves raw value is not used as the
            # photon coefficient in the audited call site.
            "0.99d0" not in driver
            and "1.0d-6" not in driver
            and "eps_sink" not in budget_region
            and "raw_epsilon" not in budget_region
            and "epsilon_eff" in budget_region
        ),
        "source_api_declares_supplied_inflow": (
            "delta_inflow_mass_code" in module_text
            and "retained BH mass" in module_text
        ),
        "source_smoke_records_supplied_mass_api": (
            "call snrt_agn_photon_budget" in source_smoke_text
            and "supplied inflow" in source_smoke_text
        ),
        "julian_year_marker_emitted": '"julian_year_days":365.25' in writer,
        "atomic_transaction_call_present": transaction_position >= 0,
        "accounting_commit_after_transaction": (
            transaction_position >= 0
            and accounting_commit_position > transaction_position
        ),
        "stable_idsink_key_map_present": (
            "idsink" in driver and "accounted_ids" in driver
        ),
        "source_loop_is_not_openmp_directive": (
            bool(source_loop_match) and "!$omp" not in source_loop.lower()
        ),
        "local_leaf_single_owner_contract_present": (
            "snrt_agn_find_local_leaf" in driver
            and "MPI owner" in driver
            and "local" in driver.lower()
        ),
        "rt_enable_latched_once": (
            "enabled_latched" in driver
            and "Runtime control is latched once per process" in driver
            and driver.count("get_environment_variable('SNRT_RT_ENABLE'") == 1
        ),
        "transaction_module_present": (
            "subroutine snrt_agn_deposit_transaction" in module_text
            and "Commit is deliberately a separate phase" in module_text
        ),
        "transaction_validates_finite_inputs": (
            "ieee_is_finite" in module_text
            and "any(.not. ieee_is_finite(emitted_photons))" in module_text
        ),
        "production_makefile_links_transaction_module": (
            "snrt_agn_source.o" in makefile_text
            and "snrt_ramses_driver.o" in makefile_text
            and "snrt_agn_source.o" in makefile_text[makefile_text.find("snrt_ramses_driver.o:") :]
        ),
        "production_makefile_links_shared_efficiency_helper": (
            "snrt_agn_efficiency.o" in makefile_text
            and "snrt_agn_efficiency.o" in makefile_text[makefile_text.find("MODOBJ") :]
            and "sink_particle.kjhan.o: snrt_agn_efficiency.o" in makefile_text
            and "snrt_ramses_driver.o: amr_commons.kjhan.o pm_commons.o snrt_agn_efficiency.o" in makefile_text
        ),
        "helper_has_amr_parameter_dependency": (
            "snrt_agn_efficiency.o: amr_parameters.jaehyun.o" in makefile_text
        ),
        "driver_has_direct_amr_pm_dependencies": (
            "snrt_ramses_driver.o: amr_commons.kjhan.o pm_commons.o" in makefile_text
        ),
        "helper_path_is_explicit_audit_input": helper_path.is_file(),
        "production_makefile_links_transport_graph": all(
            object_name in makefile_text
            for object_name in (
                "mpi_mod.o",
                "snrt_state.o",
                "snrt_amr_topology.o",
                "snrt_transport_step.o",
                "snrt_cuda_sparse_transport_interface.o",
                "snrt_cuda_multigroup_interface.o",
            )
        ),
        "production_makefile_keeps_runtime_cuda_gate": (
            "SNRT=1 requires USE_CUDA=1" in makefile_text
        ),
    }
    return static


def _audit(
    input_path: Path,
    source_path: Path,
    driver_path: Path,
    source_module_path: Path,
    makefile_path: Path,
    helper_path: Path,
    source_smoke_path: Path,
) -> dict[str, Any]:
    blockers: list[str] = []
    try:
        records, duplicate_count = _canonicalize_agn_records(input_path)
        input_audit: dict[str, Any] = {
            "passed": True,
            "record_count": len(records),
            "duplicate_count_collapsed": duplicate_count,
            "coarse_steps": sorted({int(record["nstep_coarse"]) for record in records}),
            "sink_ids": sorted({int(record["sink_id"]) for record in records}),
            "all_records_pre_reset_instantaneous": True,
            "source_algebra_validated": True,
        }
    except (OSError, ValueError, KeyError, TypeError) as error:
        input_audit = {"passed": False, "error": str(error)}
        blockers.append(f"AGN coarse-state input audit failed: {error}")

    try:
        static = _static_audit(
            source_path,
            driver_path,
            source_module_path,
            makefile_path,
            helper_path,
            source_smoke_path,
        )
    except OSError as error:
        static = {}
        blockers.append(f"production source audit could not read a required file: {error}")
    for criterion, passed in static.items():
        if not passed:
            blockers.append(f"static criterion failed: {criterion}")

    passed = not blockers
    return {
        "schema": "snrt_agn_coarse_ledger_audit_v2",
        "passed": passed,
        "bundle": "F-P1.5 AGN ledger transaction",
        "criteria": static,
        "input": input_audit,
        "blockers": blockers,
        "physical_closure_claim": False,
        "limitations": [
            "No physical AGN SED, obscuration, escape-fraction, or hydro-closure claim is made.",
            "A run UUID/dump counter is not present; conflicting rewind payloads fail closed.",
            "Cross-coarse-step deferred re-emission and a durable crash journal remain open follow-ups.",
            "The shared helper convention is now checked here; this remains an arithmetic/engineering audit, not a physical AGN SED validation.",
            "Input evidence is arithmetic/transactional and must be refreshed from the same production run.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="AGN coarse-state JSONL ledger")
    parser.add_argument(
        "--source",
        type=Path,
        default=REPO_ROOT / "patch/lagRamses/sink_particle.kjhan.f90",
    )
    parser.add_argument(
        "--driver",
        type=Path,
        default=REPO_ROOT / "patch/lagRamses/snrt_ramses_driver.f90",
    )
    parser.add_argument(
        "--source-module",
        type=Path,
        default=REPO_ROOT / "patch/lagRamses/snrt_agn_source.f90",
    )
    parser.add_argument(
        "--makefile",
        type=Path,
        default=REPO_ROOT / "bin/Makefile",
    )
    parser.add_argument(
        "--helper",
        type=Path,
        required=True,
        help="shared AGN efficiency helper source; its SHA256 is recorded",
    )
    parser.add_argument(
        "--source-smoke",
        type=Path,
        default=REPO_ROOT / "patch/lagRamses/snrt_agn_source_smoke.f90",
        help="positional photon-budget API smoke source",
    )
    parser.add_argument("--output", type=Path, required=True, help="JSON audit report")
    args = parser.parse_args()

    report = _audit(
        args.input,
        args.source,
        args.driver,
        args.source_module,
        args.makefile,
        args.helper,
        args.source_smoke,
    )
    report["provenance"] = {
        "input": str(args.input),
        "source": str(args.source),
        "driver": str(args.driver),
        "source_module": str(args.source_module),
        "makefile": str(args.makefile),
        "helper": str(args.helper),
        "source_smoke": str(args.source_smoke),
        "sha256": {
            "input": _sha256(args.input) if args.input.exists() else None,
            "source": _sha256(args.source) if args.source.exists() else None,
            "driver": _sha256(args.driver) if args.driver.exists() else None,
            "source_module": _sha256(args.source_module) if args.source_module.exists() else None,
            "makefile": _sha256(args.makefile) if args.makefile.exists() else None,
            "helper": _sha256(args.helper) if args.helper.exists() else None,
            "source_smoke": _sha256(args.source_smoke) if args.source_smoke.exists() else None,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"AGN_COARSE_LEDGER_AUDIT_{'PASS' if report['passed'] else 'FAIL'} "
        f"records={report['input'].get('record_count', 0)} output={args.output}"
    )
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
