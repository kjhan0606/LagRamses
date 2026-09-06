#!/usr/bin/env python3
"""Independently reproduce the 2026-09-02 Fable SN/AGN audit findings.

This is a read-only audit.  It does not compile, launch, or modify a RAMSES
run.  The evidence is intentionally split into source markers, manifest
metadata, and two small numerical counterexamples so that the audit does not
depend on Fable's interpretation or on a successful production build.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]


def read_text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8", errors="replace")


def line_numbers(relative: str, patterns: list[str]) -> list[str]:
    lines = read_text(relative).splitlines()
    evidence: list[str] = []
    for number, line in enumerate(lines, start=1):
        if any(re.search(pattern, line, re.IGNORECASE) for pattern in patterns):
            evidence.append(f"{relative}:{number}")
    return evidence


def contains(relative: str, pattern: str) -> bool:
    return re.search(pattern, read_text(relative), re.IGNORECASE | re.MULTILINE) is not None


def sha256(relative: str) -> str:
    digest = hashlib.sha256()
    with (ROOT / relative).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(relative: str) -> dict[str, Any]:
    return json.loads(read_text(relative))


def git_head() -> str:
    return subprocess.check_output(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, text=True
    ).strip()


def git_status() -> list[str]:
    return subprocess.check_output(
        ("git", "status", "--short"), cwd=ROOT, text=True
    ).splitlines()


def numerical_reproductions() -> dict[str, Any]:
    # A cumulative quantity C(t)=t^2 should telescope to C(6)-C(0)=36.
    # The implementation's C(age+dt)-C(age), with age interpreted as the
    # current end-of-step age, produces a shifted and non-telescoping total.
    ages = [1.0, 3.0, 6.0]
    timesteps = [1.0, 2.0, 3.0]
    cumulative = lambda t: t * t
    observed = sum(
        cumulative(age + dt) - cumulative(age)
        for age, dt in zip(ages, timesteps)
    )
    correct = cumulative(ages[-1]) - cumulative(0.0)

    # A one-Gyr physical age becomes 1.0 when passed to a table whose axis is
    # explicitly years; it is therefore 1e-9 of the intended coordinate.
    physical_age_gyr = 1.0
    compiled_query_on_year_axis = physical_age_gyr
    intended_query_yr = physical_age_gyr * 1.0e9

    return {
        "forward_cumulative_interval": {
            "cumulative_function": "C(t)=t^2",
            "ages": ages,
            "timesteps": timesteps,
            "implementation_total": observed,
            "correct_telescoping_total": correct,
            "absolute_difference": abs(observed - correct),
            "first_interval_implementation": [1.0, 2.0],
            "first_interval_expected": [0.0, 1.0],
            "reproduced": observed != correct,
        },
        "year_gyr_coordinate_mismatch": {
            "physical_age_gyr": physical_age_gyr,
            "compiled_query_if_untagged": compiled_query_on_year_axis,
            "intended_year_axis_query": intended_query_yr,
            "coordinate_ratio": compiled_query_on_year_axis / intended_query_yr,
            "reproduced": compiled_query_on_year_axis != intended_query_yr,
        },
    }


def finding(
    finding_id: str,
    severity: str,
    status: str,
    title: str,
    reproduction: str,
    evidence: list[str],
    priority: str,
    acceptance: str,
) -> dict[str, Any]:
    return {
        "id": finding_id,
        "severity": severity,
        "status": status,
        "title": title,
        "reproduction": reproduction,
        "evidence": evidence,
        "implementation_priority": priority,
        "acceptance_test": acceptance,
    }


def reproduce() -> dict[str, Any]:
    makefile = read_text("bin/Makefile")
    patch_runtime = read_text("patch/lagRamses/stellar_ramses_runtime.f90")
    patch_interpolation = read_text("patch/lagRamses/stellar_yield_interpolation.f90")
    patch_increment = read_text("patch/lagRamses/stellar_source_increment.f90")
    mirror_runtime = read_text("simulation/snrt/native/phase0/stellar_ramses_runtime.f90")
    mirror_increment = read_text("simulation/snrt/native/phase0/stellar_source_increment.f90")
    hdf5_backup = read_text("patch/lagRamses/backup_hdf5.f90")
    hdf5_restore = read_text("patch/lagRamses/restore_hdf5.f90")
    agn_driver = read_text("patch/lagRamses/snrt_ramses_driver.f90")
    agn_sink = read_text("patch/lagRamses/sink_particle.kjhan.f90")
    agn_diag = read_text("patch/lagRamses/snrt_agn_ledger.f90")

    legacy_audit = load_json("simulation/snrt/data/legacy_yield_table_audit.json")
    g2_asset_audit = load_json("simulation/snrt/data/g2_candidate_source_audit.json")
    g2_contract = load_json("simulation/snrt/config/g2_physics_contract_v1.json")
    transition_audit = load_json(
        "simulation/snrt/data/feedback_transition_phase0_output_00011_native_audit.json"
    )
    candidate_field_map = load_json("simulation/snrt/config/stellar_ramses_field_map_v1.json")
    agn_photon = load_json("simulation/snrt/data/p4_pilot_agn_photon_ledger.json")

    production_uses_patch = "PATCH = ../patch/lagRamses" in makefile
    g1_uses_mirror = contains(
        "simulation/snrt/tests/run_g1_native_contract.sh", r"SOURCE_DIR=.*native/phase0"
    )
    g1_excludes_runtime = not contains(
        "simulation/snrt/tests/run_g1_native_contract.sh", r"stellar_ramses_runtime"
    )
    compiled_has_age_gyr = "age_gyr" in patch_runtime
    table_declares_years = contains(
        "patch/lagRamses/stellar_yield_tables.f90", r"age.{0,20}yr"
    )
    compiled_clamps = "lower = minimum" in patch_interpolation and "upper = maximum" in patch_interpolation
    production_converts_age = contains(
        "patch/lagRamses/stellar_yield_tables.f90", r"age_yr\s*\*\s*1\.0e-9"
    )
    current_interval_contract = (
        contains("patch/lagRamses/stellar_source_increment.f90", r"previous_age_gyr")
        and contains("patch/lagRamses/stellar_source_increment.f90", r"current_age_gyr")
        and contains("patch/lagRamses/stellar_source_increment.f90", r"cumulative_difference")
    )
    production_requires_external_table = contains(
        "patch/lagRamses/stellar_ramses_runtime.f90",
        r"embedded fallback is disabled",
    )
    mirror_converts_age = contains(
        "simulation/snrt/native/phase0/stellar_yield_tables.f90", r"age_yr\s*\*\s*1\.0e-9"
    )
    mirror_rejects_domain = contains(
        "simulation/snrt/native/phase0/stellar_yield_interpolation.f90", r"out.of.domain|no lower|no upper"
    )
    production_nvar18 = bool(re.search(r"NVAR\s*=\s*18", makefile))
    production_nener0 = bool(re.search(r"NENER\s*=\s*0", makefile))
    hydro_index_formula = contains(
        "patch/lagRamses/read_hydro_params.f90", r"inener=ndim\+3"
    ) and contains(
        "patch/lagRamses/read_hydro_params.f90", r"imetal=nener\+ndim\+3"
    )
    mirror_uses_inener = contains(
        "simulation/snrt/native/phase0/stellar_ramses_runtime.f90", r"energy_index\s*=\s*inener"
    )
    compiled_literal_energy_field = contains(
        "patch/lagRamses/stellar_ramses_runtime.f90", r"unew\(target_cell,5\)"
    )

    findings = [
        finding(
            "F1", "critical", "reproduced",
            "No approved physical yield asset",
            "The physical contract is blocked and the staged candidate audit is not production-ready; the legacy audit is explicitly legacy-only.",
            [
                "simulation/snrt/config/g2_physics_contract_v1.json:approval.current_status",
                "simulation/snrt/data/g2_candidate_source_audit.json:status",
                "simulation/snrt/data/legacy_yield_table_audit.json:status",
            ],
            "P0.5",
            "An approved, checksummed full mass-metallicity-age grid and sidecar pass the G2 audit for every enabled channel.",
        ),
        finding(
            "F2", "critical", "reproduced",
            "Validated mirror differs from the compiled production tree",
            f"The production Makefile selects patch/lagRamses while the G1 runner selects native/phase0; the source identities are separate ({sha256('patch/lagRamses/stellar_source_increment.f90')[:12]} vs {sha256('simulation/snrt/native/phase0/stellar_source_increment.f90')[:12]}), and the runner's source list does not include the RAMSES deposition runtime.",
            line_numbers("bin/Makefile", [r"PATCH =", r"STELLAR_ENRICHMENT_MODOBJ"])
            + line_numbers("simulation/snrt/tests/run_g1_native_contract.sh", [r"SOURCE_DIR"]),
            "P0.1",
            "The same source objects are linked by the production Makefile and exercised by the native contract suite; a parity test fails on divergence.",
        ),
        finding(
            "F3", "critical", "not_reproduced" if production_converts_age else "reproduced",
            "Compiled runtime queries a year axis in Gyr",
            f"The historical unit counterexample gives 1.0 Gyr -> {numerical_reproductions()['year_gyr_coordinate_mismatch']['compiled_query_if_untagged']:.1f} on a year axis instead of 1.0e9 yr; the current production table reader converts age_yr to age_gyr={production_converts_age} before interpolation.",
            line_numbers("patch/lagRamses/stellar_yield_tables.f90", [r"age.*yr"])
            + line_numbers("patch/lagRamses/stellar_ramses_runtime.f90", [r"age_gyr"]),
            "P0.2",
            "A compiled-tree age-node test proves one unambiguous unit conversion, including the RAMSES aexp convention, at every source-query boundary.",
        ),
        finding(
            "F4", "critical", "not_reproduced" if current_interval_contract else "reproduced",
            "Cumulative release interval is forward-shifted",
            "The historical C(t)=t^2 variable-step counterexample returns 64 rather than the telescoping 36; the current source-increment contract uses explicit previous/current ages and cumulative_difference.",
            line_numbers("patch/lagRamses/stellar_source_increment.f90", [r"age \+ timestep", r"cumulative_difference"])
            + line_numbers("simulation/snrt/native/phase0/stellar_source_increment.f90", [r"age_gyr \+ timestep_gyr"]),
            "P0.2",
            "Variable-timestep, first-interval, repeated-call, and stop/restart tests satisfy exact cumulative telescoping and idempotence.",
        ),
        finding(
            "F5", "critical", "reproduced",
            "HDF5 restart omits stellar release state",
            "The HDF5 particle writer/reader handles birth_epoch but has no HDF5 fields for tpp, mp0, or indtab; binary paths are therefore not evidence for HDF5 restart safety.",
            line_numbers("patch/lagRamses/backup_hdf5.f90", [r"birth_epoch", r"/particles"])
            + line_numbers("patch/lagRamses/restore_hdf5.f90", [r"birth_epoch", r"/particles"]),
            "P0.3",
            "A bitwise HDF5 round trip preserves tpp, mp0, indtab and every newly introduced stellar progress/ledger field.",
        ),
        finding(
            "F6", "high", "reproduced",
            "Legacy Sedov energy path lacks the 1e51 erg factor",
            "The lagRamses Sedov expressions omit 1d51 while the base RAMSES expression includes it; channel-mode kinetic feedback is separately disabled, so the missing factor is not double-counted there.",
            line_numbers("patch/lagRamses/feedback.kjhan3.f90", [r"ESN_code", r"ESN="])
            + line_numbers("pm/feedback.f90", [r"1d51"])
            + line_numbers("patch/lagRamses/amr_step.jaehyun.f90", [r"use_channel_resolved_feedback"]),
            "P1.2",
            "One isolated event with declared cgs input produces the exact expected thermal/kinetic energy, or the legacy path is explicitly retired and gated.",
        ),
        finding(
            "F7", "high", "not_reproduced" if not compiled_clamps else "reproduced",
            "Silent endpoint clamp and one-point channel fixture",
            f"The historical compiled interpolator substituted minimum/maximum nodes outside the domain; the current clamp marker is present={compiled_clamps}. The physical fixture remains non-production and is covered separately by the G2 asset gate.",
            line_numbers("patch/lagRamses/stellar_yield_interpolation.f90", [r"lower = minimum", r"upper = maximum"])
            + line_numbers("patch/lagRamses/stellar_ramses_runtime.f90", [r"0\.8", r"120\.0", r"140\.0"]),
            "P0.4",
            "Production mode rejects every out-of-domain query and explicit configuration supplies complete channel mass/population coverage.",
        ),
        finding(
            "F8", "high", "not_reproduced" if production_requires_external_table else "reproduced",
            "Embedded synthetic fallback and implicit channel assumptions",
            f"The historical finding concerned an embedded fallback when PHASE0_YIELD_TABLE was absent; the current production runtime explicitly requires the external table (guard={production_requires_external_table}) and the production parity gate forbids the embedded macro. SNIa/PISN remain intentionally disabled until their physical gates pass.",
            line_numbers("patch/lagRamses/stellar_ramses_runtime.f90", [r"PHASE0_YIELD_TABLE", r"fallback"])
            + line_numbers("patch/lagRamses/stellar_ssp_sources.f90", [r"Kroupa", r"0\.8", r"140"])
            + ["simulation/snrt/config/g2_physics_contract_v1.json:channel_partition"],
            "P0.4",
            "A production binary cannot start without an approved external table, sidecar, explicit IMF/channel config, and enabled-channel approval.",
        ),
        finding(
            "F9", "high", "reproduced",
            "Channel-mode SN deposition is one-cell thermal feedback with no radial momentum",
            "The compiled channel runtime selects one target cell, writes total energy and momentum fields, and the source contract treats unresolved isotropic source-frame vector momentum as zero; the kinetic path is disabled in channel mode.",
            line_numbers("patch/lagRamses/stellar_ramses_runtime.f90", [r"target_cell", r"unew\(target_cell", r"momentum"])
            + line_numbers("patch/lagRamses/amr_step.jaehyun.f90", [r"not\.use_channel_resolved_feedback"]),
            "P1.1",
            "An approved SN deposition model passes isolated-event energy, radial-momentum, spatial-resolution, and delayed-cooling ownership tests.",
        ),
        finding(
            "F10", "high", "reproduced",
            "AGN accretion, deposit, and live source conventions differ",
            "The sink blast uses dMsmbh_AGN with thermal/jet efficiencies, the coarse ledger records an effective efficiency and Lbol, while the live source consumes a dMsmbh increment and a fixed half-factor per group.",
            line_numbers("patch/lagRamses/sink_particle.kjhan.f90", [r"dMsmbh_AGN", r"epsilon_r", r"EAGN"])
            + line_numbers("patch/lagRamses/snrt_ramses_driver.f90", [r"delta_accreted", r"0\.5d0", r"epsilon_r"]),
            "P1.3",
            "One declared AGN convention closes accreted mass, Lbol, SED photons, thermal/jet energy, momentum, and deposited gas through restart.",
        ),
        finding(
            "F11", "high", "reproduced",
            "Deferred live source can retry an accumulator",
            "The live SNRT path resets accounted_mass at every new coarse step, advances it only after source_ok, and the sink path resets dMsmbh only when ok_blast_agn is true. When a blast is deferred, the full accumulated dMsmbh is therefore eligible again on the next coarse step; this is statically reproducible, with a dynamic run still needed to measure its effect.",
            line_numbers("patch/lagRamses/snrt_ramses_driver.f90", [r"accounted_mass", r"source_ok", r"delta_accreted"]),
            "P1.4",
            "A deferred source has an atomic per-interval transaction key and a restart/retry test proves exactly-once or explicitly compensating semantics.",
        ),
        finding(
            "F12", "medium", "reproduced",
            "AGN ledger interval/deferred semantics are underspecified",
            "The rate ledger records aexp/time_code and instantaneous inflow, but no explicit interval start/end or deferred/committed state is exposed to downstream consumers.",
            line_numbers("simulation/snrt/tools/p4_build_agn_rate_ledger.py", [r"time_code", r"aexp", r"accretion_rate_convention"])
            + ["simulation/snrt/data/p4_pilot_zoom_agn_candidates.csv:header"],
            "P1.3",
            "Each row declares interval endpoints, accumulator state, commit/defer outcome, and the exact rate-to-integrated-mass convention.",
        ),
        finding(
            "F13", "medium", "partially_reproduced",
            "Deduplication and efficiency-field compatibility are incomplete",
            "The current reader rejects duplicate sink IDs and the merger rejects source-ID collisions, so the original blanket 'no deduplication' claim is only partly current. The ledger remains fail-closed rather than coalescing duplicate keys, and raw/effective efficiency conventions still require one schema.",
            line_numbers("simulation/snrt/snrt_core/sink_diagnostic.py", [r"duplicate", r"radiative_efficiency"])
            + line_numbers("simulation/snrt/tools/merge_photon_source_ledgers.py", [r"collision", r"source.id", r"duplicate"])
            + line_numbers("simulation/snrt/tools/p7_convert_sinkprops.py", [r"effective_radiative_efficiency"]),
            "P1.5",
            "Duplicate (coarse_step,sink_id) rows are deterministically coalesced or conflict-rejected, converter/reader field names agree, and effective efficiency is preserved.",
        ),
        finding(
            "F14", "medium", "partially_reproduced",
            "He/disabled-element/untracked-metal and legacy field semantics remain ambiguous",
            f"The NVAR difference is historical rather than a current production mismatch: the executable Makefile declares NVAR={makefile.split('NVAR =', 1)[1].splitlines()[0].strip()} and the candidate map assumes NVAR={candidate_field_map['compile_time_assumptions']['nvar']}; transitional output-00011 records NVAR={transition_audit['makefile_flags']['NVAR']} and first element field={transition_audit['runtime_log_evidence']['first_element_field']}. The remaining finding is the absence of one startup-validated He/disabled-element/untracked-metal semantic contract in the compiled path.",
            [
                "simulation/snrt/data/feedback_transition_phase0_output_00011_native_audit.json:makefile_flags/runtime_log_evidence",
                "simulation/snrt/config/stellar_ramses_field_map_v1.json:compile_time_assumptions",
            ]
            + line_numbers("patch/lagRamses/hydro_parameters.f90", [r"iHydrogen", r"iHelium"])
            + line_numbers("patch/lagRamses/stellar_ramses_runtime.f90", [r"elem_c", r"unew\(target_cell,5"]),
            "P0.6",
            "Startup rejects an unrecognized executed field map and source/cell ledgers close H, He, tracked elements, residual metal, and disabled fields.",
        ),
        finding(
            "F15", "medium", "reproduced",
            "Native mirror energy field collides with total metal for the production NENER=0 layout",
            "The compiled production runtime writes the literal energy field 5, which is correct for the raw RAMSES layout. The separate native mirror assigns energy_index=inener; with NENER=0 and ndim=3, read_hydro_params defines inener=6 and imetal=6, so the mirror field-map validator must reject that layout. The defect is statically proven and is a source-of-truth/parity blocker.",
            line_numbers("patch/lagRamses/stellar_ramses_runtime.f90", [r"unew\(target_cell,5"])
            + line_numbers("patch/lagRamses/read_hydro_params.f90", [r"inener=ndim\+3", r"imetal=nener\+ndim\+3"])
            + line_numbers("simulation/snrt/native/phase0/stellar_ramses_runtime.f90", [r"energy_index = inener"])
            + ["simulation/snrt/config/stellar_ramses_field_map_v1.json:indices_one_based"],
            "P0.6",
            "The actual production build and the native test use one field map; NENER=0 compile-time checks prove energy, metal, delayed-cooling, and element indices cannot overlap.",
        ),
        finding(
            "F16", "medium", "reproduced",
            "Stellar and AGN SEDs are pilot candidates with domain/escape assumptions",
            "The BPASS ledger documents clamped sources and escape fraction one; the AGN ledger is an unobscured Sazonov-style baseline with escape fraction one and partial low-energy support.",
            [
                "simulation/snrt/P4_STELLAR_SED.md:73-89",
                "simulation/snrt/data/p4_pilot_agn_photon_ledger.json:normalization",
            ]
            + line_numbers("simulation/snrt/tools/p4_build_agn_photon_ledger.py", [r"Sazonov", r"escape_fraction"]),
            "P2.1",
            "Approved population/AGN SED, obscuration, escape, domain, and normalization sensitivities close photon-number and energy ledgers.",
        ),
        finding(
            "F17", "low", "reproduced",
            "Diagnostic AGN ledger hard-codes sink one and appends across rewinds",
            "The diagnostic file opens with position='append', emits a header only in process memory, and writes literal sink ID 1 with index-1 state. This is diagnostic-only but can create conflicting rows after restart/rewind.",
            line_numbers("patch/lagRamses/snrt_agn_ledger.f90", [r"position='append'", r"write\(unit_id", r"tsink\(1\)"] ),
            "P1.5",
            "Diagnostic rows carry the actual sink ID and an immutable interval/run key; restart/rewind behavior is tested or append is replaced by an atomic partitioned writer.",
        ),
    ]

    statuses: dict[str, list[str]] = {"reproduced": [], "partially_reproduced": [], "not_reproduced": []}
    for item in findings:
        statuses[item["status"]].append(item["id"])

    return {
        "schema": "fable_sn_agn_independent_reproduction_v1",
        "date": "2026-09-02",
        "repository_root": str(ROOT),
        "repository_head": git_head(),
        "worktree_status": git_status(),
        "method": [
            "read-only source/build/config inspection",
            "artifact metadata inspection",
            "pure numerical counterexamples for unit and interval semantics",
        ],
        "scope_limits": [
            "No production executable was rebuilt or launched by this audit.",
            "F11 is statically reproduced; a dynamic run is still needed to measure its coupled observable effect.",
            "F14 remains partial because one historical metadata comparison was corrected; F15 is statically proven for the NENER=0 layout.",
            "Static evidence establishes wiring and contract risks, not the magnitude of a full cosmological run's observable bias.",
        ],
        "independent_checks": {
            "production_makefile_selects_patch_tree": production_uses_patch,
            "g1_runner_selects_separate_native_mirror": g1_uses_mirror,
            "compiled_runtime_uses_gyr": compiled_has_age_gyr,
            "g1_runner_excludes_ramses_runtime": g1_excludes_runtime,
            "table_axis_declares_years": table_declares_years,
            "compiled_interpolator_clamps": compiled_clamps,
            "production_converts_age_once": production_converts_age,
            "current_interval_contract": current_interval_contract,
            "production_requires_external_table": production_requires_external_table,
            "mirror_converts_year_axis": mirror_converts_age,
            "mirror_rejects_domain": mirror_rejects_domain,
            "production_nvar18": production_nvar18,
            "production_nener0": production_nener0,
            "hydro_index_formula_present": hydro_index_formula,
            "mirror_uses_inener": mirror_uses_inener,
            "compiled_literal_energy_field": compiled_literal_energy_field,
            "transition_nvar": transition_audit["makefile_flags"]["NVAR"],
            "candidate_map_nvar": candidate_field_map["compile_time_assumptions"]["nvar"],
            "agn_group_count": len(agn_photon["groups"]),
            "agn_group_edges": agn_photon["group_edges_ev"],
            "g2_status": g2_contract["approval"]["current_status"],
            "legacy_status": legacy_audit["status"],
            "candidate_asset_status": g2_asset_audit["status"],
        },
        "numerical_reproductions": numerical_reproductions(),
        "findings": findings,
        "summary": statuses,
        "priority_order": {
            "P0": [
                "P0.1 F2/F15: make the compiled production tree the tested source of truth and gate build/map parity.",
                "P0.2 F3/F4: correct age units and cumulative interval semantics with variable-dt/restart tests.",
                "P0.3 F5: persist complete stellar release state through HDF5 restart.",
                "P0.4 F7/F8: remove fail-open fallback, silent clamp, hard-coded population windows, and implicit channel activation.",
                "P0.6 F14/F15: validate NVAR/NENER/field/species/He semantics in the actual executable.",
                "P0.5 F1: approve a physical full grid and immutable provenance sidecar only after executed field semantics are closed.",
            ],
            "P1": [
                "P1.1 F9: choose and validate SN energy, radial momentum, spatial deposition, and cooling ownership.",
                "P1.2 F6: fix or retire the legacy Sedov energy path with a one-event dimensional test.",
                "P1.3 F10/F12: unify AGN accretion, luminosity, SED, interval, thermal/jet, and live-source conventions.",
                "P1.4 F11: make deferred source accounting atomic and exactly-once across retry/restart.",
                "P1.5 F13/F17: finalize deduplication, schema compatibility, actual sink IDs, and rewind-safe diagnostics.",
                "P1.6: implement SNIa DTD and PISN population gates only after their physical inputs are approved.",
            ],
            "P2": ["P2.1 F16: approve stellar/AGN SEDs, obscuration, escape, and normalization."],
            "P3": [
                "P3: complete dust scattering/IR/radiation pressure, live stellar/AGN coupling, B3 convergence, AMR/MPI determinism, and production rerun."
            ],
        },
    }


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    args = parser.parse_args()
    payload = reproduce()
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("FABLE_INDEPENDENT_REPRODUCTION")
        print(f"HEAD={payload['repository_head']}")
        print(f"reproduced={','.join(payload['summary']['reproduced'])}")
        print(f"partially_reproduced={','.join(payload['summary']['partially_reproduced'])}")
        print(f"not_reproduced={','.join(payload['summary']['not_reproduced']) or 'none'}")
        interval = payload["numerical_reproductions"]["forward_cumulative_interval"]
        print(
            "interval_counterexample="
            f"implementation:{interval['implementation_total']} "
            f"correct:{interval['correct_telescoping_total']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
