#!/usr/bin/env python3
"""Audit mass, metallicity, rotation, and age coverage of staged G2 sources."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

from audit_g2_huscher2025_candidate import (
    HuscherAuditError,
    audit_huscher2025_candidate,
)
from audit_g2_boccioli_roberti2026_candidate import (
    BoccioliRobertiAuditError,
    audit_boccioli_roberti2026_candidate,
)
from audit_g2_doherty2014_sagb_candidate import (
    DohertyAuditError,
    audit_doherty2014_candidate,
)
from audit_g2_stockinger2020_candidate import (
    StockingerAuditError,
    audit_stockinger2020_candidate,
)
from audit_g2_sukhbold2016_candidate import (
    SukhboldAuditError,
    audit_sukhbold2016_candidate,
)
from audit_g2_limongi2024_transition_fates import (
    LimongiTransitionAuditError,
    audit_limongi2024_transition_fates,
)
from audit_g2_roberti2024_ultralowz_candidate import (
    RobertiUltraLowZAuditError,
    audit_roberti2024_ultralowz_candidate,
)
from audit_g2_heger_woosley2010_popiii_candidate import (
    HegerWoosleyPopIIIAuditError,
    audit_heger_woosley2010_popiii_candidate,
)
from audit_g2_baseline_metallicity_demand import (
    BaselineMetallicityDemandError,
    audit_g2_baseline_metallicity_demand,
)

from adapt_g2_candidate_sources import (
    DEFAULT_ROOT,
    LIMONGI_ID,
    NUGRID_ID,
    SourceAdapterError,
    adapt_candidate,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_candidate_grid_coverage_contract_v1.json"


class CoverageAuditError(ValueError):
    """Candidate grid coverage cannot be audited without an explicit policy."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise CoverageAuditError(f"cannot read {path}: {exc}") from exc
    return digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CoverageAuditError(f"cannot read coverage contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-candidate-grid-coverage-contract"
        or contract.get("schema_version") != 1
    ):
        raise CoverageAuditError("unsupported G2 candidate-grid coverage contract")
    rules = contract.get("coverage_rules", {})
    if rules.get("cross_source_interpolation_allowed") is not False:
        raise CoverageAuditError("review contract must forbid cross-source interpolation")
    if rules.get("out_of_domain_extrapolation_allowed") is not False:
        raise CoverageAuditError("review contract must forbid extrapolation")
    if contract.get("approval", {}).get("canonical_conversion_allowed") is not False:
        raise CoverageAuditError("coverage review unexpectedly permits conversion")
    if rules.get("flattened_branch_union_is_not_interpolable_coverage") is not True:
        raise CoverageAuditError("coverage contract must reject flattened branch interpolation")
    if rules.get("terminal_candidate_nodes_require_source_node_fate_and_remnant_records") is not True:
        raise CoverageAuditError("coverage contract must require terminal source-node records")
    return contract


def _uncovered_edges(required: list[float], nodes: list[float]) -> list[list[float]]:
    if len(required) != 2 or required[1] <= required[0]:
        raise CoverageAuditError("invalid required mass range")
    if not nodes:
        return [[required[0], required[1]]]
    ordered = sorted(set(float(value) for value in nodes))
    if any(not math.isfinite(value) for value in ordered):
        raise CoverageAuditError("non-finite mass node")
    gaps: list[list[float]] = []
    if ordered[0] > required[0]:
        gaps.append([required[0], min(ordered[0], required[1])])
    if ordered[-1] < required[1]:
        gaps.append([max(ordered[-1], required[0]), required[1]])
    return [gap for gap in gaps if gap[1] > gap[0]]


def _uncovered_source_hulls(
    required: list[float], source_hulls: dict[str, list[float]]
) -> list[list[float]]:
    """Return edge and internal gaps after merging only source-owned hulls."""
    lower, upper = required
    intervals = sorted(
        (max(lower, values[0]), min(upper, values[1]))
        for values in source_hulls.values()
        if values[1] > lower and values[0] < upper
    )
    if not intervals:
        return [[lower, upper]]
    merged: list[list[float]] = []
    for left, right in intervals:
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        else:
            merged[-1][1] = max(merged[-1][1], right)
    gaps: list[list[float]] = []
    cursor = lower
    for left, right in merged:
        if left > cursor:
            gaps.append([cursor, left])
        cursor = max(cursor, right)
    if cursor < upper:
        gaps.append([cursor, upper])
    return gaps


def _nearest_log_offsets(source: list[float], target: list[float]) -> list[dict[str, float]]:
    result: list[dict[str, float]] = []
    for value in source:
        nearest = min(target, key=lambda candidate: abs(math.log10(value / candidate)))
        result.append(
            {
                "source_metallicity_mass_fraction": value,
                "nearest_other_source_metallicity_mass_fraction": nearest,
                "absolute_log10_offset_dex": abs(math.log10(value / nearest)),
            }
        )
    return result


def audit_candidate_grid_coverage(
    *, root: Path = DEFAULT_ROOT, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    limongi = adapt_candidate(LIMONGI_ID, root=Path(root), include_records=False)
    nugrid = adapt_candidate(NUGRID_ID, root=Path(root), include_records=False)
    huscher = audit_huscher2025_candidate(root=Path(root))
    boccioli = audit_boccioli_roberti2026_candidate(root=Path(root))
    doherty = audit_doherty2014_candidate(root=Path(root))
    stockinger = audit_stockinger2020_candidate(root=Path(root))
    sukhbold = audit_sukhbold2016_candidate(root=Path(root))
    transition_fate = audit_limongi2024_transition_fates(root=Path(root))
    roberti_ultralowz = audit_roberti2024_ultralowz_candidate(root=Path(root))
    heger_woosley_popiii = audit_heger_woosley2010_popiii_candidate(root=Path(root))
    baseline_demand = audit_g2_baseline_metallicity_demand()
    configured = contract["candidate_channel_nodes"]
    channel_3_branches = contract.get("channel_3_branch_nodes")
    if not isinstance(channel_3_branches, dict):
        raise CoverageAuditError("channel-3 branch-node inventory is missing")
    observed_nugrid = nugrid["source_axes"]["mass_msun"]
    if configured[NUGRID_ID]["2"] != [value for value in observed_nugrid if value <= 7.0]:
        raise CoverageAuditError("NuGrid AGB nodes drifted from coverage contract")
    if configured[NUGRID_ID]["1"] != [value for value in observed_nugrid if value >= 12.0]:
        raise CoverageAuditError("NuGrid massive-wind nodes drifted from coverage contract")
    limongi_axes = limongi["source_axes"]
    if configured[LIMONGI_ID]["1"] != limongi_axes["recommended_mass_msun"]:
        raise CoverageAuditError("Limongi wind/outcome nodes drifted from coverage contract")
    if configured[LIMONGI_ID]["3"] != limongi_axes["recommended_mass_msun"]:
        raise CoverageAuditError("Limongi SNII outcome nodes drifted from coverage contract")
    if channel_3_branches.get(LIMONGI_ID) != {
        "recommended_set_R": limongi_axes["recommended_mass_msun"]
    }:
        raise CoverageAuditError("Limongi channel-3 branch inventory drifted")
    if configured["huscher2025_agb"]["2"] != huscher["single_star_grid"]["mass_msun"]:
        raise CoverageAuditError("Huscher AGB nodes drifted from coverage contract")
    f23_single_masses = [float(value) for value in boccioli["grids"]["F23_single"]["mass_msun"]]
    if configured["boccioli_roberti2026_neutrino_ccsn"]["1"] != f23_single_masses:
        raise CoverageAuditError("Boccioli-Roberti F23 wind nodes drifted from coverage contract")
    boccioli_branch_masses = {
        branch: [float(value) for value in boccioli["grids"][branch]["mass_msun"]]
        for branch in ("F23_single", "F23_binary", "LC18", "WH07")
    }
    if channel_3_branches.get("boccioli_roberti2026_neutrino_ccsn") != boccioli_branch_masses:
        raise CoverageAuditError("Boccioli-Roberti channel-3 branch inventory drifted")
    boccioli_union = sorted(
        {value for values in boccioli_branch_masses.values() for value in values}
    )
    if configured["boccioli_roberti2026_neutrino_ccsn"]["3"] != boccioli_union:
        raise CoverageAuditError("Boccioli-Roberti F23 SNII nodes drifted from coverage contract")
    doherty_masses = [float(value) for value in doherty["primary_grid"]["mass_msun"]]
    if configured["doherty2014_sagb"]["2"] != doherty_masses:
        raise CoverageAuditError("Doherty SAGB nodes drifted from coverage contract")
    stockinger_masses = [float(value) for value in stockinger["model_grid"]["zams_mass_msun"]]
    if configured["stockinger2020_low_mass_ccsn"]["3"] != stockinger_masses:
        raise CoverageAuditError("Stockinger low-mass CCSN nodes drifted from coverage contract")
    sukhbold_masses = [float(value) for value in sukhbold["z96_grid"]["zams_mass_msun"]]
    high_mass_engines = sukhbold["high_mass_engine_evidence"]["engines"]
    sukhbold_branch_masses = {
        "Z9.6": sukhbold_masses,
        "N20_high_mass": sorted(
            float(value) for value in high_mass_engines["N20"]["high_mass_results"]
        ),
        "W18_high_mass": sorted(
            float(value) for value in high_mass_engines["W18"]["high_mass_results"]
        ),
    }
    if channel_3_branches.get("sukhbold2016_ccsn") != sukhbold_branch_masses:
        raise CoverageAuditError("Sukhbold channel-3 branch inventory drifted")
    sukhbold_union = sorted(
        {value for values in sukhbold_branch_masses.values() for value in values}
    )
    if configured["sukhbold2016_ccsn"]["3"] != sukhbold_union:
        raise CoverageAuditError("Sukhbold candidate-node union drifted from coverage contract")
    roberti_masses = [float(value) for value in roberti_ultralowz["source_grid"]["masses_msun"]]
    if configured["roberti2024_ultralowz_ccsn"]["3"] != roberti_masses:
        raise CoverageAuditError("Roberti ultra-low-Z CCSN nodes drifted from coverage contract")
    heger_woosley_masses = [
        float(value)
        for value in heger_woosley_popiii["source_grid"]["zams_masses_msun"]
    ]
    heger_woosley_runtime_masses = list(heger_woosley_masses)
    if configured["heger_woosley2010_popiii"]["3"] != heger_woosley_runtime_masses:
        raise CoverageAuditError("Heger--Woosley Pop III CCSN nodes drifted from coverage contract")
    if channel_3_branches.get("heger_woosley2010_popiii") != {
        "Z0_source_grid": heger_woosley_masses
    }:
        raise CoverageAuditError("Heger--Woosley channel-3 branch inventory drifted")
    transition_reference = contract["transition_fate_reference"]
    transition_policy = transition_fate["project_transition_policy"]
    if transition_reference["candidate_channel_nodes_contributed"] != 0:
        raise CoverageAuditError("fate reference unexpectedly contributes yield nodes")
    if transition_reference["runtime_edge_interval_msun"] != transition_policy["unresolved_runtime_edge_interval_msun"]:
        raise CoverageAuditError("transition-fate runtime edge drifted")
    if transition_reference["continuous_fate_interpolation_allowed"] is not False:
        raise CoverageAuditError("transition-fate reference unexpectedly permits interpolation")
    if transition_reference["stockinger_e8p8_anchor_may_define_population_fate_law"] is not False:
        raise CoverageAuditError("transition-fate reference promotes the e8.8 anchor")
    demand_reference = contract["runtime_demand_reference"]
    if demand_reference["baseline_id"] != baseline_demand["baseline_identity"]["baseline_id"]:
        raise CoverageAuditError("baseline metallicity-demand identity drifted")
    if demand_reference["baseline_role"] != baseline_demand["baseline_identity"]["role"]:
        raise CoverageAuditError("baseline metallicity-demand role drifted")
    if any(
        demand_reference[key] is not False
        for key in (
            "comparison_population_defines_production_domain",
            "metallicity_floor_or_clamp_allowed",
            "solar_source_extrapolation_to_ultra_low_z_allowed",
        )
    ):
        raise CoverageAuditError("baseline metallicity-demand reference is not fail closed")
    sparse_reference = contract["ultralowz_sparse_reference"]
    if sparse_reference["candidate_id"] != roberti_ultralowz["source_identity"]["candidate_id"]:
        raise CoverageAuditError("Roberti ultra-low-Z reference identity drifted")
    if sparse_reference["mass_nodes_msun"] != roberti_masses:
        raise CoverageAuditError("Roberti ultra-low-Z mass nodes drifted")
    if sparse_reference["metallicity_nodes_mass_fraction"] != roberti_ultralowz["source_grid"]["metallicity_mass_fraction"]:
        raise CoverageAuditError("Roberti ultra-low-Z metallicity nodes drifted")
    if any(
        sparse_reference[key] is not False
        for key in (
            "mass_interpolation_allowed",
            "metallicity_interpolation_allowed",
            "rotation_population_selected",
            "model_025z600_allowed",
            "source_only_mrt_merge_allowed",
        )
    ):
        raise CoverageAuditError("Roberti ultra-low-Z reference is not fail closed")
    popiii_reference = contract["popiii_mass_grid_reference"]
    if popiii_reference["candidate_id"] != heger_woosley_popiii["source_identity"]["candidate_id"]:
        raise CoverageAuditError("Heger--Woosley Pop III reference identity drifted")
    if popiii_reference["full_source_mass_hull_msun"] != [
        min(heger_woosley_masses), max(heger_woosley_masses)
    ]:
        raise CoverageAuditError("Heger--Woosley full mass hull drifted")
    if popiii_reference["runtime_channel_3_mass_hull_msun"] != [
        min(heger_woosley_runtime_masses), max(heger_woosley_runtime_masses)
    ]:
        raise CoverageAuditError("Heger--Woosley runtime mass hull drifted")
    if popiii_reference["metallicity_mass_fraction"] != 0.0:
        raise CoverageAuditError("Heger--Woosley metallicity scope drifted")
    if any(
        popiii_reference[key] is not False
        for key in (
            "mass_interpolation_allowed",
            "metallicity_extrapolation_allowed",
            "explosion_energy_distribution_selected",
            "piston_distribution_selected",
            "mixing_distribution_selected",
            "rotation_included",
        )
    ):
        raise CoverageAuditError("Heger--Woosley Pop III reference is not fail closed")

    mass_report: dict[str, Any] = {}
    for channel, required in contract["runtime_mass_ranges_msun"].items():
        by_source = {
            candidate: [float(value) for value in channels[channel]]
            for candidate, channels in configured.items()
            if channels[channel]
        }
        union = sorted({value for values in by_source.values() for value in values})
        source_hulls = {
            candidate: [min(values), max(values)] for candidate, values in by_source.items()
        }
        uncovered = _uncovered_edges([float(value) for value in required], union)
        uncovered_hulls = _uncovered_source_hulls(
            [float(value) for value in required], source_hulls
        )
        mass_report[channel] = {
            "runtime_required_range_msun": required,
            "candidate_nodes_by_source_msun": by_source,
            "candidate_union_nodes_msun": union,
            "source_hulls_msun": source_hulls,
            "union_node_hull_msun": [min(union), max(union)] if union else None,
            "uncovered_runtime_edge_intervals_msun": uncovered,
            "uncovered_runtime_source_hull_intervals_msun": uncovered_hulls,
            "runtime_edges_covered": not uncovered,
            "runtime_source_hulls_cover_domain": not uncovered_hulls,
            "cross_source_seams": (
                "present_and_not_interpolable_without_cross_calibration"
                if len(by_source) > 1
                else "none"
            ),
            "production_coverage_approved": False,
        }
        if channel == "3":
            mass_report[channel]["candidate_branch_nodes_by_source_msun"] = channel_3_branches
            mass_report[channel]["flattened_branch_union_is_interpolable"] = False
            mass_report[channel]["source_node_fate_and_remnant_records_required"] = True

    limongi_z_map = limongi_axes["metallicity_mass_fraction_from_source_article"]
    limongi_z = [float(limongi_z_map[str(value)]) for value in (-3, -2, -1, 0)]
    nugrid_z = [float(value) for value in nugrid["source_axes"]["metallicity_mass_fraction"]]
    huscher_z = [float(value) for value in huscher["single_star_grid"]["metallicity_mass_fraction"]]
    doherty_z = [float(value) for value in doherty["primary_grid"]["metallicity_mass_fraction"]]
    roberti_z = [float(value) for value in roberti_ultralowz["source_grid"]["metallicity_mass_fraction"]]
    exact_common_z = sorted(set(limongi_z) & set(nugrid_z) & set(huscher_z) & set(doherty_z))
    pairwise_common_z = {
        "limongi_nugrid": sorted(set(limongi_z) & set(nugrid_z)),
        "limongi_huscher": sorted(set(limongi_z) & set(huscher_z)),
        "nugrid_huscher": sorted(set(nugrid_z) & set(huscher_z)),
        "nugrid_doherty": sorted(set(nugrid_z) & set(doherty_z)),
        "huscher_doherty": sorted(set(huscher_z) & set(doherty_z)),
        "limongi_doherty": sorted(set(limongi_z) & set(doherty_z)),
    }
    metallicity_report = {
        "limongi_source_defined_mass_fraction": limongi_z,
        "limongi_mapping_status": limongi_axes["metallicity_mapping_status"],
        "nugrid_mass_fraction": nugrid_z,
        "huscher_mass_fraction": huscher_z,
        "doherty_mass_fraction": doherty_z,
        "boccioli_roberti_f23_metallicity": "source-labelled solar only; exact total-metallicity mass fraction not asserted by the staged release",
        "stockinger_metallicity": "e8.8 and s9.0 are source-labelled solar; z9.6 is zero metallicity; no continuous metallicity axis",
        "sukhbold_z96_metallicity": "source-labelled solar only; no continuous metallicity axis",
        "limongi2024_transition_fate_metallicity": "source-solar only; fate constraints may not be extrapolated across metallicity",
        "roberti2024_ultralowz_sparse_mass_fraction": roberti_z,
        "roberti2024_ultralowz_mass_nodes_msun": roberti_masses,
        "roberti2024_baseline_values_inside_zero_to_first_positive_coordinate_interval": (
            baseline_demand["stellar_population"]["birth_metallicity_mass_fraction"]["minimum"] >= roberti_z[0]
            and baseline_demand["stellar_population"]["birth_metallicity_mass_fraction"]["maximum"] <= roberti_z[1]
        ),
        "roberti2024_metallicity_interpolation_allowed": sparse_reference["metallicity_interpolation_allowed"],
        "roberti2024_production_domain_covered": False,
        "heger_woosley2010_popiii_exact_mass_fraction": [0.0],
        "heger_woosley2010_positive_metallicity_coverage": False,
        "heger_woosley2010_metallicity_extrapolation_allowed": popiii_reference[
            "metallicity_extrapolation_allowed"
        ],
        "pairwise_exact_common_nodes": pairwise_common_z,
        "exact_common_nodes": exact_common_z,
        "exact_common_node_count": len(exact_common_z),
        "limongi_to_nearest_nugrid_log_offsets": _nearest_log_offsets(
            limongi_z, nugrid_z
        ),
        "limongi_to_nearest_huscher_log_offsets": _nearest_log_offsets(
            limongi_z, huscher_z
        ),
        "required_runtime_domain": contract["coverage_rules"][
            "required_birth_metallicity_domain"
        ],
        "required_runtime_domain_status": contract["coverage_rules"][
            "required_birth_metallicity_domain_status"
        ],
        "production_coverage_approved": False,
    }
    duplicate_report = {
        "limongi_duplicate_coordinate_count": limongi["source_components"]
        ["evolutionary_properties"]["duplicate_model_phase_coordinate_count"],
        "limongi_all_duplicates_exactly_identical": limongi["source_components"]
        ["evolutionary_properties"]["all_duplicate_rows_physically_identical"],
        "nugrid_duplicate_policy": "three component copies at (5,0.01) were independently verified exact and collapse once in the partial projection",
        "non_identical_duplicate_count": 0,
    }
    blockers = [
        "Channel 1 has no candidate coverage below 11 Msun although its runtime range begins at 0.8 Msun.",
        "Doherty reaches the channel-2 upper edge at Z >= 0.004, but its Z <= 0.001 baseline grids stop at 7.5 Msun and do not establish full metallicity coverage.",
        "Sukhbold overlaps the Stockinger and F23 source hulls from 9--12 Msun, closing the prior internal 9.6--11 Msun hull gap, but channel 3 still has no yield-candidate hull from 8--8.8 Msun; Limongi 2024 reclassifies this edge as an unresolved, non-interpolable terminal-fate policy seam.",
        "The required production birth-metallicity domain has not been selected; all 42,342 stars in the inherited comparison catalogue lie below the lowest staged positive-Z full-grid candidate node, and no floor, solar extrapolation, or discrete zero-Z event substitution is allowed.",
        "No exact total-metallicity node is common to Limongi, NuGrid, Huscher, and Doherty, so the massive-star/AGB source seam is not cross-calibrated.",
        "The rotation distribution/marginalization for Limongi is not selected and NuGrid has no matching rotation axis.",
        "No source provides the required per-star age-resolved cumulative wind/AGB composition history.",
        "The Huscher IMF-weighted population Mdot table fails its stated unit-mass normalization check and cannot substitute for a per-star history.",
        "Boccioli-Roberti branch nodes reach 120 Msun, but F23 single/binary, LC18, and WH07 have incompatible population, metallicity, rotation, and engine coordinates; the flattened union is not interpolable and provides no canonical momentum.",
        "Doherty lacks calcium, an age-resolved history, and verified redistribution terms.",
        "Stockinger lacks nitrogen and a stable-decay/tracer projection, has no canonical event momentum, and has unresolved redistribution terms.",
        "Sukhbold is solar-only; its separated Z9.6 and W18/N20 high-mass branches do not form a continuous fate law, failed high-mass outcomes lack complete remnant/wind records, the radioactive inventory and neutrino-wind nucleosynthesis are incomplete, canonical momentum is absent, and redistribution permission is unresolved.",
        "Roberti 2024 reaches Z=0 and ultra-low positive Z but has only 15 and 25 Msun nodes, an unselected rotation distribution, four zero-Z MRT omissions, unresolved wind/terminal ownership, and a quarantined 025z600 mass-budget outlier; being inside its sparse metallicity coordinate interval is not production coverage.",
        "Heger & Woosley 2010 densely covers Z=0 from 10--100 Msun, but leaves both the runtime 8--10 and 100--120 Msun Pop III edges uncovered and has no positive-Z axis; explosion energy, piston, artificial mixing, rotation, and neutrino-wind policies are not production-selected.",
    ]
    return {
        "schema": "snrt-g2-candidate-grid-coverage-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "blocked_incomplete_mass_metallicity_rotation_age_coverage",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "contract_path": str(contract_path),
        "contract_sha256": _sha256(contract_path),
        "audit_code_sha256": _sha256(TOOL_PATH),
        "source_adapter_hashes": {
            LIMONGI_ID: limongi["adapter_code_sha256"],
            NUGRID_ID: nugrid["adapter_code_sha256"],
            "huscher2025_agb": huscher["audit_code_sha256"],
            "boccioli_roberti2026_neutrino_ccsn": boccioli["audit_code_sha256"],
            "doherty2014_sagb": doherty["audit_code_sha256"],
            "stockinger2020_low_mass_ccsn": stockinger["audit_code_sha256"],
            "sukhbold2016_ccsn": sukhbold["audit_code_sha256"],
            "limongi2024_transition_fates": transition_fate["audit_code_sha256"],
            "roberti2024_ultralowz_ccsn": roberti_ultralowz["audit_code_sha256"],
            "heger_woosley2010_popiii": heger_woosley_popiii["audit_code_sha256"],
            "baseline_metallicity_demand": baseline_demand["audit_code_sha256"],
        },
        "mass_coverage": mass_report,
        "channel_3_branch_inventory": {
            "nodes_by_source_and_branch_msun": channel_3_branches,
            "flattened_union_is_interpolable": False,
            "source_node_fate_and_remnant_records_required": True,
            "production_coverage_approved": False,
        },
        "metallicity_coverage": metallicity_report,
        "rotation_coverage": {
            LIMONGI_ID: limongi_axes["rotation_velocity_km_s"],
            NUGRID_ID: None,
            "huscher2025_agb": None,
            "boccioli_roberti2026_neutrino_ccsn": {
                "F23_single": None,
                "F23_binary": "binary-stripped branch present but population weighting not selected",
                "LC18": [0, 150, 300],
            },
            "doherty2014_sagb": None,
            "stockinger2020_low_mass_ccsn": None,
            "sukhbold2016_ccsn": None,
            "roberti2024_ultralowz_ccsn": roberti_ultralowz["source_grid"]["rotation_km_s_by_mass_metallicity"],
            "heger_woosley2010_popiii": None,
            "population_selection_status": "not_selected",
            "production_coverage_approved": False,
        },
        "release_history_coverage": {
            LIMONGI_ID: "integrated terminal/wind snapshots only",
            NUGRID_ID: "integrated total/wind/pre-explosion snapshots plus lifetime only",
            "huscher2025_agb": "integrated per-star yield snapshots plus an IMF-weighted Mdot/CNO population table with unresolved normalization",
            "boccioli_roberti2026_neutrino_ccsn": "integrated Post/Wind/Presn snapshots only; no age-resolved cumulative history",
            "doherty2014_sagb": "integrated gross/net SAGB wind yields with synthetic-pulse sensitivity columns; no per-star age-resolved history",
            "stockinger2020_low_mass_ccsn": "three post-bounce event snapshots per model; not a stellar-age cumulative history",
            "sukhbold2016_ccsn": "integrated terminal explosion and presupernova-wind yields; not an age-resolved cumulative history",
            "limongi2024_transition_fates": "thermal-pulse structural history only; no age-resolved ejecta composition or terminal-event yield",
            "roberti2024_ultralowz_ccsn": "integrated explosive yields only; no separate wind/terminal partition or stellar-age cumulative history",
            "heger_woosley2010_popiii": "integrated post-supernova ejecta at Z=0 only; no age-resolved release history and neutrino-wind nucleosynthesis is omitted",
            "age_resolved_required": True,
            "production_coverage_approved": False,
        },
        "duplicate_resolution": duplicate_report,
        "transition_fate_coverage": {
            "status": transition_fate["status"],
            "source_scope": transition_reference["source_scope"],
            "runtime_edge_interval_msun": transition_policy["unresolved_runtime_edge_interval_msun"],
            "runtime_edge_classification": transition_reference["runtime_edge_classification"],
            "potential_ecsn_model_interval_msun": transition_fate["source_reported_fate_statements"]["potential_ecsn_model_interval_msun"],
            "ordinary_core_collapse_lower_model_mass_msun": transition_fate["source_reported_fate_statements"]["ordinary_core_collapse_lower_model_mass_msun"],
            "continuous_fate_interpolation_allowed": transition_policy["continuous_fate_interpolation_allowed"],
            "stockinger_e8p8_anchor_may_define_population_fate_law": transition_policy["stockinger_e8p8_anchor_may_define_population_fate_law"],
            "candidate_channel_nodes_contributed": transition_reference["candidate_channel_nodes_contributed"],
            "production_fate_policy_approved": False,
        },
        "baseline_metallicity_demand": {
            "status": baseline_demand["status"],
            "baseline_id": baseline_demand["baseline_identity"]["baseline_id"],
            "baseline_role": baseline_demand["baseline_identity"]["role"],
            "comparison_population_defines_production_domain": baseline_demand["baseline_identity"]["defines_production_domain"],
            "star_count": baseline_demand["stellar_population"]["star_count"],
            "observed_birth_metallicity_mass_fraction": [
                baseline_demand["stellar_population"]["birth_metallicity_mass_fraction"]["minimum"],
                baseline_demand["stellar_population"]["birth_metallicity_mass_fraction"]["maximum"],
            ],
            "lowest_positive_full_grid_candidate_metallicity_mass_fraction": baseline_demand["candidate_domain_comparison"]["lowest_positive_full_grid_candidate_metallicity_mass_fraction"],
            "fraction_below_lowest_positive_full_grid_candidate": baseline_demand["candidate_domain_comparison"]["fraction_below_lowest_positive_full_grid_candidate"],
            "maximum_baseline_z_to_candidate_lower_edge_offset_dex": baseline_demand["candidate_domain_comparison"]["maximum_baseline_z_to_candidate_lower_edge_offset_dex"],
            "metallicity_floor_or_clamp_allowed": baseline_demand["policy"]["metallicity_floor_or_clamp_allowed"],
            "solar_source_extrapolation_to_ultra_low_z_allowed": baseline_demand["policy"]["solar_source_extrapolation_to_ultra_low_z_allowed"],
            "production_domain_selected": baseline_demand["policy"]["required_production_birth_metallicity_domain_selected"],
        },
        "ultralowz_sparse_candidate": {
            "status": roberti_ultralowz["status"],
            "mass_nodes_msun": roberti_masses,
            "metallicity_nodes_mass_fraction": roberti_z,
            "baseline_values_inside_zero_to_first_positive_coordinate_interval": metallicity_report["roberti2024_baseline_values_inside_zero_to_first_positive_coordinate_interval"],
            "official_mrt_model_count": roberti_ultralowz["yield_model_inventory"]["official_mrt_model_count"],
            "source_only_models_missing_from_official_mrt": roberti_ultralowz["yield_model_inventory"]["source_only_models_missing_from_official_mrt"],
            "mass_budget_outlier_models": roberti_ultralowz["mass_budget_review"]["outlier_models"],
            "mass_interpolation_allowed": sparse_reference["mass_interpolation_allowed"],
            "metallicity_interpolation_allowed": sparse_reference["metallicity_interpolation_allowed"],
            "rotation_population_selected": sparse_reference["rotation_population_selected"],
            "production_coverage_approved": False,
        },
        "popiii_mass_grid_candidate": {
            "status": heger_woosley_popiii["status"],
            "metallicity_mass_fraction": 0.0,
            "full_source_mass_hull_msun": [
                min(heger_woosley_masses), max(heger_woosley_masses)
            ],
            "runtime_channel_3_mass_hull_msun": [
                min(heger_woosley_runtime_masses), max(heger_woosley_runtime_masses)
            ],
            "runtime_channel_3_uncovered_edge_intervals_msun": _uncovered_edges(
                [float(value) for value in contract["runtime_mass_ranges_msun"]["3"]],
                heger_woosley_runtime_masses,
            ),
            "source_mass_node_count": len(heger_woosley_masses),
            "source_coordinate_count": heger_woosley_popiii["source_grid"]["coordinate_count"],
            "explosion_energy_distribution_selected": popiii_reference[
                "explosion_energy_distribution_selected"
            ],
            "piston_distribution_selected": popiii_reference["piston_distribution_selected"],
            "mixing_distribution_selected": popiii_reference["mixing_distribution_selected"],
            "mass_interpolation_allowed": popiii_reference["mass_interpolation_allowed"],
            "metallicity_extrapolation_allowed": popiii_reference[
                "metallicity_extrapolation_allowed"
            ],
            "production_coverage_approved": False,
        },
        "blockers": blockers,
        "interpretation": (
            "The source-defined Limongi [Fe/H]-to-Z mapping and exact duplicate "
            "classification are resolved. Candidate union hulls still do not "
            "constitute an interpolable production grid because edge gaps, source "
            "seams, rotation assumptions, the runtime Z domain, age histories, and the "
            "Huscher population-table normalization remain open. Sukhbold closes the "
            "candidate source-hull gap from 9.6--11 Msun through overlaps with Stockinger "
            "and F23. Channel-3 candidate coordinates now retain source branch identity "
            "through 120 Msun; their flattened union is explicitly non-interpolable. "
            "Limongi 2024 shows that the 8--8.8 Msun edge is a terminal-fate "
            "policy seam rather than an interval that may be populated by yield interpolation. "
            "The inherited comparison population additionally lies more than 4.43 dex below the "
            "lowest staged positive-Z full-grid lower edge at its maximum Z. The fate seam, "
            "Roberti 2024 now supplies exact ultra-low-Z candidate coordinates, but only at "
            "15 and 25 Msun and without an approved rotation/interpolation/component policy. "
            "Heger & Woosley 2010 adds dense exact-Z=0 mass coverage from 10--100 Msun, "
            "but no positive-Z axis or approved explosion, piston, mixing, rotation, or "
            "8--10 or 100--120 Msun edge policy. "
            "The fate seam, ultra-low-Z production domain, and all approval constraints remain open."
        ),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_candidate_grid_coverage(
            root=args.root, contract_path=args.contract
        )
    except (
        CoverageAuditError,
        SourceAdapterError,
        HuscherAuditError,
        BoccioliRobertiAuditError,
        DohertyAuditError,
        StockingerAuditError,
        SukhboldAuditError,
        LimongiTransitionAuditError,
        RobertiUltraLowZAuditError,
        HegerWoosleyPopIIIAuditError,
        BaselineMetallicityDemandError,
    ) as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, indent=2), file=sys.stderr)
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
