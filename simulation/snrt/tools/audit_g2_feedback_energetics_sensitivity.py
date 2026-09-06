#!/usr/bin/env python3
"""Audit review-only energy and scalar-momentum overlays for G2 feedback.

The source projection deliberately keeps unavailable energy and momentum null.
This tool evaluates physically distinct sensitivity budgets without mutating
those nulls, choosing a canonical model, or emitting canonical yield rows.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

from build_g2_nugrid_channel_projection import (
    ChannelProjectionError,
    build_channel_projection,
)
from audit_g2_stockinger2020_candidate import (
    StockingerAuditError,
    audit_stockinger2020_candidate,
)
from audit_g2_sukhbold2016_candidate import (
    SukhboldAuditError,
    audit_sukhbold2016_candidate,
)
from audit_g2_roberti2024_ultralowz_candidate import (
    RobertiUltraLowZAuditError,
    audit_roberti2024_ultralowz_candidate,
)
from audit_g2_heger_woosley2010_popiii_candidate import (
    HegerWoosleyPopIIIAuditError,
    audit_heger_woosley2010_popiii_candidate,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = (
    SNRT_ROOT / "config" / "g2_feedback_energetics_sensitivity_contract_v1.json"
)
DEFAULT_CANDIDATE_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"


class EnergeticsAuditError(ValueError):
    """The review overlay violates its fail-closed sensitivity contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise EnergeticsAuditError(f"cannot read {path}: {exc}") from exc
    return digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EnergeticsAuditError(f"cannot read energetics contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema")
        != "snrt-g2-feedback-energetics-sensitivity-contract"
        or contract.get("schema_version") != 1
    ):
        raise EnergeticsAuditError("unsupported G2 energetics sensitivity contract")
    approval = contract.get("approval", {})
    if approval.get("canonical_conversion_allowed") is not False:
        raise EnergeticsAuditError("sensitivity overlay unexpectedly permits conversion")
    if approval.get("runtime_deposition_allowed") is not False:
        raise EnergeticsAuditError("sensitivity overlay unexpectedly permits deposition")
    vector = contract.get("source_frame_momentum", {}).get(
        "unresolved_isotropic_vector_g_cm_s_per_star"
    )
    if vector != [0.0, 0.0, 0.0]:
        raise EnergeticsAuditError("isotropic source-frame vector must be exactly zero")
    return contract


def _terminal_models(projection: dict[str, Any]) -> dict[int, list[dict[str, Any]]]:
    rows = projection.get("rows")
    if not isinstance(rows, list) or not rows:
        raise EnergeticsAuditError("NuGrid projection rows are unavailable")
    selected: dict[tuple[int, float, float], dict[str, Any]] = {}
    for row in rows:
        if row["energy_erg_per_star"] is not None:
            raise EnergeticsAuditError("source energy null was replaced before sensitivity audit")
        if row["momentum_g_cm_s_per_star"] is not None:
            raise EnergeticsAuditError("source momentum null was replaced before sensitivity audit")
        key = (
            int(row["channel"]),
            float(row["initial_mass_msun_per_star"]),
            float(row["birth_metallicity_mass_fraction"]),
        )
        if key not in selected or row["age_yr"] > selected[key]["age_yr"]:
            selected[key] = row
    models: dict[int, list[dict[str, Any]]] = {1: [], 2: [], 3: []}
    for key, row in sorted(selected.items()):
        channel = key[0]
        if channel not in models:
            raise EnergeticsAuditError(f"unexpected channel {channel}")
        if row["returned_mass_msun_per_star"] < 0.0:
            raise EnergeticsAuditError(f"negative terminal returned mass at {key}")
        models[channel].append(row)
    if {channel: len(values) for channel, values in models.items()} != {
        1: 20,
        2: 40,
        3: 20,
    }:
        raise EnergeticsAuditError("unexpected NuGrid terminal model counts")
    return models


def _range_summary(values: list[float]) -> dict[str, float]:
    if not values or any(not math.isfinite(value) or value < 0.0 for value in values):
        raise EnergeticsAuditError("non-finite or negative sensitivity value")
    return {
        "minimum": min(values),
        "maximum": max(values),
        "unweighted_grid_sum": math.fsum(values),
    }


def _wind_scenarios(
    models: list[dict[str, Any]], velocities_km_s: list[float], solar_mass_g: float,
    km_s_to_cm_s: float,
) -> list[dict[str, Any]]:
    scenarios: list[dict[str, Any]] = []
    for velocity_km_s in velocities_km_s:
        if not math.isfinite(velocity_km_s) or velocity_km_s <= 0.0:
            raise EnergeticsAuditError("wind velocity nodes must be finite and positive")
        velocity_cgs = velocity_km_s * km_s_to_cm_s
        masses_g = [row["returned_mass_msun_per_star"] * solar_mass_g for row in models]
        energies = [0.5 * mass_g * velocity_cgs**2 for mass_g in masses_g]
        radial_momenta = [mass_g * velocity_cgs for mass_g in masses_g]
        scenarios.append(
            {
                "velocity_km_s": velocity_km_s,
                "kinetic_energy_erg_per_model": _range_summary(energies),
                "scalar_radial_momentum_g_cm_s_per_model": _range_summary(
                    radial_momenta
                ),
            }
        )
    return scenarios


def _supernova_scenarios(
    models: list[dict[str, Any]], energies_erg: list[float], densities_cm3: list[float],
    solar_mass_g: float, km_s_to_cm_s: float,
) -> dict[str, Any]:
    ejecta_masses_g = [
        row["returned_mass_msun_per_star"] * solar_mass_g for row in models
    ]
    launch: list[dict[str, Any]] = []
    for energy_erg in energies_erg:
        if not math.isfinite(energy_erg) or energy_erg <= 0.0:
            raise EnergeticsAuditError("SN energy nodes must be finite and positive")
        launch_momenta = [
            math.sqrt(2.0 * mass_g * energy_erg) for mass_g in ejecta_masses_g
        ]
        launch.append(
            {
                "explosion_energy_erg_per_event": energy_erg,
                "scalar_ejecta_launch_momentum_g_cm_s_per_model": _range_summary(
                    launch_momenta
                ),
                "unweighted_grid_energy_sum_erg": len(models) * energy_erg,
            }
        )
    terminal: list[dict[str, Any]] = []
    for density in densities_cm3:
        if not math.isfinite(density) or density <= 0.0:
            raise EnergeticsAuditError("ambient-density nodes must be finite and positive")
        momentum_msun_km_s = 2.8e5 * density ** (-0.17)
        terminal.append(
            {
                "ambient_hydrogen_number_density_cm3": density,
                "scalar_terminal_shell_momentum_msun_km_s_per_event": momentum_msun_km_s,
                "scalar_terminal_shell_momentum_g_cm_s_per_event": (
                    momentum_msun_km_s * solar_mass_g * km_s_to_cm_s
                ),
                "unweighted_grid_scalar_terminal_shell_momentum_g_cm_s": (
                    len(models) * momentum_msun_km_s * solar_mass_g * km_s_to_cm_s
                ),
            }
        )
    return {
        "launch_energy_sensitivity": launch,
        "terminal_shell_density_sensitivity_at_published_1e51_erg_calibration": terminal,
    }


def _stockinger_energy_anchors(
    report: dict[str, Any], solar_mass_g: float
) -> dict[str, Any]:
    if report["canonical_rows_emitted"] != 0 or report["production_ready"] is not False:
        raise EnergeticsAuditError("Stockinger review unexpectedly permits production use")
    if report["model_grid"]["cross_model_interpolation_allowed"] is not False:
        raise EnergeticsAuditError("Stockinger event anchors must remain non-interpolable")
    energy_models = report["diagnostic_explosion_energy"]["models"]
    grid_models = report["model_grid"]["models"]
    records: list[dict[str, Any]] = []
    for model in ("e8.8", "s9.0", "z9.6"):
        energy = float(energy_models[model]["last_finite_diagnostic_energy_erg"])
        ejecta_mass_msun = float(grid_models[model]["final_total_ejecta_msun"])
        if energy <= 0.0 or ejecta_mass_msun <= 0.0:
            raise EnergeticsAuditError("invalid Stockinger energy anchor")
        records.append(
            {
                "model": model,
                "zams_mass_msun": grid_models[model]["zams_mass_msun"],
                "metallicity": grid_models[model]["metallicity"],
                "collapse_type": grid_models[model]["collapse_type"],
                "diagnostic_energy_sample_time_postbounce_s": energy_models[model]["last_finite_time_postbounce_s"],
                "diagnostic_explosion_energy_erg": energy,
                "final_event_ejecta_mass_msun": ejecta_mass_msun,
                "derived_scalar_ejecta_launch_momentum_g_cm_s": math.sqrt(
                    2.0 * ejecta_mass_msun * solar_mass_g * energy
                ),
                "canonical_energy_selected": False,
                "canonical_momentum_selected": False,
            }
        )
    return {
        "source_candidate_id": report["source_identity"]["candidate_id"],
        "records": records,
        "cross_model_interpolation_allowed": False,
        "vsh_dataset_used": False,
        "vsh_metadata_quarantined": report["diagnostic_explosion_energy"]["vsh_quarantined"],
        "source_audit_code_sha256": report["audit_code_sha256"],
    }


def _sukhbold_energy_anchors(report: dict[str, Any]) -> dict[str, Any]:
    if report["canonical_rows_emitted"] != 0 or report["production_ready"] is not False:
        raise EnergeticsAuditError("Sukhbold review unexpectedly permits production use")
    grid = report["z96_grid"]
    if grid["cross_engine_interpolation_allowed"] is not False:
        raise EnergeticsAuditError("Sukhbold energy grid must remain engine-specific")
    records: list[dict[str, Any]] = []
    for mass in grid["zams_mass_msun"]:
        source = grid["models"][str(mass)]
        energy = float(source["final_kinetic_energy_erg"])
        if energy <= 0.0:
            raise EnergeticsAuditError("invalid Sukhbold final kinetic-energy anchor")
        records.append(
            {
                "zams_mass_msun": mass,
                "metallicity": grid["metallicity"],
                "engine": grid["engine"],
                "final_kinetic_energy_erg": energy,
                "baryonic_mass_cut_after_fallback_msun": source[
                    "baryonic_mass_cut_after_fallback_msun"
                ],
                "fallback_mass_msun": source["fallback_mass_msun"],
                "stable_event_ejecta_sum_msun": source["stable_ejecta_sum_msun"],
                "exact_total_event_ejecta_mass_available": False,
                "derived_scalar_ejecta_launch_momentum_g_cm_s": None,
                "canonical_energy_selected": False,
                "canonical_momentum_selected": False,
            }
        )
    return {
        "source_candidate_id": report["source_identity"]["candidate_id"],
        "records": records,
        "cross_engine_interpolation_allowed": False,
        "cross_source_interpolation_allowed": False,
        "exact_ejecta_launch_momentum_derived": False,
        "source_audit_code_sha256": report["audit_code_sha256"],
    }


def _roberti_ultralowz_energy_anchors(
    report: dict[str, Any], contract: dict[str, Any]
) -> dict[str, Any]:
    if report["canonical_rows_emitted"] != 0 or report["production_ready"] is not False:
        raise EnergeticsAuditError("Roberti review unexpectedly permits production use")
    expected = contract["source_table_ultralowz_energy_anchors"]
    if report["source_identity"]["candidate_id"] != expected["source_candidate_id"]:
        raise EnergeticsAuditError("Roberti energetics source identity drifted")
    if report["source_grid"]["model_count"] != expected["model_count"]:
        raise EnergeticsAuditError("Roberti energetics model count drifted")
    if report["source_grid"]["masses_msun"] != expected["zams_mass_msun"]:
        raise EnergeticsAuditError("Roberti energetics mass grid drifted")
    if report["source_grid"]["metallicity_mass_fraction"] != expected["metallicity_mass_fraction"]:
        raise EnergeticsAuditError("Roberti energetics metallicity grid drifted")
    records: list[dict[str, Any]] = []
    metallicity_by_code = {"z": 0.0, "f": 3.236e-7, "e": 3.236e-6}
    outliers = report["mass_budget_review"]["outlier_models"]
    for model, source in sorted(report["supernova_properties"]["records"].items()):
        records.append(
            {
                "model": model,
                "zams_mass_msun": float(int(model[:3])),
                "metallicity_mass_fraction": metallicity_by_code[model[3]],
                "rotation_km_s": int(model[4:]),
                "iron_core_mass_msun": source["iron_core_mass_msun"],
                "remnant_mass_msun": source["remnant_mass_msun"],
                "thermal_bomb_kinetic_energy_erg": source["explosion_kinetic_energy_erg"],
                "derived_scalar_ejecta_launch_momentum_g_cm_s": None,
                "mass_budget_quarantined": model in outliers,
                "canonical_energy_selected": False,
                "canonical_momentum_selected": False,
            }
        )
    return {
        "source_candidate_id": report["source_identity"]["candidate_id"],
        "records": records,
        "energy_semantics": report["supernova_properties"]["semantics"],
        "mass_interpolation_allowed": False,
        "metallicity_interpolation_allowed": False,
        "rotation_interpolation_or_marginalization_allowed": False,
        "exact_ejecta_launch_momentum_derived": False,
        "mass_budget_quarantined_models": outliers,
        "source_audit_code_sha256": report["audit_code_sha256"],
    }


def _heger_woosley_popiii_energy_grid(
    report: dict[str, Any], contract: dict[str, Any]
) -> dict[str, Any]:
    if report["canonical_rows_emitted"] != 0 or report["production_ready"] is not False:
        raise EnergeticsAuditError("Heger--Woosley review unexpectedly permits production use")
    expected = contract["source_table_popiii_energy_grid"]
    if report["source_identity"]["candidate_id"] != expected["source_candidate_id"]:
        raise EnergeticsAuditError("Heger--Woosley energetics source identity drifted")
    grid = report["source_grid"]
    checks = {
        "metallicity_mass_fraction": grid["metallicity_mass_fraction"],
        "zams_mass_hull_msun": [
            grid["zams_mass_msun_minimum"], grid["zams_mass_msun_maximum"]
        ],
        "zams_mass_node_count": grid["zams_mass_count"],
        "coordinate_count": grid["coordinate_count"],
        "s4_kinetic_energy_bethe": grid["s4_kinetic_energy_bethe"],
        "ye_kinetic_energy_bethe": grid["ye_kinetic_energy_bethe"],
        "mixing_normalized_to_he_core": grid["mixing_normalized_to_he_core"],
    }
    for key, observed in checks.items():
        if observed != expected[key]:
            raise EnergeticsAuditError(f"Heger--Woosley energetics {key} drifted")
    return {
        "source_candidate_id": report["source_identity"]["candidate_id"],
        **checks,
        "kinetic_energy_at_infinity_erg_range": [3.0e50, 1.0e52],
        "energy_semantics": report["physical_semantics"]["energy_quantity"],
        "piston_locations": report["physical_semantics"]["piston_locations"],
        "fallback_included": report["physical_semantics"]["fallback_included"],
        "explicit_remnant_mass_available": report["physical_semantics"][
            "explicit_remnant_mass_column_available"
        ],
        "inferred_remnant_mass_promoted": False,
        "explosion_energy_distribution_selected": False,
        "piston_distribution_selected": False,
        "mixing_distribution_selected": False,
        "mass_interpolation_allowed": False,
        "metallicity_extrapolation_allowed": False,
        "exact_ejecta_launch_momentum_derived": False,
        "source_audit_code_sha256": report["audit_code_sha256"],
    }


def audit_energetics_sensitivity(
    *, candidate_root: Path = DEFAULT_CANDIDATE_ROOT,
    contract_path: Path = DEFAULT_CONTRACT,
) -> dict[str, Any]:
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    projection = build_channel_projection(
        root=Path(candidate_root), include_rows=True
    )
    stockinger = audit_stockinger2020_candidate(root=Path(candidate_root))
    sukhbold = audit_sukhbold2016_candidate(root=Path(candidate_root))
    roberti_ultralowz = audit_roberti2024_ultralowz_candidate(root=Path(candidate_root))
    heger_woosley_popiii = audit_heger_woosley2010_popiii_candidate(
        root=Path(candidate_root)
    )
    if projection["canonical_conversion_allowed"] is not False:
        raise EnergeticsAuditError("input projection unexpectedly permits conversion")
    models = _terminal_models(projection)
    constants = contract["physical_constants"]
    solar_mass_g = float(constants["solar_mass_g"])
    km_s_to_cm_s = float(constants["km_s_to_cm_s"])
    if solar_mass_g <= 0.0 or km_s_to_cm_s <= 0.0:
        raise EnergeticsAuditError("invalid physical constants")
    channels = contract["channel_sensitivity_models"]
    massive_wind = _wind_scenarios(
        models[1],
        [float(value) for value in channels["1"]["terminal_velocity_km_s"]],
        solar_mass_g,
        km_s_to_cm_s,
    )
    agb_wind = _wind_scenarios(
        models[2],
        [float(value) for value in channels["2"]["terminal_velocity_km_s"]],
        solar_mass_g,
        km_s_to_cm_s,
    )
    supernova = _supernova_scenarios(
        models[3],
        [float(value) for value in channels["3"]["explosion_energy_erg_per_event"]],
        [
            float(value)
            for value in channels["3"]["ambient_hydrogen_number_density_cm3"]
        ],
        solar_mass_g,
        km_s_to_cm_s,
    )
    stockinger_anchors = _stockinger_energy_anchors(stockinger, solar_mass_g)
    sukhbold_anchors = _sukhbold_energy_anchors(sukhbold)
    roberti_ultralowz_anchors = _roberti_ultralowz_energy_anchors(
        roberti_ultralowz, contract
    )
    heger_woosley_popiii_grid = _heger_woosley_popiii_energy_grid(
        heger_woosley_popiii, contract
    )
    total_model_count = sum(len(values) for values in models.values())
    return {
        "schema": "snrt-g2-feedback-energetics-sensitivity-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "review_only_blocked_model_selection_required",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "runtime_deposition_allowed": False,
        "contract_path": str(contract_path),
        "contract_sha256": _sha256(contract_path),
        "audit_code_sha256": _sha256(TOOL_PATH),
        "input_projection": {
            "status": projection["status"],
            "source_candidate_id": projection["source_candidate_id"],
            "model_count_by_channel": {
                str(channel): len(values) for channel, values in models.items()
            },
            "source_energy_and_momentum_nulls_preserved": True,
            "projection_contract_sha256": projection["contract_sha256"],
            "projection_builder_code_sha256": projection["builder_code_sha256"],
        },
        "source_frame_vector_momentum": {
            "semantics": contract["source_frame_momentum"]["canonical_quantity"],
            "isotropic_vector_g_cm_s_per_star": [0.0, 0.0, 0.0],
            "candidate_isotropic_model_count": total_model_count,
            "bulk_advective_momentum_separate": True,
            "scalar_radial_momentum_separate": True,
        },
        "machine_readable_event_energy_anchors": stockinger_anchors,
        "machine_readable_mass_grid_energy_anchors": sukhbold_anchors,
        "source_table_ultralowz_energy_anchors": roberti_ultralowz_anchors,
        "source_table_popiii_energy_grid": heger_woosley_popiii_grid,
        "channel_sensitivity": {
            "1_massive_star_wind": {
                "model_count": len(models[1]),
                "returned_mass_unweighted_grid_sum_msun": math.fsum(
                    row["returned_mass_msun_per_star"] for row in models[1]
                ),
                "scenarios": massive_wind,
            },
            "2_agb_wind": {
                "model_count": len(models[2]),
                "returned_mass_unweighted_grid_sum_msun": math.fsum(
                    row["returned_mass_msun_per_star"] for row in models[2]
                ),
                "scenarios": agb_wind,
            },
            "3_core_collapse_supernova": {
                "model_count": len(models[3]),
                "returned_mass_unweighted_grid_sum_msun": math.fsum(
                    row["returned_mass_msun_per_star"] for row in models[3]
                ),
                **supernova,
            },
        },
        "interpretation": (
            "The overlay quantifies sensitivity but does not fill source nulls. "
            "Wind/AGB radial budgets are velocity-dependent, SN launch momentum is "
            "ejecta-mass and energy dependent, and terminal shell momentum is an "
            "environment-dependent scalar rather than a canonical net vector."
        ),
        "blockers": contract["approval"]["required_before_approval"],
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_energetics_sensitivity(
            candidate_root=args.candidate_root, contract_path=args.contract
        )
    except (
        EnergeticsAuditError,
        ChannelProjectionError,
        StockingerAuditError,
        SukhboldAuditError,
        RobertiUltraLowZAuditError,
        HegerWoosleyPopIIIAuditError,
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
