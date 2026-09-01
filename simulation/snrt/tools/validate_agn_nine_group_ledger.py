#!/usr/bin/env python3
"""Validate the canonical AGN ledger against the closed P0 nine-group table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.primordial import (
    H_I_FIT,
    HE_I_FIT,
    HE_II_FIT,
    group_spectral_closure_from_metadata,
)
from snrt_core.source_ledger import read_photon_source_ledger_csv


DEFAULT_LEDGER = ROOT / "data" / "p4_pilot_agn_photon_ledger.csv"
DEFAULT_METADATA = ROOT / "data" / "p4_pilot_agn_photon_ledger.json"
DEFAULT_GROUP_EDGES = ROOT / "config" / "p0_photon_group_edges_ev.txt"
GENERATOR = ROOT / "tools" / "p4_build_agn_photon_ledger.py"
SOURCE_REBIND_TOOL = ROOT / "tools" / "p4_attach_pilot_sources.py"
P4_RUNNER = ROOT / "tools" / "p4_run_transport_pilot.py"
DEFAULT_STATIC_INPUT = ROOT / "data" / "p4_coeval_static_rt_input_agn9.h5"
DEFAULT_STATIC_METADATA = ROOT / "data" / "p4_coeval_static_rt_input_agn9.json"
DEFAULT_TRANSPORT_CONTROL = (
    ROOT / "data" / "p4_validation" / "p4_agn9_stage4_0p001myr.h5"
)
DEFAULT_FAILED_FULL_CFL_PROBE = (
    ROOT / "data" / "p4_validation" / "p4_agn9_stage4_one_step.h5"
)
DEFAULT_EXTERNAL_ASSET_MANIFEST = ROOT / "data" / "agn_nine_group_external_assets.json"
EXPECTED_INTERVAL_CONVENTION = "left_closed_right_open_except_final_closed"
SED_MINIMUM_EV = 10.0
EV_TO_ERG = 1.602176634e-12


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_head() -> str:
    return subprocess.check_output(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, text=True
    ).strip()


def unique_interval_index(intervals: np.ndarray, target: tuple[float, float]) -> int:
    matches = np.flatnonzero(
        np.all(intervals == np.asarray(target, dtype=np.float64), axis=1)
    )
    if matches.size != 1:
        raise ValueError(f"expected exactly one configured interval {target}, found {matches.size}")
    return int(matches[0])


def read_candidate_contract(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("candidate CSV requires a header")
        required = {
            "source_id",
            "source_kind",
            "x_code",
            "y_code",
            "z_code",
            "bolometric_luminosity_erg_s",
        }
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"candidate CSV missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("candidate CSV is empty")
    return {
        "source_id": np.asarray([int(row["source_id"]) for row in rows], dtype=np.int64),
        "source_kind": np.asarray([row["source_kind"] for row in rows], dtype=str),
        "position_code": np.asarray(
            [[float(row[name]) for name in ("x_code", "y_code", "z_code")] for row in rows],
            dtype=np.float64,
        ),
        "bolometric_luminosity_erg_s": np.asarray(
            [float(row["bolometric_luminosity_erg_s"]) for row in rows],
            dtype=np.float64,
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--group-edges", type=Path, default=DEFAULT_GROUP_EDGES)
    parser.add_argument("--static-input", type=Path, default=DEFAULT_STATIC_INPUT)
    parser.add_argument(
        "--static-metadata", type=Path, default=DEFAULT_STATIC_METADATA
    )
    parser.add_argument(
        "--transport-control", type=Path, default=DEFAULT_TRANSPORT_CONTROL
    )
    parser.add_argument(
        "--external-asset-manifest",
        type=Path,
        default=DEFAULT_EXTERNAL_ASSET_MANIFEST,
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    ledger = read_photon_source_ledger_csv(args.ledger)
    closure = group_spectral_closure_from_metadata(metadata)
    expected_edges = np.loadtxt(args.group_edges, comments="#", dtype=np.float64)
    metadata_edges = np.asarray(metadata["group_edges_ev"], dtype=np.float64)
    intervals = np.asarray(
        [group["energy_interval_ev"] for group in metadata["groups"]],
        dtype=np.float64,
    )
    expected_intervals = np.column_stack((expected_edges[:-1], expected_edges[1:]))
    group_indices = np.asarray(
        [group["index"] for group in metadata["groups"]], dtype=np.int64
    )
    group_means = np.asarray(
        [group["photon_weighted_mean_energy_ev"] for group in metadata["groups"]],
        dtype=np.float64,
    )
    declared_totals = np.asarray(
        metadata["group_photon_rate_total_s"], dtype=np.float64
    )
    grouped_totals = np.asarray(
        [group["total_photon_rate_s"] for group in metadata["groups"]],
        dtype=np.float64,
    )
    csv_totals = ledger.photon_luminosity_s.sum(axis=0, dtype=np.float64)
    cross_sections = np.asarray(
        (
            closure.cross_sections.hydrogen_i,
            closure.cross_sections.helium_i,
            closure.cross_sections.helium_ii,
        ),
        dtype=np.float64,
    )
    excess_energy = np.asarray(
        closure.photoelectron_excess_energy_ev, dtype=np.float64
    )
    provenance = metadata.get("provenance", {})
    candidates = Path(metadata["candidates"])
    candidate_contract = read_candidate_contract(candidates)
    candidate_lbol = candidate_contract["bolometric_luminosity_erg_s"]
    static_metadata = json.loads(args.static_metadata.read_text(encoding="utf-8"))
    external_assets = json.loads(
        args.external_asset_manifest.read_text(encoding="utf-8")
    )
    with h5py.File(args.static_input, "r") as handle:
        static_source_luminosity = np.asarray(
            handle["sources/photon_luminosity_s"], dtype=np.float64
        )
        static_source_cell_index = np.asarray(
            handle["sources/cell_index"], dtype=np.int64
        )
    with h5py.File(args.transport_control, "r") as handle:
        transport_attributes = dict(handle.attrs)
        transported_group_energy = np.asarray(
            handle["group_energy_ev"], dtype=np.float64
        )
        transport_absorbed_shape = handle[
            "diagnostics/cumulative_absorbed_photons_cm3"
        ].shape
        transport_photon_ledger_shape = handle[
            "diagnostics/photon_ledger/aggregate_residual_photons"
        ].shape

    hard_group = unique_interval_index(expected_intervals, (2000.0, 10000.0))
    soft_xray_group = unique_interval_index(expected_intervals, (500.0, 2000.0))
    species_thresholds = np.asarray(
        (H_I_FIT.threshold_ev, HE_I_FIT.threshold_ev, HE_II_FIT.threshold_ev),
        dtype=np.float64,
    )
    wholly_below_species_threshold = (
        expected_edges[1:][None, :] <= species_thresholds[:, None]
    )
    wholly_below_sed_support = expected_edges[1:] <= SED_MINIMUM_EV
    partially_supported_by_sed = (expected_edges[:-1] < SED_MINIMUM_EV) & (
        expected_edges[1:] > SED_MINIMUM_EV
    )
    fully_supported_by_sed = expected_edges[:-1] >= SED_MINIMUM_EV
    expected_group_status = np.full(
        expected_intervals.shape[0], "agn_sed_fully_supported", dtype=object
    )
    expected_group_status[wholly_below_sed_support] = (
        "agn_sed_below_support_zero_photons"
    )
    expected_group_status[partially_supported_by_sed] = (
        "agn_sed_partially_supported_10ev_to_upper"
    )
    metadata_group_status = np.asarray(
        [group["closure_status"] for group in metadata["groups"]], dtype=object
    )
    metadata_supported_intervals = [
        group.get("sed_supported_interval_ev") for group in metadata["groups"]
    ]
    hard_photon_rate = float(csv_totals[hard_group])
    soft_xray_photon_rate = float(csv_totals[soft_xray_group])
    hard_power_ev_s = hard_photon_rate * group_means[hard_group]
    soft_xray_power_ev_s = soft_xray_photon_rate * group_means[soft_xray_group]
    nonzero_photon_rate = float(csv_totals.sum(dtype=np.float64))
    supported_sed_power_ev_s = float(
        np.sum(csv_totals * group_means, dtype=np.float64)
    )
    total_bolometric_luminosity_erg_s = float(
        candidate_lbol.sum(dtype=np.float64)
    )
    hard_power_erg_s = hard_power_ev_s * EV_TO_ERG
    luminous_source_count = int(np.count_nonzero(candidate_lbol > 0.0))

    gas_input = Path(static_metadata["gas_input"])
    gas_metadata_value = static_metadata.get("gas_metadata")
    gas_metadata = None if gas_metadata_value is None else Path(gas_metadata_value)
    zoom_manifest = Path(static_metadata["zoom_manifest"])
    static_provenance = static_metadata["provenance"]
    external_asset_by_id = {
        asset["id"]: asset for asset in external_assets.get("assets", [])
    }
    expected_external_assets = {
        "agn9_static_rt_input": args.static_input,
        "agn9_short_transport_control": args.transport_control,
        "agn9_full_cfl_failed_probe": DEFAULT_FAILED_FULL_CFL_PROBE,
        "coeval_gas_input": gas_input,
    }
    external_asset_contract_closes = (
        external_assets.get("schema") == "snrt_agn_nine_group_external_assets_v1"
        and external_assets.get("publication_deposit", {}).get("status")
        == "pending_final_publication_archive"
        and all(
            asset_id in external_asset_by_id
            and Path(external_asset_by_id[asset_id]["path"]).resolve()
            == path.resolve()
            and external_asset_by_id[asset_id]["sha256"] == sha256(path)
            and external_asset_by_id[asset_id]["size_bytes"] == path.stat().st_size
            for asset_id, path in expected_external_assets.items()
        )
    )

    criteria = {
        "metadata_schema_v2": metadata.get("schema")
        == "snrt_agn_photon_ledger_v2",
        "p0_default_group_table": metadata.get("group_table_mode") == "p0_default",
        "exactly_nine_groups": metadata_edges.shape == (10,)
        and ledger.photon_luminosity_s.shape[1] == 9,
        "metadata_edges_exactly_match_config": np.array_equal(
            metadata_edges, expected_edges
        ),
        "metadata_intervals_exactly_match_config": np.array_equal(
            intervals, expected_intervals
        ),
        "group_indices_contiguous": np.array_equal(group_indices, np.arange(9)),
        "half_open_interval_policy_declared": metadata.get(
            "group_interval_convention"
        )
        == EXPECTED_INTERVAL_CONVENTION,
        "group_means_inside_intervals": bool(
            np.all(group_means >= expected_edges[:-1])
            and np.all(group_means <= expected_edges[1:])
        ),
        "csv_metadata_totals_close": bool(
            np.allclose(csv_totals, declared_totals, rtol=2.0e-15, atol=0.0)
            and np.allclose(csv_totals, grouped_totals, rtol=2.0e-15, atol=0.0)
        ),
        "candidate_rows_match_ledger_identity_kind_position": bool(
            int(metadata["source_count"]) == len(ledger.source_id)
            and np.array_equal(candidate_contract["source_id"], ledger.source_id)
            and np.array_equal(candidate_contract["source_kind"], ledger.source_kind)
            and np.array_equal(candidate_contract["position_code"], ledger.position_code)
            and len(np.unique(ledger.source_id)) == len(ledger.source_id)
        ),
        "luminous_source_accounting_closes": bool(
            luminous_source_count == 5
            and metadata.get("luminous_source_count") == luminous_source_count
            and metadata.get("total_bolometric_luminosity_erg_s")
            == total_bolometric_luminosity_erg_s
        ),
        "wholly_below_sed_support_groups_zero": bool(
            np.all(csv_totals[wholly_below_sed_support] == 0.0)
            and np.all(
                ledger.photon_luminosity_s[:, wholly_below_sed_support] == 0.0
            )
        ),
        "partial_and_full_sed_support_groups_positive": bool(
            np.all(csv_totals[partially_supported_by_sed | fully_supported_by_sed] > 0.0)
        ),
        "sed_support_status_and_intervals_explicit": bool(
            np.array_equal(metadata_group_status, expected_group_status)
            and all(
                supported is None
                if below
                else np.array_equal(
                    supported,
                    (max(float(low), SED_MINIMUM_EV), float(high)),
                )
                for low, high, below, supported in zip(
                    expected_edges[:-1],
                    expected_edges[1:],
                    wholly_below_sed_support,
                    metadata_supported_intervals,
                    strict=True,
                )
            )
        ),
        "two_to_ten_kev_group_present_and_positive": bool(
            np.array_equal(intervals[hard_group], (2000.0, 10000.0))
            and hard_photon_rate > 0.0
        ),
        "subthreshold_cross_sections_exactly_zero": bool(
            all(
                np.all(cross_sections[species, mask] == 0.0)
                for species, mask in enumerate(wholly_below_species_threshold)
            )
        ),
        "hard_group_microphysics_positive": bool(
            np.all(cross_sections[:, hard_group] > 0.0)
            and np.all(excess_energy[:, hard_group] > 0.0)
        ),
        "all_values_finite_nonnegative": bool(
            np.isfinite(ledger.photon_luminosity_s).all()
            and np.all(ledger.photon_luminosity_s >= 0.0)
            and np.isfinite(group_means).all()
            and np.isfinite(cross_sections).all()
            and np.all(cross_sections >= 0.0)
            and np.isfinite(excess_energy).all()
            and np.all(excess_energy >= 0.0)
        ),
        "metadata_generator_hash_current": provenance.get("generator_sha256")
        == sha256(GENERATOR),
        "metadata_candidates_hash_current": candidates.is_file()
        and provenance.get("candidates_sha256") == sha256(candidates),
        "metadata_group_config_hash_current": provenance.get(
            "group_edges_file_sha256"
        )
        == sha256(args.group_edges),
        "static_input_has_exact_nine_group_source_matrix": static_source_luminosity.shape
        == ledger.photon_luminosity_s.shape
        and static_source_luminosity.shape[1] == 9
        and np.array_equal(static_source_luminosity, ledger.photon_luminosity_s)
        and static_source_cell_index.shape == (len(ledger.source_id), 3),
        "static_rebind_metadata_closes": static_metadata.get("schema")
        == "snrt_static_source_rebind_v2"
        and static_metadata.get("group_count") == 9
        and np.array_equal(static_metadata.get("group_edges_ev"), expected_edges)
        and static_metadata.get("coeval_within_1e-12") is True
        and static_metadata["provenance"]["source_rebind_tool_sha256"]
        == sha256(SOURCE_REBIND_TOOL)
        and static_metadata["provenance"]["photon_ledger_sha256"]
        == sha256(args.ledger)
        and static_metadata["provenance"]["photon_metadata_sha256"]
        == sha256(args.metadata)
        and gas_input.is_file()
        and static_provenance["gas_input_sha256"] == sha256(gas_input)
        and gas_metadata is not None
        and gas_metadata.is_file()
        and static_provenance["gas_metadata_sha256"] == sha256(gas_metadata)
        and zoom_manifest.is_file()
        and static_provenance["zoom_manifest_sha256"] == sha256(zoom_manifest)
        and static_provenance["group_edges_sha256"] == sha256(args.group_edges),
        "external_hdf5_asset_manifest_closes": external_asset_contract_closes,
        "transport_control_uses_nine_groups": transported_group_energy.shape
        == (9,)
        and np.array_equal(transported_group_energy, group_means)
        and transport_absorbed_shape[0] == 9
        and transport_photon_ledger_shape == (9,),
        "transport_control_internal_gates_pass": bool(
            transport_attributes.get("validation_passed", False)
            and transport_attributes["photon_ledger_relative_error"] < 1.0e-12
            and transport_attributes["hydrogen_ledger_l1_relative_error"]
            < 1.0e-3
            and transport_attributes["helium_i_ledger_l1_relative_error"]
            < 1.0e-3
            and transport_attributes["helium_ii_ledger_l1_relative_error"]
            < 1.0e-3
            and transport_attributes["maximum_fixed_point_residual"] < 1.0e-4
            and transport_attributes[
                "photoelectron_energy_ledger_l1_relative_error"
            ]
            < 1.0e-12
            and transport_attributes["electron_root_bracket_failure_count"] == 0
        ),
    }

    payload = {
        "schema": "snrt_agn_nine_group_validation_v1",
        "passed": all(criteria.values()),
        "criteria": criteria,
        "configuration": {
            "group_edges_ev": expected_edges.tolist(),
            "number_of_groups": 9,
            "interval_convention": EXPECTED_INTERVAL_CONVENTION,
            "sed_support_minimum_ev": SED_MINIMUM_EV,
            "species_threshold_ev": species_thresholds.tolist(),
        },
        "hard_xray_diagnostics": {
            "energy_interval_ev": intervals[hard_group].tolist(),
            "photon_rate_s": hard_photon_rate,
            "photon_number_ratio_to_0p5_2kev": hard_photon_rate
            / soft_xray_photon_rate,
            "energy_power_ratio_to_0p5_2kev": hard_power_ev_s
            / soft_xray_power_ev_s,
            "fraction_of_total_photon_rate": hard_photon_rate
            / nonzero_photon_rate,
            "fraction_of_supported_sed_energy_power": hard_power_ev_s
            / supported_sed_power_ev_s,
            "fraction_of_candidate_bolometric_luminosity": hard_power_erg_s
            / total_bolometric_luminosity_erg_s,
            "emitted_power_erg_s": hard_power_erg_s,
            "photon_weighted_mean_energy_ev": float(group_means[hard_group]),
            "cross_sections_cm2": cross_sections[:, hard_group].tolist(),
            "photoelectron_excess_energy_ev": excess_energy[:, hard_group].tolist(),
        },
        "ledger": {
            "source_count": len(ledger.source_id),
            "luminous_source_count": luminous_source_count,
            "total_candidate_bolometric_luminosity_erg_s": total_bolometric_luminosity_erg_s,
            "group_photon_rate_total_s": csv_totals.tolist(),
            "total_photon_rate_s": nonzero_photon_rate,
            "total_supported_sed_energy_power_ev_s": supported_sed_power_ev_s,
        },
        "transport_control": {
            "duration_myr": float(transport_attributes["elapsed_time_s"])
            / (365.25 * 86400.0 * 1.0e6),
            "number_of_groups": len(transported_group_energy),
            "validation_passed": bool(transport_attributes["validation_passed"]),
            "photon_ledger_relative_error": float(
                transport_attributes["photon_ledger_relative_error"]
            ),
            "hydrogen_ledger_l1_relative_error": float(
                transport_attributes["hydrogen_ledger_l1_relative_error"]
            ),
            "maximum_fixed_point_residual": float(
                transport_attributes["maximum_fixed_point_residual"]
            ),
            "photoelectron_energy_ledger_l1_relative_error": float(
                transport_attributes[
                    "photoelectron_energy_ledger_l1_relative_error"
                ]
            ),
            "electron_root_bracket_failure_count": int(
                transport_attributes["electron_root_bracket_failure_count"]
            ),
        },
        "provenance": {
            "git_head": git_head(),
            "validator_sha256": sha256(Path(__file__).resolve()),
            "generator_sha256": sha256(GENERATOR),
            "group_edges_sha256": sha256(args.group_edges),
            "ledger_sha256": sha256(args.ledger),
            "metadata_sha256": sha256(args.metadata),
            "candidates_sha256": sha256(candidates),
            "primordial_sha256": sha256(ROOT / "snrt_core" / "primordial.py"),
            "source_ledger_sha256": sha256(ROOT / "snrt_core" / "source_ledger.py"),
            "source_rebind_tool_sha256": sha256(SOURCE_REBIND_TOOL),
            "p4_runner_sha256": sha256(P4_RUNNER),
            "static_input_sha256": sha256(args.static_input),
            "static_metadata_sha256": sha256(args.static_metadata),
            "transport_control_sha256": sha256(args.transport_control),
            "external_asset_manifest_sha256": sha256(
                args.external_asset_manifest
            ),
        },
        "scope": (
            "AGN source-ledger group topology, SED integration, and hard-X-ray "
            "wiring; not an obscuration calibration or a transported-field promotion"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"AGN_NINE_GROUP_{'PASS' if payload['passed'] else 'FAIL'} "
        f"hard_q={hard_photon_rate:.6g} "
        f"hard_to_soft_q={payload['hard_xray_diagnostics']['photon_number_ratio_to_0p5_2kev']:.6g} "
        f"hard_supported_sed_fraction={payload['hard_xray_diagnostics']['fraction_of_supported_sed_energy_power']:.6g} "
        f"hard_bolometric_fraction={payload['hard_xray_diagnostics']['fraction_of_candidate_bolometric_luminosity']:.6g} "
        f"output={args.output}"
    )
    if not payload["passed"]:
        failed = ", ".join(name for name, passed in criteria.items() if not passed)
        raise RuntimeError(f"AGN nine-group validation failed: {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
