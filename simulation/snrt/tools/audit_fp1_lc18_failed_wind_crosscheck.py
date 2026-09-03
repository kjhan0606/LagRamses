#!/usr/bin/env python3
"""Cross-check the unresolved LC18 failed-model wind release against CDS."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable
from zipfile import BadZipFile, ZipFile

from adapt_g2_candidate_sources import (
    LIMONGI_ID,
    SourceAdapterError,
    adapt_candidate,
)
from audit_g2_boccioli_roberti2026_candidate import (
    BoccioliRobertiAuditError,
    _summary_rows,
    _yield_table,
    audit_boccioli_roberti2026_candidate,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
REPOSITORY_ROOT = SNRT_ROOT.parents[1]
DEFAULT_ROOT = REPOSITORY_ROOT / "external" / "g2_candidates"
DEFAULT_BOCCIOLI_CONTRACT = (
    SNRT_ROOT / "config" / "g2_boccioli_roberti2026_candidate_contract_v1.json"
)
DEFAULT_PHASE_CONTRACT = (
    SNRT_ROOT / "config" / "g2_limongi_phase_mass_history_contract_v1.json"
)
DEFAULT_PHYSICAL_PACKAGE_CONTRACT = (
    SNRT_ROOT / "config" / "fp1_physical_package_admission_contract_v1.json"
)
BOCCIOLI_ADMISSION_ID = "boccioli_roberti2026_lc18"
EXPECTED_ADMISSION_BLOCKERS = [
    "failed_model_wind_summary_table_anomaly_requires_author_or_corrected_release",
    "age_resolved_wind_missing",
    "per_node_injected_energy_mapping_missing",
    "canonical_momentum_and_deposition_missing",
]


class Lc18FailedWindCrosscheckError(ValueError):
    """The staged LC18/BR26 cross-check evidence is inconsistent."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise Lc18FailedWindCrosscheckError(
            f"cannot read {label} {path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise Lc18FailedWindCrosscheckError(f"{label} must be a JSON object")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise Lc18FailedWindCrosscheckError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def _coordinate(record: dict[str, Any]) -> tuple[int, int, float]:
    source = record["source_coordinate"]
    return (
        int(source["rotation_velocity_km_s"]),
        int(source["metallicity_feh"]),
        float(source["initial_mass_msun"]),
    )


def _build_phase_histories(
    records: list[dict[str, Any]], phase_order: list[str]
) -> tuple[dict[tuple[int, int, float], dict[str, Any]], dict[str, Any]]:
    phase_rank = {phase: index for index, phase in enumerate(phase_order)}
    grouped: dict[tuple[int, int, float], list[dict[str, Any]]] = defaultdict(list)
    collapsed_duplicates: list[dict[str, Any]] = []
    for record in records:
        source = record["source_coordinate"]
        occurrence = int(source["phase_occurrence"])
        if occurrence > 1:
            collapsed_duplicates.append(
                {
                    "coordinate": {
                        "rotation_velocity_km_s": int(
                            source["rotation_velocity_km_s"]
                        ),
                        "metallicity_feh": int(source["metallicity_feh"]),
                        "initial_mass_msun": float(source["initial_mass_msun"]),
                    },
                    "phase": source["phase"],
                    "occurrence": occurrence,
                    "source_line": record["source_line"],
                }
            )
            continue
        grouped[_coordinate(record)].append(record)
    if len(grouped) != 108:
        raise Lc18FailedWindCrosscheckError(
            f"expected 108 CDS phase-history models, found {len(grouped)}"
        )

    age_violations: list[dict[str, Any]] = []
    mass_violations: list[dict[str, Any]] = []
    missing_terminal_phase: list[dict[str, Any]] = []
    histories: dict[tuple[int, int, float], dict[str, Any]] = {}
    unique_phase_counts: list[int] = []
    for coordinate, model_records in grouped.items():
        try:
            ordered = sorted(
                model_records,
                key=lambda value: phase_rank[value["source_coordinate"]["phase"]],
            )
        except KeyError as exc:
            raise Lc18FailedWindCrosscheckError(
                f"unknown CDS phase {exc.args[0]!r}"
            ) from exc
        phases = [record["source_coordinate"]["phase"] for record in ordered]
        if len(phases) != len(set(phases)):
            raise Lc18FailedWindCrosscheckError(
                f"non-collapsed duplicate phase at {coordinate}"
            )
        if not phases or phases[-1] != "PSN":
            missing_terminal_phase.append(
                {
                    "rotation_velocity_km_s": coordinate[0],
                    "metallicity_feh": coordinate[1],
                    "initial_mass_msun": coordinate[2],
                }
            )

        cumulative_age = 0.0
        previous_age = 0.0
        previous_mass = coordinate[2]
        nodes: list[dict[str, Any]] = []
        for record in ordered:
            source = record["source_coordinate"]
            duration = float(record["phase_duration_yr"])
            total_mass = float(record["total_mass_msun"])
            cumulative_age += duration
            node = {
                "phase": source["phase"],
                "source_line": record["source_line"],
                "phase_duration_yr": duration,
                "cumulative_age_yr": cumulative_age,
                "total_mass_msun": total_mass,
                "cumulative_wind_mass_msun": coordinate[2] - total_mass,
            }
            nodes.append(node)
            if (
                not math.isfinite(duration)
                or duration <= 0.0
                or not math.isfinite(cumulative_age)
                or cumulative_age <= previous_age
            ):
                age_violations.append(
                    {"coordinate": list(coordinate), "node": dict(node)}
                )
            if (
                not math.isfinite(total_mass)
                or total_mass < 0.0
                or total_mass > previous_mass + 1.0e-12
            ):
                mass_violations.append(
                    {
                        "coordinate": list(coordinate),
                        "previous_mass_msun": previous_mass,
                        "node": dict(node),
                    }
                )
            previous_age = cumulative_age
            previous_mass = total_mass

        terminal = nodes[-1]
        histories[coordinate] = {
            "unique_phase_count": len(nodes),
            "terminal_phase": terminal["phase"],
            "terminal_age_yr": terminal["cumulative_age_yr"],
            "terminal_total_mass_msun": terminal["total_mass_msun"],
            "terminal_cumulative_wind_mass_msun": terminal[
                "cumulative_wind_mass_msun"
            ],
            "nodes": nodes,
        }
        unique_phase_counts.append(len(nodes))

    diagnostics = {
        "model_count": len(histories),
        "unique_phase_row_count": sum(unique_phase_counts),
        "minimum_unique_phase_count_per_model": min(unique_phase_counts),
        "maximum_unique_phase_count_per_model": max(unique_phase_counts),
        "collapsed_duplicate_row_count": len(collapsed_duplicates),
        "collapsed_duplicate_rows": collapsed_duplicates,
        "strictly_increasing_cumulative_age_violation_count": len(age_violations),
        "strictly_increasing_cumulative_age_violations": age_violations,
        "nonincreasing_total_mass_violation_count": len(mass_violations),
        "nonincreasing_total_mass_violations": mass_violations,
        "missing_psn_terminal_phase_count": len(missing_terminal_phase),
        "missing_psn_terminal_phase": missing_terminal_phase,
    }
    return histories, diagnostics


def _structure_map(
    records: list[dict[str, Any]], expected_coordinates: set[tuple[int, int, float]]
) -> tuple[dict[tuple[int, int, float], dict[str, Any]], list[list[float]]]:
    structures: dict[tuple[int, int, float], dict[str, Any]] = {}
    for record in records:
        coordinate = _coordinate(record)
        if coordinate in structures:
            raise Lc18FailedWindCrosscheckError(
                f"duplicate CDS table7 coordinate {coordinate}"
            )
        if coordinate not in expected_coordinates:
            raise Lc18FailedWindCrosscheckError(
                f"CDS table7 coordinate has no LC18 summary model: {coordinate}"
            )
        structures[coordinate] = {
            key: value
            for key, value in record.items()
            if key not in {"source_coordinate"}
        }
    missing = sorted(expected_coordinates - set(structures))
    return structures, [list(coordinate) for coordinate in missing]


def _release_rows(
    archive_path: Path, grid: dict[str, Any]
) -> tuple[dict[tuple[int, int, float], dict[str, Any]], dict[str, Any]]:
    try:
        archive = ZipFile(archive_path)
    except (OSError, BadZipFile) as exc:
        raise Lc18FailedWindCrosscheckError(
            f"cannot open LC18 archive {archive_path}: {exc}"
        ) from exc
    rows: dict[tuple[int, int, float], dict[str, Any]] = {}
    try:
        summary = _summary_rows(
            archive, "LC18/Summary_table_LC18.txt", lc18=True
        )
        metallicity_feh = grid["metallicity_feh"]
        for label in grid["metallicity_labels"]:
            for rotation in grid["rotation_km_s"]:
                member = (
                    f"LC18/Elements/Met_{label}/"
                    f"Yields_Rot_{rotation:03d}_Eles_LC18_Wind.txt"
                )
                table = _yield_table(archive, member)
                by_mass = {
                    int(row["mass_msun"]): row
                    for row in summary
                    if row["metallicity_label"] == label
                    and row["rotation_km_s"] == rotation
                }
                if set(by_mass) != set(table["mass_msun"]):
                    raise Lc18FailedWindCrosscheckError(
                        f"summary/Wind-table mass axis mismatch: {label}:{rotation}"
                    )
                for column, mass in enumerate(table["mass_msun"]):
                    summary_row = by_mass[mass]
                    coordinate = (
                        int(rotation),
                        int(metallicity_feh[label]),
                        float(mass),
                    )
                    if coordinate in rows:
                        raise Lc18FailedWindCrosscheckError(
                            f"duplicate LC18 release coordinate {coordinate}"
                        )
                    wind_sum = math.fsum(
                        table["values"][species][column]
                        for species in table["species"]
                    )
                    rows[coordinate] = {
                        "metallicity_label": label,
                        "archive_member_path": member,
                        "summary": dict(summary_row),
                        "wind_table_element_sum_msun": wind_sum,
                    }
    except BoccioliRobertiAuditError as exc:
        raise Lc18FailedWindCrosscheckError(str(exc)) from exc
    finally:
        archive.close()
    return rows, {
        "summary_row_count": len(summary),
        "joined_release_row_count": len(rows),
    }


def audit_lc18_failed_wind_crosscheck(
    *,
    root: Path = DEFAULT_ROOT,
    boccioli_contract_path: Path = DEFAULT_BOCCIOLI_CONTRACT,
    phase_contract_path: Path = DEFAULT_PHASE_CONTRACT,
    physical_package_contract_path: Path = DEFAULT_PHYSICAL_PACKAGE_CONTRACT,
) -> dict[str, Any]:
    root = Path(root).resolve()
    boccioli_contract_path = Path(boccioli_contract_path).resolve()
    phase_contract_path = Path(phase_contract_path).resolve()
    physical_package_contract_path = Path(physical_package_contract_path).resolve()

    boccioli_audit = audit_boccioli_roberti2026_candidate(
        root=root, contract_path=boccioli_contract_path
    )
    limongi = adapt_candidate(LIMONGI_ID, root=root, include_records=True)
    boccioli_contract = _read_json(
        boccioli_contract_path, "Boccioli-Roberti candidate contract"
    )
    phase_contract = _read_json(phase_contract_path, "Limongi phase contract")
    physical_contract = _read_json(
        physical_package_contract_path, "physical-package contract"
    )

    if (
        phase_contract.get("schema")
        != "snrt-g2-limongi-phase-mass-history-contract"
        or phase_contract.get("schema_version") != 1
        or phase_contract.get("limitations", {}).get(
            "canonical_conversion_allowed"
        )
        is not False
        or phase_contract.get("limitations", {}).get("runtime_deposition_allowed")
        is not False
    ):
        raise Lc18FailedWindCrosscheckError(
            "Limongi phase-history review contract is not fail-closed"
        )

    if boccioli_audit["canonical_rows_emitted"] != 0:
        raise Lc18FailedWindCrosscheckError(
            "Boccioli-Roberti review unexpectedly emitted canonical rows"
        )
    if limongi["canonical_rows_emitted"] != 0:
        raise Lc18FailedWindCrosscheckError(
            "Limongi review unexpectedly emitted canonical rows"
        )
    admission = physical_contract.get("candidate_qualification", {}).get(
        BOCCIOLI_ADMISSION_ID
    )
    approval = physical_contract.get("approval")
    if (
        not isinstance(admission, dict)
        or admission.get("hard_blockers") != EXPECTED_ADMISSION_BLOCKERS
        or admission.get("production_qualified") is not False
        or physical_contract.get("physical_node_inventory") != []
        or not isinstance(approval, dict)
        or any(
            approval.get(name) is not False
            for name in (
                "physical_package_selected",
                "canonical_conversion_allowed",
                "runtime_deposition_allowed",
                "production_ready",
                "publication_ready",
            )
        )
    ):
        raise Lc18FailedWindCrosscheckError(
            "Boccioli-Roberti admission blockers or qualification drifted"
        )

    release_rows, release_diagnostics = _release_rows(
        root / "boccioli_roberti2026_ccsn" / "LC18.zip",
        boccioli_contract["grids"]["LC18"],
    )
    components = limongi["source_components"]
    properties = components["evolutionary_properties"]
    if properties["all_duplicate_rows_physically_identical"] is not True:
        raise Lc18FailedWindCrosscheckError(
            "non-identical CDS phase duplicates cannot be collapsed"
        )
    histories, phase_diagnostics = _build_phase_histories(
        properties["records"],
        phase_contract["phase_order"],
    )
    release_coordinates = set(release_rows)
    history_coordinates = set(histories)
    unmatched_release = sorted(release_coordinates - history_coordinates)
    unmatched_cds = sorted(history_coordinates - release_coordinates)
    if unmatched_release or unmatched_cds or len(release_coordinates) != 108:
        raise Lc18FailedWindCrosscheckError(
            "LC18 release/CDS model-coordinate join is not one-to-one"
        )
    structures, missing_structures = _structure_map(
        components["presupernova_properties"]["records"], release_coordinates
    )

    rows: list[dict[str, Any]] = []
    for coordinate in sorted(release_coordinates):
        release = release_rows[coordinate]
        summary = release["summary"]
        history = histories[coordinate]
        summary_wind = float(summary["wind_mass_msun"])
        release_wind = float(release["wind_table_element_sum_msun"])
        cds_wind = float(history["terminal_cumulative_wind_mass_msun"])
        failed = not bool(summary["exploded"])
        rows.append(
            {
                "coordinate": {
                    "rotation_velocity_km_s": coordinate[0],
                    "metallicity_feh": coordinate[1],
                    "metallicity_label": release["metallicity_label"],
                    "initial_mass_msun": coordinate[2],
                },
                "archive_member_path": release["archive_member_path"],
                "exploded": not failed,
                "summary_wind_mass_msun": summary_wind,
                "release_wind_table_element_sum_msun": release_wind,
                "cds_terminal_cumulative_wind_mass_msun": cds_wind,
                "cds_terminal_age_yr": history["terminal_age_yr"],
                "differences_msun": {
                    "summary_minus_release_wind_table": summary_wind
                    - release_wind,
                    "summary_minus_cds_terminal_wind": summary_wind - cds_wind,
                    "release_wind_table_minus_cds_terminal_wind": release_wind
                    - cds_wind,
                },
                "cds_phase_history": history,
                "cds_presupernova_structure": structures.get(coordinate),
                "resolution": (
                    "unresolved_failed_wind_anomaly"
                    if failed
                    else "successful_model_release_control_review_only"
                ),
            }
        )

    successful = [row for row in rows if row["exploded"]]
    failed = [row for row in rows if not row["exploded"]]
    successful_internal_residuals = [
        abs(row["differences_msun"]["summary_minus_release_wind_table"])
        for row in successful
    ]
    successful_cds_residuals = [
        abs(row["differences_msun"]["summary_minus_cds_terminal_wind"])
        for row in successful
    ]
    all_cds_residuals = [
        abs(row["differences_msun"]["summary_minus_cds_terminal_wind"])
        for row in rows
    ]
    cds_half_bin = 0.5 * float(
        phase_contract["limitations"]["phase_endpoint_total_mass_precision_msun"]
    )
    successful_control = {
        "model_count": len(successful),
        "summary_wind_positive_count": sum(
            row["summary_wind_mass_msun"] > 0.0 for row in successful
        ),
        "release_wind_table_nonzero_count": sum(
            row["release_wind_table_element_sum_msun"] > 0.0
            for row in successful
        ),
        "maximum_absolute_summary_minus_release_wind_table_msun": max(
            successful_internal_residuals
        ),
        "maximum_absolute_summary_minus_cds_terminal_wind_msun": max(
            successful_cds_residuals
        ),
        "summary_minus_cds_above_nominal_half_bin_count": sum(
            value > cds_half_bin for value in successful_cds_residuals
        ),
        "interpretation": (
            "Successful models verify that the BR26 release normally carries "
            "nonzero Wind tables. Their CDS endpoint differences are measured "
            "cross-source discrepancies, not a promotion tolerance."
        ),
    }
    failed_anomaly = {
        "model_count": len(failed),
        "summary_wind_positive_count": sum(
            row["summary_wind_mass_msun"] > 0.0 for row in failed
        ),
        "release_wind_table_exact_zero_count": sum(
            row["release_wind_table_element_sum_msun"] == 0.0 for row in failed
        ),
        "cds_terminal_wind_positive_count": sum(
            row["cds_terminal_cumulative_wind_mass_msun"] > 0.0 for row in failed
        ),
        "cds_terminal_wind_zero_count": sum(
            row["cds_terminal_cumulative_wind_mass_msun"] == 0.0 for row in failed
        ),
        "unresolved_count": sum(
            row["resolution"] == "unresolved_failed_wind_anomaly"
            for row in failed
        ),
    }
    cross_source = {
        "comparison": "BR26 summary wind minus LC18 CDS initial-minus-PSN mass",
        "nominal_cds_total_mass_half_bin_msun": cds_half_bin,
        "model_count": len(rows),
        "above_nominal_half_bin_count": sum(
            value > cds_half_bin for value in all_cds_residuals
        ),
        "maximum_absolute_difference_msun": max(all_cds_residuals),
        "agreement_required_for_this_review": False,
        "interpretation": (
            "The two staged sources do not agree at nominal table5 precision. "
            "The discrepancy is retained as an author/source question and is "
            "not silently reconciled."
        ),
    }

    expected_failed = boccioli_audit["quality_findings"][
        "lc18_failed_models_with_reported_wind_but_zero_wind_table_count"
    ]
    if (
        len(successful) != 52
        or len(failed) != expected_failed
        or failed_anomaly["summary_wind_positive_count"] != expected_failed
        or failed_anomaly["release_wind_table_exact_zero_count"] != expected_failed
    ):
        raise Lc18FailedWindCrosscheckError(
            "LC18 successful-control or failed-wind anomaly counts drifted"
        )
    if len(structures) != 96 or len(missing_structures) != 12:
        raise Lc18FailedWindCrosscheckError(
            "CDS table7 structure coverage drifted from 96 present / 12 absent"
        )

    source_terms = limongi["source_semantics_evidence"][
        "source_use_terms_evidence"
    ]["candidate_record"]
    return {
        "schema": "snrt-fp1-lc18-failed-wind-crosscheck",
        "schema_version": 1,
        "gate": "F-P1H-E-review",
        "status": "failed_wind_anomaly_independently_crosschecked_unresolved",
        "production_ready": False,
        "publication_ready": False,
        "canonical_conversion_allowed": False,
        "runtime_deposition_allowed": False,
        "canonical_rows_emitted": 0,
        "physical_nodes_emitted": 0,
        "inquiry_sent": False,
        "source_identity": {
            "boccioli_roberti": boccioli_audit["source_identity"],
            "limongi_chieffi_cds": limongi["verified_acquisition"],
            "limongi_chieffi_cds_redistribution_status": source_terms[
                "redistribution_status"
            ],
            "limongi_chieffi_cds_production_license_status": source_terms[
                "production_license_status"
            ],
            "review_use_only": True,
        },
        "join": {
            **release_diagnostics,
            "cds_model_count": len(histories),
            "joined_model_count": len(rows),
            "unmatched_release_coordinates": [list(value) for value in unmatched_release],
            "unmatched_cds_coordinates": [list(value) for value in unmatched_cds],
            "one_to_one": True,
        },
        "successful_release_control": successful_control,
        "failed_wind_anomaly": failed_anomaly,
        "quality_findings": {
            "lc18_readme_consistency_pass": boccioli_audit["quality_findings"][
                "lc18_readme_consistency_pass"
            ],
            "failed_wind_anomaly_resolved": False,
            "cross_source_difference_silently_reconciled": False,
        },
        "cross_source_wind_comparison": cross_source,
        "phase_history": phase_diagnostics,
        "presupernova_structure": {
            "available_model_count": len(structures),
            "explicit_null_model_count": len(missing_structures),
            "missing_coordinates": missing_structures,
            "binding_energy_is_not_injected_explosion_energy": True,
        },
        "rows": rows,
        "admission": {
            "candidate_id": BOCCIOLI_ADMISSION_ID,
            "hard_blockers_unchanged": True,
            "hard_blockers": list(admission["hard_blockers"]),
            "production_qualified": False,
            "physical_package_contract_path": str(physical_package_contract_path),
            "physical_package_contract_sha256": _sha256(
                physical_package_contract_path
            ),
        },
        "scientific_limitations": [
            "failed-model BR26 Wind tables remain unresolved and are not reconstructed",
            "CDS table5 supplies phase-resolved total mass but not phase-resolved isotopic wind composition",
            "CDS table7 is absent for 12 models and missing values remain null",
            "CDS table7 binding energy is not an injected explosion energy",
            "machine-readable per-model injected energy and canonical momentum are absent",
            "CDS catalogue redistribution has no explicit identified license and this artifact is review-only",
        ],
        "tool_sha256": _sha256(TOOL_PATH),
        "boccioli_contract_sha256": _sha256(boccioli_contract_path),
        "phase_contract_sha256": _sha256(phase_contract_path),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument(
        "--boccioli-contract", type=Path, default=DEFAULT_BOCCIOLI_CONTRACT
    )
    parser.add_argument("--phase-contract", type=Path, default=DEFAULT_PHASE_CONTRACT)
    parser.add_argument(
        "--physical-package-contract",
        type=Path,
        default=DEFAULT_PHYSICAL_PACKAGE_CONTRACT,
    )
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_lc18_failed_wind_crosscheck(
            root=args.root,
            boccioli_contract_path=args.boccioli_contract,
            phase_contract_path=args.phase_contract,
            physical_package_contract_path=args.physical_package_contract,
        )
    except (
        Lc18FailedWindCrosscheckError,
        BoccioliRobertiAuditError,
        SourceAdapterError,
    ) as exc:
        print(
            json.dumps({"status": "error", "error": str(exc)}, indent=2),
            file=sys.stderr,
        )
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
