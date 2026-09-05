#!/usr/bin/env python3
"""Cross-check the unresolved LC18 failed-model wind release against CDS."""

from __future__ import annotations

import argparse
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
from fp1_publication_rights import (
    PublicationRightsError,
    evaluate_derived_artifact_publication,
)
from fp1_limongi_phase_history import (
    PhaseHistoryInvariantError,
    build_phase_histories,
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
PHASE_HISTORY_TOOL_PATH = Path(__file__).with_name("fp1_limongi_phase_history.py")
BOCCIOLI_ADMISSION_ID = "boccioli_roberti2026_lc18"
EXPECTED_ADMISSION_BLOCKERS = [
    "failed_model_wind_summary_table_anomaly_requires_author_or_corrected_release",
    "age_resolved_wind_missing",
    "per_node_injected_energy_mapping_missing",
    "canonical_momentum_and_deposition_missing",
]


class Lc18FailedWindCrosscheckError(ValueError):
    """The staged LC18/BR26 cross-check evidence is inconsistent."""

    def __init__(self, message: str, *, diagnostics: dict[str, Any] | None = None):
        super().__init__(message)
        self.diagnostics = diagnostics or {}


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
    try:
        canonical_histories, diagnostics = build_phase_histories(
            records, phase_order
        )
    except PhaseHistoryInvariantError as exc:
        raise Lc18FailedWindCrosscheckError(
            str(exc), diagnostics=exc.diagnostics
        ) from exc
    histories = {
        coordinate: {
            "unique_phase_count": canonical["unique_phase_count"],
            "terminal_phase": canonical["terminal_phase"],
            "terminal_age_yr": canonical["terminal_age_yr"],
            "terminal_total_mass_msun": canonical["terminal_total_mass_msun"],
            "terminal_cumulative_wind_mass_msun": canonical[
                "terminal_cumulative_wind_mass_msun"
            ],
            "nodes": canonical["nodes"],
        }
        for coordinate, canonical in canonical_histories.items()
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


def _residual_statistics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    signed = [
        float(row["differences_msun"]["summary_minus_cds_terminal_wind"])
        for row in rows
    ]
    relative: list[float] = []
    relative_null_count = 0
    for row, difference in zip(rows, signed):
        denominator = float(row["summary_wind_mass_msun"])
        if denominator == 0.0:
            relative_null_count += 1
        else:
            relative.append(difference / denominator)
    return {
        "model_count": len(rows),
        "comparison": "summary_wind_mass_msun minus cds_terminal_cumulative_wind_msun",
        "signed_difference_unit": "Msun",
        "positive_signed_difference_count": sum(value > 0.0 for value in signed),
        "negative_signed_difference_count": sum(value < 0.0 for value in signed),
        "zero_signed_difference_count": sum(value == 0.0 for value in signed),
        "minimum_signed_difference_msun": min(signed) if signed else None,
        "maximum_signed_difference_msun": max(signed) if signed else None,
        "maximum_absolute_difference_msun": max(
            (abs(value) for value in signed), default=None
        ),
        "relative_difference_definition": (
            "(summary_wind_mass_msun - cds_terminal_cumulative_wind_msun) / "
            "summary_wind_mass_msun"
        ),
        "relative_difference_denominator": "summary_wind_mass_msun",
        "relative_null_count_zero_denominator": relative_null_count,
        "minimum_relative_difference": min(relative) if relative else None,
        "maximum_relative_difference": max(relative) if relative else None,
    }


def _parsed_cds_zero_metadata(precision_msun: float) -> dict[str, Any]:
    """Describe parsed Table 5 zeros without inferring a physical zero wind."""

    if not math.isfinite(precision_msun) or precision_msun <= 0.0:
        raise Lc18FailedWindCrosscheckError(
            "CDS phase-endpoint mass precision must be a positive finite value"
        )
    return {
        "cds_phase_endpoint_total_mass_precision_msun": precision_msun,
        "cds_phase_endpoint_total_mass_half_bin_msun": 0.5 * precision_msun,
        "physical_zero_inferred": False,
        "parsed_exact_zero_definition": (
            "cds_terminal_cumulative_wind_mass_msun == 0.0 after parsing "
            "M_initial - M_total(PSN) from CDS table5 at the declared source "
            "endpoint precision"
        ),
        "parsed_exact_zero_interpretation": (
            "An exact parsed zero is compatible with physical mass loss below "
            "the source table half-bin; it is not evidence that physical wind "
            "is zero. The pipeline does not round the wind values."
        ),
    }


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
    hard_blockers_unchanged = isinstance(admission, dict) and admission.get(
        "hard_blockers"
    ) == EXPECTED_ADMISSION_BLOCKERS
    if (
        not isinstance(admission, dict)
        or not hard_blockers_unchanged
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
    one_to_one = (
        not unmatched_release
        and not unmatched_cds
        and len(release_coordinates) == 108
        and len(history_coordinates) == 108
    )
    if not one_to_one:
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
    cds_zero_metadata = _parsed_cds_zero_metadata(
        float(phase_contract["limitations"]["phase_endpoint_total_mass_precision_msun"])
    )
    cds_half_bin = cds_zero_metadata["cds_phase_endpoint_total_mass_half_bin_msun"]
    successful_control = {
        "model_count": len(successful),
        "summary_wind_positive_count": sum(
            row["summary_wind_mass_msun"] > 0.0 for row in successful
        ),
        "release_wind_table_nonzero_count": sum(
            row["release_wind_table_element_sum_msun"] > 0.0
            for row in successful
        ),
        "cds_terminal_wind_parsed_positive_count": sum(
            row["cds_terminal_cumulative_wind_mass_msun"] > 0.0
            for row in successful
        ),
        "cds_terminal_wind_parsed_exact_zero_count": sum(
            row["cds_terminal_cumulative_wind_mass_msun"] == 0.0
            for row in successful
        ),
        "maximum_absolute_summary_minus_release_wind_table_msun": max(
            successful_internal_residuals
        ),
        "maximum_absolute_summary_minus_cds_terminal_wind_msun": max(
            successful_cds_residuals
        ),
        "summary_minus_cds_above_phase_endpoint_half_bin_count": sum(
            value > cds_half_bin for value in successful_cds_residuals
        ),
        **cds_zero_metadata,
        "interpretation": (
            "Successful models verify that the BR26 release normally carries "
            "nonzero Wind tables. Their CDS endpoint differences are measured "
            "cross-source discrepancies, not a promotion tolerance. The four "
            "parsed exact-zero CDS endpoints are successful controls and do not "
            "establish physical zero wind."
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
        "cds_terminal_wind_parsed_positive_count": sum(
            row["cds_terminal_cumulative_wind_mass_msun"] > 0.0 for row in failed
        ),
        "cds_terminal_wind_parsed_exact_zero_count": sum(
            row["cds_terminal_cumulative_wind_mass_msun"] == 0.0 for row in failed
        ),
        "unresolved_count": sum(
            row["resolution"] == "unresolved_failed_wind_anomaly"
            for row in failed
        ),
        **cds_zero_metadata,
        "interpretation": (
            "The three parsed exact-zero CDS endpoints are inside the failed "
            "BR26 zero-Wind release anomaly; they do not define or resolve "
            "that anomaly, and parsed zero is not inferred physical zero wind."
        ),
    }
    cross_source = {
        "comparison": "BR26 summary wind minus LC18 CDS initial-minus-PSN mass",
        **cds_zero_metadata,
        "model_count": len(rows),
        "cds_terminal_wind_parsed_positive_count": sum(
            row["cds_terminal_cumulative_wind_mass_msun"] > 0.0
            for row in rows
        ),
        "cds_terminal_wind_parsed_exact_zero_count": sum(
            row["cds_terminal_cumulative_wind_mass_msun"] == 0.0
            for row in rows
        ),
        "cds_terminal_wind_parsed_exact_zero_count_by_outcome": {
            "successful": sum(
                row["exploded"]
                and row["cds_terminal_cumulative_wind_mass_msun"] == 0.0
                for row in rows
            ),
            "failed": sum(
                not row["exploded"]
                and row["cds_terminal_cumulative_wind_mass_msun"] == 0.0
                for row in rows
            ),
        },
        "above_cds_phase_endpoint_half_bin_count": sum(
            value > cds_half_bin for value in all_cds_residuals
        ),
        "maximum_absolute_difference_msun": max(all_cds_residuals),
        "agreement_required_for_this_review": False,
        "residual_statistics": {
            "all_models": _residual_statistics(rows),
            "successful_release_control": _residual_statistics(successful),
            "failed_wind_anomaly": _residual_statistics(failed),
        },
        "interpretation": (
            "The two staged sources do not agree at the declared CDS table5 "
            "endpoint precision. The discrepancy is retained as an author/"
            "source question and is not silently reconciled. Parsed exact-zero "
            "endpoints are source-precision observations, not physical-zero "
            "inferences."
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
        or successful_control["cds_terminal_wind_parsed_exact_zero_count"] != 4
        or successful_control["cds_terminal_wind_parsed_positive_count"] != 48
        or failed_anomaly["cds_terminal_wind_parsed_exact_zero_count"] != 3
        or failed_anomaly["cds_terminal_wind_parsed_positive_count"] != 53
        or cross_source["cds_terminal_wind_parsed_exact_zero_count"] != 7
        or cross_source["cds_terminal_wind_parsed_positive_count"] != 101
        or cross_source["cds_terminal_wind_parsed_exact_zero_count_by_outcome"]
        != {"successful": 4, "failed": 3}
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
    source_terms_evidence = limongi["source_semantics_evidence"][
        "source_use_terms_evidence"
    ]
    terms_path_value = source_terms_evidence.get("path")
    if not isinstance(terms_path_value, str) or not terms_path_value:
        raise Lc18FailedWindCrosscheckError(
            "Limongi source-use-terms evidence path is missing"
        )
    try:
        publication_gate = evaluate_derived_artifact_publication(
            candidate_id="limongi_chieffi_2018_cds",
            terms_path=Path(terms_path_value),
            approval_record={},
            review_use_only=True,
            derived_artifact_kind="fp1_lc18_failed_wind_crosscheck",
        )
    except PublicationRightsError as exc:
        raise Lc18FailedWindCrosscheckError(
            f"Limongi derived-artifact publication gate failed: {exc}"
        ) from exc
    return {
        "schema": "snrt-fp1-lc18-failed-wind-crosscheck",
        "schema_version": 1,
        "gate": "F-P1H-E-review",
        "status": "failed_wind_anomaly_independently_crosschecked_unresolved",
        "production_ready": False,
        "publication_ready": publication_gate["publication_ready"],
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
            "limongi_chieffi_cds_rights": {
                "redistribution_status": source_terms["redistribution_status"],
                "production_license_status": source_terms["production_license_status"],
                "authoritative_for_verdict": False,
            },
            "review_use_only": publication_gate["review_use_only"],
        },
        "join": {
            **release_diagnostics,
            "cds_model_count": len(histories),
            "joined_model_count": len(rows),
            "unmatched_release_coordinates": [list(value) for value in unmatched_release],
            "unmatched_cds_coordinates": [list(value) for value in unmatched_cds],
            "one_to_one": one_to_one,
        },
        "successful_release_control": successful_control,
        "failed_wind_anomaly": failed_anomaly,
        "quality_findings": {
            "lc18_readme_consistency_pass": boccioli_audit["quality_findings"][
                "lc18_readme_consistency_pass"
            ],
            "failed_wind_anomaly_resolved": False,
            "cross_source_difference_silently_reconciled": False,
            "derived_artifact_publication_gate_pass": publication_gate["allowed"],
        },
        "publication_gate": publication_gate,
        "cross_source_wind_comparison": cross_source,
        "phase_history": {
            **phase_diagnostics,
            "phase_order": list(phase_contract["phase_order"]),
            "phase_order_provenance": {
                "source": "g2_limongi_phase_mass_history_contract_v1",
                "source_attested_for_intermediate_burning_order": False,
                "interpretation": (
                    "MS/H/He/PSN labels are source-attested; the C/Ne/O/Si "
                    "ordering is a project contract assumption."
                ),
            },
        },
        "presupernova_structure": {
            "available_model_count": len(structures),
            "explicit_null_model_count": len(missing_structures),
            "missing_coordinates": missing_structures,
            "binding_energy_is_not_injected_explosion_energy": True,
        },
        "rows": rows,
        "admission": {
            "candidate_id": BOCCIOLI_ADMISSION_ID,
            "hard_blockers_unchanged": hard_blockers_unchanged,
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
        "phase_history_shared_code_sha256": _sha256(PHASE_HISTORY_TOOL_PATH),
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
        error_payload: dict[str, Any] = {"status": "error", "error": str(exc)}
        if isinstance(exc, Lc18FailedWindCrosscheckError) and exc.diagnostics:
            error_payload["diagnostics"] = exc.diagnostics
        print(json.dumps(error_payload, indent=2), file=sys.stderr)
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
