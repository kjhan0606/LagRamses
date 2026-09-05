"""Shared, fail-closed Limongi phase-history aggregation.

The G2 phase-history audit and the F-P1 LC18 cross-check must derive their
model/phase counts and cumulative wind histories from one implementation.
This module performs only source-record aggregation; it does not approve
physical source data or emit canonical/runtime payloads.
"""

from __future__ import annotations

from collections import defaultdict
import math
from typing import Any


class PhaseHistoryInvariantError(ValueError):
    """The source phase history cannot be aggregated without ambiguity."""

    def __init__(
        self, message: str, *, diagnostics: dict[str, Any] | None = None
    ) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics or {}


def _coordinate(record: dict[str, Any]) -> tuple[int, int, float]:
    source = record["source_coordinate"]
    return (
        int(source["rotation_velocity_km_s"]),
        int(source["metallicity_feh"]),
        float(source["initial_mass_msun"]),
    )


def build_phase_histories(
    records: list[dict[str, Any]], phase_order: list[str]
) -> tuple[dict[tuple[int, int, float], dict[str, Any]], dict[str, Any]]:
    """Collapse exact duplicate occurrences and build validated histories."""
    phase_rank = {phase: index for index, phase in enumerate(phase_order)}
    grouped: dict[tuple[int, int, float], list[dict[str, Any]]] = defaultdict(list)
    collapsed_duplicates: list[dict[str, Any]] = []
    for record in records:
        try:
            source = record["source_coordinate"]
            occurrence = int(source["phase_occurrence"])
            coordinate = _coordinate(record)
        except (KeyError, TypeError, ValueError) as exc:
            raise PhaseHistoryInvariantError(
                "phase-history record has an invalid source coordinate"
            ) from exc
        if occurrence > 1:
            collapsed_duplicates.append(
                {
                    "coordinate": {
                        "rotation_velocity_km_s": coordinate[0],
                        "metallicity_feh": coordinate[1],
                        "initial_mass_msun": coordinate[2],
                    },
                    "phase": source.get("phase"),
                    "occurrence": occurrence,
                    "source_line": record.get("source_line"),
                }
            )
            continue
        grouped[coordinate].append(record)

    if len(grouped) != 108:
        raise PhaseHistoryInvariantError(
            f"expected 108 CDS phase-history models, found {len(grouped)}",
            diagnostics={
                "model_count": len(grouped),
                "collapsed_duplicate_row_count": len(collapsed_duplicates),
            },
        )

    age_violations: list[dict[str, Any]] = []
    mass_violations: list[dict[str, Any]] = []
    negative_cumulative_mass: list[dict[str, Any]] = []
    missing_terminal_phase: list[dict[str, Any]] = []
    phase_order_violations: list[dict[str, Any]] = []
    histories: dict[tuple[int, int, float], dict[str, Any]] = {}
    unique_phase_counts: list[int] = []
    terminal_ages: list[float] = []

    for coordinate, model_records in grouped.items():
        try:
            source_ranks = [
                phase_rank[record["source_coordinate"]["phase"]]
                for record in model_records
            ]
            ordered = sorted(
                model_records,
                key=lambda value: phase_rank[value["source_coordinate"]["phase"]],
            )
        except KeyError as exc:
            raise PhaseHistoryInvariantError(
                f"unknown CDS phase {exc.args[0]!r}"
            ) from exc
        if source_ranks != sorted(source_ranks):
            phase_order_violations.append(
                {
                    "coordinate": list(coordinate),
                    "observed_phase_ranks": source_ranks,
                }
            )
        phases = [record["source_coordinate"]["phase"] for record in ordered]
        if len(phases) != len(set(phases)):
            raise PhaseHistoryInvariantError(
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
        previous_total_mass = coordinate[2]
        nodes: list[dict[str, Any]] = []
        for record in ordered:
            source = record["source_coordinate"]
            duration = float(record["phase_duration_yr"])
            total_mass = float(record["total_mass_msun"])
            cumulative_age += duration
            cumulative_wind = coordinate[2] - total_mass
            node = {
                "phase": source["phase"],
                "source_line": record["source_line"],
                "phase_duration_yr": duration,
                "cumulative_age_yr": cumulative_age,
                "total_mass_msun": total_mass,
                "cumulative_wind_mass_msun": cumulative_wind,
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
                or total_mass > previous_total_mass + 1.0e-12
            ):
                mass_violations.append(
                    {
                        "coordinate": list(coordinate),
                        "previous_mass_msun": previous_total_mass,
                        "node": dict(node),
                    }
                )
            if cumulative_wind < -1.0e-12:
                negative_cumulative_mass.append(
                    {"coordinate": list(coordinate), "node": dict(node)}
                )
            previous_age = cumulative_age
            previous_total_mass = total_mass

        if not nodes:
            raise PhaseHistoryInvariantError(
                f"phase-history model has no non-duplicate phases at {coordinate}"
            )
        terminal_ages.append(cumulative_age)
        unique_phase_counts.append(len(nodes))
        histories[coordinate] = {
            "unique_phase_count": len(nodes),
            "terminal_phase": nodes[-1]["phase"],
            "terminal_age_yr": nodes[-1]["cumulative_age_yr"],
            "terminal_total_mass_msun": nodes[-1]["total_mass_msun"],
            "terminal_cumulative_wind_mass_msun": nodes[-1][
                "cumulative_wind_mass_msun"
            ],
            "nodes": nodes,
        }

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
        "negative_cumulative_mass_count": len(negative_cumulative_mass),
        "negative_cumulative_mass": negative_cumulative_mass,
        "missing_psn_terminal_phase_count": len(missing_terminal_phase),
        "missing_psn_terminal_phase": missing_terminal_phase,
        "observed_source_order_matches_contract_rank": not phase_order_violations,
        "phase_order_violation_count": len(phase_order_violations),
        "phase_order_violations": phase_order_violations,
        "minimum_terminal_age_yr": min(terminal_ages),
        "maximum_terminal_age_yr": max(terminal_ages),
    }
    if (
        age_violations
        or mass_violations
        or negative_cumulative_mass
        or missing_terminal_phase
        or phase_order_violations
    ):
        raise PhaseHistoryInvariantError(
            "phase-history invariants violated", diagnostics=diagnostics
        )
    return histories, diagnostics
