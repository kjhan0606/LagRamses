"""Validate one canonical cumulative row against its F-P1 physical source node."""

from __future__ import annotations

import math
from typing import Any


class SourceNodeProjectionError(ValueError):
    """A canonical row is not the declared source node's cumulative state."""


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SourceNodeProjectionError(f"source node lacks numeric {field}")
    number = float(value)
    if not math.isfinite(number):
        raise SourceNodeProjectionError(f"source node has nonfinite {field}")
    return number


def _vector(value: Any, length: int, field: str) -> list[float]:
    if not isinstance(value, list) or len(value) != length:
        raise SourceNodeProjectionError(f"source node lacks explicit {field}")
    return [_number(item, f"{field}[{index}]") for index, item in enumerate(value)]


def _interpolate(age: float, ages: list[float], values: list[float]) -> float:
    if age < ages[0]:
        return 0.0
    if age >= ages[-1]:
        return values[-1]
    for left in range(len(ages) - 1):
        if ages[left] <= age <= ages[left + 1]:
            fraction = (age - ages[left]) / (ages[left + 1] - ages[left])
            return values[left] + fraction * (values[left + 1] - values[left])
    raise SourceNodeProjectionError("wind interpolation interval was not found")


def _assert_close(observed: float, expected: float, field: str) -> None:
    if not math.isclose(observed, expected, rel_tol=1.0e-10, abs_tol=1.0e-12):
        raise SourceNodeProjectionError(
            f"canonical {field} disagrees with source node: {observed} != {expected}"
        )


def validate_canonical_row_against_source_node(
    row: dict[str, Any], node: dict[str, Any], energy_semantics: str
) -> None:
    """Check mass, time, channel, and cumulative physical payload exactly."""
    channel = int(row["channel"])
    age = _number(row["age_yr"], "canonical age_yr")
    _assert_close(
        _number(row["initial_mass_msun_per_star"], "canonical initial mass"),
        _number(node.get("zams_mass_msun"), "zams_mass_msun"),
        "initial mass",
    )
    _assert_close(
        _number(row["birth_metallicity_mass_fraction"], "canonical metallicity"),
        _number(node.get("birth_metallicity_value"), "birth_metallicity_value"),
        "birth metallicity",
    )

    expected_returned = 0.0
    expected_remnant = 0.0
    expected_energy = 0.0
    expected_momentum = [0.0, 0.0, 0.0]
    expected_ejecta = [0.0] * 11

    if channel == 1:
        raw_ages = node.get("wind_release_age_yr")
        raw_masses = node.get("cumulative_wind_mass_msun")
        if not isinstance(raw_ages, list) or not isinstance(raw_masses, list):
            raise SourceNodeProjectionError("source node lacks an age-resolved wind history")
        ages = [_number(value, "wind_release_age_yr") for value in raw_ages]
        masses = [
            _number(value, "cumulative_wind_mass_msun")
            for value in raw_masses
        ]
        tracked = node.get("cumulative_wind_tracked_elements_msun")
        if not ages or len(ages) != len(masses) or not isinstance(tracked, list) or len(tracked) != len(ages):
            raise SourceNodeProjectionError("source node lacks an age-resolved wind history")
        expected_returned = _interpolate(age, ages, masses)
        for element in range(11):
            history = [
                _number(_vector(record, 11, "cumulative wind tracked elements")[element], "wind element")
                for record in tracked
            ]
            expected_ejecta[element] = _interpolate(age, ages, history)
    elif channel == 3:
        outcome = node.get("outcome")
        if outcome == "pisn_complete_disruption":
            raise SourceNodeProjectionError("PISN complete disruption cannot be mapped to channel 3")
        lifetime = _number(
            node.get("lifetime_yr_or_declared_no_terminal_horizon"),
            "terminal lifetime",
        )
        if age >= lifetime and outcome != "not_terminal_within_horizon":
            if node.get("terminal_remnant_owner_channel") != 3:
                raise SourceNodeProjectionError("source node does not assign terminal ownership to channel 3")
            expected_returned = _number(
                node.get("terminal_ejecta_mass_msun"), "terminal_ejecta_mass_msun"
            )
            expected_ejecta = _vector(
                node.get("terminal_ejecta_tracked_elements_msun"),
                11,
                "terminal_ejecta_tracked_elements_msun",
            )
            expected_remnant = _number(
                node.get("baryonic_remnant_mass_msun"),
                "baryonic_remnant_mass_msun",
            )
            if energy_semantics == "cumulative_physical_erg_per_initial_star":
                expected_energy = _number(
                    node.get("final_kinetic_energy_erg"), "final_kinetic_energy_erg"
                )
            elif energy_semantics == "cumulative_injected_erg_per_initial_star":
                expected_energy = _number(
                    node.get("injected_energy_erg_or_null"),
                    "injected_energy_erg_or_null",
                )
            else:
                raise SourceNodeProjectionError("unsupported canonical energy semantics")
            expected_momentum = _vector(
                node.get("source_frame_vector_momentum_g_cm_s"),
                3,
                "source_frame_vector_momentum_g_cm_s",
            )
    else:
        raise SourceNodeProjectionError(
            "F-P1 high-mass source nodes may map only to wind channel 1 or SNII channel 3"
        )

    observed_momentum = _vector(row.get("momentum_g_cm_s_per_star"), 3, "canonical momentum")
    observed_ejecta = _vector(row.get("ejecta_msun_per_star"), 11, "canonical ejecta")
    _assert_close(_number(row["returned_mass_msun_per_star"], "canonical returned mass"), expected_returned, "returned mass")
    _assert_close(_number(row["remnant_mass_msun_per_star"], "canonical remnant mass"), expected_remnant, "remnant mass")
    _assert_close(_number(row["energy_erg_per_star"], "canonical energy"), expected_energy, "energy")
    for index, (observed, expected) in enumerate(zip(observed_momentum, expected_momentum)):
        _assert_close(observed, expected, f"momentum[{index}]")
    for index, (observed, expected) in enumerate(zip(observed_ejecta, expected_ejecta)):
        _assert_close(observed, expected, f"ejecta[{index}]")
