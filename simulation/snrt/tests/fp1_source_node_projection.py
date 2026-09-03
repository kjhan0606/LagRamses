#!/usr/bin/env python3
"""Regression tests for canonical-row/source-node physical binding."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from fp1_source_node_fixture import approved_source_node_contract  # noqa: E402
from fp1_source_node_projection import (  # noqa: E402
    SourceNodeProjectionError,
    validate_canonical_row_against_source_node,
)


def _row(channel: int, age: float) -> dict:
    return {
        "channel": channel,
        "initial_mass_msun_per_star": 60.0,
        "birth_metallicity_mass_fraction": 0.001,
        "age_yr": age,
        "returned_mass_msun_per_star": 0.0,
        "remnant_mass_msun_per_star": 0.0,
        "energy_erg_per_star": 0.0,
        "momentum_g_cm_s_per_star": [0.0, 0.0, 0.0],
        "ejecta_msun_per_star": [0.0] * 11,
    }


def _expect_error(row: dict, node: dict, fragment: str) -> None:
    try:
        validate_canonical_row_against_source_node(
            row, node, "cumulative_physical_erg_per_initial_star"
        )
    except SourceNodeProjectionError as exc:
        assert fragment in str(exc), str(exc)
    else:
        raise AssertionError(f"expected source-node projection error containing {fragment!r}")


def main() -> int:
    node = approved_source_node_contract()["physical_nodes"][0]

    wind = _row(1, 0.5)
    wind["returned_mass_msun_per_star"] = 3.0
    wind["ejecta_msun_per_star"][0] = 3.0
    validate_canonical_row_against_source_node(
        wind, node, "cumulative_physical_erg_per_initial_star"
    )

    terminal = _row(3, 1.0)
    terminal["remnant_mass_msun_per_star"] = 54.0
    validate_canonical_row_against_source_node(
        terminal, node, "cumulative_physical_erg_per_initial_star"
    )

    inconsistent = _row(3, 5.0e6)
    inconsistent["returned_mass_msun_per_star"] = 50.0
    inconsistent["ejecta_msun_per_star"][0] = 50.0
    inconsistent["energy_erg_per_star"] = 1.0e51
    _expect_error(inconsistent, node, "returned mass disagrees")

    wrong_channel = _row(2, 1.0)
    _expect_error(wrong_channel, node, "may map only")

    wrong_metallicity = _row(1, 0.0)
    wrong_metallicity["birth_metallicity_mass_fraction"] = 0.01
    _expect_error(wrong_metallicity, node, "birth metallicity disagrees")

    print("FP1_SOURCE_NODE_PROJECTION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
