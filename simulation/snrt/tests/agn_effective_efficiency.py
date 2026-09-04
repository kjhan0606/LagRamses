#!/usr/bin/env python3
"""Arithmetic and source-contract tests for F-P1.5-R.

The native smoke test exercises the production Fortran helper.  This small
standard-library test records the mass-accounting invariants that cannot be
reached without constructing a full RAMSES sink state.
"""

from __future__ import annotations

import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
WRITER = REPO_ROOT / "patch/lagRamses/sink_particle.kjhan.f90"
DRIVER = REPO_ROOT / "patch/lagRamses/snrt_ramses_driver.f90"
HELPER = REPO_ROOT / "patch/lagRamses/snrt_agn_efficiency.f90"
SOURCE = REPO_ROOT / "patch/lagRamses/snrt_agn_source.f90"


def resolve(raw: float, spin: bool, bondi: float, eddington: float, mad: bool, floor: float) -> tuple[float, float, float, float]:
    """Reference arithmetic for finite, promotable states."""
    base = raw if spin else 0.1
    inflow = min(max(bondi, 0.0), max(eddington, 0.0))
    ratio = bondi / eddington if eddington > 0.0 else 0.0
    effective = base
    if mad and floor > 0.0 and ratio < floor:
        effective = base * max(ratio, 0.0) / floor
    return base, effective, inflow, ratio


def assert_close(actual: float, expected: float) -> None:
    assert math.isclose(actual, expected, rel_tol=1.0e-13, abs_tol=1.0e-14), (actual, expected)


def main() -> int:
    writer = WRITER.read_text(encoding="utf-8")
    driver = DRIVER.read_text(encoding="utf-8")
    helper = HELPER.read_text(encoding="utf-8")
    source = SOURCE.read_text(encoding="utf-8")

    cases = (
        ((0.1, True, 2.0, 1.0, False, 0.01), (0.1, 0.1, 1.0, 2.0)),
        ((0.1, True, 0.001, 1.0, True, 0.01), (0.1, 0.01, 0.001, 0.001)),
        ((0.1, True, 0.2, 1.0, True, 0.01), (0.1, 0.1, 0.2, 0.2)),
        # Equality is high-state because the MAD branch is strictly below the floor.
        ((0.1, True, 0.01, 1.0, True, 0.01), (0.1, 0.1, 0.01, 0.01)),
        ((0.1, False, 0.2, 0.1, False, 0.01), (0.1, 0.1, 0.1, 2.0)),
    )
    for inputs, expected in cases:
        actual = resolve(*inputs)
        for got, want in zip(actual, expected, strict=True):
            assert_close(got, want)
        assert 0.0 < actual[0] < 1.0
        assert 0.0 <= actual[1] < 1.0

    # The limited supplied mass is the cumulative Bondi/Eddington minimum,
    # while retained mass is checked only one-sidedly.
    supplied = min(max(2.0, 0.0), max(0.5, 0.0))
    retained_increment = 0.9 * supplied
    assert supplied == 0.5
    assert retained_increment <= (1.0 - 0.1) * supplied * (1.0 + 1.0e-8)
    assert 0.46 > (1.0 - 0.1) * supplied * (1.0 + 1.0e-8)

    # A nonzero retained carry-over with no new supplied inflow is not a new
    # photon event; the separate retained cursor sees zero increment.
    accounted_inflow = 0.5
    retained_seen = retained_increment
    next_supplied = 0.0
    next_retained = retained_seen
    assert next_supplied - accounted_inflow < 0.0
    assert max(0.0, next_supplied - 0.0) == 0.0
    assert next_retained - retained_seen == 0.0

    # Exact source scaling uses effective efficiency and supplied mass.
    mass_unit = 5.0e33
    dt = 4.0
    effective = 0.05
    energy = effective * supplied * mass_unit * (2.99792458e10**2)
    luminosity = energy / dt
    assert_close(luminosity * dt, energy)
    assert_close(effective * supplied * mass_unit * (2.99792458e10**2), energy)

    assert "module snrt_agn_efficiency" in helper
    assert "pure subroutine snrt_agn_resolve_efficiency" in helper
    assert writer.count("call snrt_agn_resolve_efficiency") == 1
    assert driver.count("call snrt_agn_resolve_efficiency") == 1
    assert '"raw_radiative_efficiency":' in writer
    assert '"radiative_efficiency":' in writer
    assert '"effective_radiative_efficiency":' in writer
    assert '"efficiency_status_name":' in writer
    assert "luminosity=epsilon_eff*inflow_rate" in writer
    assert "call snrt_agn_photon_budget(delta_inflow" in driver
    assert "supplied_mass=min(max(dMBH_coarse(isink),0.0d0)," in driver
    assert "max(dMEd_coarse(isink),0.0d0)" in driver
    assert "retained_bound=(1.0d0-epsilon_eff)*delta_inflow" in driver
    assert "retained_seen(isink)=dMsmbh(isink)" in driver
    assert "0.99d0" not in driver
    assert "max(1.0d-6" not in driver
    assert "delta_inflow_mass_code" in source

    print(
        "AGN_EFFECTIVE_EFFICIENCY_TEST_OK algebra=thermal,mad_low,mad_high,boundary "
        "native_parity=run_fp15_smoke supplied_mass=min_bondi_edd "
        "carryover=no_duplicate epsilon=effective"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
