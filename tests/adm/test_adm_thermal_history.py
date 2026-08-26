#!/usr/bin/env python3
"""Unit tests for the explicit post-decoupling ADM temperature prescription."""

from __future__ import annotations

from adm_thermal_history import (
    T_CMB0_K,
    dark_cmb_temperature,
    floor_equivalent_z_kd,
    post_decoupling_temperature,
)


def test_coupled_limit() -> None:
    xi = 0.5
    z_init = 100.0
    assert dark_cmb_temperature(xi, z_init) == post_decoupling_temperature(
        xi, z_init, z_init
    )


def test_post_decoupling_scaling() -> None:
    xi = 0.5
    z_init = 100.0
    z_kd = 1000.0
    expected = xi * T_CMB0_K * (1.0 + z_init) ** 2 / (1.0 + z_kd)
    assert abs(post_decoupling_temperature(xi, z_init, z_kd) / expected - 1.0) < 1.0e-14


def test_floor_equivalent_decoupling() -> None:
    z_kd = floor_equivalent_z_kd(0.5, 100.0, 1.0)
    assert abs(post_decoupling_temperature(0.5, 100.0, z_kd) - 1.0) < 1.0e-13


def test_rejects_coupled_start() -> None:
    try:
        post_decoupling_temperature(0.5, 100.0, 99.0)
    except ValueError:
        return
    raise AssertionError("z_kd below z_init must be rejected")


if __name__ == "__main__":
    test_coupled_limit()
    test_post_decoupling_scaling()
    test_floor_equivalent_decoupling()
    test_rejects_coupled_start()
    print("ADM thermal-history regression passed")
