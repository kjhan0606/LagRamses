#!/usr/bin/env python3
"""Regression checks for the ADM leaf-density normalisation."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    ROOT / "patch" / "lagRamses" / "dark_cooling_fine.f90",
    ROOT / "patch" / "cuRamses" / "dark_cooling_fine.f90",
)


def leaf_density(masses: list[float], volume: float) -> float:
    return sum(masses) / volume


def test_leaf_mass_normalisation() -> None:
    volume = 8.0
    one_particle = leaf_density([2.0], volume)
    two_particles = leaf_density([2.0, 2.0], volume)
    unequal_particles = leaf_density([1.0, 3.0], volume)

    assert one_particle == 0.25
    assert two_particles == 0.50
    assert unequal_particles == 0.50

    other_leaf_density = leaf_density([7.0], volume)
    assert other_leaf_density == 0.875
    assert two_particles != leaf_density([2.0, 2.0, 7.0], volume)


def test_production_uses_leaf_mass_and_proper_volume() -> None:
    for source in SOURCES:
        text = source.read_text()
        assert "dm_mass_cell(ind) = dm_mass_cell(ind) +" in text
        assert "mp(ipart) * scale_d * scale_l**3" in text
        assert "rho_D = dm_mass_cell(ind) / vol_phys" in text
        assert "vol_phys = (dx_loc*scale_l)**3" in text
        assert "vol_phys = (dx_loc*scale_l/aexp)**3" not in text


if __name__ == "__main__":
    test_leaf_mass_normalisation()
    test_production_uses_leaf_mass_and_proper_volume()
    print("ADM leaf-density regression passed")
