#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "paper"))
import analyse_lowz_sidm_profiles as lowz  # noqa: E402


def make_profile() -> dict[str, np.ndarray]:
    radius = np.sqrt(lowz.EDGES_PKPC_H[:-1] * lowz.EDGES_PKPC_H[1:])
    return {
        "radius_pkpc_h": radius,
        "density": np.full_like(radius, 2.0),
        "enclosed_mass_hinv_msun": np.linspace(1.0, 34.0, radius.size),
        "circular_velocity_kms": np.full_like(radius, 200.0),
        "sigma_1d_kms": np.full_like(radius, 100.0),
        "count": np.full(radius.size, 80),
    }


def copy_profile(profile: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {key: np.array(value, copy=True) for key, value in profile.items()}


def main() -> None:
    cdm = make_profile()
    sidm = copy_profile(cdm)
    sidm["density"] *= 0.8
    sidm["enclosed_mass_hinv_msun"] *= 0.9
    sidm["circular_velocity_kms"] *= np.sqrt(0.9)
    sidm["sigma_1d_kms"] *= 1.1

    result = lowz.aperture_diagnostics({"cdm": cdm, "sidm3": sidm})
    samples = result["sidm3"]
    assert len(samples) == len(lowz.SCIENCE_APERTURES_PKPC_H)
    assert all(sample["valid"] for sample in samples)
    assert all(
        abs(sample["enclosed_mass_deficit_fraction"] - 0.1) < 1.0e-12
        for sample in samples
    )
    assert all(
        abs(sample["circular_velocity_ratio_to_cdm"] - np.sqrt(0.9)) < 1.0e-12
        for sample in samples
    )

    first_index = samples[0]["bin_index"]
    sidm["sigma_1d_kms"][first_index] = np.nan
    invalid = lowz.aperture_diagnostics({"cdm": cdm, "sidm3": sidm})
    assert not invalid["sidm3"][0]["valid"]
    assert invalid["sidm3"][0]["density_ratio_to_cdm"] is None
    print("low-z aperture diagnostics regression passed")


if __name__ == "__main__":
    main()
