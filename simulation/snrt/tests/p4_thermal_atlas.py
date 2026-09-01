"""P4 offline-atlas interpolation and HDF5 round-trip checks."""

from __future__ import annotations

from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.thermal_atlas import ThermalAtlas, read_thermal_atlas, write_thermal_atlas


def main() -> None:
    axes = (np.array([0.1, 0.2]), np.array([-4.0, 0.0]), np.array([2.0, 6.0]))
    shape = tuple(len(axis) for axis in axes)
    mean_mu = np.empty(shape)
    rate = np.empty(shape)
    for ia, aexp in enumerate(axes[0]):
        mean_mu[ia] = 1.0 + aexp
        rate[ia] = np.array([-1.0, 1.0])
    equilibrium = np.full(shape[:-1], 4.0)
    atlas = ThermalAtlas(
        *axes,
        net_rate_erg_s_cm3=rate,
        mean_molecular_weight=mean_mu,
        equilibrium_log_temperature_k=equilibrium,
        provenance={
            "thermal_component": "metal_only",
            "primordial_rates_included": "false",
            "uv_background_included": "false",
            "photoheating_included": "false",
            "metallicity_scaling": "linear_z_solar",
            "metallicity_application": "analytic_runtime_multiplier",
            "rate_sign_convention": "heating_positive_cooling_negative",
            "source_data_name": "synthetic_no_uvb.h5",
            "source_data_sha256": "1" * 64,
            "source_cooling_dataset": "CoolingRates/Metals/Cooling",
            "source_repository": "synthetic-test",
            "source_repository_revision": "synthetic-test-v1",
            "generator_name": "p4_thermal_atlas.py",
            "generator_version": "test",
            "generator_sha256": "2" * 64,
            "grackle_version": "test",
            "grackle_revision": "3" * 40,
            "generated_utc": "2026-09-01T00:00:00+00:00",
            "cmb_metal_floor": "continuous_subtraction_at_tcmb",
        },
    )
    assert np.isclose(atlas.mean_mu(0.15, 1.0e4, 1.0e-2, 1.0e-3), 1.15)
    assert np.isclose(atlas.net_rate(0.15, 1.0e4, 1.0e-2, 1.0e-3), 0.0)
    assert np.isclose(
        atlas.net_rate(0.1, 1.0e2, 1.0e-4, 3.7e-3),
        3.7e-3 * atlas.net_rate(0.1, 1.0e2, 1.0e-4, 1.0),
    )
    assert np.isclose(atlas.equilibrium_temperature(0.15, 1.0e-2, 1.0e-3), 1.0e4)
    with TemporaryDirectory() as temporary_directory:
        path = Path(temporary_directory) / "atlas.h5"
        write_thermal_atlas(path, atlas)
        restored = read_thermal_atlas(path)
    assert np.array_equal(restored.mean_molecular_weight, atlas.mean_molecular_weight)
    print("P4_THERMAL_ATLAS_OK slices=2")


if __name__ == "__main__":
    main()
