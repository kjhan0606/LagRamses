"""Synthetic conservative HDF5/AMR staging test.

Run with: ``python tests/p4_hdf5_staging.py``
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import h5py
import numpy as np

from snrt_core.snapshot import read_static_rt_input
from snrt_core.thermal_atlas import ThermalAtlas, write_thermal_atlas


PROJECT_ROOT = Path(__file__).resolve().parents[1]
STAGE_PATH = PROJECT_ROOT / "tools" / "p4_stage_hdf5_level15.py"


def _load_stage_module():
    spec = importlib.util.spec_from_file_location("p4_stage_hdf5_level15", STAGE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("unable to load HDF5 staging module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    stage = _load_stage_module()
    shape = (2, 2, 2)
    density_code = np.arange(1.0, 9.0).reshape(-1)
    velocity_x = np.arange(11.0, 19.0)
    velocity_y = np.arange(21.0, 29.0)
    velocity_z = np.arange(31.0, 39.0)
    pressure_code = density_code * 6.9e-3
    refined_density_code = np.arange(10.0, 18.0)
    refined_velocity_x = np.arange(101.0, 109.0)
    refined_velocity_y = np.arange(201.0, 209.0)
    refined_velocity_z = np.arange(301.0, 309.0)
    gamma = 1.4
    metal_mass_fraction = 0.02 * 1.0e-2
    refined_pressure_code = refined_density_code * 6.9e-3

    def total_energy(density, vx, vy, vz, pressure):
        return pressure / (gamma - 1.0) + 0.5 * density * (vx**2 + vy**2 + vz**2)

    total_energy_code = total_energy(density_code, velocity_x, velocity_y, velocity_z, pressure_code)
    refined_total_energy_code = total_energy(
        refined_density_code,
        refined_velocity_x,
        refined_velocity_y,
        refined_velocity_z,
        refined_pressure_code,
    )
    with TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        snapshot_path = root / "data_00001.h5"
        info_path = root / "info_00001.txt"
        manifest_path = root / "zoom.json"
        field_map_path = root / "field_map.json"
        atlas_path = root / "atlas.h5"
        source_path = root / "sources.csv"
        output_path = root / "staged.h5"
        metadata_path = root / "staged.json"

        info_path.write_text(
            "aexp = 0.500000\n"
            "unit_l = 1.000000e5\n"
            "unit_d = 1.000000e-24\n"
            "unit_t = 1.000000\n"
        )
        manifest_path.write_text(
            json.dumps({"final": {"left_edge_code": [0.0, 0.0, 0.0], "width_code": 1.0, "shape": list(shape)}})
        )
        field_map_path.write_text(
            json.dumps(
                {
                    "schema": "snrt-hdf5-field-map",
                    "schema_version": 1,
                    "fields": {
                        "density": {"dataset": "uold_1", "unit": "code_density", "averaging": "volume"},
                        "velocity": [
                            {
                                "dataset": "uold_2",
                                "unit": "code_momentum_density",
                                "averaging": "volume",
                                "quantity": "momentum_density",
                            },
                            {
                                "dataset": "uold_3",
                                "unit": "code_momentum_density",
                                "averaging": "volume",
                                "quantity": "momentum_density",
                            },
                            {
                                "dataset": "uold_4",
                                "unit": "code_momentum_density",
                                "averaging": "volume",
                                "quantity": "momentum_density",
                            },
                        ],
                        "thermal_pressure": {
                            "derive": "ideal_gas_pressure_from_conservative",
                            "depends_on": ["density", "velocity", "total_energy_density"],
                            "gamma": gamma,
                        },
                        "total_energy_density": {
                            "dataset": "uold_5",
                            "unit": "code_energy_density",
                            "averaging": "volume",
                        },
                        "metal_density": {
                            "dataset": "uold_6",
                            "unit": "code_density",
                            "averaging": "volume",
                        },
                        "metallicity_solar": {
                            "derive": "metallicity_from_density_and_metal_density",
                            "depends_on": ["density", "metal_density"],
                            "solar_mass_fraction": 0.02,
                        },
                        "dust_to_metal": {
                            "constant": 0.5,
                            "unit": "dimensionless",
                            "averaging": "volume",
                            "reason": "synthetic test dust ratio",
                        },
                        "x_hii": {
                            "constant": 0.1,
                            "unit": "dimensionless",
                            "averaging": "volume",
                            "reason": "synthetic test ionization",
                        },
                        "x_heii": {
                            "constant": 0.2,
                            "unit": "dimensionless",
                            "averaging": "volume",
                            "reason": "synthetic test ionization",
                        },
                        "x_heiii": {
                            "constant": 0.3,
                            "unit": "dimensionless",
                            "averaging": "volume",
                            "reason": "synthetic test ionization",
                        },
                        "x_h2": {
                            "constant": 0.4,
                            "unit": "dimensionless",
                            "averaging": "volume",
                            "reason": "synthetic test molecular fraction",
                        },
                    },
                }
            )
        )
        source_path.write_text(
            "source_id,source_kind,x_code,y_code,z_code,q_group_0_s\n"
            "1,agn,0.25,0.25,0.25,1.0e49\n"
        )
        atlas_shape = (2, 2, 2)
        write_thermal_atlas(
            atlas_path,
            ThermalAtlas(
                scale_factor=np.array([0.4, 0.6]),
                log_hydrogen_number_density_cm3=np.array([-2.0, 2.0]),
                log_temperature_k=np.array([1.0, 5.0]),
                net_rate_erg_s_cm3=np.zeros(atlas_shape),
                mean_molecular_weight=np.full(atlas_shape, 1.2),
                equilibrium_log_temperature_k=np.full((2, 2), 3.0),
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
                    "generator_name": "p4_hdf5_staging.py",
                    "generator_version": "test",
                    "generator_sha256": "2" * 64,
                    "grackle_version": "test",
                    "grackle_revision": "3" * 40,
                    "generated_utc": "2026-09-01T00:00:00+00:00",
                    "cmb_metal_floor": "continuous_subtraction_at_tcmb",
                },
            ),
        )
        with h5py.File(snapshot_path, "w") as handle:
            header = handle.create_group("header")
            header.attrs["gamma"] = gamma
            amr = handle.create_group("amr/level_1")
            amr.create_dataset("xg_1", data=[0.5])
            amr.create_dataset("xg_2", data=[0.5])
            amr.create_dataset("xg_3", data=[0.5])
            son_flags = np.zeros(8)
            son_flags[0] = 1.0
            amr.create_dataset("son_flag", data=son_flags)
            hydro = handle.create_group("hydro/level_1")
            hydro.create_dataset("uold_1", data=density_code)
            hydro.create_dataset("uold_2", data=density_code * velocity_x)
            hydro.create_dataset("uold_3", data=density_code * velocity_y)
            hydro.create_dataset("uold_4", data=density_code * velocity_z)
            hydro.create_dataset("uold_5", data=total_energy_code)
            hydro.create_dataset("uold_6", data=density_code * metal_mass_fraction)
            refined_amr = handle.create_group("amr/level_2")
            refined_amr.create_dataset("xg_1", data=[0.25])
            refined_amr.create_dataset("xg_2", data=[0.25])
            refined_amr.create_dataset("xg_3", data=[0.25])
            refined_amr.create_dataset("son_flag", data=np.zeros(8))
            refined_hydro = handle.create_group("hydro/level_2")
            refined_hydro.create_dataset("uold_1", data=refined_density_code)
            refined_hydro.create_dataset("uold_2", data=refined_density_code * refined_velocity_x)
            refined_hydro.create_dataset("uold_3", data=refined_density_code * refined_velocity_y)
            refined_hydro.create_dataset("uold_4", data=refined_density_code * refined_velocity_z)
            refined_hydro.create_dataset("uold_5", data=refined_total_energy_code)
            refined_hydro.create_dataset("uold_6", data=refined_density_code * metal_mass_fraction)

        original_argv = sys.argv
        try:
            sys.argv = [
                str(STAGE_PATH),
                "--snapshot",
                str(snapshot_path),
                "--info",
                str(info_path),
                "--zoom-manifest",
                str(manifest_path),
                "--field-map",
                str(field_map_path),
                "--thermal-atlas",
                str(atlas_path),
                "--source-ledger",
                str(source_path),
                "--output",
                str(output_path),
                "--metadata-output",
                str(metadata_path),
                "--analysis-level",
                "1",
            ]
            stage.main()
            sys.argv = [
                str(STAGE_PATH),
                "--snapshot",
                str(snapshot_path),
                "--info",
                str(info_path),
                "--zoom-manifest",
                str(manifest_path),
                "--field-map",
                str(field_map_path),
                "--thermal-atlas",
                str(atlas_path),
                "--source-ledger",
                str(source_path),
                "--output",
                str(output_path),
                "--metadata-output",
                str(metadata_path),
                "--analysis-level",
                "1",
                "--preflight-only",
            ]
            stage.main()
        finally:
            sys.argv = original_argv

        restored = read_static_rt_input(output_path)
        metadata = json.loads(metadata_path.read_text())
        assert stage._production_contract_complete(json.loads(field_map_path.read_text()), sources_present=True) is False
        expected_density = density_code.reshape(shape, order="F")
        expected_density[0, 0, 0] = refined_density_code.mean()
        assert np.allclose(
            restored.hydrogen_number_density_cm3 / (0.76 / 1.67262192369e-24),
            expected_density * 1.0e-24,
        )
        expected_velocity_x = velocity_x.reshape(shape, order="F")
        expected_velocity_y = velocity_y.reshape(shape, order="F")
        expected_velocity_z = velocity_z.reshape(shape, order="F")
        expected_velocity_x[0, 0, 0] = np.average(refined_velocity_x, weights=refined_density_code)
        expected_velocity_y[0, 0, 0] = np.average(refined_velocity_y, weights=refined_density_code)
        expected_velocity_z[0, 0, 0] = np.average(refined_velocity_z, weights=refined_density_code)
        assert np.allclose(restored.velocity_cm_s[0] / 1.0e5, expected_velocity_x)
        assert np.allclose(restored.velocity_cm_s[1] / 1.0e5, expected_velocity_y)
        assert np.allclose(restored.velocity_cm_s[2] / 1.0e5, expected_velocity_z)
        assert np.allclose(restored.metallicity_solar, 1.0e-2)
        assert np.allclose(restored.dust_to_metal, 0.5)
        assert np.allclose(restored.dust_relative_abundance, 5.0e-3)
        assert np.allclose(restored.x_hii, 0.1)
        assert np.allclose(restored.x_heii, 0.2)
        assert np.allclose(restored.x_heiii, 0.3)
        assert np.allclose(restored.x_h2, 0.4)
        assert np.all(restored.cell_level == 1)
        assert restored.sources is not None and restored.sources.cell_index.tolist() == [[0, 0, 0]]
        assert metadata["coverage_min"] == 1.0 and metadata["coverage_max"] == 1.0
        assert metadata["mass_relative_error"] < 1.0e-14
        assert len(metadata["resampled_fields"]) == 11
        assert metadata["level_summary"][0]["refined_cells_in_cube"] == 1
        assert metadata["level_summary"][1]["deposited_leaf_cells"] == 8
        assert metadata["production_contract_complete"] is False
        assert metadata["field_contract"]["thermal_pressure"]["status"] == "derived"
        assert metadata["snapshot_header_gamma"] == gamma
    print("P4_HDF5_STAGING_OK raw_conservative_fields=11 coverage=1 production_contract=False")


if __name__ == "__main__":
    main()
