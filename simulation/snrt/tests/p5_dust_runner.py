#!/usr/bin/env python3
"""Exercise the physical Draine dust sidecar through the P5 thermal runner."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.snapshot import GridSpec, SourceCatalog, neutral_primordial_input, write_static_rt_input
from tools import p4_build_agn_photon_ledger as agn_ledger
from tools.build_draine_dust_thermal import build_thermal_metadata


def main() -> int:
    dust_metadata_path = ROOT.parents[1] / "external" / "draine_wd01_rv31" / "p0_dust_opacity_rv31_photon_index1.json"
    atlas_path = ROOT / "data" / "production_metal_thermal_atlas_v2.h5"
    if not dust_metadata_path.is_file():
        raise FileNotFoundError(dust_metadata_path)
    if not atlas_path.is_file():
        raise FileNotFoundError(atlas_path)

    edges = agn_ledger._read_group_edges(agn_ledger.DEFAULT_P0_GROUP_EDGES)
    photon_per_lbol, mean_energy, closure = agn_ledger._group_conversion(0.1, edges)
    source_luminosity = 1.0e35 * photon_per_lbol
    photon_metadata = {
        "schema": "agn_photon_source_ledger_v1",
        "source_kind": "agn",
        "group_edges_ev": edges.tolist(),
        "groups": [
            {
                "index": index,
                "energy_interval_ev": [float(edges[index]), float(edges[index + 1])],
                "photon_weighted_mean_energy_ev": float(mean_energy[index]),
                "total_photon_rate_s": float(source_luminosity[index]),
            }
            for index in range(edges.size - 1)
        ],
        "group_spectral_closure": {
            "cross_sections_cm2": {
                "hydrogen_i": np.asarray(closure.cross_sections.hydrogen_i).tolist(),
                "helium_i": np.asarray(closure.cross_sections.helium_i).tolist(),
                "helium_ii": np.asarray(closure.cross_sections.helium_ii).tolist(),
            },
            "photoelectron_excess_energy_ev": {
                "hydrogen_i": np.asarray(closure.photoelectron_excess_energy_ev[0]).tolist(),
                "helium_i": np.asarray(closure.photoelectron_excess_energy_ev[1]).tolist(),
                "helium_ii": np.asarray(closure.photoelectron_excess_energy_ev[2]).tolist(),
            },
        },
    }

    with TemporaryDirectory(prefix="p5-dust-runner-test-") as directory:
        work = Path(directory)
        input_path = work / "input.h5"
        photon_path = work / "photon.json"
        output_path = work / "output.h5"
        shape = (2, 2, 2)
        snapshot = neutral_primordial_input(
            GridSpec(cell_width_cm=1.0e18, left_edge_cm=np.zeros(3)),
            np.full(shape, 1.0e-24),
            np.full(shape, 1.0e4),
            dust_relative_abundance=1.0,
            sources=SourceCatalog(
                cell_index=np.asarray([[0, 0, 0]]),
                photon_luminosity_s=source_luminosity[None, :],
            ),
        )
        write_static_rt_input(input_path, snapshot)
        photon_path.write_text(json.dumps(photon_metadata, indent=2) + "\n", encoding="utf-8")
        environment = os.environ.copy()
        environment["PYTHONPATH"] = str(ROOT)
        environment["JAX_PLATFORMS"] = "cpu"
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "p5_run_thermochemical_pilot.py"),
                "--input",
                str(input_path),
                "--photon-metadata",
                str(photon_path),
                "--dust-opacity-metadata",
                str(dust_metadata_path),
                "--thermal-atlas",
                str(atlas_path),
                "--scale-factor",
                "0.208497764676753",
                "--metallicity-solar",
                "1.0e-6",
                "--output",
                str(output_path),
                "--duration-myr",
                "1.0e-5",
                "--sn-order",
                "4",
                "--thermal-subcycles",
                "4",
                "--thermal-implicit-iterations",
                "24",
                "--time-averaged-absorption-iterations",
                "20",
            ],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        assert "P5_THERMOCHEMICAL_PILOT_OK" in result.stdout
        with h5py.File(output_path, "r") as handle:
            assert handle.attrs["dust_model"] == "metadata"
            assert handle.attrs["dust_opacity_schema"] == "snrt_dust_opacity_v1"
            assert handle.attrs["dust_binding_status"] == "reference_control"
            assert handle.attrs["dust_opacity_metadata_sha256"]
            assert handle.attrs["dust_payload_sha256"] == ""
            assert handle.attrs["dust_source_table_sha256"] == ""
            assert handle.attrs["dust_builder_sha256"] == ""
            assert handle.attrs["source_sed_identity"] == ""
            assert np.asarray(handle["diagnostics/cumulative_dust_absorbed_photons_cm3"]).shape == (9, *shape)
            assert np.asarray(handle["diagnostics/cumulative_dust_scattered_photons_cm3"]).shape == (9, *shape)
            assert np.max(np.asarray(handle["diagnostics/cumulative_dust_absorbed_photons_cm3"])) > 0.0
            assert np.allclose(handle["diagnostics/cumulative_dust_scattered_photons_cm3"], 0.0)
            assert np.max(np.asarray(handle["thermal/cumulative_dust_heating_energy_erg_cm3"])) > 0.0
            assert handle.attrs["primary_absorption_closure_relative_error"] <= 1.0e-5

        scattering_metadata = json.loads(dust_metadata_path.read_text(encoding="utf-8"))
        for key in (
            "group_edges_path",
            "group_edges_sha256",
            "source_table",
            "builder",
            "closure_code_manifest",
            "payload_hash_scheme",
            "payload_sha256",
        ):
            scattering_metadata.pop(key, None)
        scattering_metadata.update(
            {
                "schema": "snrt_dust_opacity_v3",
                "schema_version": 3,
                "status": "reference_scattering_control",
                "phase_function": "phase_isotropic_candidate",
                "scattering_cross_section_per_h_cm2": [5.0e-21] * (edges.size - 1),
                "scattering_weighted_energy_ev": np.sqrt(edges[:-1] * edges[1:]).tolist(),
                "scattering_angle_cosine": [0.0] * (edges.size - 1),
                "scattering_angle_cosine_squared": [0.5] * (edges.size - 1),
                "transport_corrected_scattering_cross_section_per_h_cm2": [
                    5.0e-21
                ]
                * (edges.size - 1),
                "isotropic_candidate_momentum_overestimate_factor": [1.0] * (edges.size - 1),
                "isotropic_candidate_momentum_bound_unbounded": [False] * (edges.size - 1),
            }
        )
        scattering_metadata_path = work / "dust-scattering.json"
        scattering_output_path = work / "scattering-output.h5"
        scattering_metadata_path.write_text(
            json.dumps(scattering_metadata, indent=2) + "\n", encoding="utf-8"
        )
        scattering_command = list(result.args)
        scattering_command[scattering_command.index(str(dust_metadata_path))] = str(
            scattering_metadata_path
        )
        scattering_command[scattering_command.index(str(output_path))] = str(scattering_output_path)
        scattering_command.extend(("--dust-scattering", "isotropic"))
        scattering_result = subprocess.run(
            scattering_command,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        if scattering_result.returncode != 0:
            raise RuntimeError(scattering_result.stdout + scattering_result.stderr)
        with h5py.File(scattering_output_path, "r") as handle:
            assert handle.attrs["dust_opacity_schema"] == "snrt_dust_opacity_v3"
            assert handle.attrs["dust_scattering"] == "isotropic"
            scattered = np.asarray(handle["diagnostics/cumulative_dust_scattered_photons_cm3"])
            assert np.max(scattered) > 0.0
            assert handle.attrs["dust_momentum_rate_semantics"] == "total_absorption_plus_scattering"
            total = np.asarray(handle["rates/dust_total_momentum_rate_dyn_cm3"])
            absorption = np.asarray(handle["rates/dust_absorption_momentum_rate_dyn_cm3"])
            scattering = np.asarray(handle["rates/dust_scattering_momentum_rate_dyn_cm3"])
            assert np.allclose(total, absorption + scattering)
            assert handle.attrs["primary_absorption_closure_relative_error"] <= 1.0e-5

        physical_opacity_path = (
            ROOT.parents[1]
            / "external"
            / "draine_wd01_rv31"
            / "p0_dust_opacity_rv31_photon_index1_scattering.json"
        )
        source_table_path = (
            ROOT.parents[1]
            / "external"
            / "draine_wd01_rv31"
            / "kext_albedo_WD_MW_3.1_60_D03.all"
        )
        thermal_metadata_path = work / "dust-thermal.json"
        thermal_metadata = build_thermal_metadata(
            source_table_path,
            agn_ledger.DEFAULT_P0_GROUP_EDGES,
            temperature_grid_k=np.geomspace(5.0, 300.0, 64),
        )
        thermal_metadata_path.write_text(
            json.dumps(thermal_metadata, indent=2) + "\n", encoding="utf-8"
        )
        thermal_output_path = work / "thermal-output.h5"
        thermal_command = list(scattering_command)
        thermal_command[thermal_command.index(str(scattering_metadata_path))] = str(physical_opacity_path)
        thermal_command[thermal_command.index(str(scattering_output_path))] = str(thermal_output_path)
        physical_control_output_path = work / "physical-control-output.h5"
        physical_control_command = list(thermal_command)
        physical_control_command[physical_control_command.index(str(thermal_output_path))] = str(
            physical_control_output_path
        )
        physical_control_result = subprocess.run(
            physical_control_command,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        if physical_control_result.returncode != 0:
            raise RuntimeError(physical_control_result.stdout + physical_control_result.stderr)
        thermal_command.extend(("--dust-thermal-metadata", str(thermal_metadata_path)))
        thermal_result = subprocess.run(
            thermal_command,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        if thermal_result.returncode != 0:
            raise RuntimeError(thermal_result.stdout + thermal_result.stderr)
        assert "dust_ir=snrt_dust_thermal_v1" in thermal_result.stdout
        with h5py.File(physical_control_output_path, "r") as baseline, h5py.File(thermal_output_path, "r") as handle:
            assert handle.attrs["dust_thermal_schema"] == "snrt_dust_thermal_v1"
            assert handle.attrs["dust_thermal_status"] == "candidate_kirchhoff_equilibrium"
            assert handle.attrs["dust_ir_transport_semantics"] == "recorded_not_transport_reemitted"
            assert handle.attrs["dust_ir_out_of_range_count"] == 0
            assert handle.attrs["dust_ir_energy_closure_relative_error"] <= 1.0e-5
            assert handle.attrs["dust_ir_power_residual_relative_error"] <= 1.0e-5
            assert np.asarray(handle["thermal/dust_grain_temperature_k"]).shape == shape
            assert np.asarray(handle["thermal/dust_ir_reemitted_energy_erg_cm3"]).shape == shape
            assert np.asarray(handle["thermal/dust_ir_untracked_energy_erg_cm3"]).shape == shape
            assert np.asarray(handle["sources/dust_ir_photon_rate_cm3_s"]).shape == (9, *shape)
            assert np.max(np.asarray(handle["thermal/dust_ir_reemitted_energy_erg_cm3"])) > 0.0
            assert np.max(np.asarray(handle["thermal/dust_ir_untracked_energy_erg_cm3"])) > 0.0
            assert np.allclose(
                baseline["ionization/x_hii"],
                handle["ionization/x_hii"],
                rtol=0.0,
                atol=0.0,
            )
            assert np.allclose(
                baseline["thermal/internal_energy_density_erg_cm3"],
                handle["thermal/internal_energy_density_erg_cm3"],
                rtol=0.0,
                atol=0.0,
            )
            for dataset in (
                "diagnostics/cumulative_dust_absorbed_photons_cm3",
                "diagnostics/cumulative_dust_scattered_photons_cm3",
                "thermal/cumulative_dust_heating_energy_erg_cm3",
            ):
                assert np.allclose(baseline[dataset], handle[dataset], rtol=0.0, atol=0.0)

        zero_input_path = work / "zero-dust-input.h5"
        zero_output_path = work / "zero-dust-thermal-output.h5"
        zero_snapshot = neutral_primordial_input(
            GridSpec(cell_width_cm=1.0e18, left_edge_cm=np.zeros(3)),
            np.full(shape, 1.0e-24),
            np.full(shape, 1.0e4),
            dust_relative_abundance=0.0,
            sources=SourceCatalog(
                cell_index=np.asarray([[0, 0, 0]]),
                photon_luminosity_s=source_luminosity[None, :],
            ),
        )
        write_static_rt_input(zero_input_path, zero_snapshot)
        zero_command = list(thermal_command)
        zero_command[zero_command.index(str(input_path))] = str(zero_input_path)
        zero_command[zero_command.index(str(thermal_output_path))] = str(zero_output_path)
        zero_result = subprocess.run(
            zero_command, check=False, capture_output=True, text=True, env=environment
        )
        if zero_result.returncode != 0:
            raise RuntimeError(zero_result.stdout + zero_result.stderr)
        with h5py.File(zero_output_path, "r") as handle:
            assert handle.attrs["dust_ir_out_of_range_count"] == 0
            assert np.allclose(handle["thermal/cumulative_dust_heating_energy_erg_cm3"], 0.0)
            assert np.allclose(handle["thermal/dust_grain_temperature_k"], 0.0)
            assert np.allclose(handle["thermal/dust_ir_reemitted_energy_erg_cm3"], 0.0)
            assert np.allclose(handle["thermal/dust_ir_untracked_energy_erg_cm3"], 0.0)
            assert np.allclose(handle["sources/dust_ir_photon_rate_cm3_s"], 0.0)

        non_v3_output_path = work / "non-v3-thermal-output.h5"
        non_v3_command = list(thermal_command)
        non_v3_command[non_v3_command.index(str(physical_opacity_path))] = str(dust_metadata_path)
        non_v3_command[non_v3_command.index(str(thermal_output_path))] = str(non_v3_output_path)
        non_v3_result = subprocess.run(
            non_v3_command, check=False, capture_output=True, text=True, env=environment
        )
        assert non_v3_result.returncode != 0
        non_v3_message = non_v3_result.stdout + non_v3_result.stderr
        assert (
            "unbound reference dust closure cannot be used" in non_v3_message
            or "requires snrt_dust_opacity_v3" in non_v3_message
        )

    print("P5_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
