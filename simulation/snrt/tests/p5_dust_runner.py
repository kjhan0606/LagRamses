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
            assert np.max(np.asarray(handle["diagnostics/cumulative_dust_absorbed_photons_cm3"])) > 0.0
            assert np.max(np.asarray(handle["thermal/cumulative_dust_heating_energy_erg_cm3"])) > 0.0
            assert handle.attrs["primary_absorption_closure_relative_error"] <= 1.0e-5

    print("P5_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
