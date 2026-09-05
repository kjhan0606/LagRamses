#!/usr/bin/env python3
"""Exercise the non-zero-dust sidecar through the complete P4 runner."""

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
    edges = agn_ledger._read_group_edges(agn_ledger.DEFAULT_P0_GROUP_EDGES)
    photon_per_lbol, mean_energy, closure = agn_ledger._group_conversion(0.1, edges)
    # Keep this one-step fixture in a modest-ionization regime so the runner's
    # coupled chemistry ledger remains a useful integration check.
    source_luminosity = 1.0e35 * photon_per_lbol
    group_metadata = {
        "schema": "agn_photon_source_ledger_v1",
        "source_kind": "agn",
        "groups": [
            {
                "index": index,
                "energy_interval_ev": [float(edges[index]), float(edges[index + 1])],
                "photon_weighted_mean_energy_ev": float(mean_energy[index]),
                "total_photon_rate_s": float(source_luminosity[index]),
            }
            for index in range(edges.size - 1)
        ],
        "group_edges_ev": edges.tolist(),
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
    dust_metadata = {
        "schema": "snrt_dust_opacity_v1",
        "group_edges_ev": edges.tolist(),
        "absorption_cross_section_per_h_cm2": [1.0e-21] * (edges.size - 1),
        "absorption_weighted_energy_ev": np.sqrt(edges[:-1] * edges[1:]).tolist(),
        "reference_mixture": "synthetic runner fixture",
        "opacity_source": "synthetic runner fixture",
        "spectral_weighting": "group geometric-mean fixture",
    }

    with TemporaryDirectory(prefix="p4-dust-runner-test-") as directory:
        work = Path(directory)
        input_path = work / "input.h5"
        photon_metadata_path = work / "agn.json"
        dust_metadata_path = work / "dust.json"
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
        photon_metadata_path.write_text(json.dumps(group_metadata, indent=2) + "\n", encoding="utf-8")
        dust_metadata_path.write_text(json.dumps(dust_metadata, indent=2) + "\n", encoding="utf-8")
        environment = os.environ.copy()
        environment["PYTHONPATH"] = str(ROOT)
        environment["JAX_PLATFORMS"] = "cpu"
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "p4_run_transport_pilot.py"),
                "--input",
                str(input_path),
                "--photon-metadata",
                str(photon_metadata_path),
                "--dust-opacity-metadata",
                str(dust_metadata_path),
                "--steps",
                "1",
                "--sn-order",
                "4",
                "--output",
                str(output_path),
            ],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        assert "P4_TRANSPORT_PILOT_OK" in result.stdout
        with h5py.File(output_path, "r") as handle:
            dust_heating = np.asarray(handle["rates/dust_heating_erg_cm3_s"])
            dust_momentum = np.asarray(handle["rates/dust_momentum_rate_dyn_cm3"])
            dust_absorbed = np.asarray(handle["diagnostics/cumulative_dust_absorbed_photons_cm3"])
            assert handle.attrs["dust_model"] == "metadata"
            assert handle.attrs["dust_opacity_schema"] == "snrt_dust_opacity_v1"
            assert handle.attrs["dust_binding_status"] == "reference_control"
            assert handle.attrs["dust_opacity_metadata_sha256"]
            assert handle.attrs["dust_payload_sha256"] == ""
            assert handle.attrs["dust_source_table_sha256"] == ""
            assert handle.attrs["dust_builder_sha256"] == ""
            assert handle.attrs["source_sed_identity"] == ""
            assert np.asarray(handle["group_energy_ev"]).shape == (9,)
            assert dust_momentum.shape == (3, *shape)
            assert dust_absorbed.shape == (9, *shape)
            assert np.all(np.isfinite(dust_heating)) and np.max(dust_heating) > 0.0
            assert np.all(np.isfinite(dust_momentum))

    print("P4_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
