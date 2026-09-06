#!/usr/bin/env python3
"""Check the DUST-2 thermal sidecar, operator, and fail-closed paths."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import jax
import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.dust import (  # noqa: E402
    EV_ERG,
    dust_thermal_model_from_metadata,
    evaluate_dust_thermal,
    read_dust_opacity_metadata,
    read_dust_thermal_metadata,
)
from tools.build_draine_dust_thermal import build_thermal_metadata  # noqa: E402


EDGES = np.asarray((0.01, 1.0, 5.6, 11.2, 13.6, 24.59, 54.42, 500.0, 2000.0, 10000.0))


def main() -> int:
    jax.config.update("jax_enable_x64", True)
    source = ROOT.parents[1] / "external" / "draine_wd01_rv31" / "kext_albedo_WD_MW_3.1_60_D03.all"
    edges_path = ROOT / "config" / "p0_photon_group_edges_ev.txt"
    opacity_path = (
        ROOT.parents[1]
        / "external"
        / "draine_wd01_rv31"
        / "p0_dust_opacity_rv31_photon_index1_scattering.json"
    )
    opacity = read_dust_opacity_metadata(opacity_path, expected_group_edges_ev=EDGES)
    metadata = build_thermal_metadata(
        source,
        edges_path,
        temperature_grid_k=np.geomspace(5.0, 300.0, 64),
    )

    with TemporaryDirectory(prefix="dust-thermal-test-") as directory:
        path = Path(directory) / "thermal.json"
        path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        closure = read_dust_thermal_metadata(
            path,
            expected_group_edges_ev=EDGES,
            expected_group_edges_sha256=opacity.group_edges_sha256,
            expected_source_table_sha256=opacity.source_table_sha256,
            expected_dust_mass_per_h_g=opacity.dust_mass_per_h_g,
        )
        assert closure.ir_group_indices.tolist() == [0]
        assert np.allclose(
            closure.ir_energy_fraction.sum(axis=1) + closure.untracked_energy_fraction,
            1.0,
            rtol=0.0,
            atol=1.0e-10,
        )
        assert np.all(np.diff(closure.emitted_power_per_h_erg_s) > 0.0)
        assert np.all(closure.ir_mean_photon_energy_ev[:, 0] >= EDGES[0])
        assert np.all(closure.ir_mean_photon_energy_ev[:, 0] <= EDGES[1])

        model = dust_thermal_model_from_metadata(path, dtype=jnp.float64, expected_group_edges_ev=EDGES)
        heating = jnp.asarray(((1.0e-25, 2.0e-25), (0.0, 1.0e-25)), dtype=jnp.float64)
        result = evaluate_dust_thermal(
            model,
            heating,
            jnp.ones_like(heating),
            jnp.ones_like(heating),
            13.1,
        )
        tracked = np.asarray(result.reemitted_energy_rate)
        untracked = np.asarray(result.untracked_energy_rate)
        assert not np.asarray(result.out_of_range).any()
        assert np.allclose(tracked + untracked, np.asarray(heating), rtol=0.0, atol=3.0e-8 * 1.0e-25)
        assert np.asarray(result.grain_temperature_k)[1, 0] == 0.0
        assert np.asarray(result.reemitted_energy_rate)[1, 0] == 0.0
        assert np.all(np.diff(np.asarray(result.grain_temperature_k)[0]) >= 0.0)
        assert np.all(np.isfinite(np.asarray(result.ir_photon_rate)))

        multi_metadata = build_thermal_metadata(
            source,
            edges_path,
            ir_group_indices=[0, 1],
            temperature_grid_k=np.geomspace(5.0, 300.0, 64),
        )
        multi_path = Path(directory) / "thermal-multi.json"
        multi_path.write_text(json.dumps(multi_metadata) + "\n", encoding="utf-8")
        multi_model = dust_thermal_model_from_metadata(
            multi_path, dtype=jnp.float64, expected_group_edges_ev=EDGES
        )
        multi_result = evaluate_dust_thermal(
            multi_model,
            jnp.asarray([[1.0e-18]], dtype=jnp.float64),
            jnp.ones((1, 1), dtype=jnp.float64),
            jnp.ones((1, 1), dtype=jnp.float64),
            13.1,
        )
        assert not np.asarray(multi_result.out_of_range).any()
        multi_temperature = float(np.asarray(multi_result.grain_temperature_k)[0, 0])
        multi_photon_rate = np.asarray(multi_result.ir_photon_rate)[:, 0, 0]
        expected_multi_energy = 0.0
        for column, group in enumerate((0, 1)):
            mean_energy = np.interp(
                np.log(multi_temperature),
                np.log(np.asarray(multi_metadata["temperature_k"])),
                np.asarray(multi_metadata["ir_mean_photon_energy_ev"])[:, column],
            )
            expected_multi_energy += multi_photon_rate[group] * mean_energy * EV_ERG
            assert multi_photon_rate[group] > 0.0
        assert np.allclose(
            expected_multi_energy,
            np.asarray(multi_result.reemitted_energy_rate)[0, 0],
            rtol=0.0,
            atol=3.0e-26,
        )
        assert np.allclose(
            expected_multi_energy + np.asarray(multi_result.untracked_energy_rate)[0, 0],
            1.0e-18,
            rtol=0.0,
            atol=3.0e-26,
        )
        assert np.allclose(multi_photon_rate[2:], 0.0)

        too_hot = evaluate_dust_thermal(
            model,
            jnp.asarray([[1.0e-8]], dtype=jnp.float64),
            jnp.ones((1, 1), dtype=jnp.float64),
            jnp.ones((1, 1), dtype=jnp.float64),
            13.1,
        )
        assert bool(np.asarray(too_hot.out_of_range)[0, 0])
        cold_background = evaluate_dust_thermal(
            model,
            jnp.asarray([[1.0e-25]], dtype=jnp.float64),
            jnp.ones((1, 1), dtype=jnp.float64),
            jnp.ones((1, 1), dtype=jnp.float64),
            2.0,
        )
        assert bool(np.asarray(cold_background.out_of_range)[0, 0])

        malformed = dict(metadata)
        malformed["untracked_energy_fraction"] = list(metadata["untracked_energy_fraction"])
        malformed["untracked_energy_fraction"][0] += 0.1
        malformed_path = Path(directory) / "malformed.json"
        malformed_path.write_text(json.dumps(malformed) + "\n", encoding="utf-8")
        try:
            read_dust_thermal_metadata(malformed_path, expected_group_edges_ev=EDGES)
        except ValueError:
            pass
        else:
            raise AssertionError("malformed thermal fractions were accepted")

        try:
            read_dust_thermal_metadata(path, expected_group_edges_ev=EDGES * (1.0 + 1.0e-6))
        except ValueError:
            pass
        else:
            raise AssertionError("mismatched thermal group edges were accepted")

        wrong_source = json.loads(json.dumps(metadata))
        wrong_source["source_table"]["sha256"] = "0" * 64
        wrong_source_path = Path(directory) / "wrong-source.json"
        wrong_source_path.write_text(json.dumps(wrong_source) + "\n", encoding="utf-8")
        try:
            read_dust_thermal_metadata(wrong_source_path, expected_group_edges_ev=EDGES)
        except ValueError:
            pass
        else:
            raise AssertionError("thermal source-table hash mismatch was accepted")

        wrong_mass = json.loads(json.dumps(metadata))
        wrong_mass["source_table"]["dust_mass_per_h_g"] *= 2.0
        wrong_mass_path = Path(directory) / "wrong-mass.json"
        wrong_mass_path.write_text(json.dumps(wrong_mass) + "\n", encoding="utf-8")
        try:
            read_dust_thermal_metadata(
                wrong_mass_path,
                expected_group_edges_ev=EDGES,
                expected_dust_mass_per_h_g=opacity.dust_mass_per_h_g,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("thermal dust-mass mismatch was accepted")

    print(
        "DUST_THERMAL_TEST_OK "
        f"temperatures={closure.temperature_k.size} ir_groups={closure.ir_group_indices.size} "
        f"source_sha256={closure.source_table_sha256}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
