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
from tools.build_draine_dust_thermal import build_thermal_metadata, build_native_reference_namelist  # noqa: E402
from tools.build_dl01_dust_material import build_material, KB, AMU, CARBON_AMU, SILICATE_AMU  # noqa: E402


EDGES = np.asarray((0.01, 1.0, 5.6, 11.2, 13.6, 24.59, 54.42, 500.0, 2000.0, 10000.0))


def main() -> int:
    jax.config.update("jax_enable_x64", True)
    knots = np.geomspace(5., 300., 160)
    bulk = build_material(knots, .3)
    bulk64 = build_material(knots, .3, order=64)
    np.testing.assert_allclose(bulk["internal_energy_erg_g"], bulk64["internal_energy_erg_g"], rtol=1e-10)
    for fraction, atom_mass in ((1., CARBON_AMU), (0., SILICATE_AMU)):
        hot = build_material(np.array([1e8, 2e8]), fraction)
        c = np.diff(hot["internal_energy_erg_g"])[0] / 1e8
        np.testing.assert_allclose(c, 3 * KB / (atom_mass * AMU), rtol=1e-8)
    print("DL01_MATERIAL_PASS monotonic=1 quadrature=1 high_T_three_modes=1")
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
        ledger = ROOT / "data" / "p4_pilot_agn_photon_ledger.json"
        constant_native = build_native_reference_namelist(
            source, ledger, edges_path, np.geomspace(5., 300., 64), 10., 1e-24)
        assert "contract_version=3" in constant_native
        material_path = Path(directory) / "material.json"
        material = dict(schema="snrt_dust_material_energy_v1", source_id="synthetic_cubic_control",
                        source_url="urn:snrt:synthetic-test", composition="synthetic; NOT physical dust",
                        energy_zero="U(0)=0; no zero-point energy",
                        temperature_k=[5., 10., 20., 100., 300.],
                        internal_energy_erg_g=[125., 1000., 8000., 1e6, 27e6])
        material_path.write_text(json.dumps(material))
        native = build_native_reference_namelist(
            source, ledger, edges_path, np.geomspace(5., 300., 64), 10., None,
            material_path=material_path)
        assert "contract_version=4" in native and "internal_energy_per_h_erg_input=" in native
        assert "material_sha256=" in native and "approval_id=''" in native
        assert "ntemperature_input=67" in native  # preserves material knots 10,20,100
        scattered = build_native_reference_namelist(
            source, ledger, edges_path, np.geomspace(5., 300., 64), 10., None,
            material_path=material_path, scattering="isotropic_elastic")
        assert "scattering_input=" not in native
        assert "scattering_model='isotropic_elastic'" in scattered and "scattering_input=" in scattered
        assert "NOT the measured Draine phase function" in scattered
        exchanged = build_native_reference_namelist(
            source, ledger, edges_path, np.geomspace(5.,300.,64), 10., None,
            material_path=material_path, gas_exchange="hydrogen_accommodation",
            collision_area=3.495e-22, accommodation=.5)
        assert "gas_exchange_model=" not in native
        assert "gas_exchange_model='hydrogen_accommodation'" in exchanged
        assert "collision_area_per_h_input=" in exchanged
        for area, alpha in ((None,.5),(0.,.5),(1e-21,None),(1e-21,1.01),(1e-21,-.1)):
            try:
                build_native_reference_namelist(source, ledger, edges_path,
                    np.geomspace(5.,300.,64),10.,None,material_path=material_path,
                    gas_exchange="hydrogen_accommodation",collision_area=area,accommodation=alpha)
            except ValueError:
                pass
            else:
                raise AssertionError("invalid gas exchange parameters accepted")
        for option in ("isotropic_elastic", "unrecognized"):
            try:
                build_native_reference_namelist(source, ledger, edges_path,
                    np.geomspace(5., 300., 64), 10., 1e-24, scattering=option)
            except ValueError:
                pass
            else:
                raise AssertionError("invalid v3/scattering selection accepted")
        for key, value in (("internal_energy_erg_g", [1., 1., 2., 3., 4.]),
                           ("energy_zero", "unknown"), ("source_id", "bad'quote"),
                           ("temperature_k", [6., 10., 20., 100., 300.])):
            material_path.write_text(json.dumps(dict(material, **{key: value})))
            try:
                build_native_reference_namelist(source, ledger, edges_path,
                                               np.geomspace(5., 300., 64), 10., None,
                                               material_path=material_path)
            except ValueError:
                pass
            else:
                raise AssertionError(f"invalid material accepted: {key}")
        print("DUST_MATERIAL_EXPORT_PASS v3=1 v4=1 knots_preserved=1 failures=4 reference_only=1")
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
