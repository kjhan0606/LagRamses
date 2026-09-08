#!/usr/bin/env python3
"""Validate the staged Draine table and its P0 dust-opacity closure."""

from __future__ import annotations

import importlib.util
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
STAGED_SIDECAR = ROOT.parents[1] / "external" / "draine_wd01_rv31" / "p0_dust_opacity_rv31_photon_index1_scattering.json"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

def main() -> int:
    # The independent grain/Mie data check does not need the historical JAX
    # P0 runner. Load that dependency only for the original sidecar test.
    from snrt_core.dust import read_dust_opacity_metadata
    tool = ROOT / "tools/build_draine_dust_opacity.py"
    spec = importlib.util.spec_from_file_location("build_draine_dust_opacity", tool)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load Draine opacity builder")
    MODULE = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(MODULE)
    source = ROOT.parents[1] / "external" / "draine_wd01_rv31" / "kext_albedo_WD_MW_3.1_60_D03.all"
    if not source.is_file():
        raise FileNotFoundError(source)
    edges_path = ROOT / "config" / "p0_photon_group_edges_ev.txt"
    metadata = MODULE.build_opacity_metadata(source, edges_path, photon_index=1.0)
    edges = np.asarray(metadata["group_edges_ev"], dtype=np.float64)
    opacity = np.asarray(metadata["absorption_cross_section_per_h_cm2"], dtype=np.float64)
    weighted_energy = np.asarray(metadata["absorption_weighted_energy_ev"], dtype=np.float64)
    assert edges.size == 10
    assert opacity.shape == (9,)
    assert np.isfinite(opacity).all() and np.all(opacity > 0.0)
    assert np.isfinite(weighted_energy).all()
    assert np.all(weighted_energy >= edges[:-1])
    assert np.all(weighted_energy <= edges[1:])
    assert metadata["source_table"]["row_count"] == 812
    assert metadata["source_table"]["absorption_consistency_max_relative_error"] < 1.0e-2

    with __import__("tempfile").TemporaryDirectory(prefix="draine-opacity-test-") as directory:
        output = Path(directory) / "dust.json"
        output.write_text(json.dumps(metadata) + "\n", encoding="utf-8")
        closure = read_dust_opacity_metadata(output, expected_group_edges_ev=edges)
    assert np.allclose(closure.absorption_cross_section_per_h_cm2, opacity)

    scattering_metadata = MODULE.build_opacity_metadata(
        source,
        edges_path,
        photon_index=1.0,
        include_scattering=True,
    )
    scattering = np.asarray(
        scattering_metadata["scattering_cross_section_per_h_cm2"], dtype=np.float64
    )
    cosine = np.asarray(scattering_metadata["scattering_angle_cosine"], dtype=np.float64)
    cosine_squared = np.asarray(
        scattering_metadata["scattering_angle_cosine_squared"], dtype=np.float64
    )
    assert scattering_metadata["schema"] == "snrt_dust_opacity_v3"
    assert scattering_metadata["status"] == "reference_scattering_control"
    assert scattering.shape == (9,) and np.all(scattering >= 0.0)
    assert np.all((cosine >= -1.0) & (cosine <= 1.0))
    assert np.all(cosine_squared + 1.0e-4 >= cosine**2)
    assert "moment_inequality_max_violation" in scattering_metadata["source_table"]
    assert scattering_metadata["source_table"]["moment_inequality_max_violation"] <= 1.0e-4
    factors = scattering_metadata["isotropic_candidate_momentum_overestimate_factor"]
    unbounded = scattering_metadata["isotropic_candidate_momentum_bound_unbounded"]
    assert len(factors) == len(unbounded) == 9
    assert any(unbounded)
    assert all(value is None for value, flag in zip(factors, unbounded, strict=True) if flag)
    assert all(float(value) >= 1.0 for value, flag in zip(factors, unbounded, strict=True) if not flag)
    with __import__("tempfile").TemporaryDirectory(prefix="draine-scattering-test-") as directory:
        output = Path(directory) / "dust-v3.json"
        output.write_text(json.dumps(scattering_metadata) + "\n", encoding="utf-8")
        scattering_closure = read_dust_opacity_metadata(
            output,
            expected_group_edges_ev=edges,
        )
    assert scattering_closure.scattering_phase_function == "phase_isotropic_candidate"
    assert np.allclose(scattering_closure.scattering_cross_section_per_h_cm2, scattering)

    expected_staged_sha256 = "7521ef988a47b590f375f49cdedf375109f5ee306968e54749b38f5e43a1faa8"
    digest = hashlib.sha256(STAGED_SIDECAR.read_bytes()).hexdigest()
    assert digest == expected_staged_sha256
    staged = json.loads(STAGED_SIDECAR.read_text(encoding="utf-8"))
    staged_edges = np.asarray(staged["group_edges_ev"], dtype=np.float64)
    staged_closure = read_dust_opacity_metadata(STAGED_SIDECAR, expected_group_edges_ev=staged_edges)
    assert staged_closure.schema == "snrt_dust_opacity_v3"
    assert staged_closure.scattering_phase_function == "phase_isotropic_candidate"
    rebuilt_staged = MODULE.build_opacity_metadata(
        source,
        edges_path,
        photon_index=1.0,
        include_scattering=True,
    )
    assert staged == rebuilt_staged
    print(
        "DRAINE_DUST_OPACITY_TEST_OK "
        f"rows={metadata['source_table']['row_count']} groups={opacity.size} "
        f"scattering_groups={np.count_nonzero(scattering > 0.0)} "
        f"max_consistency={metadata['source_table']['absorption_consistency_max_relative_error']:.3e}"
    )
    return 0


def check_d03() -> int:
    """Offline Mie data checks, reusing this existing optical test entry point."""
    tool = ROOT / "tools/build_d03_grain_optics.py"
    spec = importlib.util.spec_from_file_location("d03", tool)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    import miepython
    # Independent published Mie amplitude case (Wiscombe/Prahl documentation).
    an, bn = miepython.an_bn(4/3, 50)
    assert abs(an[0] - (.531105889295 - .499031485631j)) < 2e-11
    assert abs(bn[0] - (.791924475935 - .405931152229j)) < 2e-11
    m, x = 1.5-.01j, 1e-5
    ext, sca, _, _ = miepython.efficiencies_mx(m, x)
    polar = (m*m-1)/(m*m+2)
    assert abs((ext-sca)/(-4*x*polar.imag)-1) < 1e-8
    assert abs(sca/(8/3*x**4*abs(polar)**2)-1) < 1e-8
    manifest = json.loads((ROOT / "data/dust_d03_optics_generation_v1.json").read_text())
    tables = {}
    for name, identity in zip(module.FILES, manifest["sources"], strict=True):
        raw = (ROOT / "data/draine_d03" / name).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == identity["sha256"]
        data = module.read_index(raw, name)
        module.index_at(data, np.array([1.3e-4, 2000., 10000.]))
        for bad in (np.array([20000.]), np.array([1e-8]), np.array([np.nan])):
            try:
                module.index_at(data, bad)
            except ValueError:
                pass
            else:
                raise AssertionError("out-of-range spectral input was accepted")
        tables[name] = data
    energies = np.array([.1, 17.662918001204492, 4023.594574013186, 10000.])
    q = module.grain_efficiencies(tables, energies)
    assert np.all(q[:2] > 0) and np.all(abs(q[2]) <= 1)
    # Graphite orientation averaging must preserve scattering-weighted moments.
    a = module.efficiencies(tables["callindex.out_CpaD03_0.01"], energies, .01)
    b = module.efficiencies(tables["callindex.out_CpeD03_0.01"], energies, .01)
    np.testing.assert_allclose(q[1, :, 0]*q[2, :, 0], (a[1]*a[2]+2*b[1]*b[2])/3, rtol=2e-15)
    assert np.all(q[2, -1] > .999)
    assert np.all(q[0, -1] > 0)  # never drop the hard-X group
    compiled = ROOT.parents[1] / "patch/lagRamses/dust_d03_optics_data.inc"
    assert hashlib.sha256(compiled.read_bytes()).hexdigest() == manifest["compiled_sha256"]
    assert hashlib.sha256(tool.read_bytes()).hexdigest() == manifest["generator_sha256"]
    print("D03_MIE_AMPLITUDES_RAYLEIGH_SPECTRAL_SUPPORT_ANGULAR_MOMENT_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(check_d03() if sys.argv[1:] == ["--d03"] else main())
