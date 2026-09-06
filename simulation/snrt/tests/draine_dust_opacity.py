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

from snrt_core.dust import read_dust_opacity_metadata


TOOL = ROOT / "tools" / "build_draine_dust_opacity.py"
SPEC = importlib.util.spec_from_file_location("build_draine_dust_opacity", TOOL)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load Draine opacity builder")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def main() -> int:
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

    expected_staged_sha256 = "61350545eea164c8db94ff830abd7f9e57cd7efc6f1e36f389d61627d364b9da"
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


if __name__ == "__main__":
    raise SystemExit(main())
