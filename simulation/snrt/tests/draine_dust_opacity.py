#!/usr/bin/env python3
"""Validate the staged Draine table and its P0 dust-opacity closure."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
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
    print(
        "DRAINE_DUST_OPACITY_TEST_OK "
        f"rows={metadata['source_table']['row_count']} groups={opacity.size} "
        f"max_consistency={metadata['source_table']['absorption_consistency_max_relative_error']:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
