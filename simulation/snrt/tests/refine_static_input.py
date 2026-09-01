#!/usr/bin/env python3
"""Validate conservative factor-of-two static-input refinement."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from snrt_core.snapshot import GridSpec, SourceCatalog, StaticRTInput


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "refine_static_rt_input.py"


def _load_tool():
    spec = importlib.util.spec_from_file_location("refine_static_rt_input", TOOL)
    if spec is None or spec.loader is None:
        raise RuntimeError("unable to load refinement tool")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    refine = _load_tool().refine_static_rt_input
    shape = (2, 3, 1)
    scalar = np.arange(np.prod(shape), dtype=np.float64).reshape(shape) + 1.0
    snapshot = StaticRTInput(
        grid=GridSpec(8.0, np.asarray((1.0, 2.0, 3.0))),
        hydrogen_number_density_cm3=scalar,
        helium_number_density_cm3=0.1 * scalar,
        temperature_k=1.0e4 + scalar,
        dust_relative_abundance=np.zeros(shape),
        x_hii=np.zeros(shape),
        x_heii=np.zeros(shape),
        x_heiii=np.zeros(shape),
        sources=SourceCatalog(
            np.asarray(((0, 0, 0), (1, 2, 0))),
            np.asarray(((1.0e49, 2.0e49), (3.0e49, 4.0e49))),
        ),
    )
    refined = refine(snapshot, 2)
    assert refined.shape == (4, 6, 2)
    assert refined.grid.cell_width_cm == 4.0
    assert np.array_equal(refined.grid.left_edge_cm, snapshot.grid.left_edge_cm)
    assert np.array_equal(refined.hydrogen_number_density_cm3[::2, ::2, ::2], scalar)
    assert np.all(refined.hydrogen_number_density_cm3[0:2, 0:2, :] == scalar[0, 0, 0])
    assert len(refined.sources.cell_index) == 16
    assert np.allclose(
        refined.sources.photon_luminosity_s.sum(axis=0),
        snapshot.sources.photon_luminosity_s.sum(axis=0),
        rtol=1.0e-13,
        atol=0.0,
    )
    assert np.all(refined.sources.cell_index[:8] >= 0)
    assert np.all(refined.sources.cell_index[:8] < np.asarray(refined.shape))
    assert np.allclose(refined.sources.photon_luminosity_s[:8], snapshot.sources.photon_luminosity_s[0] / 8.0)
    print("REFINE_STATIC_INPUT_TEST_OK factor=2 source_luminosity_conserved=true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
