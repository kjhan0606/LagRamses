#!/usr/bin/env python3
"""Check the short factor-2 refined-mesh source-cell control pair."""

from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "p5_validation"
PATHS = (
    DATA / "p5_refined2_point_courant0p1_sub1_0p5myr_source_limit0p25_f64.h5",
    DATA / "p5_refined2_point_courant0p05_sub1_0p5myr_source_limit0p25_f64.h5",
)


def main() -> int:
    fields = []
    source_subcycles = []
    for path in PATHS:
        with h5py.File(path, "r") as handle:
            assert bool(handle.attrs["validation_passed"])
            assert handle.attrs["source_deposition_mode"] == "point"
            assert float(handle.attrs["source_cell_photons_per_neutral_target"]) == 0.25
            source_subcycles.append(int(handle.attrs["source_cell_subcycles"]))
            field = np.asarray(handle["ionization/x_hii"][...], dtype=np.float64)
            assert field.shape == (64, 64, 64)
            fields.append(field)

    assert source_subcycles == [25, 13]
    difference = np.abs(fields[0] - fields[1])
    assert np.isfinite(difference).all()
    assert float(difference.mean()) < 2.0e-6
    assert float(difference.max()) < 1.0e-2
    assert tuple(np.unravel_index(np.argmax(difference), difference.shape)) == (28, 34, 38)

    coarse0 = fields[0].reshape(32, 2, 32, 2, 32, 2).mean(axis=(1, 3, 5))
    coarse1 = fields[1].reshape(32, 2, 32, 2, 32, 2).mean(axis=(1, 3, 5))
    coarse_difference = np.abs(coarse0 - coarse1)
    assert float(coarse_difference.max()) < 1.0e-2
    assert tuple(np.unravel_index(np.argmax(coarse_difference), coarse_difference.shape)) == (14, 17, 19)
    print(
        "P5_REFINED_MESH_CONVERGENCE_OK "
        f"mean_abs={difference.mean():.6g} max={difference.max():.6g} "
        f"coarse_block_max={coarse_difference.max():.6g} "
        "short_control=true science_gate=not_promoted"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
