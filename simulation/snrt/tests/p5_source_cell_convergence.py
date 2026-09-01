"""Check the recorded full-duration source-cell convergence controls."""

from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "p5_validation"
PATHS = (
    DATA / "p5_coeval_s8_courant0p1_sub1_6p37myr_source_limit0p25_f64.h5",
    DATA / "p5_coeval_s8_courant0p05_sub1_6p37myr_source_limit0p25_f64.h5",
)


def main() -> None:
    fields = []
    source_subcycles = []
    for path in PATHS:
        with h5py.File(path, "r") as handle:
            assert bool(handle.attrs["validation_passed"])
            assert float(handle.attrs["source_cell_photons_per_neutral_target"]) == 0.25
            source_subcycles.append(int(handle.attrs["source_cell_subcycles"]))
            fields.append(np.asarray(handle["ionization/x_hii"][...], dtype=np.float64))

    assert source_subcycles == [50, 25]
    difference = np.abs(fields[0] - fields[1])
    assert difference.shape == (32, 32, 32)
    assert float(difference.mean()) < 1.0e-5
    assert float(difference.max()) < 3.0e-2
    assert tuple(np.unravel_index(np.argmax(difference), difference.shape)) == (14, 17, 19)
    assert float(difference[14, 17, 19]) > 1.0e-2
    print(
        "P5_SOURCE_CELL_CONVERGENCE_OK "
        f"mean_abs={difference.mean():.6g} max={difference.max():.6g} "
        "source_cell_gate=open"
    )


if __name__ == "__main__":
    main()
