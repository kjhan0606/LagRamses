"""P4 canonical snapshot staging checks.

Run with: JAX_PLATFORMS=cpu .venv/bin/python tests/p4_ingestion.py
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from snrt_core.snapshot import GridSpec, SourceCatalog, neutral_primordial_input, read_static_rt_input, write_static_rt_input


def main() -> None:
    shape = (4, 3, 2)
    density = np.full(shape, 2.0e-24)
    temperature = np.linspace(1.0e2, 2.0e4, np.prod(shape)).reshape(shape)
    sources = SourceCatalog(
        cell_index=np.array([[0, 0, 0], [3, 2, 1]]),
        photon_luminosity_s=np.array([[1.0e49, 2.0e48], [3.0e47, 0.0]]),
    )
    snapshot = neutral_primordial_input(
        GridSpec(cell_width_cm=3.085677581e18, left_edge_cm=np.array([1.0, 2.0, 3.0])),
        density,
        temperature,
        dust_relative_abundance=0.3,
        sources=sources,
    )
    with TemporaryDirectory() as temporary_directory:
        path = Path(temporary_directory) / "static_rt_input.h5"
        write_static_rt_input(path, snapshot)
        restored = read_static_rt_input(path)
    assert restored.shape == shape
    assert np.array_equal(restored.temperature_k, temperature)
    assert np.allclose(restored.hydrogen_number_density_cm3, density * 0.76 / 1.67262192369e-24)
    assert np.allclose(restored.helium_number_density_cm3, density * 0.24 / (4.0 * 1.67262192369e-24))
    assert restored.sources is not None
    assert np.array_equal(restored.sources.cell_index, sources.cell_index)
    assert np.array_equal(restored.sources.photon_luminosity_s, sources.photon_luminosity_s)
    print("P4_INGESTION_OK shape=4x3x2 sources=2 groups=2")


if __name__ == "__main__":
    main()
