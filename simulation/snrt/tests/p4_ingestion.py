"""P4 canonical snapshot staging checks.

Run with: JAX_PLATFORMS=cpu .venv/bin/python tests/p4_ingestion.py
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np

from snrt_core.snapshot import (
    GridSpec,
    RamsesFieldMap,
    SourceCatalog,
    neutral_primordial_input,
    read_static_rt_input,
    resolve_dust_abundance,
    write_static_rt_input,
)


def main() -> None:
    shape = (4, 3, 2)
    density = np.full(shape, 2.0e-24)
    temperature = np.linspace(1.0e2, 2.0e4, np.prod(shape)).reshape(shape)
    sources = SourceCatalog(
        cell_index=np.array([[0, 0, 0], [3, 2, 1]]),
        photon_luminosity_s=np.array([[1.0e49, 2.0e48], [3.0e47, 0.0]]),
    )
    velocity = np.stack([np.full(shape, 1.0e5), np.full(shape, -2.0e5), np.full(shape, 3.0e5)])
    metallicity = np.full(shape, 0.02)
    dust_to_metal = np.full(shape, 0.5)
    x_h2 = np.full(shape, 0.1)
    cell_level = np.full(shape, 15, dtype=np.int16)
    snapshot = neutral_primordial_input(
        GridSpec(cell_width_cm=3.085677581e18, left_edge_cm=np.array([1.0, 2.0, 3.0])),
        density,
        temperature,
        dust_relative_abundance=0.3,
        sources=sources,
        velocity_cm_s=velocity,
        metallicity_solar=metallicity,
        dust_to_metal=dust_to_metal,
        x_h2=x_h2,
        cell_level=cell_level,
        x_hii=0.1,
        x_heii=0.2,
        x_heiii=0.3,
    )
    with TemporaryDirectory() as temporary_directory:
        path = Path(temporary_directory) / "static_rt_input.h5"
        write_static_rt_input(path, snapshot)
        restored = read_static_rt_input(path)
        with h5py.File(path, "r") as handle:
            assert int(handle.attrs["format_version"]) == 3
        with h5py.File(path, "a") as handle:
            del handle["gas"].attrs["dust_relative_abundance_origin"]
        try:
            read_static_rt_input(path)
        except ValueError as error:
            assert "lacks the required dust_relative_abundance_origin" in str(error)
        else:
            raise AssertionError("legacy non-zero-dust input without origin was accepted")
    assert restored.shape == shape
    assert np.array_equal(restored.temperature_k, temperature)
    assert np.allclose(restored.hydrogen_number_density_cm3, density * 0.76 / 1.67262192369e-24)
    assert np.allclose(restored.helium_number_density_cm3, density * 0.24 / (4.0 * 1.67262192369e-24))
    assert restored.sources is not None
    assert np.array_equal(restored.sources.cell_index, sources.cell_index)
    assert np.array_equal(restored.sources.photon_luminosity_s, sources.photon_luminosity_s)
    assert np.array_equal(restored.velocity_cm_s, velocity)
    assert np.array_equal(restored.metallicity_solar, metallicity)
    assert np.array_equal(restored.dust_to_metal, dust_to_metal)
    assert restored.dust_relative_abundance_origin == "direct"
    assert np.array_equal(restored.x_h2, x_h2)
    assert np.array_equal(restored.cell_level, cell_level)
    assert np.allclose(restored.x_hii, 0.1)
    assert np.allclose(restored.x_heii, 0.2)
    assert np.allclose(restored.x_heiii, 0.3)
    restored.validate_production_contract(require_sources=True)

    resolved, origin = resolve_dust_abundance(
        shape,
        dust_relative_abundance=None,
        metallicity_solar=metallicity,
        dust_to_metal=dust_to_metal,
    )
    assert origin == "metallicity_solar_times_dust_to_metal"
    assert np.array_equal(resolved, metallicity * dust_to_metal)
    for keyword in ("metallicity_solar", "dust_to_metal"):
        partial = {
            "dust_relative_abundance": None,
            "metallicity_solar": metallicity if keyword == "metallicity_solar" else None,
            "dust_to_metal": dust_to_metal if keyword == "dust_to_metal" else None,
        }
        try:
            resolve_dust_abundance(shape, **partial)
        except ValueError as error:
            assert "supplied together" in str(error)
        else:
            raise AssertionError(f"partial dust mapping accepted for {keyword}")
    for keyword in ("metallicity_solar", "dust_to_metal"):
        try:
            RamsesFieldMap(density="density", temperature="temperature", **{keyword: "field"})
        except ValueError as error:
            assert "supplied together" in str(error)
        else:
            raise AssertionError(f"partial RAMSES dust field map accepted for {keyword}")

    derived_snapshot = neutral_primordial_input(
        GridSpec(cell_width_cm=3.085677581e18, left_edge_cm=np.array([1.0, 2.0, 3.0])),
        density,
        temperature,
        dust_relative_abundance=0.01,
        metallicity_solar=metallicity,
        dust_to_metal=dust_to_metal,
        dust_relative_abundance_origin="metallicity_solar_times_dust_to_metal",
    )
    with TemporaryDirectory() as temporary_directory:
        derived_path = Path(temporary_directory) / "derived_static_rt_input.h5"
        write_static_rt_input(derived_path, derived_snapshot)
        derived_restored = read_static_rt_input(derived_path)
    assert derived_restored.dust_relative_abundance_origin == "metallicity_solar_times_dust_to_metal"
    assert np.allclose(
        derived_restored.dust_relative_abundance,
        derived_restored.metallicity_solar * derived_restored.dust_to_metal,
    )
    try:
        neutral_primordial_input(
            GridSpec(cell_width_cm=3.085677581e18, left_edge_cm=np.array([1.0, 2.0, 3.0])),
            density,
            temperature,
            dust_relative_abundance=0.011,
            metallicity_solar=metallicity,
            dust_to_metal=dust_to_metal,
            dust_relative_abundance_origin="metallicity_solar_times_dust_to_metal",
        )
    except ValueError as error:
        assert "disagrees" in str(error)
    else:
        raise AssertionError("inconsistent derived dust abundance was accepted")
    print("P4_INGESTION_OK format=v3 shape=4x3x2 sources=2 groups=2")


if __name__ == "__main__":
    main()
