"""B1 non-equilibrium primordial cooling and metal-atlas acceptance tests."""

from __future__ import annotations

from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory

import h5py
import jax
import jax.numpy as jnp
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.implicit import implicit_atomic_chemistry_with_transitions
from snrt_core.jax_thermal_atlas import from_numpy_atlas, net_rate as jax_metal_net_rate
from snrt_core.primordial import PrimordialState
from snrt_core.primordial_cooling import (
    collisional_ionization_coefficients,
    primordial_cooling_components,
    primordial_net_rate,
)
from snrt_core.thermal_atlas import read_thermal_atlas
from tools.build_metal_thermal_atlas import read_uvb_free_metal_cooling


ATLAS = PROJECT_ROOT / "data" / "production_metal_thermal_atlas_v2.h5"
LEGACY_ATLAS = PROJECT_ROOT / "data" / "p4_validation_thermal_atlas.h5"
NO_UVB_DATA = PROJECT_ROOT.parents[1] / "external" / "grackle" / "CloudyData_noUVB.h5"
HM2012_DATA = Path(
    "/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5"
)


def state(x_hii: float, x_heii: float = 0.01, x_heiii: float = 0.0) -> PrimordialState:
    return PrimordialState(
        n_hydrogen=jnp.asarray(1.0, dtype=jnp.float64),
        n_helium=jnp.asarray(0.079, dtype=jnp.float64),
        x_hydrogen_ii=jnp.asarray(x_hii, dtype=jnp.float64),
        x_helium_ii=jnp.asarray(x_heii, dtype=jnp.float64),
        x_helium_iii=jnp.asarray(x_heiii, dtype=jnp.float64),
    )


def main() -> None:
    jax.config.update("jax_enable_x64", True)
    coefficient = collisional_ionization_coefficients(jnp.asarray(1.0e5, dtype=jnp.float64))
    assert np.isclose(float(coefficient.hydrogen_i), 3.824258925061223e-9, rtol=2.0e-13)
    assert np.isclose(float(coefficient.helium_i), 4.788994257071327e-10, rtol=2.0e-13)
    assert np.isclose(float(coefficient.helium_ii), 4.508883088873742e-12, rtol=2.0e-13)

    mostly_neutral = state(0.01)
    mostly_ionized = state(0.99, x_heii=0.01)
    neutral_rate = float(primordial_net_rate(mostly_neutral, jnp.asarray(1.0e5), 0.2085))
    ionized_rate = float(primordial_net_rate(mostly_ionized, jnp.asarray(1.0e5), 0.2085))
    assert neutral_rate < 0.0 and ionized_rate < 0.0
    assert not np.isclose(neutral_rate, ionized_rate, rtol=0.1, atol=0.0)

    cmb_temperature = 2.73 / 0.2085
    below_cmb = primordial_cooling_components(mostly_ionized, jnp.asarray(0.8 * cmb_temperature), 0.2085)
    above_cmb = primordial_cooling_components(mostly_ionized, jnp.asarray(2.0 * cmb_temperature), 0.2085)
    assert float(below_cmb.compton_net_heating) > 0.0
    assert float(above_cmb.compton_net_heating) < 0.0

    initial = state(0.1)
    updated = implicit_atomic_chemistry_with_transitions(
        initial,
        jnp.asarray(1.0e5, dtype=jnp.float64),
        1.0e12,
        48,
    )
    next_state = updated[0]
    recombination_h, recombination_heii, recombination_heiii = map(float, updated[1:4])
    ionization_h, ionization_hei, ionization_heii = map(float, updated[4:7])
    hydrogen_change = float(next_state.x_hydrogen_ii - initial.x_hydrogen_ii)
    helium_ii_change = 0.079 * float(next_state.x_helium_ii - initial.x_helium_ii)
    helium_iii_change = 0.079 * float(next_state.x_helium_iii - initial.x_helium_iii)
    assert abs(ionization_h - recombination_h - hydrogen_change) < 2.0e-12
    assert abs(ionization_hei - recombination_heii - helium_ii_change - helium_iii_change) < 2.0e-12
    assert abs(ionization_heii - recombination_heiii - helium_iii_change) < 2.0e-12

    atlas = read_thermal_atlas(ATLAS)
    assert atlas.provenance["thermal_component"] == "metal_only"
    assert atlas.provenance["uv_background_included"] == "false"
    density = 10.0 ** atlas.log_hydrogen_number_density_cm3[12]
    temperature = 10.0 ** atlas.log_temperature_k[100]
    solar_rate = float(atlas.net_rate(atlas.scale_factor[0], temperature, density, 1.0))
    metallicities = np.asarray([0.0, 10.0**-2.35, 0.37, np.sqrt(10.0)], dtype=np.float64)
    host_rates = np.asarray(
        atlas.net_rate(atlas.scale_factor[0], temperature, density, metallicities),
        dtype=np.float64,
    )
    assert solar_rate < 0.0
    assert np.allclose(host_rates, solar_rate * metallicities, rtol=2.0e-12, atol=0.0)
    device_rates = np.asarray(
        jax.device_get(
            jax_metal_net_rate(
                from_numpy_atlas(atlas, dtype=jnp.float64),
                atlas.scale_factor[0],
                jnp.full(metallicities.shape, temperature, dtype=jnp.float64),
                jnp.full(metallicities.shape, density, dtype=jnp.float64),
                jnp.asarray(metallicities, dtype=jnp.float64),
            )
        )
    )
    assert np.allclose(device_rates, solar_rate * metallicities, rtol=2.0e-12, atol=0.0)
    try:
        atlas.net_rate(atlas.scale_factor[0], temperature, density, -1.0e-3)
    except ValueError as error:
        assert "non-negative" in str(error)
    else:
        raise AssertionError("negative metallicity was accepted by the host atlas")
    invalid_device_rate = jax_metal_net_rate(
        from_numpy_atlas(atlas, dtype=jnp.float64),
        atlas.scale_factor[0],
        jnp.asarray(temperature, dtype=jnp.float64),
        jnp.asarray(density, dtype=jnp.float64),
        jnp.asarray(-1.0e-3, dtype=jnp.float64),
    )
    assert np.isnan(float(invalid_device_rate))

    source_log_density, source_log_temperature, source_coefficient = read_uvb_free_metal_cooling(NO_UVB_DATA)
    density_index = 12
    high_temperature_index = -8
    cmb_log_temperature = np.log10(2.73 / atlas.scale_factor[0])
    cmb_coefficient = 10.0 ** np.interp(
        cmb_log_temperature,
        source_log_temperature,
        np.log10(source_coefficient[density_index]),
    )
    expected_continuous_rate = (
        -source_coefficient[density_index, high_temperature_index] + cmb_coefficient
    ) * (10.0 ** source_log_density[density_index]) ** 2
    assert source_log_temperature[high_temperature_index] - cmb_log_temperature > 2.0
    assert np.isclose(
        atlas.net_rate_erg_s_cm3[0, density_index, high_temperature_index],
        expected_continuous_rate,
        rtol=2.0e-12,
    )
    assert atlas.provenance["cmb_metal_floor"] == "continuous_subtraction_at_tcmb"

    try:
        read_thermal_atlas(LEGACY_ATLAS)
    except ValueError as error:
        assert "format v3" in str(error)
    else:
        raise AssertionError("legacy full-equilibrium atlas was accepted")

    read_uvb_free_metal_cooling(NO_UVB_DATA)
    try:
        read_uvb_free_metal_cooling(HM2012_DATA)
    except ValueError as error:
        assert "UV background" in str(error)
    else:
        raise AssertionError("HM2012 UVB cooling data were accepted as no-UVB metal cooling")

    with TemporaryDirectory() as temporary_directory:
        tampered = Path(temporary_directory) / "tampered.h5"
        shutil.copyfile(ATLAS, tampered)
        with h5py.File(tampered, "r+") as handle:
            handle["provenance"].attrs["uv_background_included"] = "true"
        try:
            read_thermal_atlas(tampered)
        except ValueError as error:
            assert "uv_background_included" in str(error)
        else:
            raise AssertionError("tampered thermal-atlas provenance was accepted")

    print(
        "B1_THERMAL_COUPLING_OK "
        f"neutral_net={neutral_rate:.9e} ionized_net={ionized_rate:.9e} "
        f"metal_z_ratio={host_rates[-1] / solar_rate:.9e} chemistry_ledger=conservative "
        "legacy=rejected uvb_source=rejected"
    )


if __name__ == "__main__":
    main()
