"""Non-equilibrium atomic H/He cooling evaluated at the SNRT ion state.

The rate coefficients follow Grackle 3.4.2-dev, revision
``f93091ff8456962d7017a5bff7472945a30e3dad``.  Grackle in turn identifies
the atomic excitation and ionization fits with Black (1981), Cen (1992), and
Abel et al. (1997), the recombination fits with Hui & Gnedin (1997), and the
free-free fit with Black (1981) and Spitzer & Hart (1979).

All public cooling functions return volumetric rates in erg cm^-3 s^-1.
``primordial_net_rate`` uses the SNRT thermal sign convention: heating is
positive and cooling is negative.
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp

from .primordial import (
    PrimordialState,
    electron_number_density,
    hui_gnedin_case_b_helium_ii_radiative,
    neutral_number_densities,
)


BOLTZMANN_ERG_K = 1.380649e-16
CMB_TEMPERATURE_K = 2.73


class CollisionalIonizationCoefficients(NamedTuple):
    """Electron-impact ionization coefficients in cm^3 s^-1."""

    hydrogen_i: jnp.ndarray
    helium_i: jnp.ndarray
    helium_ii: jnp.ndarray


class PrimordialCoolingComponents(NamedTuple):
    """Positive atomic cooling components plus signed CMB Compton exchange."""

    collisional_excitation: jnp.ndarray
    collisional_ionization: jnp.ndarray
    recombination: jnp.ndarray
    bremsstrahlung: jnp.ndarray
    compton_net_heating: jnp.ndarray


def _exp_polynomial(log_temperature_ev: jnp.ndarray, coefficients: tuple[float, ...]) -> jnp.ndarray:
    """Evaluate ``exp(sum(c_i log(T_eV)^i))`` with Horner's rule."""

    value = jnp.zeros_like(log_temperature_ev) + coefficients[-1]
    for coefficient in reversed(coefficients[:-1]):
        value = value * log_temperature_ev + coefficient
    return jnp.exp(value)


def collisional_ionization_coefficients(
    temperature_k: jnp.ndarray,
) -> CollisionalIonizationCoefficients:
    """Return Abel et al. electron-impact H I, He I, and He II rates."""

    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    temperature_ev = temperature / 11605.0
    log_temperature_ev = jnp.log(temperature_ev)
    hydrogen_i = _exp_polynomial(
        log_temperature_ev,
        (
            -32.71396786375,
            13.53655609057,
            -5.739328757388,
            1.563154982022,
            -0.2877056004391,
            0.03482559773736999,
            -0.00263197617559,
            0.0001119543953861,
            -2.039149852002e-6,
        ),
    )
    helium_i = _exp_polynomial(
        log_temperature_ev,
        (
            -44.09864886561001,
            23.91596563469,
            -10.75323019821,
            3.058038757198,
            -0.5685118909884001,
            0.06795391233790001,
            -0.005009056101857001,
            0.0002067236157507,
            -3.649161410833e-6,
        ),
    )
    helium_ii = _exp_polynomial(
        log_temperature_ev,
        (
            -68.71040990212001,
            43.93347632635,
            -18.48066993568,
            4.701626486759002,
            -0.7692466334492,
            0.08113042097303,
            -0.005324020628287001,
            0.0001975705312221,
            -3.165581065665e-6,
        ),
    )
    active_helium = temperature_ev > 0.8
    return CollisionalIonizationCoefficients(
        hydrogen_i=jnp.maximum(hydrogen_i, jnp.finfo(temperature.dtype).tiny),
        helium_i=jnp.where(active_helium, helium_i, 0.0),
        helium_ii=jnp.where(active_helium, helium_ii, 0.0),
    )


def _hui_gnedin_case_b_hydrogen_recombination_cooling(
    temperature_k: jnp.ndarray,
) -> jnp.ndarray:
    """Return the Hui--Gnedin H II case-B cooling coefficient."""

    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    lambda_hi = 2.0 * 157807.0 / temperature
    return (
        3.435e-30
        * temperature
        * lambda_hi**1.970
        / (1.0 + (lambda_hi / 2.25) ** 0.376) ** 3.720
    )


def _case_b_recombination_cooling_coefficients(
    temperature_k: jnp.ndarray,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Return H II, He II radiative/dielectronic, and He III cooling fits."""

    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    hydrogen_ii = _hui_gnedin_case_b_hydrogen_recombination_cooling(temperature)
    helium_ii_radiative = (
        hui_gnedin_case_b_helium_ii_radiative(temperature)
        * BOLTZMANN_ERG_K
        * temperature
    )
    helium_ii_dielectronic = (
        1.24e-13
        * temperature**-1.5
        * jnp.exp(-470000.0 / temperature)
        * (1.0 + 0.3 * jnp.exp(-94000.0 / temperature))
    )
    helium_iii = 8.0 * _hui_gnedin_case_b_hydrogen_recombination_cooling(
        temperature / 4.0
    )
    return hydrogen_ii, helium_ii_radiative, helium_ii_dielectronic, helium_iii


def primordial_cooling_components(
    state: PrimordialState,
    temperature_k: jnp.ndarray,
    scale_factor: float | jnp.ndarray,
) -> PrimordialCoolingComponents:
    """Evaluate atomic cooling at the network's non-equilibrium ion fractions."""

    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    electron = electron_number_density(state)
    hydrogen_i, helium_i, helium_ii = neutral_number_densities(state)
    hydrogen_ii = state.n_hydrogen * state.x_hydrogen_ii
    helium_iii = state.n_helium * state.x_helium_iii

    denominator = 1.0 + jnp.sqrt(temperature / 1.0e5)
    excitation_hydrogen_i = 7.5e-19 * jnp.exp(-118348.0 / temperature) / denominator
    excitation_helium_i = 9.1e-27 * temperature**-0.1687 * jnp.exp(-13179.0 / temperature) / denominator
    excitation_helium_ii = 5.54e-17 * temperature**-0.397 * jnp.exp(-473638.0 / temperature) / denominator
    collisional_excitation = (
        excitation_hydrogen_i * hydrogen_i * electron
        + excitation_helium_i * helium_ii * electron**2
        + excitation_helium_ii * helium_ii * electron
    )

    ionization = collisional_ionization_coefficients(temperature)
    ionization_helium_i_excited = (
        5.01e-27
        * temperature**-0.1687
        * jnp.exp(-55338.0 / temperature)
        / denominator
    )
    collisional_ionization = (
        2.18e-11 * ionization.hydrogen_i * hydrogen_i * electron
        + 3.94e-11 * ionization.helium_i * helium_i * electron
        + 8.72e-11 * ionization.helium_ii * helium_ii * electron
        + ionization_helium_i_excited * helium_ii * electron**2
    )

    recombination_coefficients = _case_b_recombination_cooling_coefficients(temperature)
    recombination = electron * (
        recombination_coefficients[0] * hydrogen_ii
        + (recombination_coefficients[1] + recombination_coefficients[2]) * helium_ii
        + recombination_coefficients[3] * helium_iii
    )

    gaunt_factor = 1.1 + 0.34 * jnp.exp(-(5.5 - jnp.log10(temperature)) ** 2 / 3.0)
    bremsstrahlung = (
        1.43e-27
        * jnp.sqrt(temperature)
        * gaunt_factor
        * electron
        * (hydrogen_ii + helium_ii + 4.0 * helium_iii)
    )

    expansion = jnp.maximum(jnp.asarray(scale_factor, dtype=temperature.dtype), jnp.finfo(temperature.dtype).tiny)
    redshift_factor = 1.0 / expansion
    cmb_temperature = CMB_TEMPERATURE_K * redshift_factor
    compton_net_heating = -5.65e-36 * redshift_factor**4 * (temperature - cmb_temperature) * electron
    return PrimordialCoolingComponents(
        collisional_excitation,
        collisional_ionization,
        recombination,
        bremsstrahlung,
        compton_net_heating,
    )


def primordial_net_rate(
    state: PrimordialState,
    temperature_k: jnp.ndarray,
    scale_factor: float | jnp.ndarray,
) -> jnp.ndarray:
    """Return signed non-equilibrium primordial thermal rate."""

    components = primordial_cooling_components(state, temperature_k, scale_factor)
    atomic_cooling = (
        components.collisional_excitation
        + components.collisional_ionization
        + components.recombination
        + components.bremsstrahlung
    )
    return components.compton_net_heating - atomic_cooling
