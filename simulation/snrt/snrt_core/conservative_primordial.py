"""Conservative H/He finite-volume S_N chemistry closure."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .implicit import coupled_photo_collisional_hhe_update
from .primordial import (
    EV_ERG,
    PhotoCrossSections,
    PrimordialState,
    case_b_helium_recombination,
    default_photoelectron_excess_energy,
    hui_gnedin_case_b_hydrogen,
)
from .secondary import (
    HELIUM_I_IONIZATION_ENERGY_EV,
    HELIUM_II_IONIZATION_ENERGY_EV,
    HYDROGEN_I_IONIZATION_ENERGY_EV,
    furlanetto_stoever_2010,
)
from .sharding import XShardings
from .transport import TransportConfig, advance_with_absorption, angular_integral
from .primordial_cooling import collisional_ionization_coefficients


class _OpacityState(NamedTuple):
    mean_x_hydrogen_i: jnp.ndarray
    x_helium_ii: jnp.ndarray
    x_helium_iii: jnp.ndarray


class ConservativePrimordialStepResult(NamedTuple):
    """A conservative primordial-chemistry S_N update and species ledgers."""

    intensity: jnp.ndarray
    x_hydrogen_ii: jnp.ndarray
    x_hydrogen_i: jnp.ndarray
    x_helium_ii: jnp.ndarray
    x_helium_iii: jnp.ndarray
    mean_x_hydrogen_ii: jnp.ndarray
    mean_x_hydrogen_i: jnp.ndarray
    absorbed_photons: jnp.ndarray
    gas_photoheating_rate: jnp.ndarray
    photoelectron_excitation_rate: jnp.ndarray
    photoelectron_energy: jnp.ndarray
    photoelectron_energy_ledger_residual: jnp.ndarray
    hydrogen_photoionizations: jnp.ndarray
    helium_i_photoionizations: jnp.ndarray
    helium_ii_photoionizations: jnp.ndarray
    secondary_hydrogen_ionizations: jnp.ndarray
    secondary_helium_i_ionizations: jnp.ndarray
    secondary_helium_ii_ionizations: jnp.ndarray
    hydrogen_collisional_ionizations: jnp.ndarray
    helium_i_collisional_ionizations: jnp.ndarray
    helium_ii_collisional_ionizations: jnp.ndarray
    hydrogen_recombinations: jnp.ndarray
    helium_ii_recombinations: jnp.ndarray
    helium_iii_recombinations: jnp.ndarray
    hydrogen_ledger_residual: jnp.ndarray
    helium_i_ledger_residual: jnp.ndarray
    helium_ii_ledger_residual: jnp.ndarray
    fixed_point_residual: jnp.ndarray
    electron_root_bracket_found: jnp.ndarray


def build_conservative_primordial_step(
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    transport: TransportConfig,
    cross_sections: PhotoCrossSections,
    group_energy_ev: jnp.ndarray,
    *,
    photoelectron_excess_energy_ev: jnp.ndarray | None = None,
    fixed_point_iterations: int = 16,
    fixed_point_relaxation: float = 0.5,
    use_secondary_ionization: bool = False,
    use_tensor_core_angular_reduction: bool = False,
    in_shardings=None,
    out_shardings=None,
):
    """Build a H/He photon-conserving local-implicit S_N step.

    H uses the analytic relaxation and time-averaged H I opacity validated by
    the Strömgren test. Helium opacity is evaluated at the implicit end state;
    its three-state backward-Euler equations make each He I/He II absorbed
    photon equal to the corresponding discrete photoionization term. Optional
    high-energy secondary ionization is charged to a separate photoelectron
    energy ledger. A global fixed point synchronizes those local states with
    directional attenuation.
    """

    if fixed_point_iterations < 1:
        raise ValueError("fixed_point_iterations must be positive")
    if not 0.0 < fixed_point_relaxation <= 1.0:
        raise ValueError("fixed_point_relaxation must lie in (0, 1]")
    group_energy_ev = jnp.asarray(group_energy_ev)
    if photoelectron_excess_energy_ev is None:
        photoelectron_excess_energy_ev = default_photoelectron_excess_energy(group_energy_ev)
    else:
        photoelectron_excess_energy_ev = jnp.asarray(photoelectron_excess_energy_ev)
        if photoelectron_excess_energy_ev.shape != (3, group_energy_ev.shape[0]):
            raise ValueError("photoelectron_excess_energy_ev must have shape (3, n_group)")

    def step(
        intensity: jnp.ndarray,
        emissivity: jnp.ndarray,
        n_hydrogen: jnp.ndarray,
        n_helium: jnp.ndarray,
        x_hydrogen_ii: jnp.ndarray,
        x_helium_ii: jnp.ndarray,
        x_helium_iii: jnp.ndarray,
        temperature_k: jnp.ndarray,
        x_hydrogen_i: jnp.ndarray | None = None,
    ) -> ConservativePrimordialStepResult:
        if x_hydrogen_i is None:
            x_hydrogen_i = 1.0 - x_hydrogen_ii
        extra_axes = (1,) * n_hydrogen.ndim
        sigma_hi = cross_sections.hydrogen_i.reshape((-1,) + extra_axes)
        sigma_hei = cross_sections.helium_i.reshape((-1,) + extra_axes)
        sigma_heii = cross_sections.helium_ii.reshape((-1,) + extra_axes)
        alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
        alpha_heii, alpha_heiii = case_b_helium_recombination(temperature_k)
        collisional = collisional_ionization_coefficients(temperature_k)
        tiny = jnp.finfo(n_hydrogen.dtype).tiny
        minimum_n_hi = jnp.maximum(1.0e-12 * n_hydrogen, tiny)
        minimum_n_he = jnp.maximum(1.0e-12 * n_helium, tiny)
        # Target availability is a start-of-step property. Keeping it fixed
        # avoids switching a secondary channel on or off inside the global
        # opacity iteration when a species crosses the numerical floor.
        target_n_hi = n_hydrogen * x_hydrogen_i
        target_n_hei = n_helium * (1.0 - x_helium_ii - x_helium_iii)
        target_n_heii = n_helium * x_helium_ii

        def solve_at_opacity(opacity_state: _OpacityState) -> ConservativePrimordialStepResult:
            mean_n_hi = n_hydrogen * opacity_state.mean_x_hydrogen_i
            n_hei = n_helium * (1.0 - opacity_state.x_helium_ii - opacity_state.x_helium_iii)
            n_heii = n_helium * opacity_state.x_helium_ii
            kappa_hi = sigma_hi * mean_n_hi[None, ...]
            kappa_hei = sigma_hei * n_hei[None, ...]
            kappa_heii = sigma_heii * n_heii[None, ...]
            total_kappa = kappa_hi + kappa_hei + kappa_heii
            next_intensity, absorbed_directional = advance_with_absorption(
                transport,
                directions,
                intensity,
                emissivity,
                total_kappa,
            )
            absorbed_photons = angular_integral(
                absorbed_directional,
                weights,
                use_tensor_core_reduction=use_tensor_core_angular_reduction,
            )
            safe_kappa = jnp.maximum(total_kappa, jnp.finfo(total_kappa.dtype).tiny)
            hydrogen_photoionizations = jnp.sum(absorbed_photons * kappa_hi / safe_kappa, axis=0)
            helium_i_photoionizations = jnp.sum(absorbed_photons * kappa_hei / safe_kappa, axis=0)
            helium_ii_photoionizations = jnp.sum(absorbed_photons * kappa_heii / safe_kappa, axis=0)
            secondary_hydrogen_ionizations = jnp.zeros_like(n_hydrogen)
            secondary_helium_i_ionizations = jnp.zeros_like(n_hydrogen)
            secondary_helium_ii_ionizations = jnp.zeros_like(n_hydrogen)
            gas_photoheating_energy = jnp.zeros_like(n_hydrogen)
            photoelectron_excitation_energy = jnp.zeros_like(n_hydrogen)
            photoelectron_energy = jnp.zeros_like(n_hydrogen)
            hydrogen_ionized_fraction = 1.0 - opacity_state.mean_x_hydrogen_i
            for absorbed, electron_energy in zip(
                (
                    absorbed_photons * kappa_hi / safe_kappa,
                    absorbed_photons * kappa_hei / safe_kappa,
                    absorbed_photons * kappa_heii / safe_kappa,
                ),
                photoelectron_excess_energy_ev,
                strict=True,
            ):
                deposited_energy = absorbed * electron_energy.reshape((-1,) + extra_axes)
                photoelectron_energy = photoelectron_energy + jnp.sum(deposited_energy, axis=0)
                if use_secondary_ionization:
                    fractions = furlanetto_stoever_2010(
                        electron_energy,
                        hydrogen_ionized_fraction,
                    )
                    hydrogen_fraction = jnp.where(
                        (target_n_hi > minimum_n_hi)[None, ...],
                        fractions.hydrogen_i_ionization,
                        0.0,
                    )
                    helium_i_fraction = jnp.where(
                        (target_n_hei > minimum_n_he)[None, ...],
                        fractions.helium_i_ionization,
                        0.0,
                    )
                    helium_ii_fraction = jnp.where(
                        (target_n_heii > minimum_n_he)[None, ...],
                        fractions.helium_ii_ionization,
                        0.0,
                    )
                    unavailable_ionization_fraction = (
                        fractions.hydrogen_i_ionization
                        + fractions.helium_i_ionization
                        + fractions.helium_ii_ionization
                        - hydrogen_fraction
                        - helium_i_fraction
                        - helium_ii_fraction
                    )
                    secondary_hydrogen_ionizations = secondary_hydrogen_ionizations + jnp.sum(
                        deposited_energy
                        * hydrogen_fraction
                        / HYDROGEN_I_IONIZATION_ENERGY_EV,
                        axis=0,
                    )
                    secondary_helium_i_ionizations = secondary_helium_i_ionizations + jnp.sum(
                        deposited_energy
                        * helium_i_fraction
                        / HELIUM_I_IONIZATION_ENERGY_EV,
                        axis=0,
                    )
                    secondary_helium_ii_ionizations = (
                        secondary_helium_ii_ionizations
                        + jnp.sum(
                            deposited_energy
                            * helium_ii_fraction
                            / HELIUM_II_IONIZATION_ENERGY_EV,
                            axis=0,
                        )
                    )
                    gas_photoheating_energy = gas_photoheating_energy + jnp.sum(
                        deposited_energy
                        * (fractions.heating + unavailable_ionization_fraction),
                        axis=0,
                    )
                    photoelectron_excitation_energy = photoelectron_excitation_energy + jnp.sum(
                        deposited_energy * fractions.excitation, axis=0
                    )
                else:
                    gas_photoheating_energy = gas_photoheating_energy + jnp.sum(deposited_energy, axis=0)
            photoelectron_energy_ledger_residual = photoelectron_energy - (
                HYDROGEN_I_IONIZATION_ENERGY_EV * secondary_hydrogen_ionizations
                + HELIUM_I_IONIZATION_ENERGY_EV * secondary_helium_i_ionizations
                + HELIUM_II_IONIZATION_ENERGY_EV * secondary_helium_ii_ionizations
                + gas_photoheating_energy
                + photoelectron_excitation_energy
            )
            gas_photoheating_rate = gas_photoheating_energy * EV_ERG / transport.dt
            photoelectron_excitation_rate = photoelectron_excitation_energy * EV_ERG / transport.dt
            hydrogen_photoionization_rate = (
                (hydrogen_photoionizations + secondary_hydrogen_ionizations) / transport.dt
            ) / jnp.maximum(mean_n_hi, minimum_n_hi)
            helium_i_photoionization_rate = (
                (helium_i_photoionizations + secondary_helium_i_ionizations) / transport.dt
            ) / jnp.maximum(n_hei, minimum_n_he)
            helium_ii_photoionization_rate = (
                (helium_ii_photoionizations + secondary_helium_ii_ionizations)
                / transport.dt
            ) / jnp.maximum(n_heii, minimum_n_he)
            chemistry_state = PrimordialState(
                n_hydrogen,
                n_helium,
                1.0 - x_hydrogen_i,
                x_helium_ii,
                x_helium_iii,
            )
            (
                solved_state,
                solved_mean_x_hii,
                solved_mean_x_hi,
                electron_density,
                electron_root_bracket_found,
            ) = (
                coupled_photo_collisional_hhe_update(
                    chemistry_state,
                    hydrogen_photoionization_rate,
                    helium_i_photoionization_rate,
                    helium_ii_photoionization_rate,
                    temperature_k,
                    transport.dt,
                )
            )
            next_x_hii = solved_state.x_hydrogen_ii
            next_x_hi = 1.0 - next_x_hii
            next_x_heii = solved_state.x_helium_ii
            next_x_heiii = solved_state.x_helium_iii
            hydrogen_collisional_ionizations = (
                transport.dt
                * collisional.hydrogen_i
                * electron_density
                * n_hydrogen
                * solved_mean_x_hi
            )
            helium_i_collisional_ionizations = (
                transport.dt
                * collisional.helium_i
                * electron_density
                * n_helium
                * (1.0 - next_x_heii - next_x_heiii)
            )
            helium_ii_collisional_ionizations = (
                transport.dt
                * collisional.helium_ii
                * electron_density
                * n_helium
                * next_x_heii
            )
            hydrogen_recombinations = (
                transport.dt
                * alpha_hii
                * electron_density
                * n_hydrogen
                * solved_mean_x_hii
            )
            helium_ii_recombinations = (
                transport.dt * alpha_heii * electron_density * n_helium * next_x_heii
            )
            helium_iii_recombinations = (
                transport.dt * alpha_heiii * electron_density * n_helium * next_x_heiii
            )
            hydrogen_change = n_hydrogen * (x_hydrogen_i - next_x_hi)
            helium_ii_change = n_helium * (next_x_heii - x_helium_ii)
            helium_iii_change = n_helium * (next_x_heiii - x_helium_iii)
            fixed_point_residual = jnp.maximum(
                jnp.abs(solved_mean_x_hi - opacity_state.mean_x_hydrogen_i),
                jnp.maximum(
                    jnp.abs(next_x_heii - opacity_state.x_helium_ii),
                    jnp.abs(next_x_heiii - opacity_state.x_helium_iii),
                ),
            )
            return ConservativePrimordialStepResult(
                intensity=next_intensity,
                x_hydrogen_ii=next_x_hii,
                x_hydrogen_i=next_x_hi,
                x_helium_ii=next_x_heii,
                x_helium_iii=next_x_heiii,
                mean_x_hydrogen_ii=solved_mean_x_hii,
                mean_x_hydrogen_i=solved_mean_x_hi,
                absorbed_photons=absorbed_photons,
                gas_photoheating_rate=gas_photoheating_rate,
                photoelectron_excitation_rate=photoelectron_excitation_rate,
                photoelectron_energy=photoelectron_energy,
                photoelectron_energy_ledger_residual=photoelectron_energy_ledger_residual,
                hydrogen_photoionizations=hydrogen_photoionizations,
                helium_i_photoionizations=helium_i_photoionizations,
                helium_ii_photoionizations=helium_ii_photoionizations,
                secondary_hydrogen_ionizations=secondary_hydrogen_ionizations,
                secondary_helium_i_ionizations=secondary_helium_i_ionizations,
                secondary_helium_ii_ionizations=secondary_helium_ii_ionizations,
                hydrogen_collisional_ionizations=hydrogen_collisional_ionizations,
                helium_i_collisional_ionizations=helium_i_collisional_ionizations,
                helium_ii_collisional_ionizations=helium_ii_collisional_ionizations,
                hydrogen_recombinations=hydrogen_recombinations,
                helium_ii_recombinations=helium_ii_recombinations,
                helium_iii_recombinations=helium_iii_recombinations,
                hydrogen_ledger_residual=hydrogen_photoionizations
                + secondary_hydrogen_ionizations
                + hydrogen_collisional_ionizations
                - hydrogen_change
                - hydrogen_recombinations,
                helium_i_ledger_residual=helium_i_photoionizations
                + secondary_helium_i_ionizations
                + helium_i_collisional_ionizations
                - helium_ii_change
                - helium_iii_change
                - helium_ii_recombinations,
                helium_ii_ledger_residual=helium_ii_photoionizations
                + secondary_helium_ii_ionizations
                + helium_ii_collisional_ionizations
                - helium_iii_change
                - helium_iii_recombinations,
                fixed_point_residual=fixed_point_residual,
                electron_root_bracket_found=electron_root_bracket_found,
            )

        def refine(_: int, opacity_state: _OpacityState) -> _OpacityState:
            solved = solve_at_opacity(opacity_state)
            return _OpacityState(
                (1.0 - fixed_point_relaxation) * opacity_state.mean_x_hydrogen_i
                + fixed_point_relaxation * solved.mean_x_hydrogen_i,
                (1.0 - fixed_point_relaxation) * opacity_state.x_helium_ii
                + fixed_point_relaxation * solved.x_helium_ii,
                (1.0 - fixed_point_relaxation) * opacity_state.x_helium_iii
                + fixed_point_relaxation * solved.x_helium_iii,
            )

        initial_opacity = _OpacityState(x_hydrogen_i, x_helium_ii, x_helium_iii)
        opacity_state = jax.lax.fori_loop(0, fixed_point_iterations, refine, initial_opacity)
        return solve_at_opacity(opacity_state)

    return jax.jit(step, in_shardings=in_shardings, out_shardings=out_shardings)


def build_x_sharded_conservative_primordial_step(
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    transport: TransportConfig,
    cross_sections: PhotoCrossSections,
    group_energy_ev: jnp.ndarray,
    shardings: XShardings,
    *,
    photoelectron_excess_energy_ev: jnp.ndarray | None = None,
    fixed_point_iterations: int = 16,
    fixed_point_relaxation: float = 0.5,
    use_secondary_ionization: bool = False,
    use_tensor_core_angular_reduction: bool = False,
):
    """Build a full conservative S_N step sharded over the x cell dimension.

    Angular groups and directions are replicated.  Every cell-centered state
    and diagnostic follows the same x partition, so the fixed-point opacity,
    photon ledger, and chemistry closure remain local to the transport shard.
    XLA inserts the directional x-face exchanges required by the upwind stencil.
    """

    scalar = shardings.scalar_field
    group = shardings.group_field
    result_shardings = ConservativePrimordialStepResult(
        intensity=shardings.intensity,
        x_hydrogen_ii=scalar,
        x_hydrogen_i=scalar,
        x_helium_ii=scalar,
        x_helium_iii=scalar,
        mean_x_hydrogen_ii=scalar,
        mean_x_hydrogen_i=scalar,
        absorbed_photons=group,
        gas_photoheating_rate=scalar,
        photoelectron_excitation_rate=scalar,
        photoelectron_energy=scalar,
        photoelectron_energy_ledger_residual=scalar,
        hydrogen_photoionizations=scalar,
        helium_i_photoionizations=scalar,
        helium_ii_photoionizations=scalar,
        secondary_hydrogen_ionizations=scalar,
        secondary_helium_i_ionizations=scalar,
        secondary_helium_ii_ionizations=scalar,
        hydrogen_collisional_ionizations=scalar,
        helium_i_collisional_ionizations=scalar,
        helium_ii_collisional_ionizations=scalar,
        hydrogen_recombinations=scalar,
        helium_ii_recombinations=scalar,
        helium_iii_recombinations=scalar,
        hydrogen_ledger_residual=scalar,
        helium_i_ledger_residual=scalar,
        helium_ii_ledger_residual=scalar,
        fixed_point_residual=scalar,
        electron_root_bracket_found=scalar,
    )
    return build_conservative_primordial_step(
        directions,
        weights,
        transport,
        cross_sections,
        group_energy_ev,
        photoelectron_excess_energy_ev=photoelectron_excess_energy_ev,
        fixed_point_iterations=fixed_point_iterations,
        fixed_point_relaxation=fixed_point_relaxation,
        use_secondary_ionization=use_secondary_ionization,
        use_tensor_core_angular_reduction=use_tensor_core_angular_reduction,
        in_shardings=(
            shardings.intensity,
            group,
            scalar,
            scalar,
            scalar,
            scalar,
            scalar,
            scalar,
            scalar,
        ),
        out_shardings=result_shardings,
    )
