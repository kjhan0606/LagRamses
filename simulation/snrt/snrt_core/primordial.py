"""Primordial H/He photo-chemistry primitives for the static S_N core.

Cross sections use the analytic form of Verner et al. (1996). The H II,
He II, and He III radiative case-B coefficients follow Hui & Gnedin (1997).
He II dielectronic recombination is added separately to its radiative case-B
coefficient.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import NamedTuple

import jax.numpy as jnp
import numpy as np

from snrt_core.provenance import (
    require_sha256,
    sha256_file,
    validate_code_manifest,
    validate_payload_hash,
)
from snrt_core.sed import (
    QuadratureDiagnostics,
    SOURCE_SED_CANDIDATE_STATUS,
    SOURCE_SED_CONTRACT_STATUSES,
    SED_INTERPOLATION_CONVENTION,
    SED_QUADRATURE_BASE_SUBDIVISIONS,
    SED_QUADRATURE_RELATIVE_TOLERANCE,
    SED_QUADRATURE_REFINED_SUBDIVISIONS,
    SED_QUADRATURE_SCHEME,
)


EV_ERG = 1.602176634e-12
SNRT_ROOT = Path(__file__).resolve().parents[1]
AGN_PHOTON_CLOSURE_CODE_MANIFEST = {
    "agn_ledger_builder": SNRT_ROOT / "tools" / "p4_build_agn_photon_ledger.py",
    "source_sed": SNRT_ROOT / "snrt_core" / "sed.py",
    "primordial_closure": SNRT_ROOT / "snrt_core" / "primordial.py",
    "integrity_helper": SNRT_ROOT / "snrt_core" / "provenance.py",
}


class VernerFit(NamedTuple):
    threshold_ev: float
    maximum_ev: float
    energy_scale_ev: float
    sigma_scale_mb: float
    y_a: float
    power: float
    y_w: float
    y_0: float
    y_1: float


# Table 1 of Verner et al. (1996), ground-state H I, He I, and He II.
H_I_FIT = VernerFit(13.60, 5.0e4, 4.298e-1, 5.475e4, 3.288e1, 2.963, 0.0, 0.0, 0.0)
HE_I_FIT = VernerFit(24.59, 5.0e4, 1.361e1, 9.492e2, 1.469, 3.188, 2.039, 4.434e-1, 2.136)
HE_II_FIT = VernerFit(54.42, 5.0e4, 1.720, 1.369e4, 3.288e1, 2.963, 0.0, 0.0, 0.0)


class PhotoCrossSections(NamedTuple):
    """Group-averaged cross sections in cm^2, indexed as ``(group,)``."""

    hydrogen_i: jnp.ndarray
    helium_i: jnp.ndarray
    helium_ii: jnp.ndarray


class GroupSpectralClosure(NamedTuple):
    """SED closure shared by source conversion and the RT microphysics.

    ``photoelectron_excess_energy_ev`` has shape ``(3, n_group)`` and uses
    the species order H I, He I, He II.  Its value is weighted by the same
    absorber cross section used to construct ``cross_sections``; this keeps
    photoheating consistent with a group-integrated photon budget.
    """

    cross_sections: PhotoCrossSections
    photon_weighted_energy_ev: jnp.ndarray
    photoelectron_excess_energy_ev: jnp.ndarray


class PrimordialState(NamedTuple):
    """Number densities and ion fractions; all densities are in cm^-3."""

    n_hydrogen: jnp.ndarray
    n_helium: jnp.ndarray
    x_hydrogen_ii: jnp.ndarray
    x_helium_ii: jnp.ndarray
    x_helium_iii: jnp.ndarray


class PhotoRates(NamedTuple):
    """Per-absorber photoionization rates in s^-1."""

    hydrogen_i: jnp.ndarray
    helium_i: jnp.ndarray
    helium_ii: jnp.ndarray


def _verner_cross_section_numpy(energy_ev: np.ndarray, fit: VernerFit) -> np.ndarray:
    """Evaluate a Verner fit in NumPy for offline SED quadrature."""

    energy = np.asarray(energy_ev, dtype=np.float64)
    x = energy / fit.energy_scale_ev - fit.y_0
    y = np.sqrt(x**2 + fit.y_1**2)
    profile = ((x - 1.0) ** 2 + fit.y_w**2) * y ** (0.5 * fit.power - 5.5)
    profile *= (1.0 + np.sqrt(y / fit.y_a)) ** (-fit.power)
    sigma = fit.sigma_scale_mb * profile * 1.0e-18
    return np.where((energy >= fit.threshold_ev) & (energy <= fit.maximum_ev), sigma, 0.0)


def verner_cross_section(energy_ev: jnp.ndarray, fit: VernerFit) -> jnp.ndarray:
    """Evaluate a Verner et al. ground-state photoionization fit in cm^2."""
    energy = jnp.asarray(energy_ev)
    x = energy / fit.energy_scale_ev - fit.y_0
    y = jnp.sqrt(x**2 + fit.y_1**2)
    profile = ((x - 1.0) ** 2 + fit.y_w**2) * y ** (0.5 * fit.power - 5.5)
    profile *= (1.0 + jnp.sqrt(y / fit.y_a)) ** (-fit.power)
    sigma = fit.sigma_scale_mb * profile * 1.0e-18
    return jnp.where((energy >= fit.threshold_ev) & (energy <= fit.maximum_ev), sigma, 0.0)


def primordial_cross_sections(group_energy_ev: jnp.ndarray) -> PhotoCrossSections:
    """Return compatibility cross sections evaluated at group representative energies.

    Production source metadata should use :func:`sed_weighted_group_closure`.
    This centre-energy helper remains for analytic benchmarks whose groups are
    deliberately monochromatic.
    """
    energies = jnp.asarray(group_energy_ev)
    return PhotoCrossSections(
        hydrogen_i=verner_cross_section(energies, H_I_FIT),
        helium_i=verner_cross_section(energies, HE_I_FIT),
        helium_ii=verner_cross_section(energies, HE_II_FIT),
    )


def default_photoelectron_excess_energy(group_energy_ev: jnp.ndarray) -> jnp.ndarray:
    """Return the legacy representative-energy excess closure ``(3, group)``."""

    energies = jnp.asarray(group_energy_ev)
    thresholds = jnp.asarray((H_I_FIT.threshold_ev, HE_I_FIT.threshold_ev, HE_II_FIT.threshold_ev), dtype=energies.dtype)
    return jnp.maximum(energies[None, :] - thresholds[:, None], 0.0)


def sed_weighted_group_closure(
    group_edges_ev: np.ndarray | jnp.ndarray,
    energy_ev: np.ndarray | jnp.ndarray,
    photon_number_spectrum_per_ev: np.ndarray | jnp.ndarray,
    *,
    allow_empty_groups: bool = False,
) -> GroupSpectralClosure:
    """Integrate Verner cross sections over a photon-number SED.

    For each group ``g`` and absorber ``s`` this computes

    ``sigma_bar[g,s] = ∫ N_E sigma_s(E) dE / ∫ N_E dE``

    and the corresponding absorption-weighted photoelectron excess energy.
    The SED is an arbitrary non-negative shape; its normalization cancels from
    the closure.  This is an offline operation so that the resulting arrays
    are static inputs to JAX/XLA.  When ``allow_empty_groups`` is true, a group
    with zero source photons receives zero absorber closure and a geometric-
    mean representative energy; this is valid only for an inactive source
    group and is recorded by the caller.
    """

    closure, _ = sed_weighted_group_closure_with_diagnostics(
        group_edges_ev,
        energy_ev,
        photon_number_spectrum_per_ev,
        allow_empty_groups=allow_empty_groups,
    )
    return closure


def _sed_group_grids(
    energies: np.ndarray,
    edges: np.ndarray,
    extra_energies: np.ndarray,
    *,
    subdivisions: int,
) -> tuple[np.ndarray, ...]:
    if not isinstance(subdivisions, int) or subdivisions < 1:
        raise ValueError("quadrature subdivisions must be a positive integer")
    grids: list[np.ndarray] = []
    for lower, upper in zip(edges[:-1], edges[1:], strict=True):
        refinement = np.geomspace(lower, upper, subdivisions + 1, dtype=np.float64)
        source_nodes = energies[(energies > lower) & (energies < upper)]
        extra_nodes = extra_energies[
            (extra_energies > lower) & (extra_energies < upper)
        ]
        grids.append(np.unique(np.concatenate((refinement, source_nodes, extra_nodes))))
    return tuple(grids)


def _closure_once(
    edges: np.ndarray,
    energies: np.ndarray,
    spectrum: np.ndarray,
    *,
    energy_fraction_per_ev: np.ndarray | None,
    allow_empty_groups: bool,
    subdivisions: int,
) -> GroupSpectralClosure:
    """Evaluate one deterministic quadrature level of the Verner closure."""

    species_fits = (H_I_FIT, HE_I_FIT, HE_II_FIT)
    thresholds = np.asarray([fit.threshold_ev for fit in species_fits], dtype=np.float64)
    grids = _sed_group_grids(
        energies,
        edges,
        thresholds,
        subdivisions=subdivisions,
    )
    averaged_sigma = np.zeros((3, len(edges) - 1), dtype=np.float64)
    excess_energy = np.zeros_like(averaged_sigma)
    photon_mean_energy = np.sqrt(edges[:-1] * edges[1:])

    for group, group_energy in enumerate(grids):
        if energy_fraction_per_ev is None:
            group_spectrum = np.interp(group_energy, energies, spectrum)
        else:
            group_fraction = np.interp(
                group_energy, energies, energy_fraction_per_ev
            )
            group_spectrum = group_fraction / (group_energy * EV_ERG)
        photon_count = float(np.trapezoid(group_spectrum, group_energy))
        if not np.isfinite(photon_count) or photon_count < 0.0:
            raise ValueError(f"SED has an invalid photon integral in group {group}")
        if photon_count == 0.0:
            if not allow_empty_groups:
                raise ValueError(f"SED has no photons in group {group}")
            continue
        photon_mean_energy[group] = (
            np.trapezoid(group_spectrum * group_energy, group_energy) / photon_count
        )
        for species, fit in enumerate(species_fits):
            absorbing_lower = max(float(edges[group]), fit.threshold_ev)
            absorbing_upper = min(float(edges[group + 1]), fit.maximum_ev)
            if absorbing_lower >= absorbing_upper:
                continue
            absorbing_mask = (group_energy >= absorbing_lower) & (
                group_energy <= absorbing_upper
            )
            absorbing_energy = group_energy[absorbing_mask]
            absorbing_spectrum = group_spectrum[absorbing_mask]
            if absorbing_energy[0] != absorbing_lower:
                absorbing_energy = np.insert(absorbing_energy, 0, absorbing_lower)
                if energy_fraction_per_ev is None:
                    absorbing_spectrum = np.insert(
                        absorbing_spectrum,
                        0,
                        np.interp(absorbing_lower, energies, spectrum),
                    )
                else:
                    absorbing_spectrum = np.insert(
                        absorbing_spectrum,
                        0,
                        np.interp(absorbing_lower, energies, energy_fraction_per_ev)
                        / (absorbing_lower * EV_ERG),
                    )
            if absorbing_energy[-1] != absorbing_upper:
                absorbing_energy = np.append(absorbing_energy, absorbing_upper)
                if energy_fraction_per_ev is None:
                    absorbing_spectrum = np.append(
                        absorbing_spectrum,
                        np.interp(absorbing_upper, energies, spectrum),
                    )
                else:
                    absorbing_spectrum = np.append(
                        absorbing_spectrum,
                        np.interp(absorbing_upper, energies, energy_fraction_per_ev)
                        / (absorbing_upper * EV_ERG),
                    )
            sigma = _verner_cross_section_numpy(absorbing_energy, fit)
            weighted_sigma = float(
                np.trapezoid(absorbing_spectrum * sigma, absorbing_energy)
            )
            averaged_sigma[species, group] = weighted_sigma / photon_count
            if weighted_sigma > 0.0:
                excess_energy[species, group] = (
                    np.trapezoid(
                        absorbing_spectrum
                        * sigma
                        * (absorbing_energy - fit.threshold_ev),
                        absorbing_energy,
                    )
                    / weighted_sigma
                )

    return GroupSpectralClosure(
        cross_sections=PhotoCrossSections(
            hydrogen_i=averaged_sigma[0],
            helium_i=averaged_sigma[1],
            helium_ii=averaged_sigma[2],
        ),
        photon_weighted_energy_ev=photon_mean_energy,
        photoelectron_excess_energy_ev=excess_energy,
    )


def _closure_relative_error_by_group(
    reference: GroupSpectralClosure,
    refined: GroupSpectralClosure,
) -> np.ndarray:
    reference_values = np.concatenate(
        (
            np.asarray(reference.cross_sections, dtype=np.float64),
            np.asarray(reference.photon_weighted_energy_ev, dtype=np.float64)[None, :],
            np.asarray(reference.photoelectron_excess_energy_ev, dtype=np.float64),
        ),
        axis=0,
    )
    refined_values = np.concatenate(
        (
            np.asarray(refined.cross_sections, dtype=np.float64),
            np.asarray(refined.photon_weighted_energy_ev, dtype=np.float64)[None, :],
            np.asarray(refined.photoelectron_excess_energy_ev, dtype=np.float64),
        ),
        axis=0,
    )
    scale = np.maximum(
        np.maximum(np.abs(reference_values), np.abs(refined_values)), 1.0e-300
    )
    return np.max(np.abs(refined_values - reference_values) / scale, axis=0)


def sed_weighted_group_closure_with_diagnostics(
    group_edges_ev: np.ndarray | jnp.ndarray,
    energy_ev: np.ndarray | jnp.ndarray,
    photon_number_spectrum_per_ev: np.ndarray | jnp.ndarray,
    *,
    energy_fraction_per_ev: np.ndarray | None = None,
    interpolation_convention: str = "piecewise_linear_photon_number_per_ev",
    allow_empty_groups: bool = False,
) -> tuple[GroupSpectralClosure, QuadratureDiagnostics]:
    """Return a Verner closure plus a base-versus-refined convergence record.

    If ``energy_fraction_per_ev`` is supplied, it is the authoritative
    tabulated quantity and the photon spectrum is derived as
    ``f_E/(E*EV_ERG)``.  This is the explicit-source convention.  Without it,
    callers retain the historical piecewise-linear photon-number spectrum
    interface used by stellar and analytic controls.
    """

    edges = np.asarray(group_edges_ev, dtype=np.float64)
    energies = np.asarray(energy_ev, dtype=np.float64)
    spectrum = np.asarray(photon_number_spectrum_per_ev, dtype=np.float64)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or np.any(edges <= 0.0):
        raise ValueError("group_edges_ev must be a finite, positive one-dimensional edge array")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("group_edges_ev must be strictly increasing")
    if energies.ndim != 1 or spectrum.shape != energies.shape or len(energies) < 2:
        raise ValueError("energy and photon SED arrays must be one-dimensional with at least two samples")
    if not np.isfinite(energies).all() or np.any(energies <= 0.0):
        raise ValueError("SED energies must be finite and positive")
    if not np.isfinite(spectrum).all() or np.any(spectrum < 0.0):
        raise ValueError("photon-number SED must be finite and non-negative")
    if energy_fraction_per_ev is not None:
        energy_fraction = np.asarray(energy_fraction_per_ev, dtype=np.float64)
        if energy_fraction.shape != energies.shape:
            raise ValueError("energy-fraction and photon SED arrays must have identical shapes")
        if not np.isfinite(energy_fraction).all() or np.any(energy_fraction < 0.0):
            raise ValueError("energy-fraction SED must be finite and non-negative")
        if interpolation_convention != SED_INTERPOLATION_CONVENTION:
            raise ValueError(
                "energy-fraction SEDs require the declared "
                f"{SED_INTERPOLATION_CONVENTION} convention"
            )
    elif interpolation_convention != "piecewise_linear_photon_number_per_ev":
        raise ValueError(
            "photon-number SEDs require piecewise_linear_photon_number_per_ev "
            "unless an energy-fraction table is supplied"
        )
    else:
        energy_fraction = None

    order = np.argsort(energies)
    energies = energies[order]
    spectrum = spectrum[order]
    if energy_fraction is not None:
        energy_fraction = energy_fraction[order]
    if np.any(np.diff(energies) <= 0.0):
        unique = np.r_[True, np.diff(energies) > 0.0]
        energies = energies[unique]
        spectrum = spectrum[unique]
        if energy_fraction is not None:
            energy_fraction = energy_fraction[unique]
    if edges[0] < energies[0] or edges[-1] > energies[-1]:
        raise ValueError("SED energy support must cover every requested group edge")

    base = _closure_once(
        edges,
        energies,
        spectrum,
        energy_fraction_per_ev=energy_fraction,
        allow_empty_groups=allow_empty_groups,
        subdivisions=SED_QUADRATURE_BASE_SUBDIVISIONS,
    )
    refined = _closure_once(
        edges,
        energies,
        spectrum,
        energy_fraction_per_ev=energy_fraction,
        allow_empty_groups=allow_empty_groups,
        subdivisions=SED_QUADRATURE_REFINED_SUBDIVISIONS,
    )
    group_error = _closure_relative_error_by_group(base, refined)
    maximum_error = float(np.max(group_error))
    diagnostics = QuadratureDiagnostics(
        scheme=SED_QUADRATURE_SCHEME,
        interpolation_convention=interpolation_convention,
        base_subdivisions=SED_QUADRATURE_BASE_SUBDIVISIONS,
        refined_subdivisions=SED_QUADRATURE_REFINED_SUBDIVISIONS,
        relative_tolerance=SED_QUADRATURE_RELATIVE_TOLERANCE,
        group_max_relative_error=tuple(float(value) for value in group_error),
        maximum_relative_error=maximum_error,
        converged=maximum_error <= SED_QUADRATURE_RELATIVE_TOLERANCE,
    )
    if not diagnostics.converged:
        raise ValueError(
            "Verner SED quadrature did not converge: "
            f"maximum relative error {maximum_error:.6g} exceeds "
            f"{SED_QUADRATURE_RELATIVE_TOLERANCE:.6g}"
        )
    return refined, diagnostics


def group_spectral_closure_from_metadata(
    metadata: Mapping[str, object],
    *,
    require_code_manifest: bool = False,
) -> GroupSpectralClosure:
    """Load and validate the serialized SED closure in photon metadata."""

    schema = metadata.get("schema")
    source_identity = metadata.get("source_sed_identity")
    source_hash = metadata.get("source_sed_sha256")
    if source_identity is not None:
        source_identity = require_sha256(
            source_identity, "source_sed_identity", "serialized photon metadata"
        )
        contract = metadata.get("source_sed_contract")
        if not isinstance(contract, Mapping) or contract.get("identity") != source_identity:
            raise ValueError("serialized source SED contract does not match source_sed_identity")
    if source_hash is not None:
        source_hash = require_sha256(
            source_hash, "source_sed_sha256", "serialized photon metadata"
        )
    if isinstance(source_identity, str) and source_hash is not None:
        contract = metadata.get("source_sed_contract")
        if isinstance(contract, Mapping) and contract.get("input_sha256") != source_hash:
            raise ValueError("serialized source SED contract does not match source_sed_sha256")
    if schema == "snrt_agn_photon_ledger_v2" and source_identity is not None:
        if not isinstance(source_hash, str):
            raise ValueError("source-bound AGN photon metadata must include source_sed_sha256")
        contract = metadata["source_sed_contract"]
        assert isinstance(contract, Mapping)
        if contract.get("status") not in SOURCE_SED_CONTRACT_STATUSES:
            raise ValueError("source-bound AGN photon metadata has an invalid source SED status")
        if contract.get("status") != SOURCE_SED_CANDIDATE_STATUS:
            raise ValueError(
                "source-bound AGN photon metadata must use the candidate source SED status"
            )
        if contract.get("status") == SOURCE_SED_CANDIDATE_STATUS:
            if contract.get("interpolation_convention") != SED_INTERPOLATION_CONVENTION:
                raise ValueError(
                    "source-bound AGN photon metadata has an unsupported SED interpolation convention"
                )
            quadrature = contract.get("quadrature")
            if not isinstance(quadrature, Mapping) or quadrature.get("scheme") != SED_QUADRATURE_SCHEME:
                raise ValueError(
                    "source-bound AGN photon metadata lacks the declared SED quadrature scheme"
                )
        source_input_path = contract.get("input_path")
        if not isinstance(source_input_path, str) or not source_input_path.strip():
            raise ValueError("source-bound AGN photon metadata lacks source SED input_path")
        try:
            actual_source_hash = sha256_file(source_input_path)
        except (FileNotFoundError, OSError) as error:
            raise ValueError("source-bound AGN photon SED input file is unavailable") from error
        if actual_source_hash != source_hash:
            raise ValueError("source-bound AGN photon SED input hash does not match its file")
        validate_code_manifest(
            metadata,
            AGN_PHOTON_CLOSURE_CODE_MANIFEST,
            context="source-bound AGN photon metadata",
        )
        validate_payload_hash(metadata, context="source-bound AGN photon metadata")
    elif require_code_manifest and source_identity is not None:
        raise ValueError(
            "source-bound photon metadata does not declare the supported AGN closure manifest"
        )

    try:
        groups = metadata["groups"]
        group_edges = np.asarray(metadata["group_edges_ev"], dtype=np.float64)
        intervals = np.asarray(
            [group["energy_interval_ev"] for group in groups], dtype=np.float64  # type: ignore[index]
        )
        closure = metadata["group_spectral_closure"]
        group_energy = np.asarray(
            [group["photon_weighted_mean_energy_ev"] for group in groups], dtype=np.float64  # type: ignore[index]
        )
        cross_sections = closure["cross_sections_cm2"]  # type: ignore[index]
        excess = closure["photoelectron_excess_energy_ev"]  # type: ignore[index]
        sigma = np.asarray(
            [cross_sections[name] for name in ("hydrogen_i", "helium_i", "helium_ii")], dtype=np.float64
        )
        excess_array = np.asarray(
            [excess[name] for name in ("hydrogen_i", "helium_i", "helium_ii")], dtype=np.float64
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("photon metadata lacks a validated group_spectral_closure") from error
    if (
        group_edges.ndim != 1
        or len(group_edges) < 2
        or not np.isfinite(group_edges).all()
        or np.any(group_edges <= 0.0)
        or np.any(np.diff(group_edges) <= 0.0)
    ):
        raise ValueError("serialized photon group edges must be finite, positive, and increasing")
    if intervals.shape != (len(group_edges) - 1, 2) or not np.array_equal(
        intervals, np.column_stack((group_edges[:-1], group_edges[1:]))
    ):
        raise ValueError(
            "serialized photon group intervals must exactly match group_edges_ev"
        )
    if (
        group_energy.ndim != 1
        or len(group_energy) != len(group_edges) - 1
        or not np.isfinite(group_energy).all()
        or np.any(group_energy <= 0.0)
        or np.any(np.diff(group_energy) <= 0.0)
    ):
        raise ValueError("serialized group mean energies must be finite and strictly increasing")
    if sigma.shape != (3, len(group_energy)) or excess_array.shape != sigma.shape:
        raise ValueError("serialized group spectral closure has inconsistent group dimensions")
    if not np.isfinite(sigma).all() or np.any(sigma < 0.0) or not np.isfinite(excess_array).all() or np.any(excess_array < 0.0):
        raise ValueError("serialized group spectral closure contains invalid values")
    return GroupSpectralClosure(
        PhotoCrossSections(
            hydrogen_i=sigma[0],
            helium_i=sigma[1],
            helium_ii=sigma[2],
        ),
        group_energy,
        excess_array,
    )


def electron_number_density(state: PrimordialState) -> jnp.ndarray:
    """Return n_e implied by the H/He ion fractions."""
    return state.n_hydrogen * state.x_hydrogen_ii + state.n_helium * (
        state.x_helium_ii + 2.0 * state.x_helium_iii
    )


def neutral_number_densities(state: PrimordialState) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Return n_HI, n_HeI, n_HeII in cm^-3."""
    n_hi = state.n_hydrogen * (1.0 - state.x_hydrogen_ii)
    n_hei = state.n_helium * (1.0 - state.x_helium_ii - state.x_helium_iii)
    n_heii = state.n_helium * state.x_helium_ii
    return n_hi, n_hei, n_heii


def total_absorption_coefficient(state: PrimordialState, cross_sections: PhotoCrossSections) -> jnp.ndarray:
    """Return total number absorption coefficient [group, cell] in cm^-1."""
    n_hi, n_hei, n_heii = neutral_number_densities(state)
    extra_axes = (1,) * state.n_hydrogen.ndim
    sigma_hi = cross_sections.hydrogen_i.reshape((-1,) + extra_axes)
    sigma_hei = cross_sections.helium_i.reshape((-1,) + extra_axes)
    sigma_heii = cross_sections.helium_ii.reshape((-1,) + extra_axes)
    return sigma_hi * n_hi[None, ...] + sigma_hei * n_hei[None, ...] + sigma_heii * n_heii[None, ...]


def photoionization_rates(
    photon_number_density: jnp.ndarray,
    reduced_light_speed: float,
    cross_sections: PhotoCrossSections,
) -> PhotoRates:
    """Convert group photon densities [group, cell] into photo-rates [s^-1]."""
    photon_flux = reduced_light_speed * photon_number_density
    return PhotoRates(
        hydrogen_i=jnp.tensordot(cross_sections.hydrogen_i, photon_flux, axes=((0,), (0,))),
        helium_i=jnp.tensordot(cross_sections.helium_i, photon_flux, axes=((0,), (0,))),
        helium_ii=jnp.tensordot(cross_sections.helium_ii, photon_flux, axes=((0,), (0,))),
    )


def photoheating_rate(
    state: PrimordialState,
    photon_number_density: jnp.ndarray,
    group_energy_ev: jnp.ndarray,
    reduced_light_speed: float,
    cross_sections: PhotoCrossSections,
) -> jnp.ndarray:
    """Return primary photoelectron heating in erg cm^-3 s^-1.

    Secondary X-ray ionizations are intentionally excluded and will be supplied
    by the dedicated X-ray closure rather than folded into UV chemistry.
    """
    n_hi, n_hei, n_heii = neutral_number_densities(state)
    photon_flux = reduced_light_speed * photon_number_density
    energies = jnp.asarray(group_energy_ev)
    extra_axes = (1,) * state.n_hydrogen.ndim

    def species_heating(number_density: jnp.ndarray, sigma: jnp.ndarray, threshold_ev: float) -> jnp.ndarray:
        excess = jnp.maximum(energies - threshold_ev, 0.0).reshape((-1,) + extra_axes)
        opacity = sigma.reshape((-1,) + extra_axes) * number_density[None, ...]
        return jnp.sum(photon_flux * opacity * excess * EV_ERG, axis=0)

    return (
        species_heating(n_hi, cross_sections.hydrogen_i, H_I_FIT.threshold_ev)
        + species_heating(n_hei, cross_sections.helium_i, HE_I_FIT.threshold_ev)
        + species_heating(n_heii, cross_sections.helium_ii, HE_II_FIT.threshold_ev)
    )


def hui_gnedin_case_b_hydrogen(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """H II case-B recombination coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    lam = 315614.0 / temperature
    return 2.753e-14 * lam**1.5 / (1.0 + (lam / 2.740) ** 0.407) ** 2.242


def hui_gnedin_case_a_helium_ii_radiative(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """He II -> He I radiative case-A coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    lambda_helium_i = 2.0 * 285335.0 / temperature
    return 3.0e-14 * lambda_helium_i**0.654


def hui_gnedin_case_b_helium_ii_radiative(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """He II -> He I radiative case-B coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    lambda_helium_i = 2.0 * 285335.0 / temperature
    return 1.26e-14 * lambda_helium_i**0.75


def helium_ii_dielectronic_recombination(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """He II -> He I dielectronic coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    coefficient = 1.9e-3 / temperature**1.5 * jnp.exp(-4.7e5 / temperature)
    return coefficient * (1.0 + 0.3 * jnp.exp(-9.4e4 / temperature))


def hui_gnedin_case_b_helium_iii(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """He III -> He II hydrogenic case-B coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    return 2.0 * hui_gnedin_case_b_hydrogen(temperature / 4.0)


def case_b_helium_recombination(temperature_k: jnp.ndarray) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return total He II and radiative He III case-B coefficients.

    The He II coefficient is the Hui--Gnedin radiative case-B rate plus the
    separate dielectronic contribution. He III uses the hydrogenic scaling
    ``alpha_HeIII,B(T) = 2 alpha_HII,B(T/4)``.
    """
    return (
        hui_gnedin_case_b_helium_ii_radiative(temperature_k)
        + helium_ii_dielectronic_recombination(temperature_k),
        hui_gnedin_case_b_helium_iii(temperature_k),
    )


def evolve_primordial_fractions(
    state: PrimordialState,
    photo_rates: PhotoRates,
    temperature_k: jnp.ndarray,
    dt: float,
) -> PrimordialState:
    """Advance H/He ion fractions with conservative finite-time transitions.

    Electron density is lagged by one chemistry substep. This explicit closure is
    for P1 benchmarks; dense gas will use the planned local implicit Newton step.
    """
    n_electron = electron_number_density(state)
    alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
    alpha_heii, alpha_heiii = case_b_helium_recombination(temperature_k)

    probability_hi = -jnp.expm1(-photo_rates.hydrogen_i * dt)
    probability_hii = -jnp.expm1(-alpha_hii * n_electron * dt)
    x_hii = state.x_hydrogen_ii * (1.0 - probability_hii) + (1.0 - state.x_hydrogen_ii) * probability_hi

    x_hei = 1.0 - state.x_helium_ii - state.x_helium_iii
    rate_heii_out = photo_rates.helium_ii + alpha_heii * n_electron
    probability_heii_out = -jnp.expm1(-rate_heii_out * dt)
    heii_to_heiii = (
        state.x_helium_ii
        * probability_heii_out
        * photo_rates.helium_ii
        / jnp.maximum(rate_heii_out, jnp.finfo(rate_heii_out.dtype).tiny)
    )
    heii_to_hei = state.x_helium_ii * probability_heii_out - heii_to_heiii
    hei_to_heii = x_hei * (-jnp.expm1(-photo_rates.helium_i * dt))
    heiii_to_heii = state.x_helium_iii * (-jnp.expm1(-alpha_heiii * n_electron * dt))

    x_heii = state.x_helium_ii - heii_to_hei - heii_to_heiii + hei_to_heii + heiii_to_heii
    x_heiii = state.x_helium_iii - heiii_to_heii + heii_to_heiii
    return PrimordialState(
        n_hydrogen=state.n_hydrogen,
        n_helium=state.n_helium,
        x_hydrogen_ii=jnp.clip(x_hii, 0.0, 1.0),
        x_helium_ii=jnp.clip(x_heii, 0.0, 1.0),
        x_helium_iii=jnp.clip(x_heiii, 0.0, 1.0),
    )
