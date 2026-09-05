"""Dust absorption primitives and audited group-opacity loading."""

from __future__ import annotations

import json
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
    SOURCE_SED_CANDIDATE_STATUS,
    SOURCE_SED_CONTRACT_STATUSES,
    SED_INTERPOLATION_CONVENTION,
    SED_QUADRATURE_SCHEME,
)


EV_ERG = 1.602176634e-12
LIGHT_SPEED_CM_S = 2.99792458e10
DUST_BINDING_STATUSES = frozenset(("reference_control", "candidate_source_sed_matched"))
SNRT_ROOT = Path(__file__).resolve().parents[1]
DUST_CLOSURE_CODE_MANIFEST = {
    "dust_builder": SNRT_ROOT / "tools" / "build_draine_dust_opacity.py",
    "source_sed": SNRT_ROOT / "snrt_core" / "sed.py",
    "dust_loader": SNRT_ROOT / "snrt_core" / "dust.py",
    "integrity_helper": SNRT_ROOT / "snrt_core" / "provenance.py",
}


class DustModel(NamedTuple):
    """Absorbing dust scaled to a supplied reference cross section per H.

    ``absorption_cross_section_per_h`` is [group] in cm^2 per H nucleus for
    the reference dust mixture. ``relative_abundance`` is a non-negative
    cell field relative to that mixture; it may encode metallicity and a
    dust-to-metal prescription outside the transport kernel.
    ``absorption_weighted_energy_ev`` is the per-group photon energy used for
    dust heating and absorption momentum; it is optional for compatibility
    with legacy callers, which fall back to the supplied transport group
    energy.
    """

    absorption_cross_section_per_h: jnp.ndarray
    relative_abundance: jnp.ndarray
    absorption_weighted_energy_ev: jnp.ndarray | None = None


class DustOpacityClosure(NamedTuple):
    """Validated, group-averaged dust opacity metadata.

    The cross section is per H nucleus for the reference dust mixture.  The
    energy is weighted by the same dust absorption opacity and is therefore
    the appropriate energy per absorbed dust photon for local heating.
    """

    group_edges_ev: np.ndarray
    absorption_cross_section_per_h_cm2: np.ndarray
    absorption_weighted_energy_ev: np.ndarray
    schema: str = "snrt_dust_opacity_v1"
    source_sed_identity: str | None = None
    source_sed_sha256: str | None = None
    binding_status: str = "reference_control"
    group_edges_sha256: str | None = None
    source_table_sha256: str | None = None
    builder_sha256: str | None = None
    payload_sha256: str | None = None


def read_dust_opacity_metadata(
    path: str | Path,
    *,
    expected_group_edges_ev: np.ndarray | None = None,
    expected_group_edges_sha256: str | None = None,
    expected_source_sed_identity: str | None = None,
    expected_source_sed_sha256: str | None = None,
    require_source_match: bool = False,
) -> DustOpacityClosure:
    """Read and validate a source-SED-dependent dust opacity closure.

    The JSON schema is ``snrt_dust_opacity_v1`` or the source-bound
    ``snrt_dust_opacity_v2`` and requires
    ``group_edges_ev``, ``absorption_cross_section_per_h_cm2``, and
    ``absorption_weighted_energy_ev`` plus non-empty ``reference_mixture``,
    ``opacity_source``, and ``spectral_weighting`` provenance strings.  No
    opacity normalization or group ordering is inferred.  A v2 sidecar also
    carries a source-SED identity.  If ``require_source_match`` is true, that
    identity must match the photon-ledger identity; v1 remains available only
    as an explicitly labeled reference control.
    ``expected_group_edges_ev`` should be supplied by the photon-ledger
    metadata when the closure is attached to a run.
    """

    opacity_path = Path(path)
    with opacity_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict):
        raise ValueError(f"{opacity_path}: dust metadata root must be an object")
    schema = metadata.get("schema")
    if schema not in ("snrt_dust_opacity_v1", "snrt_dust_opacity_v2"):
        raise ValueError(f"{opacity_path}: unsupported dust opacity schema")
    binding_status = metadata.get("status", "reference_control")
    if binding_status not in DUST_BINDING_STATUSES:
        raise ValueError(
            f"{opacity_path}: dust binding status must be one of {sorted(DUST_BINDING_STATUSES)}"
        )
    if schema == "snrt_dust_opacity_v1" and binding_status != "reference_control":
        raise ValueError(f"{opacity_path}: v1 dust metadata must use reference_control status")
    if schema == "snrt_dust_opacity_v2" and binding_status != "candidate_source_sed_matched":
        raise ValueError(
            f"{opacity_path}: v2 source-bound dust metadata must use candidate_source_sed_matched status"
        )
    group_edges_sha256 = None
    source_table_sha256 = None
    builder_sha256 = None
    payload_sha256 = None
    if schema == "snrt_dust_opacity_v2":
        source_identity = metadata.get("source_sed_identity")
        source_hash = metadata.get("source_sed_sha256")
        source_identity = require_sha256(source_identity, "source_sed_identity", opacity_path)
        source_hash = require_sha256(source_hash, "source_sed_sha256", opacity_path)
        source_contract = metadata.get("source_sed_contract")
        if (
            not isinstance(source_contract, dict)
            or source_contract.get("identity") != source_identity
            or source_contract.get("input_sha256") != source_hash
        ):
            raise ValueError(f"{opacity_path}: source SED contract does not match its identity fields")
        if source_contract.get("status") not in SOURCE_SED_CONTRACT_STATUSES:
            raise ValueError(f"{opacity_path}: source SED contract has an invalid status")
        if source_contract.get("status") != SOURCE_SED_CANDIDATE_STATUS:
            raise ValueError(f"{opacity_path}: source-bound dust requires the candidate source SED status")
        if source_contract.get("status") == SOURCE_SED_CANDIDATE_STATUS:
            if source_contract.get("interpolation_convention") != SED_INTERPOLATION_CONVENTION:
                raise ValueError(f"{opacity_path}: source SED contract has an unsupported interpolation convention")
            quadrature = source_contract.get("quadrature")
            if not isinstance(quadrature, dict) or quadrature.get("scheme") != SED_QUADRATURE_SCHEME:
                raise ValueError(f"{opacity_path}: source SED contract lacks the declared quadrature scheme")
        source_input_path = source_contract.get("input_path")
        if not isinstance(source_input_path, str) or not source_input_path.strip():
            raise ValueError(f"{opacity_path}: source SED contract lacks input_path")
        try:
            actual_source_hash = sha256_file(source_input_path)
        except (FileNotFoundError, OSError) as error:
            raise ValueError(f"{opacity_path}: source SED input file is unavailable") from error
        if actual_source_hash != source_hash:
            raise ValueError(f"{opacity_path}: source SED input hash does not match its file")
        if require_source_match and expected_source_sed_identity != source_identity:
            raise ValueError(f"{opacity_path}: dust source SED identity does not match photon metadata")
        if expected_source_sed_sha256 is not None and expected_source_sed_sha256 != source_hash:
            raise ValueError(f"{opacity_path}: dust source SED input hash does not match photon metadata")

        group_edges_sha256 = require_sha256(
            metadata.get("group_edges_sha256"), "group_edges_sha256", opacity_path
        )
        group_edges_path = metadata.get("group_edges_path")
        if not isinstance(group_edges_path, str) or not group_edges_path.strip():
            raise ValueError(f"{opacity_path}: source-bound sidecar lacks group_edges_path")
        try:
            actual_edges_hash = sha256_file(group_edges_path)
        except (FileNotFoundError, OSError) as error:
            raise ValueError(f"{opacity_path}: group-edge input file is unavailable") from error
        if actual_edges_hash != group_edges_sha256:
            raise ValueError(f"{opacity_path}: group-edge hash does not match its file")
        if expected_group_edges_sha256 is not None and expected_group_edges_sha256 != group_edges_sha256:
            raise ValueError(f"{opacity_path}: dust group-edge hash does not match photon metadata")

        source_table = metadata.get("source_table")
        if not isinstance(source_table, dict):
            raise ValueError(f"{opacity_path}: source-bound sidecar lacks source_table provenance")
        source_table_path = source_table.get("path")
        source_table_sha256 = require_sha256(
            source_table.get("sha256"), "source_table.sha256", opacity_path
        )
        if not isinstance(source_table_path, str) or not source_table_path.strip():
            raise ValueError(f"{opacity_path}: source_table lacks path")
        try:
            actual_table_hash = sha256_file(source_table_path)
        except (FileNotFoundError, OSError) as error:
            raise ValueError(f"{opacity_path}: Draine source table is unavailable") from error
        if actual_table_hash != source_table_sha256:
            raise ValueError(f"{opacity_path}: Draine source-table hash does not match its file")

        builder = metadata.get("builder")
        if not isinstance(builder, dict):
            raise ValueError(f"{opacity_path}: source-bound sidecar lacks builder provenance")
        builder_path = builder.get("path")
        builder_sha256 = require_sha256(builder.get("sha256"), "builder.sha256", opacity_path)
        if not isinstance(builder_path, str) or not builder_path.strip():
            raise ValueError(f"{opacity_path}: builder provenance lacks path")
        try:
            actual_builder_hash = sha256_file(builder_path)
        except (FileNotFoundError, OSError) as error:
            raise ValueError(f"{opacity_path}: builder source is unavailable") from error
        if actual_builder_hash != builder_sha256:
            raise ValueError(f"{opacity_path}: builder hash does not match its file")
        validate_code_manifest(
            metadata,
            DUST_CLOSURE_CODE_MANIFEST,
            context=opacity_path,
        )
        payload_sha256 = validate_payload_hash(metadata, context=opacity_path)
    elif require_source_match and expected_source_sed_identity is not None:
        raise ValueError(
            f"{opacity_path}: v1 reference dust closure cannot be used with a source-bound photon ledger"
        )
    required = (
        "group_edges_ev",
        "absorption_cross_section_per_h_cm2",
        "absorption_weighted_energy_ev",
    )
    missing = [name for name in required if name not in metadata]
    if missing:
        raise ValueError(f"{opacity_path}: missing dust opacity fields {missing}")
    for name in ("reference_mixture", "opacity_source", "spectral_weighting"):
        if not isinstance(metadata.get(name), str) or not metadata[name].strip():
            raise ValueError(f"{opacity_path}: {name} must be a non-empty provenance string")

    edges = np.asarray(metadata["group_edges_ev"], dtype=np.float64)
    cross_section = np.asarray(metadata["absorption_cross_section_per_h_cm2"], dtype=np.float64)
    weighted_energy = np.asarray(metadata["absorption_weighted_energy_ev"], dtype=np.float64)
    if (
        edges.ndim != 1
        or edges.size < 2
        or not np.isfinite(edges).all()
        or np.any(edges <= 0.0)
        or np.any(np.diff(edges) <= 0.0)
    ):
        raise ValueError(f"{opacity_path}: group edges must be finite, positive, and increasing")
    number_of_groups = edges.size - 1
    if cross_section.shape != (number_of_groups,) or weighted_energy.shape != (number_of_groups,):
        raise ValueError(f"{opacity_path}: dust arrays do not match the number of groups")
    if (
        not np.isfinite(cross_section).all()
        or np.any(cross_section < 0.0)
        or not np.isfinite(weighted_energy).all()
        or np.any(weighted_energy <= 0.0)
    ):
        raise ValueError(f"{opacity_path}: dust cross sections/energies are invalid")
    tolerance = 1.0e-12 * np.maximum(1.0, edges[1:])
    if np.any(weighted_energy < edges[:-1] - tolerance) or np.any(weighted_energy > edges[1:] + tolerance):
        raise ValueError(f"{opacity_path}: absorption-weighted energies lie outside their groups")
    if schema == "snrt_dust_opacity_v2":
        serialized_fraction = np.asarray(
            metadata.get("source_sed_group_energy_fraction_of_lbol"),
            dtype=np.float64,
        )
        serialized_groups = metadata.get("groups")
        if (
            serialized_fraction.shape != (number_of_groups,)
            or not np.isfinite(serialized_fraction).all()
            or np.any(serialized_fraction < 0.0)
            or not isinstance(serialized_groups, list)
            or len(serialized_groups) != number_of_groups
        ):
            raise ValueError(f"{opacity_path}: source-bound group energy-fraction ledger is invalid")
        group_fraction = np.empty(number_of_groups, dtype=np.float64)
        for index, group in enumerate(serialized_groups):
            if not isinstance(group, dict):
                raise ValueError(f"{opacity_path}: source-bound group entries must be objects")
            interval = np.asarray(group.get("energy_interval_ev"), dtype=np.float64)
            if interval.shape != (2,) or not np.array_equal(
                interval, np.asarray((edges[index], edges[index + 1]), dtype=np.float64)
            ):
                raise ValueError(f"{opacity_path}: source-bound group intervals do not match group_edges_ev")
            try:
                group_fraction[index] = float(group["source_energy_fraction_of_lbol"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"{opacity_path}: source-bound group lacks energy fraction") from error
        if not np.isfinite(group_fraction).all() or np.any(group_fraction < 0.0):
            raise ValueError(f"{opacity_path}: source-bound group energy fractions are invalid")
        if not np.allclose(serialized_fraction, group_fraction, rtol=1.0e-14, atol=1.0e-15):
            raise ValueError(f"{opacity_path}: duplicated source energy fractions disagree")
    if expected_group_edges_ev is not None:
        expected = np.asarray(expected_group_edges_ev, dtype=np.float64)
        if expected.shape != edges.shape or not np.allclose(expected, edges, rtol=0.0, atol=1.0e-12):
            raise ValueError(f"{opacity_path}: dust groups do not match the photon-ledger groups")
    return DustOpacityClosure(
        edges,
        cross_section,
        weighted_energy,
        str(schema),
        metadata.get("source_sed_identity") if schema == "snrt_dust_opacity_v2" else None,
        metadata.get("source_sed_sha256") if schema == "snrt_dust_opacity_v2" else None,
        binding_status,
        group_edges_sha256,
        source_table_sha256,
        builder_sha256,
        payload_sha256,
    )


def dust_model_from_metadata(
    path: str | Path,
    relative_abundance: jnp.ndarray,
    *,
    dtype: jnp.dtype = jnp.float32,
    expected_group_edges_ev: np.ndarray | None = None,
    expected_group_edges_sha256: str | None = None,
    expected_source_sed_identity: str | None = None,
    expected_source_sed_sha256: str | None = None,
    require_source_match: bool = False,
) -> DustModel:
    """Build a JAX dust model from a validated opacity sidecar."""

    closure = read_dust_opacity_metadata(
        path,
        expected_group_edges_ev=expected_group_edges_ev,
        expected_group_edges_sha256=expected_group_edges_sha256,
        expected_source_sed_identity=expected_source_sed_identity,
        expected_source_sed_sha256=expected_source_sed_sha256,
        require_source_match=require_source_match,
    )
    abundance = jnp.asarray(relative_abundance, dtype=dtype)
    if not np.isfinite(np.asarray(abundance)).all() or np.any(np.asarray(abundance) < 0.0):
        raise ValueError("relative dust abundance must be finite and non-negative")
    return DustModel(
        absorption_cross_section_per_h=jnp.asarray(closure.absorption_cross_section_per_h_cm2, dtype=dtype),
        relative_abundance=abundance,
        absorption_weighted_energy_ev=jnp.asarray(closure.absorption_weighted_energy_ev, dtype=dtype),
    )


def absorption_coefficient(n_hydrogen: jnp.ndarray, dust: DustModel) -> jnp.ndarray:
    """Return dust absorption coefficient [group, cell] in cm^-1."""
    extra_axes = (1,) * n_hydrogen.ndim
    cross_section = dust.absorption_cross_section_per_h.reshape((-1,) + extra_axes)
    return cross_section * n_hydrogen[None, ...] * jnp.maximum(dust.relative_abundance[None, ...], 0.0)


def absorbed_dust_momentum_rate(
    absorbed_intensity: jnp.ndarray,
    dust_fraction: jnp.ndarray,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    absorption_weighted_energy_ev: jnp.ndarray,
    dt: float,
) -> jnp.ndarray:
    """Return dust absorption momentum deposition per volume and time.

    ``absorbed_intensity`` is the directional photon-number loss over one
    transport step, ordered ``[group, direction, x, y, z]``.  The returned
    vector is in dyn cm^-3 and uses the physical speed of light for photon
    momentum, even when transport uses a reduced light speed.  This is an
    absorption-only force diagnostic; scattering and IR re-emission are not
    included.
    """

    absorbed = jnp.asarray(absorbed_intensity)
    fraction = jnp.asarray(dust_fraction)
    direction = jnp.asarray(directions)
    angular_weight = jnp.asarray(weights)
    energy = jnp.asarray(absorption_weighted_energy_ev)
    if absorbed.ndim != 5 or fraction.shape != (absorbed.shape[0], *absorbed.shape[2:]):
        raise ValueError("absorbed intensity/dust fraction must have shapes [group,direction,x,y,z] and [group,x,y,z]")
    if direction.shape != (absorbed.shape[1], 3) or angular_weight.shape != (absorbed.shape[1],):
        raise ValueError("directions and weights must match the directional intensity axis")
    if energy.shape != (absorbed.shape[0],) or not np.isfinite(np.asarray(energy)).all() or np.any(np.asarray(energy) <= 0.0):
        raise ValueError("absorption-weighted photon energy must have shape (group,) and be positive")
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be positive and finite")
    energy_shape = (-1, 1, 1, 1, 1)
    dust_absorbed_energy = absorbed * fraction[:, None, ...] * energy.reshape(energy_shape) * EV_ERG
    return jnp.einsum(
        "d,gdxyz,di->igxyz",
        angular_weight,
        dust_absorbed_energy / (LIGHT_SPEED_CM_S * dt),
        direction,
    ).sum(axis=1)


def zero_dust(number_of_groups: int, shape: tuple[int, ...], dtype: jnp.dtype = jnp.float32) -> DustModel:
    """Return a no-dust model compatible with a static transport shape."""
    return DustModel(
        absorption_cross_section_per_h=jnp.zeros((number_of_groups,), dtype=dtype),
        relative_abundance=jnp.zeros(shape, dtype=dtype),
    )
