"""Dust absorption/scattering primitives and audited opacity loading."""

from __future__ import annotations

import json
from pathlib import Path
from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpy as np

from snrt_core.provenance import (
    PAYLOAD_HASH_SCHEME,
    canonical_payload_sha256,
    build_code_manifest,
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
DUST_EXTINCTION_ALBEDO_MAX_RELATIVE_ERROR = 1.0e-2
DUST_MOMENT_ROUNDING_ENVELOPE = 1.0e-4
DUST_BINDING_STATUSES = frozenset(
    (
        "reference_control",
        "candidate_source_sed_matched",
        "candidate_scattering_isotropic",
        "reference_scattering_control",
    )
)
SNRT_ROOT = Path(__file__).resolve().parents[1]
DUST_CLOSURE_CODE_MANIFEST = {
    "dust_builder": SNRT_ROOT / "tools" / "build_draine_dust_opacity.py",
    "source_sed": SNRT_ROOT / "snrt_core" / "sed.py",
    "dust_loader": SNRT_ROOT / "snrt_core" / "dust.py",
    "integrity_helper": SNRT_ROOT / "snrt_core" / "provenance.py",
}
DUST_THERMAL_CODE_MANIFEST = {
    "thermal_builder": SNRT_ROOT / "tools" / "build_draine_dust_thermal.py",
    "dust_loader": SNRT_ROOT / "snrt_core" / "dust.py",
    "integrity_helper": SNRT_ROOT / "snrt_core" / "provenance.py",
}


class DustModel(NamedTuple):
    """Dust scaled to supplied reference cross sections per H.

    ``absorption_cross_section_per_h`` is [group] in cm^2 per H nucleus for
    the reference dust mixture. ``relative_abundance`` is a non-negative
    cell field relative to that mixture; it may encode metallicity and a
    dust-to-metal prescription outside the transport kernel.
    ``absorption_weighted_energy_ev`` and ``scattering_weighted_energy_ev``
    are the per-group photon energies for their respective ledgers.  The
    scattering fields are optional so old absorption-only callers remain
    valid.  A non-``off`` ``scattering_phase_function`` is set only by a
    validated scattering sidecar or an explicit test fixture.
    """

    absorption_cross_section_per_h: jnp.ndarray
    relative_abundance: jnp.ndarray
    absorption_weighted_energy_ev: jnp.ndarray | None = None
    scattering_cross_section_per_h: jnp.ndarray | None = None
    scattering_weighted_energy_ev: jnp.ndarray | None = None
    scattering_phase_function: str = "off"


class DustOpacityClosure(NamedTuple):
    """Validated, group-averaged dust opacity metadata.

    The cross section is per H nucleus for the reference dust mixture.  The
    energy is weighted by the same dust absorption opacity and is therefore
    the appropriate energy per absorbed dust photon for local heating.
    """

    group_edges_ev: np.ndarray
    absorption_cross_section_per_h_cm2: np.ndarray
    absorption_weighted_energy_ev: np.ndarray
    scattering_cross_section_per_h_cm2: np.ndarray | None = None
    scattering_weighted_energy_ev: np.ndarray | None = None
    scattering_angle_cosine: np.ndarray | None = None
    scattering_angle_cosine_squared: np.ndarray | None = None
    transport_corrected_scattering_cross_section_per_h_cm2: np.ndarray | None = None
    isotropic_candidate_momentum_overestimate_factor: np.ndarray | None = None
    isotropic_candidate_momentum_bound_unbounded: np.ndarray | None = None
    scattering_phase_function: str = "off"
    schema: str = "snrt_dust_opacity_v1"
    source_sed_identity: str | None = None
    source_sed_sha256: str | None = None
    binding_status: str = "reference_control"
    group_edges_sha256: str | None = None
    source_table_sha256: str | None = None
    dust_mass_per_h_g: float | None = None
    builder_sha256: str | None = None
    payload_sha256: str | None = None


class DustThermalClosure(NamedTuple):
    """Validated equilibrium dust-emission table.

    The power curve is per reference-mixture H nucleus.  Fractions and photon
    energies are tabulated only for the configured IR groups; the complement
    is an explicit out-of-band energy ledger.
    """

    group_edges_ev: np.ndarray
    ir_group_indices: np.ndarray
    temperature_k: np.ndarray
    emitted_power_per_h_erg_s: np.ndarray
    ir_energy_fraction: np.ndarray
    ir_mean_photon_energy_ev: np.ndarray
    untracked_energy_fraction: np.ndarray
    reference_mixture: str
    thermal_source: str
    schema: str = "snrt_dust_thermal_v1"
    binding_status: str = "reference_thermal_control"
    group_edges_sha256: str | None = None
    source_table_sha256: str | None = None
    dust_mass_per_h_g: float | None = None
    builder_sha256: str | None = None
    payload_sha256: str | None = None


class DustThermalModel(NamedTuple):
    """JAX-ready equilibrium dust-emission table."""

    temperature_k: jnp.ndarray
    log_temperature_k: jnp.ndarray
    emitted_power_per_h_erg_s: jnp.ndarray
    log_emitted_power_per_h_erg_s: jnp.ndarray
    ir_energy_fraction: jnp.ndarray
    ir_mean_photon_energy_ev: jnp.ndarray
    untracked_energy_fraction: jnp.ndarray
    ir_group_indices: tuple[int, ...]


class DustThermalStep(NamedTuple):
    """One local thermal/emission closure evaluation."""

    grain_temperature_k: jnp.ndarray
    reemitted_energy_rate: jnp.ndarray
    untracked_energy_rate: jnp.ndarray
    ir_photon_rate: jnp.ndarray
    equilibrium_power_residual_rate: jnp.ndarray
    out_of_range: jnp.ndarray


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

    The JSON schema is ``snrt_dust_opacity_v1``, the source-bound
    ``snrt_dust_opacity_v2``, or the scattering-enabled
    ``snrt_dust_opacity_v3`` and requires
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
    if schema not in ("snrt_dust_opacity_v1", "snrt_dust_opacity_v2", "snrt_dust_opacity_v3"):
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
    if schema == "snrt_dust_opacity_v3" and binding_status not in (
        "candidate_scattering_isotropic",
        "reference_scattering_control",
    ):
        raise ValueError(
            f"{opacity_path}: v3 scattering metadata has an invalid scattering status"
        )
    group_edges_sha256 = None
    source_table_sha256 = None
    dust_mass_per_h_g = None
    builder_sha256 = None
    payload_sha256 = None
    source_bound = schema == "snrt_dust_opacity_v2" or (
        schema == "snrt_dust_opacity_v3" and "source_sed_identity" in metadata
    )
    if schema == "snrt_dust_opacity_v3":
        if binding_status == "candidate_scattering_isotropic" and not source_bound:
            raise ValueError(f"{opacity_path}: source-bound v3 candidate lacks source SED identity")
        if binding_status == "reference_scattering_control" and source_bound:
            raise ValueError(f"{opacity_path}: reference v3 control cannot carry source SED identity")
    if source_bound:
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
        try:
            dust_mass_per_h_g = float(source_table.get("dust_mass_per_h_g"))
        except (TypeError, ValueError) as error:
            raise ValueError(f"{opacity_path}: source_table lacks dust_mass_per_h_g") from error
        if not np.isfinite(dust_mass_per_h_g) or dust_mass_per_h_g <= 0.0:
            raise ValueError(f"{opacity_path}: source_table dust_mass_per_h_g is invalid")
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
            f"{opacity_path}: unbound reference dust closure cannot be used with a source-bound photon ledger"
        )
    elif schema == "snrt_dust_opacity_v3" and any(
        key in metadata
        for key in (
            "group_edges_path",
            "group_edges_sha256",
            "source_table",
            "builder",
            "closure_code_manifest",
            "payload_hash_scheme",
            "payload_sha256",
        )
    ):
        # The unbound v3 reference control is still reproducibility-pinned;
        # it simply has no source-SED identity to match.
        reference_provenance_fields = (
            "group_edges_path",
            "group_edges_sha256",
            "source_table",
            "builder",
            "closure_code_manifest",
            "payload_hash_scheme",
            "payload_sha256",
        )
        missing_reference_provenance = [
            key for key in reference_provenance_fields if key not in metadata
        ]
        if missing_reference_provenance:
            raise ValueError(
                f"{opacity_path}: v3 reference control lacks provenance fields "
                f"{missing_reference_provenance}"
            )
        group_edges_sha256 = require_sha256(
            metadata.get("group_edges_sha256"), "group_edges_sha256", opacity_path
        )
        group_edges_path = metadata.get("group_edges_path")
        if not isinstance(group_edges_path, str) or not group_edges_path.strip():
            raise ValueError(f"{opacity_path}: v3 reference control lacks group_edges_path")
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
            raise ValueError(f"{opacity_path}: v3 reference control lacks source_table provenance")
        source_table_path = source_table.get("path")
        source_table_sha256 = require_sha256(source_table.get("sha256"), "source_table.sha256", opacity_path)
        try:
            dust_mass_per_h_g = float(source_table.get("dust_mass_per_h_g"))
        except (TypeError, ValueError) as error:
            raise ValueError(f"{opacity_path}: source_table lacks dust_mass_per_h_g") from error
        if not np.isfinite(dust_mass_per_h_g) or dust_mass_per_h_g <= 0.0:
            raise ValueError(f"{opacity_path}: source_table dust_mass_per_h_g is invalid")
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
            raise ValueError(f"{opacity_path}: v3 reference control lacks builder provenance")
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
        validate_code_manifest(metadata, DUST_CLOSURE_CODE_MANIFEST, context=opacity_path)
        payload_sha256 = validate_payload_hash(metadata, context=opacity_path)
    required = (
        "group_edges_ev",
        "absorption_cross_section_per_h_cm2",
        "absorption_weighted_energy_ev",
    )
    if schema == "snrt_dust_opacity_v3":
        required += (
            "phase_function",
            "scattering_cross_section_per_h_cm2",
            "scattering_weighted_energy_ev",
            "scattering_angle_cosine",
            "scattering_angle_cosine_squared",
            "transport_corrected_scattering_cross_section_per_h_cm2",
            "isotropic_candidate_momentum_overestimate_factor",
            "isotropic_candidate_momentum_bound_unbounded",
        )
    missing = [name for name in required if name not in metadata]
    if missing:
        raise ValueError(f"{opacity_path}: missing dust opacity fields {missing}")
    for name in ("reference_mixture", "opacity_source", "spectral_weighting"):
        if not isinstance(metadata.get(name), str) or not metadata[name].strip():
            raise ValueError(f"{opacity_path}: {name} must be a non-empty provenance string")
    source_table_metadata = metadata.get("source_table")
    if isinstance(source_table_metadata, dict) and "absorption_consistency_max_relative_error" in source_table_metadata:
        try:
            consistency_error = float(source_table_metadata["absorption_consistency_max_relative_error"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"{opacity_path}: extinction/albedo residual is not numeric") from error
        if not np.isfinite(consistency_error) or consistency_error > DUST_EXTINCTION_ALBEDO_MAX_RELATIVE_ERROR:
            raise ValueError(
                f"{opacity_path}: extinction/albedo residual exceeds "
                f"{DUST_EXTINCTION_ALBEDO_MAX_RELATIVE_ERROR:.3g}"
            )
    if isinstance(source_table_metadata, dict) and "moment_inequality_max_violation" in source_table_metadata:
        try:
            moment_error = float(source_table_metadata["moment_inequality_max_violation"])
            moment_envelope = float(source_table_metadata.get("moment_rounding_envelope", DUST_MOMENT_ROUNDING_ENVELOPE))
        except (TypeError, ValueError) as error:
            raise ValueError(f"{opacity_path}: moment residual/envelope is not numeric") from error
        if (
            not np.isfinite(moment_error)
            or not np.isfinite(moment_envelope)
            or moment_error < 0.0
            or moment_envelope < 0.0
            or moment_error > moment_envelope
            or moment_envelope > DUST_MOMENT_ROUNDING_ENVELOPE
        ):
            raise ValueError(f"{opacity_path}: moment residual exceeds its declared rounding envelope")

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
    scattering_cross_section = np.zeros(number_of_groups, dtype=np.float64)
    scattering_weighted_energy = None
    scattering_angle_cosine = None
    scattering_angle_cosine_squared = None
    transport_corrected_scattering = None
    isotropic_candidate_momentum_overestimate_factor = None
    isotropic_candidate_momentum_bound_unbounded = None
    scattering_phase_function = "off"
    if schema == "snrt_dust_opacity_v3":
        scattering_phase_function = metadata.get("phase_function")
        if scattering_phase_function != "phase_isotropic_candidate":
            raise ValueError(f"{opacity_path}: v3 dust metadata has an unsupported phase function")
        scattering_field_names = (
            "scattering_cross_section_per_h_cm2",
            "scattering_weighted_energy_ev",
            "scattering_angle_cosine",
            "scattering_angle_cosine_squared",
            "transport_corrected_scattering_cross_section_per_h_cm2",
        )
        if any(metadata.get(name) is None for name in scattering_field_names):
            raise ValueError(f"{opacity_path}: v3 scattering fields cannot be null")
        scattering_cross_section = np.asarray(
            metadata.get("scattering_cross_section_per_h_cm2"), dtype=np.float64
        )
        scattering_weighted_energy = np.asarray(
            metadata.get("scattering_weighted_energy_ev"), dtype=np.float64
        )
        scattering_angle_cosine = np.asarray(
            metadata.get("scattering_angle_cosine"), dtype=np.float64
        )
        scattering_angle_cosine_squared = np.asarray(
            metadata.get("scattering_angle_cosine_squared"), dtype=np.float64
        )
        transport_corrected_scattering = np.asarray(
            metadata.get("transport_corrected_scattering_cross_section_per_h_cm2"),
            dtype=np.float64,
        )
        arrays = (
            scattering_cross_section,
            scattering_weighted_energy,
            scattering_angle_cosine,
            scattering_angle_cosine_squared,
            transport_corrected_scattering,
        )
        if any(array.shape != (number_of_groups,) for array in arrays):
            raise ValueError(f"{opacity_path}: v3 scattering arrays do not match the number of groups")
        if (
            not all(np.isfinite(array).all() for array in arrays)
            or np.any(scattering_cross_section < 0.0)
            or np.any(scattering_weighted_energy <= 0.0)
            or np.any(scattering_angle_cosine < -1.0 - 1.0e-12)
            or np.any(scattering_angle_cosine > 1.0 + 1.0e-12)
            or np.any(scattering_angle_cosine_squared < -1.0e-12)
            or np.any(scattering_angle_cosine_squared > 1.0 + 1.0e-12)
            # The source table prints rounded moments; retain the values and
            # allow only the measured raw-column rounding envelope.
            or np.any(
                scattering_angle_cosine_squared
                + DUST_MOMENT_ROUNDING_ENVELOPE
                < scattering_angle_cosine**2
            )
            or np.any(transport_corrected_scattering < 0.0)
        ):
            raise ValueError(f"{opacity_path}: v3 scattering values are invalid")
        if np.any(scattering_weighted_energy < edges[:-1] - tolerance) or np.any(
            scattering_weighted_energy > edges[1:] + tolerance
        ):
            raise ValueError(f"{opacity_path}: scattering-weighted energies lie outside their groups")
        if not np.allclose(
            transport_corrected_scattering,
            scattering_cross_section * (1.0 - scattering_angle_cosine),
            rtol=1.0e-12,
            atol=1.0e-300,
        ):
            raise ValueError(f"{opacity_path}: transport-corrected scattering values are inconsistent")
        raw_factor = metadata["isotropic_candidate_momentum_overestimate_factor"]
        raw_unbounded = metadata["isotropic_candidate_momentum_bound_unbounded"]
        if not isinstance(raw_factor, list) or not isinstance(raw_unbounded, list):
            raise ValueError(f"{opacity_path}: isotropic anisotropy bound fields must be arrays")
        if any(not isinstance(value, bool) for value in raw_unbounded):
            raise ValueError(f"{opacity_path}: isotropic anisotropy bound flags must be booleans")
        if len(raw_factor) != number_of_groups or len(raw_unbounded) != number_of_groups:
            raise ValueError(f"{opacity_path}: isotropic anisotropy bound fields have the wrong length")
        isotropic_candidate_momentum_overestimate_factor = np.full(number_of_groups, np.nan)
        isotropic_candidate_momentum_bound_unbounded = np.asarray(raw_unbounded, dtype=bool)
        for index, (value, unbounded, cosine) in enumerate(
            zip(raw_factor, isotropic_candidate_momentum_bound_unbounded, scattering_angle_cosine, strict=True)
        ):
            expected_unbounded = bool(cosine >= 1.0 - 1.0e-12)
            if bool(unbounded) != expected_unbounded:
                raise ValueError(f"{opacity_path}: isotropic anisotropy bound disagrees with g")
            if expected_unbounded:
                if value is not None:
                    raise ValueError(f"{opacity_path}: unbounded isotropic anisotropy must be null")
            else:
                try:
                    factor = float(value)
                except (TypeError, ValueError) as error:
                    raise ValueError(f"{opacity_path}: isotropic anisotropy bound is not numeric") from error
                expected_factor = 1.0 / (1.0 - cosine)
                if not np.isfinite(factor) or factor < 1.0 or not np.isclose(
                    factor, expected_factor, rtol=1.0e-12, atol=1.0e-12
                ):
                    raise ValueError(f"{opacity_path}: isotropic anisotropy bound is inconsistent with g")
                isotropic_candidate_momentum_overestimate_factor[index] = factor
    if source_bound:
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
        group_edges_ev=edges,
        absorption_cross_section_per_h_cm2=cross_section,
        absorption_weighted_energy_ev=weighted_energy,
        scattering_cross_section_per_h_cm2=scattering_cross_section,
        scattering_weighted_energy_ev=scattering_weighted_energy,
        scattering_angle_cosine=scattering_angle_cosine,
        scattering_angle_cosine_squared=scattering_angle_cosine_squared,
        transport_corrected_scattering_cross_section_per_h_cm2=transport_corrected_scattering,
        isotropic_candidate_momentum_overestimate_factor=isotropic_candidate_momentum_overestimate_factor,
        isotropic_candidate_momentum_bound_unbounded=isotropic_candidate_momentum_bound_unbounded,
        scattering_phase_function=scattering_phase_function,
        schema=str(schema),
        source_sed_identity=metadata.get("source_sed_identity")
        if source_bound
        else None,
        source_sed_sha256=metadata.get("source_sed_sha256")
        if source_bound
        else None,
        binding_status=binding_status,
        group_edges_sha256=group_edges_sha256,
        source_table_sha256=source_table_sha256,
        dust_mass_per_h_g=dust_mass_per_h_g,
        builder_sha256=builder_sha256,
        payload_sha256=payload_sha256,
    )


def read_dust_thermal_metadata(
    path: str | Path,
    *,
    expected_group_edges_ev: np.ndarray | None = None,
    expected_group_edges_sha256: str | None = None,
    expected_source_table_sha256: str | None = None,
    expected_dust_mass_per_h_g: float | None = None,
) -> DustThermalClosure:
    """Read and validate the Kirchhoff-derived equilibrium dust table."""

    thermal_path = Path(path)
    with thermal_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict) or metadata.get("schema") != "snrt_dust_thermal_v1":
        raise ValueError(f"{thermal_path}: unsupported dust thermal schema")
    status = metadata.get("status")
    if status not in ("reference_thermal_control", "candidate_kirchhoff_equilibrium"):
        raise ValueError(f"{thermal_path}: invalid dust thermal status")
    required_strings = ("reference_mixture", "thermal_source", "single_temperature_assumption")
    for name in required_strings:
        if not isinstance(metadata.get(name), str) or not metadata[name].strip():
            raise ValueError(f"{thermal_path}: {name} must be a non-empty string")

    edges = np.asarray(metadata.get("group_edges_ev"), dtype=np.float64)
    if (
        edges.ndim != 1
        or edges.size < 2
        or not np.isfinite(edges).all()
        or np.any(edges <= 0.0)
        or np.any(np.diff(edges) <= 0.0)
    ):
        raise ValueError(f"{thermal_path}: group edges are invalid")
    if expected_group_edges_ev is not None:
        expected_edges = np.asarray(expected_group_edges_ev, dtype=np.float64)
        if expected_edges.shape != edges.shape or not np.array_equal(expected_edges, edges):
            raise ValueError(f"{thermal_path}: thermal groups do not match the photon-ledger groups")
    group_edges_sha256 = require_sha256(
        metadata.get("group_edges_sha256"), "group_edges_sha256", thermal_path
    )
    group_edges_path = metadata.get("group_edges_path")
    if not isinstance(group_edges_path, str) or not group_edges_path.strip():
        raise ValueError(f"{thermal_path}: thermal sidecar lacks group_edges_path")
    try:
        if sha256_file(group_edges_path) != group_edges_sha256:
            raise ValueError(f"{thermal_path}: group-edge hash does not match its file")
    except (FileNotFoundError, OSError) as error:
        raise ValueError(f"{thermal_path}: group-edge input file is unavailable") from error
    if expected_group_edges_sha256 is not None and expected_group_edges_sha256 != group_edges_sha256:
        raise ValueError(f"{thermal_path}: thermal group-edge hash does not match the photon metadata")

    indices_raw = metadata.get("ir_group_indices")
    if not isinstance(indices_raw, list) or not indices_raw:
        raise ValueError(f"{thermal_path}: ir_group_indices must be a non-empty list")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in indices_raw):
        raise ValueError(f"{thermal_path}: IR group indices must be integers")
    ir_indices = np.asarray(indices_raw, dtype=np.int64)
    if np.any(ir_indices < 0) or np.any(ir_indices >= edges.size - 1) or np.unique(ir_indices).size != ir_indices.size:
        raise ValueError(f"{thermal_path}: IR group indices are invalid or duplicated")

    temperature = np.asarray(metadata.get("temperature_k"), dtype=np.float64)
    power = np.asarray(metadata.get("emitted_power_per_h_erg_s"), dtype=np.float64)
    fractions = np.asarray(metadata.get("ir_energy_fraction"), dtype=np.float64)
    photon_energy = np.asarray(metadata.get("ir_mean_photon_energy_ev"), dtype=np.float64)
    untracked = np.asarray(metadata.get("untracked_energy_fraction"), dtype=np.float64)
    number_of_temperatures = temperature.size
    number_of_ir_groups = ir_indices.size
    if (
        temperature.ndim != 1
        or number_of_temperatures < 2
        or power.shape != temperature.shape
        or fractions.shape != (number_of_temperatures, number_of_ir_groups)
        or photon_energy.shape != fractions.shape
        or untracked.shape != temperature.shape
    ):
        raise ValueError(f"{thermal_path}: thermal arrays have inconsistent shapes")
    if (
        not np.isfinite(temperature).all()
        or np.any(temperature <= 0.0)
        or np.any(np.diff(temperature) <= 0.0)
        or not np.isfinite(power).all()
        or np.any(power <= 0.0)
        or np.any(np.diff(power) <= 0.0)
        or not np.isfinite(fractions).all()
        or np.any(fractions < 0.0)
        or not np.isfinite(photon_energy).all()
        or np.any(photon_energy <= 0.0)
        or not np.isfinite(untracked).all()
        or np.any(untracked < 0.0)
    ):
        raise ValueError(f"{thermal_path}: thermal arrays are non-finite or out of bounds")
    fraction_tolerance = float(metadata.get("fraction_tolerance", 1.0e-10))
    if not np.isfinite(fraction_tolerance) or fraction_tolerance <= 0.0:
        raise ValueError(f"{thermal_path}: fraction_tolerance is invalid")
    fraction_sum = fractions.sum(axis=1) + untracked
    if not np.allclose(fraction_sum, 1.0, rtol=0.0, atol=fraction_tolerance):
        raise ValueError(f"{thermal_path}: tracked plus untracked energy fractions do not close")
    lower = edges[ir_indices]
    upper = edges[ir_indices + 1]
    energy_tolerance = 1.0e-12 * np.maximum(1.0, upper)
    if np.any(photon_energy < lower[None, :] - energy_tolerance[None, :]) or np.any(
        photon_energy > upper[None, :] + energy_tolerance[None, :]
    ):
        raise ValueError(f"{thermal_path}: thermal photon energies lie outside their IR groups")

    source_table = metadata.get("source_table")
    if not isinstance(source_table, dict):
        raise ValueError(f"{thermal_path}: thermal sidecar lacks source_table provenance")
    source_table_path = source_table.get("path")
    source_table_sha256 = require_sha256(
        source_table.get("sha256"), "source_table.sha256", thermal_path
    )
    if not isinstance(source_table_path, str) or not source_table_path.strip():
        raise ValueError(f"{thermal_path}: source_table lacks path")
    try:
        if sha256_file(source_table_path) != source_table_sha256:
            raise ValueError(f"{thermal_path}: source-table hash does not match its file")
    except (FileNotFoundError, OSError) as error:
        raise ValueError(f"{thermal_path}: source-table input is unavailable") from error
    try:
        dust_mass_per_h_g = float(source_table.get("dust_mass_per_h_g"))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{thermal_path}: source_table lacks dust_mass_per_h_g") from error
    if not np.isfinite(dust_mass_per_h_g) or dust_mass_per_h_g <= 0.0:
        raise ValueError(f"{thermal_path}: source_table dust_mass_per_h_g is invalid")
    if expected_source_table_sha256 is not None and source_table_sha256 != expected_source_table_sha256:
        raise ValueError(f"{thermal_path}: thermal and opacity source tables do not match")
    if expected_dust_mass_per_h_g is not None and not np.isclose(
        dust_mass_per_h_g, expected_dust_mass_per_h_g, rtol=1.0e-12, atol=0.0
    ):
        raise ValueError(f"{thermal_path}: thermal and opacity dust masses per H do not match")

    builder = metadata.get("builder")
    if not isinstance(builder, dict):
        raise ValueError(f"{thermal_path}: thermal sidecar lacks builder provenance")
    builder_path = builder.get("path")
    if not isinstance(builder_path, str) or not builder_path.strip():
        raise ValueError(f"{thermal_path}: builder provenance lacks path")
    builder_sha256 = require_sha256(builder.get("sha256"), "builder.sha256", thermal_path)
    if Path(builder_path).resolve() != DUST_THERMAL_CODE_MANIFEST["thermal_builder"].resolve():
        raise ValueError(f"{thermal_path}: thermal builder path is not the canonical builder")
    if sha256_file(builder_path) != builder_sha256:
        raise ValueError(f"{thermal_path}: thermal builder hash does not match its file")
    validate_code_manifest(metadata, DUST_THERMAL_CODE_MANIFEST, context=thermal_path)
    payload_sha256 = validate_payload_hash(metadata, context=thermal_path)
    return DustThermalClosure(
        group_edges_ev=edges,
        ir_group_indices=ir_indices,
        temperature_k=temperature,
        emitted_power_per_h_erg_s=power,
        ir_energy_fraction=fractions,
        ir_mean_photon_energy_ev=photon_energy,
        untracked_energy_fraction=untracked,
        reference_mixture=str(metadata["reference_mixture"]),
        thermal_source=str(metadata["thermal_source"]),
        binding_status=str(status),
        group_edges_sha256=group_edges_sha256,
        source_table_sha256=source_table_sha256,
        dust_mass_per_h_g=dust_mass_per_h_g,
        builder_sha256=str(metadata["builder"]["sha256"]),
        payload_sha256=payload_sha256,
    )


def dust_thermal_model_from_metadata(
    path: str | Path,
    *,
    dtype: jnp.dtype = jnp.float32,
    expected_group_edges_ev: np.ndarray | None = None,
    expected_group_edges_sha256: str | None = None,
    expected_source_table_sha256: str | None = None,
    expected_dust_mass_per_h_g: float | None = None,
) -> DustThermalModel:
    """Build a fixed-shape JAX thermal model from a validated sidecar."""

    closure = read_dust_thermal_metadata(
        path,
        expected_group_edges_ev=expected_group_edges_ev,
        expected_group_edges_sha256=expected_group_edges_sha256,
        expected_source_table_sha256=expected_source_table_sha256,
        expected_dust_mass_per_h_g=expected_dust_mass_per_h_g,
    )
    number_of_groups = closure.group_edges_ev.size - 1
    full_fraction = np.zeros((closure.temperature_k.size, number_of_groups), dtype=np.float64)
    full_photon_energy = np.ones_like(full_fraction)
    full_fraction[:, closure.ir_group_indices] = closure.ir_energy_fraction
    full_photon_energy[:, closure.ir_group_indices] = closure.ir_mean_photon_energy_ev
    return DustThermalModel(
        temperature_k=jnp.asarray(closure.temperature_k, dtype=dtype),
        log_temperature_k=jnp.log(jnp.asarray(closure.temperature_k, dtype=dtype)),
        emitted_power_per_h_erg_s=jnp.asarray(closure.emitted_power_per_h_erg_s, dtype=dtype),
        log_emitted_power_per_h_erg_s=jnp.log(
            jnp.asarray(closure.emitted_power_per_h_erg_s, dtype=dtype)
        ),
        ir_energy_fraction=jnp.asarray(full_fraction, dtype=dtype),
        ir_mean_photon_energy_ev=jnp.asarray(full_photon_energy, dtype=dtype),
        untracked_energy_fraction=jnp.asarray(closure.untracked_energy_fraction, dtype=dtype),
        ir_group_indices=tuple(int(value) for value in closure.ir_group_indices),
    )


def _log_linear_interpolate(
    query_log: jnp.ndarray,
    grid_log: jnp.ndarray,
    values: jnp.ndarray,
) -> jnp.ndarray:
    """Interpolate a one- or multi-column table in log-temperature space."""

    index = jnp.searchsorted(grid_log, query_log, side="right") - 1
    index = jnp.clip(index, 0, grid_log.size - 2)
    lower = grid_log[index]
    upper = grid_log[index + 1]
    weight = (query_log - lower) / jnp.maximum(upper - lower, jnp.finfo(query_log.dtype).tiny)
    if values.ndim == 1:
        lower_value = values[index]
        upper_value = values[index + 1]
    else:
        lower_value = values[index, ...]
        upper_value = values[index + 1, ...]
    return lower_value + weight[..., None] * (upper_value - lower_value) if values.ndim > 1 else lower_value + weight * (upper_value - lower_value)


def evaluate_dust_thermal(
    model: DustThermalModel,
    absorbed_dust_power_rate: jnp.ndarray,
    n_hydrogen: jnp.ndarray,
    relative_abundance: jnp.ndarray,
    background_temperature_k: float,
) -> DustThermalStep:
    """Evaluate local equilibrium and return a one-pass IR source ledger."""

    absorbed = jnp.asarray(absorbed_dust_power_rate)
    hydrogen = jnp.asarray(n_hydrogen)
    abundance = jnp.asarray(relative_abundance)
    if absorbed.shape != hydrogen.shape or abundance.shape != hydrogen.shape:
        raise ValueError("dust thermal fields must share the cell shape")
    if not np.isfinite(background_temperature_k) or background_temperature_k <= 0.0:
        raise ValueError("background_temperature_k must be positive and finite")
    dtype = absorbed.dtype
    background_temperature = jnp.asarray(background_temperature_k, dtype=dtype)
    background_log_temperature = jnp.log(background_temperature)
    table_min = model.temperature_k[0]
    table_max = model.temperature_k[-1]
    background_in_range = (background_temperature >= table_min) & (background_temperature <= table_max)
    background_power = _log_linear_interpolate(
        jnp.clip(background_log_temperature, model.log_temperature_k[0], model.log_temperature_k[-1]),
        model.log_temperature_k,
        model.emitted_power_per_h_erg_s,
    )
    dust_density = jnp.maximum(hydrogen * jnp.maximum(abundance, 0.0), 0.0)
    invalid_input = (
        (~jnp.isfinite(absorbed))
        | (~jnp.isfinite(hydrogen))
        | (~jnp.isfinite(abundance))
        | (absorbed < 0.0)
        | (hydrogen <= 0.0)
        | (abundance < 0.0)
    )
    local_power = jnp.maximum(jnp.where(jnp.isfinite(absorbed), absorbed, 0.0), 0.0)
    active = (dust_density > 0.0) & (local_power > 0.0)
    local_per_h = jnp.where(active, local_power / jnp.maximum(dust_density, jnp.finfo(dtype).tiny), 0.0)
    target_power = background_power + local_per_h
    target_in_range = (target_power >= model.emitted_power_per_h_erg_s[0]) & (
        target_power <= model.emitted_power_per_h_erg_s[-1]
    )
    safe_target = jnp.clip(
        target_power,
        model.emitted_power_per_h_erg_s[0],
        model.emitted_power_per_h_erg_s[-1],
    )
    lower_log_temperature = jnp.full_like(target_power, model.log_temperature_k[0])
    upper_log_temperature = jnp.full_like(target_power, model.log_temperature_k[-1])

    def bisect_temperature(
        _: int,
        bounds: tuple[jnp.ndarray, jnp.ndarray],
    ) -> tuple[jnp.ndarray, jnp.ndarray]:
        lower_bound, upper_bound = bounds
        midpoint = 0.5 * (lower_bound + upper_bound)
        midpoint_power = _log_linear_interpolate(
            midpoint,
            model.log_temperature_k,
            model.emitted_power_per_h_erg_s,
        )
        return (
            jnp.where(midpoint_power < safe_target, midpoint, lower_bound),
            jnp.where(midpoint_power < safe_target, upper_bound, midpoint),
        )

    lower_log_temperature, upper_log_temperature = jax.lax.fori_loop(
        0,
        32,
        bisect_temperature,
        (lower_log_temperature, upper_log_temperature),
    )
    solved_log_temperature = 0.5 * (lower_log_temperature + upper_log_temperature)
    solved_temperature = jnp.exp(solved_log_temperature)
    solved_temperature = jnp.where(active, solved_temperature, 0.0)
    evaluation_temperature = jnp.where(active, solved_temperature, model.temperature_k[0])
    evaluated_power = _log_linear_interpolate(
        jnp.log(jnp.maximum(evaluation_temperature, model.temperature_k[0])),
        model.log_temperature_k,
        model.emitted_power_per_h_erg_s,
    )
    fractions = _log_linear_interpolate(
        jnp.log(jnp.maximum(evaluation_temperature, model.temperature_k[0])),
        model.log_temperature_k,
        model.ir_energy_fraction,
    )
    photon_energy = _log_linear_interpolate(
        jnp.log(jnp.maximum(evaluation_temperature, model.temperature_k[0])),
        model.log_temperature_k,
        model.ir_mean_photon_energy_ev,
    )
    untracked = _log_linear_interpolate(
        jnp.log(jnp.maximum(evaluation_temperature, model.temperature_k[0])),
        model.log_temperature_k,
        model.untracked_energy_fraction,
    )
    emitted_excess_per_h = jnp.where(
        active,
        jnp.maximum(evaluated_power - background_power, 0.0),
        0.0,
    )
    emitted_excess_rate = emitted_excess_per_h * dust_density
    ir_energy_rate = jnp.moveaxis(emitted_excess_rate[..., None] * fractions, -1, 0)
    ir_photon_rate = ir_energy_rate / jnp.maximum(
        jnp.moveaxis(photon_energy, -1, 0) * EV_ERG,
        jnp.finfo(dtype).tiny,
    )
    untracked_rate = emitted_excess_rate * untracked
    residual_rate = jnp.where(
        active,
        (evaluated_power - target_power) * dust_density,
        0.0,
    )
    out_of_range = invalid_input | (active & ((~background_in_range) | (~target_in_range)))
    return DustThermalStep(
        grain_temperature_k=solved_temperature,
        reemitted_energy_rate=jnp.sum(ir_energy_rate, axis=0),
        untracked_energy_rate=untracked_rate,
        ir_photon_rate=ir_photon_rate,
        equilibrium_power_residual_rate=residual_rate,
        out_of_range=out_of_range,
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
        scattering_cross_section_per_h=(
            None
            if closure.scattering_phase_function == "off"
            else jnp.asarray(closure.scattering_cross_section_per_h_cm2, dtype=dtype)
        ),
        scattering_weighted_energy_ev=(
            None
            if closure.scattering_phase_function == "off"
            else jnp.asarray(closure.scattering_weighted_energy_ev, dtype=dtype)
        ),
        scattering_phase_function=closure.scattering_phase_function,
    )


def absorption_coefficient(n_hydrogen: jnp.ndarray, dust: DustModel) -> jnp.ndarray:
    """Return dust absorption coefficient [group, cell] in cm^-1."""
    extra_axes = (1,) * n_hydrogen.ndim
    cross_section = dust.absorption_cross_section_per_h.reshape((-1,) + extra_axes)
    return cross_section * n_hydrogen[None, ...] * jnp.maximum(dust.relative_abundance[None, ...], 0.0)


def scattering_coefficient(n_hydrogen: jnp.ndarray, dust: DustModel) -> jnp.ndarray:
    """Return the dust scattering coefficient [group, cell] in cm^-1."""

    if dust.scattering_cross_section_per_h is None:
        number_of_groups = dust.absorption_cross_section_per_h.shape[0]
        return jnp.zeros((number_of_groups, *n_hydrogen.shape), dtype=n_hydrogen.dtype)
    extra_axes = (1,) * n_hydrogen.ndim
    cross_section = dust.scattering_cross_section_per_h.reshape((-1,) + extra_axes)
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


def scattered_dust_momentum_rate(
    incoming_scattering: jnp.ndarray,
    outgoing_scattering: jnp.ndarray,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    scattering_weighted_energy_ev: jnp.ndarray,
    dt: float,
) -> jnp.ndarray:
    """Return momentum transfer from dust scattering alone."""

    incoming = jnp.asarray(incoming_scattering)
    outgoing = jnp.asarray(outgoing_scattering)
    direction = jnp.asarray(directions)
    angular_weight = jnp.asarray(weights)
    energy = jnp.asarray(scattering_weighted_energy_ev)
    if incoming.ndim != 5 or outgoing.shape != incoming.shape:
        raise ValueError("scattering event arrays must have shape [group,direction,x,y,z]")
    if direction.shape != (incoming.shape[1], 3) or angular_weight.shape != (incoming.shape[1],):
        raise ValueError("directions and weights must match the directional scattering axis")
    if energy.shape != (incoming.shape[0],) or not np.isfinite(np.asarray(energy)).all() or np.any(np.asarray(energy) <= 0.0):
        raise ValueError("scattering-weighted photon energy must have shape (group,) and be positive")
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be positive and finite")
    energy_shape = (-1, 1, 1, 1, 1)
    momentum_energy = energy.reshape(energy_shape) * EV_ERG / (LIGHT_SPEED_CM_S * dt)
    return jnp.einsum(
        "d,gdxyz,di->igxyz",
        angular_weight,
        (incoming - outgoing) * momentum_energy,
        direction,
    ).sum(axis=1)


def zero_dust(number_of_groups: int, shape: tuple[int, ...], dtype: jnp.dtype = jnp.float32) -> DustModel:
    """Return a no-dust model compatible with a static transport shape."""
    return DustModel(
        absorption_cross_section_per_h=jnp.zeros((number_of_groups,), dtype=dtype),
        relative_abundance=jnp.zeros(shape, dtype=dtype),
        scattering_cross_section_per_h=jnp.zeros((number_of_groups,), dtype=dtype),
        scattering_weighted_energy_ev=jnp.ones((number_of_groups,), dtype=dtype),
        scattering_phase_function="off",
    )
