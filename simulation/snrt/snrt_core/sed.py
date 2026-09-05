"""Validated tabulated source spectra and deterministic group moments.

The transport solver consumes photon-number groups, but the source and dust
closures must be derived from the same continuous spectrum.  This module is a
small NumPy-only boundary utility for that conversion.  It deliberately does
not choose an astrophysical SED or an escape/obscuration model.
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


EV_ERG = 1.602176634e-12
SED_SCHEMA = "snrt_source_sed_v1"
ENERGY_COLUMN = "energy_ev"
FRACTION_COLUMN = "energy_fraction_per_ev"
SED_INTERPOLATION_CONVENTION = "piecewise_linear_energy_fraction_per_ev"
SED_QUADRATURE_SCHEME = "composite_trapezoid_log_refined_union_grid_v1"
SED_QUADRATURE_BASE_SUBDIVISIONS = 2048
SED_QUADRATURE_REFINED_SUBDIVISIONS = 4096
SED_QUADRATURE_RELATIVE_TOLERANCE = 5.0e-6
SOURCE_SED_CONTRACT_STATUSES = frozenset(
    ("candidate_explicit_tabulated_sed", "reference_control_parameterized_pilot")
)
SOURCE_SED_CANDIDATE_STATUS = "candidate_explicit_tabulated_sed"


@dataclass(frozen=True)
class PhotonSED:
    """A photon-number SED normalized to a declared luminosity scale.

    ``photon_rate_per_norm_per_ev`` is photons s^-1 eV^-1 per unit of the
    declared normalization.  For the built-in CSV reader the normalization is
    ``L_bol`` in erg s^-1, so multiplying by a luminosity gives photons s^-1.
    ``energy_fraction_per_ev`` is retained to make the bolometric check
    auditable: ``q_E * E_eV * EV_ERG`` is the represented energy fraction per
    eV.
    """

    energy_ev: np.ndarray
    photon_rate_per_norm_per_ev: np.ndarray
    energy_fraction_per_ev: np.ndarray
    input_sha256: str
    identity: str
    normalization: str
    represented_bolometric_fraction: float
    path: str
    interpolation_convention: str = SED_INTERPOLATION_CONVENTION


@dataclass(frozen=True)
class SpectralGroupMoments:
    """Integrals of one validated source spectrum over configured groups."""

    group_photon_rate_per_norm: np.ndarray
    group_energy_fraction_per_norm: np.ndarray
    photon_weighted_mean_energy_ev: np.ndarray
    support_status: tuple[str, ...]
    quadrature_diagnostics: "QuadratureDiagnostics | None" = None


@dataclass(frozen=True)
class QuadratureDiagnostics:
    """Convergence evidence for a deterministic source-closure quadrature."""

    scheme: str
    interpolation_convention: str
    base_subdivisions: int
    refined_subdivisions: int
    relative_tolerance: float
    group_max_relative_error: tuple[float, ...]
    maximum_relative_error: float
    converged: bool

    def as_dict(self) -> dict[str, object]:
        """Return a JSON-serializable provenance record."""

        return {
            "scheme": self.scheme,
            "interpolation_convention": self.interpolation_convention,
            "base_subdivisions": self.base_subdivisions,
            "refined_subdivisions": self.refined_subdivisions,
            "relative_tolerance": self.relative_tolerance,
            "group_max_relative_error": list(self.group_max_relative_error),
            "maximum_relative_error": self.maximum_relative_error,
            "converged": self.converged,
        }


def sha256_file(path: str | Path) -> str:
    """Return the raw-byte SHA-256 of a source file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_identity(
    input_sha256: str,
    *,
    normalization: str,
    energy_column: str,
    spectrum_column: str,
    represented_bolometric_fraction: float,
    support_ev: tuple[float, float],
    interpolation_convention: str = SED_INTERPOLATION_CONVENTION,
    quadrature_scheme: str = SED_QUADRATURE_SCHEME,
    quadrature_base_subdivisions: int = SED_QUADRATURE_BASE_SUBDIVISIONS,
    quadrature_refined_subdivisions: int = SED_QUADRATURE_REFINED_SUBDIVISIONS,
    quadrature_relative_tolerance: float = SED_QUADRATURE_RELATIVE_TOLERANCE,
) -> str:
    """Hash only path-free validated contract fields and input bytes."""

    payload = {
        "schema": SED_SCHEMA,
        "input_sha256": input_sha256,
        "normalization": normalization,
        "energy_column": energy_column,
        "spectrum_column": spectrum_column,
        "represented_bolometric_fraction": float(represented_bolometric_fraction),
        "support_ev": [float(support_ev[0]), float(support_ev[1])],
        "interpolation_convention": interpolation_convention,
        "quadrature_scheme": quadrature_scheme,
        "quadrature_base_subdivisions": int(quadrature_base_subdivisions),
        "quadrature_refined_subdivisions": int(quadrature_refined_subdivisions),
        "quadrature_relative_tolerance": float(quadrature_relative_tolerance),
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_group_edges(group_edges_ev: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
    edges = np.asarray(group_edges_ev, dtype=np.float64)
    if (
        edges.ndim != 1
        or edges.size < 2
        or not np.isfinite(edges).all()
        or np.any(edges <= 0.0)
        or np.any(np.diff(edges) <= 0.0)
    ):
        raise ValueError("group edges must be finite, positive, and strictly increasing")
    return edges


def read_lbol_photon_sed(
    path: str | Path,
    *,
    expected_bolometric_fraction: float | None = None,
) -> PhotonSED:
    """Read an explicit ``L_bol``-normalized energy-fraction CSV.

    The CSV must contain ``energy_ev`` and ``energy_fraction_per_ev``.  The
    latter is a dimensionless fraction per eV; its integral is the fraction of
    bolometric luminosity represented by the table.  No path is included in
    the identity digest, although the metadata retains it for provenance.
    """

    sed_path = Path(path)
    if not sed_path.is_file():
        raise FileNotFoundError(sed_path)
    with sed_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{sed_path}: missing CSV header")
        required = {ENERGY_COLUMN, FRACTION_COLUMN}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{sed_path}: missing SED columns: {sorted(missing)}")
        rows: list[tuple[float, float]] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                values = (float(row[ENERGY_COLUMN]), float(row[FRACTION_COLUMN]))
            except (TypeError, ValueError) as error:
                raise ValueError(f"{sed_path}:{line_number}: non-numeric SED row") from error
            rows.append(values)
    if len(rows) < 2:
        raise ValueError(f"{sed_path}: SED requires at least two samples")
    data = np.asarray(rows, dtype=np.float64)
    energy = data[:, 0]
    energy_fraction = data[:, 1]
    if (
        not np.isfinite(data).all()
        or np.any(energy <= 0.0)
        or np.any(energy_fraction < 0.0)
        or np.any(np.diff(energy) <= 0.0)
    ):
        raise ValueError(
            f"{sed_path}: energy must be strictly increasing and positive; "
            "energy fractions must be finite and non-negative"
        )
    represented_fraction = float(np.trapezoid(energy_fraction, energy))
    if not np.isfinite(represented_fraction) or represented_fraction <= 0.0:
        raise ValueError(f"{sed_path}: represented bolometric fraction must be positive")
    if represented_fraction > 1.0 + 1.0e-10:
        raise ValueError(f"{sed_path}: represented bolometric fraction exceeds unity")
    if expected_bolometric_fraction is not None:
        if not np.isfinite(expected_bolometric_fraction) or expected_bolometric_fraction <= 0.0:
            raise ValueError("expected bolometric fraction must be finite and positive")
        if not np.isclose(represented_fraction, expected_bolometric_fraction, rtol=1.0e-8, atol=1.0e-12):
            raise ValueError(
                f"{sed_path}: represented bolometric fraction {represented_fraction:.17g} "
                f"does not match declared {expected_bolometric_fraction:.17g}"
            )
    photon_rate = energy_fraction / (energy * EV_ERG)
    input_hash = sha256_file(sed_path)
    support = (float(energy[0]), float(energy[-1]))
    identity = _canonical_identity(
        input_hash,
        normalization="L_bol_erg_s",
        energy_column=ENERGY_COLUMN,
        spectrum_column=FRACTION_COLUMN,
        represented_bolometric_fraction=represented_fraction,
        support_ev=support,
        interpolation_convention=SED_INTERPOLATION_CONVENTION,
        quadrature_scheme=SED_QUADRATURE_SCHEME,
        quadrature_base_subdivisions=SED_QUADRATURE_BASE_SUBDIVISIONS,
        quadrature_refined_subdivisions=SED_QUADRATURE_REFINED_SUBDIVISIONS,
        quadrature_relative_tolerance=SED_QUADRATURE_RELATIVE_TOLERANCE,
    )
    return PhotonSED(
        energy_ev=energy.copy(),
        photon_rate_per_norm_per_ev=photon_rate,
        energy_fraction_per_ev=energy_fraction.copy(),
        input_sha256=input_hash,
        identity=identity,
        normalization="L_bol_erg_s",
        represented_bolometric_fraction=represented_fraction,
        path=str(sed_path.resolve()),
        interpolation_convention=SED_INTERPOLATION_CONVENTION,
    )


def source_sed_metadata(sed: PhotonSED) -> dict[str, object]:
    """Return the pathful-but-hash-bound provenance block for a ``PhotonSED``."""

    return {
        "schema": SED_SCHEMA,
        "identity": sed.identity,
        "input_sha256": sed.input_sha256,
        "input_path": sed.path,
        "normalization": sed.normalization,
        "energy_column": ENERGY_COLUMN,
        "spectrum_column": FRACTION_COLUMN,
        "spectrum_units": "dimensionless energy fraction per eV per L_bol",
        "represented_bolometric_fraction": sed.represented_bolometric_fraction,
        "support_ev": [float(sed.energy_ev[0]), float(sed.energy_ev[-1])],
        "interpolation_convention": sed.interpolation_convention,
        "quadrature": {
            "scheme": SED_QUADRATURE_SCHEME,
            "base_subdivisions": SED_QUADRATURE_BASE_SUBDIVISIONS,
            "refined_subdivisions": SED_QUADRATURE_REFINED_SUBDIVISIONS,
            "relative_tolerance": SED_QUADRATURE_RELATIVE_TOLERANCE,
        },
    }


def _validate_quadrature_subdivisions(subdivisions: int) -> int:
    if not isinstance(subdivisions, int) or subdivisions < 1:
        raise ValueError("quadrature subdivisions must be a positive integer")
    return subdivisions


def source_sed_group_grid(
    sed: PhotonSED,
    group_edges_ev: np.ndarray | list[float] | tuple[float, ...],
    *,
    subdivisions: int = SED_QUADRATURE_REFINED_SUBDIVISIONS,
    extra_energy_ev: np.ndarray | list[float] | tuple[float, ...] | None = None,
) -> tuple[np.ndarray, ...]:
    """Build log-refined per-group grids for piecewise-linear SED integration."""

    subdivisions = _validate_quadrature_subdivisions(subdivisions)
    edges = _validate_group_edges(group_edges_ev)
    energy = np.asarray(sed.energy_ev, dtype=np.float64)
    if energy.ndim != 1 or len(energy) < 2 or not np.isfinite(energy).all():
        raise ValueError("PhotonSED energy samples are invalid")
    if edges[0] < energy[0] or edges[-1] > energy[-1]:
        raise ValueError("source SED support must cover every configured group edge")
    if extra_energy_ev is None:
        extra = np.empty(0, dtype=np.float64)
    else:
        extra = np.asarray(extra_energy_ev, dtype=np.float64)
        if extra.ndim != 1 or not np.isfinite(extra).all() or np.any(extra <= 0.0):
            raise ValueError("extra quadrature energies must be finite and positive")
    grids: list[np.ndarray] = []
    for lower, upper in zip(edges[:-1], edges[1:], strict=True):
        refinement = np.geomspace(lower, upper, subdivisions + 1, dtype=np.float64)
        source_nodes = energy[(energy > lower) & (energy < upper)]
        extra_nodes = extra[(extra > lower) & (extra < upper)]
        grids.append(np.unique(np.concatenate((refinement, source_nodes, extra_nodes))))
    return tuple(grids)


def _relative_error_by_group(reference: np.ndarray, refined: np.ndarray) -> np.ndarray:
    scale = np.maximum(np.maximum(np.abs(reference), np.abs(refined)), 1.0e-300)
    return np.abs(refined - reference) / scale


def _integrate_photon_sed_groups_once(
    sed: PhotonSED,
    edges: np.ndarray,
    *,
    subdivisions: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    grids = source_sed_group_grid(sed, edges, subdivisions=subdivisions)
    photon = np.zeros(edges.size - 1, dtype=np.float64)
    energy_fraction = np.zeros(edges.size - 1, dtype=np.float64)
    mean_energy = np.sqrt(edges[:-1] * edges[1:])
    statuses: list[str] = []
    for group, group_energy in enumerate(grids):
        group_fraction = np.interp(
            group_energy, sed.energy_ev, sed.energy_fraction_per_ev
        )
        group_spectrum = group_fraction / (group_energy * EV_ERG)
        photon[group] = float(np.trapezoid(group_spectrum, group_energy))
        energy_fraction[group] = float(np.trapezoid(group_fraction, group_energy))
        if not np.isfinite(photon[group]) or photon[group] < 0.0:
            raise ValueError(f"source SED produced invalid photon integral in group {group}")
        if photon[group] <= 0.0:
            statuses.append("empty_source_group_zero_photons")
            continue
        mean_energy[group] = float(
            np.trapezoid(group_fraction, group_energy) / (EV_ERG * photon[group])
        )
        statuses.append("source_sed_fully_supported")
    if not np.isfinite(energy_fraction).all() or np.any(energy_fraction < 0.0):
        raise ValueError("source SED produced invalid energy integrals")
    return photon, energy_fraction, mean_energy, tuple(statuses)


def integrate_photon_sed_groups(
    sed: PhotonSED,
    group_edges_ev: np.ndarray | list[float] | tuple[float, ...],
    *,
    allow_empty_groups: bool = True,
) -> SpectralGroupMoments:
    """Integrate a tabulated SED with a declared interpolation and convergence guard."""

    edges = _validate_group_edges(group_edges_ev)
    energy = np.asarray(sed.energy_ev, dtype=np.float64)
    fraction = np.asarray(sed.energy_fraction_per_ev, dtype=np.float64)
    if energy.ndim != 1 or fraction.shape != energy.shape or len(energy) < 2:
        raise ValueError("PhotonSED arrays have inconsistent shapes")
    if edges[0] < energy[0] or edges[-1] > energy[-1]:
        raise ValueError(
            "source SED support must cover every configured group edge: "
            f"sed=[{energy[0]:.17g},{energy[-1]:.17g}] "
            f"groups=[{edges[0]:.17g},{edges[-1]:.17g}]"
        )
    if not np.isfinite(fraction).all() or np.any(fraction < 0.0):
        raise ValueError("PhotonSED energy-fraction samples are invalid")
    base = _integrate_photon_sed_groups_once(
        sed, edges, subdivisions=SED_QUADRATURE_BASE_SUBDIVISIONS
    )
    refined = _integrate_photon_sed_groups_once(
        sed, edges, subdivisions=SED_QUADRATURE_REFINED_SUBDIVISIONS
    )
    photon_error = _relative_error_by_group(base[0], refined[0])
    fraction_error = _relative_error_by_group(base[1], refined[1])
    mean_error = _relative_error_by_group(base[2], refined[2])
    group_error = np.maximum.reduce((photon_error, fraction_error, mean_error))
    maximum_error = float(np.max(group_error))
    diagnostics = QuadratureDiagnostics(
        scheme=SED_QUADRATURE_SCHEME,
        interpolation_convention=sed.interpolation_convention,
        base_subdivisions=SED_QUADRATURE_BASE_SUBDIVISIONS,
        refined_subdivisions=SED_QUADRATURE_REFINED_SUBDIVISIONS,
        relative_tolerance=SED_QUADRATURE_RELATIVE_TOLERANCE,
        group_max_relative_error=tuple(float(value) for value in group_error),
        maximum_relative_error=maximum_error,
        converged=maximum_error <= SED_QUADRATURE_RELATIVE_TOLERANCE,
    )
    if not diagnostics.converged:
        raise ValueError(
            "source SED quadrature did not converge: "
            f"maximum relative error {maximum_error:.6g} exceeds "
            f"{SED_QUADRATURE_RELATIVE_TOLERANCE:.6g}"
        )
    photon, group_energy_fraction, mean_energy, statuses = refined
    if not allow_empty_groups and any(status.startswith("empty_") for status in statuses):
        raise ValueError("source SED has no photons in at least one requested group")
    total_fraction = float(np.sum(group_energy_fraction, dtype=np.float64))
    if not np.isclose(total_fraction, sed.represented_bolometric_fraction, rtol=2.0e-6, atol=1.0e-12):
        raise ValueError(
            "group energy integrals do not reproduce the source SED represented fraction: "
            f"groups={total_fraction:.17g} source={sed.represented_bolometric_fraction:.17g}"
        )
    return SpectralGroupMoments(
        group_photon_rate_per_norm=photon,
        group_energy_fraction_per_norm=group_energy_fraction,
        photon_weighted_mean_energy_ev=mean_energy,
        support_status=statuses,
        quadrature_diagnostics=diagnostics,
    )
