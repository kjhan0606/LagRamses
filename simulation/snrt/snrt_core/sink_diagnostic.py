"""Reader for the source-audited RAMSES ``sink_*.dat`` diagnostic."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import struct

import numpy as np


M_SUN_G = 1.98847e33
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
SPEED_OF_LIGHT_CM_S = 2.99792458e10


@dataclass(frozen=True)
class SinkDiagnostic:
    """Instantaneous accretion state written by ``sinkprops=.true.``.

    ``bondi_rate_code`` and ``eddington_rate_code`` are rates at one coarse
    step, not checkpoint accumulators. The applicable inflow rate is their
    non-negative minimum, matching the Bondi/Eddington cap in the active
    sink routine.
    """

    sink_id: np.ndarray
    mass_code: np.ndarray
    position_code: np.ndarray
    velocity_code: np.ndarray
    bondi_rate_code: np.ndarray
    eddington_rate_code: np.ndarray
    accreted_mass_code: np.ndarray
    radiative_efficiency: np.ndarray
    expansion_factor: float
    time_code: float
    length_scale_cm: float
    density_scale_g_cm3: float
    time_scale_s: float

    def __post_init__(self) -> None:
        sink_id = np.asarray(self.sink_id, dtype=np.int64)
        one_dimensional = {
            "mass_code": self.mass_code,
            "bondi_rate_code": self.bondi_rate_code,
            "eddington_rate_code": self.eddington_rate_code,
            "accreted_mass_code": self.accreted_mass_code,
            "radiative_efficiency": self.radiative_efficiency,
        }
        normalized = {name: np.asarray(value, dtype=np.float64) for name, value in one_dimensional.items()}
        position = np.asarray(self.position_code, dtype=np.float64)
        velocity = np.asarray(self.velocity_code, dtype=np.float64)
        if sink_id.ndim != 1 or any(value.shape != sink_id.shape for value in normalized.values()):
            raise ValueError("sink diagnostic scalar fields must match the sink ID array")
        if position.shape != (len(sink_id), 3) or velocity.shape != position.shape:
            raise ValueError("sink diagnostic position and velocity must have shape (n_sink, 3)")
        if not np.isfinite(position).all() or not np.isfinite(velocity).all():
            raise ValueError("sink diagnostic positions and velocities must be finite")
        if any(not np.isfinite(value).all() for value in normalized.values()):
            raise ValueError("sink diagnostic scalar fields must be finite")
        if np.any(normalized["mass_code"] <= 0.0):
            raise ValueError("sink diagnostic masses must be positive")
        if np.any(normalized["bondi_rate_code"] < 0.0) or np.any(normalized["eddington_rate_code"] < 0.0):
            raise ValueError("sink diagnostic accretion rates must be non-negative")
        if np.any((normalized["radiative_efficiency"] <= 0.0) | (normalized["radiative_efficiency"] >= 1.0)):
            raise ValueError("sink radiative efficiencies must lie in (0, 1)")
        if not 0.0 < self.expansion_factor <= 1.0:
            raise ValueError("diagnostic expansion factor must lie in (0, 1]")
        if any(value <= 0.0 or not np.isfinite(value) for value in (self.length_scale_cm, self.density_scale_g_cm3, self.time_scale_s)):
            raise ValueError("diagnostic unit scales must be positive and finite")
        object.__setattr__(self, "sink_id", sink_id)
        object.__setattr__(self, "position_code", position)
        object.__setattr__(self, "velocity_code", velocity)
        for name, value in normalized.items():
            object.__setattr__(self, name, value)

    @property
    def mass_scale_msun(self) -> float:
        return self.density_scale_g_cm3 * self.length_scale_cm**3 / M_SUN_G

    @property
    def rate_scale_msun_per_year(self) -> float:
        return self.mass_scale_msun * SECONDS_PER_YEAR / self.time_scale_s

    @property
    def bondi_rate_msun_per_year(self) -> np.ndarray:
        return self.bondi_rate_code * self.rate_scale_msun_per_year

    @property
    def eddington_rate_msun_per_year(self) -> np.ndarray:
        return self.eddington_rate_code * self.rate_scale_msun_per_year

    @property
    def inflow_rate_msun_per_year(self) -> np.ndarray:
        return np.minimum(self.bondi_rate_msun_per_year, self.eddington_rate_msun_per_year)

    @property
    def bolometric_luminosity_erg_s(self) -> np.ndarray:
        inflow_rate_g_s = self.inflow_rate_msun_per_year * M_SUN_G / SECONDS_PER_YEAR
        return self.radiative_efficiency * inflow_rate_g_s * SPEED_OF_LIGHT_CM_S**2


@dataclass(frozen=True)
class AgnCoarseState:
    """One coarse-step AGN source state from ``agn_coarse_state_v1.jsonl``."""

    sink_id: np.ndarray
    mass_msun: np.ndarray
    position_code: np.ndarray
    velocity_code: np.ndarray
    bondi_rate_msun_per_year: np.ndarray
    eddington_rate_msun_per_year: np.ndarray
    inflow_rate_msun_per_year: np.ndarray
    radiative_efficiency: np.ndarray
    bolometric_luminosity_erg_s: np.ndarray
    expansion_factor: float
    time_code: float
    nstep_coarse: int


def read_agn_coarse_state(
    path: str | Path,
    *,
    expansion_factor: float,
    expansion_factor_tolerance: float = 1.0e-10,
) -> AgnCoarseState:
    """Read one unambiguous coarse step from the active AGN JSONL diagnostic.

    The diagnostic is emitted before AGN feedback resets its mass accumulators.
    Selection is by the snapshot expansion factor, never by a reset-prone
    accreted-mass field.
    """

    if not 0.0 < expansion_factor <= 1.0 or expansion_factor_tolerance <= 0.0:
        raise ValueError("invalid AGN coarse-state expansion-factor selection")
    matches = []
    for line_number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        if not line.strip():
            continue
        record = json.loads(line)
        if record.get("record_type") != "agn_coarse_state":
            continue
        if abs(float(record["aexp"]) - expansion_factor) <= expansion_factor_tolerance:
            matches.append(record)
    if not matches:
        raise ValueError("no AGN coarse-state records match the requested expansion factor")
    steps = {int(record["nstep_coarse"]) for record in matches}
    if len(steps) != 1:
        raise ValueError("requested expansion factor matches multiple AGN coarse steps")
    matches.sort(key=lambda record: int(record["sink_id"]))
    sink_id = np.asarray([record["sink_id"] for record in matches], dtype=np.int64)
    if len(np.unique(sink_id)) != len(sink_id):
        raise ValueError("AGN coarse-state selection contains duplicate sink IDs")
    position = np.asarray([record["position_code"] for record in matches], dtype=np.float64)
    velocity = np.asarray([record["velocity_code"] for record in matches], dtype=np.float64)
    mass = np.asarray([record["mass_msun"] for record in matches], dtype=np.float64)
    rate_scale = np.asarray(
        [record["unit_mass_cgs"] / M_SUN_G * SECONDS_PER_YEAR / record["unit_time_cgs"] for record in matches],
        dtype=np.float64,
    )
    bondi = np.asarray([record["bondi_rate_code"] for record in matches], dtype=np.float64) * rate_scale
    eddington = np.asarray([record["eddington_rate_code"] for record in matches], dtype=np.float64) * rate_scale
    inflow = np.asarray([record["inflow_rate_msun_per_yr"] for record in matches], dtype=np.float64)
    efficiency = np.asarray([record["effective_radiative_efficiency"] for record in matches], dtype=np.float64)
    luminosity = np.asarray([record["bolometric_luminosity_erg_s"] for record in matches], dtype=np.float64)
    arrays = (position, velocity, mass, bondi, eddington, inflow, efficiency, luminosity)
    if position.shape != (len(sink_id), 3) or velocity.shape != position.shape or any(not np.isfinite(value).all() for value in arrays):
        raise ValueError("AGN coarse-state contains invalid source arrays")
    if np.any(mass <= 0.0) or np.any(bondi < 0.0) or np.any(eddington < 0.0) or np.any(inflow < 0.0):
        raise ValueError("AGN coarse-state contains invalid masses or accretion rates")
    if np.any((efficiency <= 0.0) | (efficiency >= 1.0)) or np.any(luminosity < 0.0):
        raise ValueError("AGN coarse-state contains invalid radiative efficiencies or luminosities")
    return AgnCoarseState(
        sink_id=sink_id,
        mass_msun=mass,
        position_code=position,
        velocity_code=velocity,
        bondi_rate_msun_per_year=bondi,
        eddington_rate_msun_per_year=eddington,
        inflow_rate_msun_per_year=inflow,
        radiative_efficiency=efficiency,
        bolometric_luminosity_erg_s=luminosity,
        expansion_factor=expansion_factor,
        time_code=float(matches[0]["t_code"]),
        nstep_coarse=steps.pop(),
    )


def _read_fortran_records(path: Path) -> list[bytes]:
    records: list[bytes] = []
    with path.open("rb") as handle:
        while marker := handle.read(4):
            if len(marker) != 4:
                raise ValueError("truncated Fortran record marker")
            record_length = struct.unpack("<I", marker)[0]
            payload = handle.read(record_length)
            if len(payload) != record_length:
                raise ValueError("truncated Fortran record payload")
            trailer = handle.read(4)
            if len(trailer) != 4 or struct.unpack("<I", trailer)[0] != record_length:
                raise ValueError("Fortran record marker mismatch")
            records.append(payload)
    return records


def _int_scalar(record: bytes, label: str) -> int:
    if len(record) != 4:
        raise ValueError(f"{label} must be one 32-bit integer record")
    return struct.unpack("<i", record)[0]


def _float_scalar(record: bytes, label: str) -> float:
    if len(record) != 8:
        raise ValueError(f"{label} must be one 64-bit real record")
    return struct.unpack("<d", record)[0]


def _array(record: bytes, dtype: str, count: int, label: str) -> np.ndarray:
    values = np.frombuffer(record, dtype=dtype)
    if values.shape != (count,):
        raise ValueError(f"{label} has {len(values)} values; expected {count}")
    return values.copy()


def read_sink_diagnostic(path: str | Path) -> SinkDiagnostic:
    """Read the exact 30-record ``sink_*.dat`` format of the active writer."""

    records = _read_fortran_records(Path(path))
    if len(records) != 30:
        raise ValueError(f"sink diagnostic has {len(records)} records; expected 30")
    n_sink = _int_scalar(records[0], "nsink")
    ndim = _int_scalar(records[1], "ndim")
    if n_sink <= 0 or ndim != 3:
        raise ValueError("sink diagnostic requires positive nsink and ndim=3")
    scalar = lambda index, label: _float_scalar(records[index], label)
    array = lambda index, label: _array(records[index], "<f8", n_sink, label)
    return SinkDiagnostic(
        sink_id=_array(records[6], "<i4", n_sink, "idsink"),
        mass_code=array(7, "msink"),
        position_code=np.column_stack((array(8, "xsink_1"), array(9, "xsink_2"), array(10, "xsink_3"))),
        velocity_code=np.column_stack((array(11, "vsink_1"), array(12, "vsink_2"), array(13, "vsink_3"))),
        bondi_rate_code=array(17, "dMBHoverdt"),
        eddington_rate_code=array(18, "dMEdoverdt"),
        accreted_mass_code=array(19, "dMsmbh"),
        radiative_efficiency=array(28, "eps_sink"),
        expansion_factor=scalar(2, "aexp"),
        length_scale_cm=scalar(3, "scale_l"),
        density_scale_g_cm3=scalar(4, "scale_d"),
        time_scale_s=scalar(5, "scale_t"),
        time_code=scalar(29, "time"),
    )
