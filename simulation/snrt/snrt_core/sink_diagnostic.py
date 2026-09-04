"""Reader for the source-audited RAMSES ``sink_*.dat`` diagnostic."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
import struct

import numpy as np


M_SUN_G = 1.98847e33
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
SPEED_OF_LIGHT_CM_S = 2.99792458e10
JULIAN_YEAR_DAYS = 365.25
AGN_EFF_STATUS_SPIN_DISABLED_DEFAULT = 1


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
    # ``raw_radiative_efficiency`` is the sink-array value.  The resolved
    # ``radiative_efficiency`` is the helper-selected base coefficient, and
    # ``effective_radiative_efficiency`` includes the active feedback-mode
    # reduction used for bolometric luminosity.  Keeping all three prevents a
    # silent convention change at the ledger boundary.
    raw_radiative_efficiency: np.ndarray
    radiative_efficiency: np.ndarray
    effective_radiative_efficiency: np.ndarray
    efficiency_status: np.ndarray
    efficiency_contract_ok: np.ndarray
    bolometric_luminosity_erg_s: np.ndarray
    expansion_factor: float
    time_code: float
    nstep_coarse: int


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant is forbidden: {value}")


def _finite_number(record: dict[str, object], name: str) -> float:
    value = record.get(name)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"AGN coarse-state field {name!r} is missing or not numeric")
    value_float = float(value)
    if not math.isfinite(value_float):
        raise ValueError(f"AGN coarse-state field {name!r} must be finite")
    return value_float


def _nullable_finite_number(record: dict[str, object], name: str) -> float | None:
    """Read a numeric diagnostic field, allowing writer-emitted JSON null."""

    if name not in record:
        raise ValueError(f"AGN coarse-state field {name!r} is missing")
    value = record[name]
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"AGN coarse-state field {name!r} is not numeric or null")
    value_float = float(value)
    if not math.isfinite(value_float):
        raise ValueError(f"AGN coarse-state field {name!r} must be finite or null")
    return value_float


def _integer_field(record: dict[str, object], name: str, *, minimum: int | None = None) -> int:
    value = record.get(name)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"AGN coarse-state field {name!r} is missing or not integral")
    value_float = float(value)
    if not math.isfinite(value_float) or not value_float.is_integer():
        raise ValueError(f"AGN coarse-state field {name!r} must be a finite integer")
    integer = int(value_float)
    if minimum is not None and integer < minimum:
        raise ValueError(f"AGN coarse-state field {name!r} is below {minimum}")
    return integer


def _boolean_field(record: dict[str, object], name: str, *, default: bool) -> bool:
    value = record.get(name, default)
    if not isinstance(value, bool):
        raise ValueError(f"AGN coarse-state field {name!r} must be boolean")
    return value


def _time_field(record: dict[str, object]) -> float:
    has_t_code = "t_code" in record and record["t_code"] is not None
    has_time_code = "time_code" in record and record["time_code"] is not None
    if not has_t_code and not has_time_code:
        raise ValueError("AGN coarse-state requires t_code (or legacy time_code)")
    t_code = _finite_number(record, "t_code") if has_t_code else None
    time_code = _finite_number(record, "time_code") if has_time_code else None
    if t_code is not None and time_code is not None and not math.isclose(
        t_code, time_code, rel_tol=0.0, abs_tol=1.0e-13
    ):
        raise ValueError("AGN coarse-state t_code and time_code disagree")
    return t_code if t_code is not None else time_code  # type: ignore[return-value]


def _vector_field(record: dict[str, object], name: str) -> list[float]:
    value = record.get(name)
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"AGN coarse-state field {name!r} must be a length-3 array")
    vector = []
    for index, component in enumerate(value):
        if isinstance(component, bool) or not isinstance(component, (int, float)):
            raise ValueError(f"AGN coarse-state field {name!r}[{index}] is not numeric")
        component_float = float(component)
        if not math.isfinite(component_float):
            raise ValueError(f"AGN coarse-state field {name!r}[{index}] must be finite")
        vector.append(component_float)
    return vector


def _require_string(record: dict[str, object], name: str, expected: str) -> None:
    if record.get(name) != expected:
        raise ValueError(f"AGN coarse-state field {name!r} must equal {expected!r}")


def _validate_agn_coarse_record(record: dict[str, object], line_number: int) -> None:
    """Validate the source-owned algebra and the pre-reset boundary contract."""

    if record.get("record_type") != "agn_coarse_state":
        raise ValueError(f"line {line_number}: not an AGN coarse-state record")
    _integer_field(record, "nstep_coarse", minimum=0)
    _integer_field(record, "sink_id", minimum=1)
    aexp = _finite_number(record, "aexp")
    if not 0.0 < aexp <= 1.0:
        raise ValueError(f"line {line_number}: aexp must lie in (0, 1]")
    _time_field(record)
    year_days = _finite_number(record, "julian_year_days")
    if not math.isclose(year_days, JULIAN_YEAR_DAYS, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"line {line_number}: unsupported year convention {year_days}")
    _require_string(record, "ledger_phase", "pre_feedback_pre_reset")
    _require_string(record, "source_interval_kind", "instantaneous_pre_reset_state")

    mass_msun = _finite_number(record, "mass_msun")
    if mass_msun <= 0.0:
        raise ValueError(f"line {line_number}: mass_msun must be positive")
    _vector_field(record, "position_code")
    _vector_field(record, "velocity_code")
    unit_mass = _finite_number(record, "unit_mass_cgs")
    unit_time = _finite_number(record, "unit_time_cgs")
    if unit_mass <= 0.0 or unit_time <= 0.0:
        raise ValueError(f"line {line_number}: source units must be positive")

    bondi = _finite_number(record, "bondi_rate_code")
    eddington = _finite_number(record, "eddington_rate_code")
    inflow = _finite_number(record, "inflow_rate_code")
    inflow_msun_per_year = _finite_number(record, "inflow_rate_msun_per_yr")
    if bondi < 0.0 or eddington < 0.0 or inflow < 0.0 or inflow_msun_per_year < 0.0:
        raise ValueError(f"line {line_number}: accretion rates must be non-negative")
    expected_inflow = min(bondi, eddington)
    if not math.isclose(inflow, expected_inflow, rel_tol=2.0e-12, abs_tol=1.0e-300):
        raise ValueError(f"line {line_number}: inflow_rate_code is not min(Bondi,Eddington)")
    expected_rate = expected_inflow * unit_mass / M_SUN_G * SECONDS_PER_YEAR / unit_time
    if not math.isclose(
        inflow_msun_per_year,
        expected_rate,
        rel_tol=2.0e-12,
        abs_tol=max(1.0e-300, abs(expected_rate) * 2.0e-12),
    ):
        raise ValueError(f"line {line_number}: inflow rate unit conversion failed")

    raw_efficiency = _nullable_finite_number(record, "raw_radiative_efficiency")
    resolved_efficiency = _nullable_finite_number(record, "radiative_efficiency")
    effective_efficiency = _nullable_finite_number(record, "effective_radiative_efficiency")
    if "efficiency_status" not in record or "efficiency_contract_ok" not in record:
        raise ValueError(
            f"line {line_number}: AGN coarse-state requires efficiency_status "
            "and efficiency_contract_ok"
        )
    efficiency_status = _integer_field(record, "efficiency_status", minimum=0)
    efficiency_contract_ok = _boolean_field(record, "efficiency_contract_ok", default=True)
    if not efficiency_contract_ok and efficiency_status == 0:
        raise ValueError(
            f"line {line_number}: false efficiency contract requires a nonzero status"
        )
    if efficiency_contract_ok:
        # The helper's raw input and resolved base are strict (0,1).  A
        # spin-disabled RAMSES branch can legitimately leave the diagnostic
        # raw sink-array value at zero while selecting the explicit .1 default;
        # accept that only with the corresponding status bit and a promotable
        # contract.  A spin-enabled zero is instead a readable,
        # non-promotable initialization divergence (handled below).
        raw_spin_disabled_default = (
            raw_efficiency == 0.0
            and (efficiency_status & AGN_EFF_STATUS_SPIN_DISABLED_DEFAULT) != 0
        )
        if (
            raw_efficiency is None
            or resolved_efficiency is None
            or effective_efficiency is None
            or not (0.0 < raw_efficiency < 1.0 or raw_spin_disabled_default)
            or not 0.0 < resolved_efficiency < 1.0
            or not 0.0 <= effective_efficiency < 1.0
        ):
            raise ValueError(f"line {line_number}: invalid raw/effective efficiency")
    luminosity = _finite_number(record, "bolometric_luminosity_erg_s")
    if luminosity < 0.0:
        raise ValueError(f"line {line_number}: bolometric luminosity must be non-negative")
    if effective_efficiency is not None:
        expected_luminosity = effective_efficiency * expected_inflow * unit_mass / unit_time * SPEED_OF_LIGHT_CM_S**2
        if not math.isclose(
            luminosity,
            expected_luminosity,
            rel_tol=5.0e-12,
            abs_tol=max(1.0e-20, abs(expected_luminosity) * 5.0e-12),
        ):
            raise ValueError(f"line {line_number}: bolometric luminosity algebra failed")

    if "mass_code" in record and record["mass_code"] is not None:
        mass_code = _finite_number(record, "mass_code")
        if mass_code <= 0.0 or not math.isclose(
            mass_code * unit_mass / M_SUN_G,
            mass_msun,
            rel_tol=5.0e-12,
            abs_tol=max(1.0e-20, abs(mass_msun) * 5.0e-12),
        ):
            raise ValueError(f"line {line_number}: mass_code/mass_msun conversion failed")


def _canonicalize_agn_records(path: str | Path) -> tuple[list[dict[str, object]], int]:
    """Parse, validate, and canonicalize AGN records by stable source key.

    An identical semantic duplicate is a harmless restart re-read and is
    counted once.  A same-key payload conflict is fail-closed: without a run
    identity or dump counter this is also how a rewind/replay ambiguity is
    surfaced rather than hidden.
    """

    records_by_key: dict[tuple[int, int], dict[str, object]] = {}
    duplicate_count = 0
    for line_number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        try:
            record = json.loads(line, parse_constant=_reject_json_constant)
        except (TypeError, ValueError, json.JSONDecodeError) as error:
            raise ValueError(f"line {line_number}: invalid JSONL record: {error}") from error
        if not isinstance(record, dict):
            raise ValueError(f"line {line_number}: JSONL record must be an object")
        record_type = record.get("record_type")
        if record_type == "agn_coarse_state_header":
            continue
        if record_type != "agn_coarse_state":
            raise ValueError(f"line {line_number}: unknown record_type {record_type!r}")
        _validate_agn_coarse_record(record, line_number)
        key = (
            _integer_field(record, "nstep_coarse", minimum=0),
            _integer_field(record, "sink_id", minimum=1),
        )
        previous = records_by_key.get(key)
        if previous is None:
            records_by_key[key] = record
        elif previous == record:
            duplicate_count += 1
        else:
            raise ValueError(
                "conflicting AGN coarse-state duplicate for "
                f"(nstep_coarse={key[0]}, sink_id={key[1]})"
            )

    records = [records_by_key[key] for key in sorted(records_by_key)]
    if not records:
        raise ValueError("AGN coarse-state ledger contains no active records")
    step_context: dict[int, tuple[float, float, str, str]] = {}
    for record in records:
        step = _integer_field(record, "nstep_coarse", minimum=0)
        context = (
            _finite_number(record, "aexp"),
            _time_field(record),
            str(record["ledger_phase"]),
            str(record["source_interval_kind"]),
        )
        previous = step_context.setdefault(step, context)
        if previous != context:
            raise ValueError(f"AGN coarse-state step {step} has inconsistent epoch metadata")
    return records, duplicate_count


def read_agn_coarse_records(path: str | Path) -> list[dict[str, object]]:
    """Return validated, stable-key canonical AGN records."""

    return _canonicalize_agn_records(path)[0]


def read_agn_coarse_state(
    path: str | Path,
    *,
    expansion_factor: float,
    expansion_factor_tolerance: float = 1.0e-10,
) -> AgnCoarseState:
    """Read one unambiguous coarse step from the active AGN JSONL diagnostic.

    The diagnostic is emitted before AGN feedback resets its mass accumulators.
    Selection is by the snapshot expansion factor, never by a reset-prone
    accreted-mass field.  Identical restart duplicates are collapsed before
    selection; conflicting same-key payloads fail closed.
    """

    if not 0.0 < expansion_factor <= 1.0 or expansion_factor_tolerance <= 0.0:
        raise ValueError("invalid AGN coarse-state expansion-factor selection")
    records = read_agn_coarse_records(path)
    matches = [
        record
        for record in records
        if abs(_finite_number(record, "aexp") - expansion_factor) <= expansion_factor_tolerance
    ]
    if not matches:
        raise ValueError("no AGN coarse-state records match the requested expansion factor")
    steps = {_integer_field(record, "nstep_coarse", minimum=0) for record in matches}
    if len(steps) != 1:
        raise ValueError("requested expansion factor matches multiple AGN coarse steps")
    non_promotable = [
        _integer_field(record, "sink_id", minimum=1)
        for record in matches
        if not bool(record["efficiency_contract_ok"])
    ]
    if non_promotable:
        raise ValueError(
            "requested AGN coarse step contains non-promotable efficiency "
            f"contract for sink IDs {non_promotable}; refusing state promotion"
        )
    matches.sort(key=lambda record: _integer_field(record, "sink_id", minimum=1))
    sink_id = np.asarray([_integer_field(record, "sink_id", minimum=1) for record in matches], dtype=np.int64)
    position = np.asarray([_vector_field(record, "position_code") for record in matches], dtype=np.float64)
    velocity = np.asarray([_vector_field(record, "velocity_code") for record in matches], dtype=np.float64)
    mass = np.asarray([_finite_number(record, "mass_msun") for record in matches], dtype=np.float64)
    rate_scale = np.asarray(
        [
            _finite_number(record, "unit_mass_cgs")
            / M_SUN_G
            * SECONDS_PER_YEAR
            / _finite_number(record, "unit_time_cgs")
            for record in matches
        ],
        dtype=np.float64,
    )
    bondi = np.asarray([_finite_number(record, "bondi_rate_code") for record in matches], dtype=np.float64) * rate_scale
    eddington = np.asarray([_finite_number(record, "eddington_rate_code") for record in matches], dtype=np.float64) * rate_scale
    inflow = np.asarray([_finite_number(record, "inflow_rate_msun_per_yr") for record in matches], dtype=np.float64)
    raw_efficiency = np.asarray(
        [_finite_number(record, "raw_radiative_efficiency") for record in matches], dtype=np.float64
    )
    resolved_efficiency = np.asarray(
        [_finite_number(record, "radiative_efficiency") for record in matches], dtype=np.float64
    )
    effective_efficiency = np.asarray(
        [_finite_number(record, "effective_radiative_efficiency") for record in matches], dtype=np.float64
    )
    efficiency_status = np.asarray(
        [_integer_field(record, "efficiency_status", minimum=0) for record in matches], dtype=np.int64
    )
    efficiency_contract_ok = np.asarray(
        [_boolean_field(record, "efficiency_contract_ok", default=True) for record in matches], dtype=bool
    )
    luminosity = np.asarray(
        [_finite_number(record, "bolometric_luminosity_erg_s") for record in matches], dtype=np.float64
    )
    arrays = (
        position,
        velocity,
        mass,
        bondi,
        eddington,
        inflow,
        raw_efficiency,
        resolved_efficiency,
        effective_efficiency,
        luminosity,
    )
    if position.shape != (len(sink_id), 3) or velocity.shape != position.shape or any(
        not np.isfinite(value).all() for value in arrays
    ):
        raise ValueError("AGN coarse-state contains invalid source arrays")
    if np.any(mass <= 0.0) or np.any(bondi < 0.0) or np.any(eddington < 0.0) or np.any(inflow < 0.0):
        raise ValueError("AGN coarse-state contains invalid masses or accretion rates")
    # The row-level validator above already admits the explicit
    # spin-disabled raw=0 fallback and rejects spin-enabled initialization
    # divergence.  Keep this vector check consistent with that policy.
    if np.any((raw_efficiency < 0.0) | (raw_efficiency >= 1.0)) or np.any(
        (resolved_efficiency <= 0.0) | (resolved_efficiency >= 1.0)
    ) or np.any(
        (effective_efficiency < 0.0) | (effective_efficiency >= 1.0)
    ) or np.any(luminosity < 0.0):
        raise ValueError("AGN coarse-state contains invalid radiative efficiencies or luminosities")
    return AgnCoarseState(
        sink_id=sink_id,
        mass_msun=mass,
        position_code=position,
        velocity_code=velocity,
        bondi_rate_msun_per_year=bondi,
        eddington_rate_msun_per_year=eddington,
        inflow_rate_msun_per_year=inflow,
        raw_radiative_efficiency=raw_efficiency,
        radiative_efficiency=resolved_efficiency,
        effective_radiative_efficiency=effective_efficiency,
        efficiency_status=efficiency_status,
        efficiency_contract_ok=efficiency_contract_ok,
        bolometric_luminosity_erg_s=luminosity,
        expansion_factor=float(matches[0]["aexp"]),
        time_code=_time_field(matches[0]),
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
