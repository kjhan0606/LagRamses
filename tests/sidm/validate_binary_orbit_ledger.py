#!/usr/bin/env python3
"""Validate the dump-free, synchronized two-sink orbit ledger."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


STATE_HEADER = (
    "record_type,step,time_code,dt_since_previous_checkpoint,"
    "dt_completed_coarse_step,aexp,requested_stop_time,level,phase,sink_id,mass,"
    "x,y,z,vx,vy,vz,overshoot"
).split(",")
TERMINAL_HEADER = (
    "terminal_state,reason,step_at_record,time_code,requested_stop_time,"
    "overshoot,phase"
).split(",")
_FLOAT_FIELDS = (
    "time_code", "dt_since_previous_checkpoint", "dt_completed_coarse_step",
    "aexp", "requested_stop_time", "mass", "x", "y", "z", "vx", "vy",
    "vz", "overshoot",
)
_COMMON_FIELDS = (
    "time_code", "dt_since_previous_checkpoint", "dt_completed_coarse_step",
    "aexp", "requested_stop_time", "level", "phase", "overshoot",
)


class BinaryOrbitLedgerError(ValueError):
    """Raised when the binary-orbit ledger violates its output contract."""


@dataclass(frozen=True)
class BinaryOrbitLedgerSummary:
    checkpoints: int
    terminal_step: int
    terminal_time: float
    requested_stop_time: float
    overshoot: float
    sink_ids: tuple[int, int]


def _read_csv(path: Path, expected_header: list[str]) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="ascii", newline="") as stream:
            reader = csv.reader(stream)
            header = next((row for row in reader if row), None)
            if header != expected_header:
                raise BinaryOrbitLedgerError(f"{path}: unexpected CSV header")
            rows = []
            for row in reader:
                if not row or (row[0].lstrip().startswith("#")):
                    continue
                if len(row) != len(expected_header):
                    raise BinaryOrbitLedgerError(f"{path}: malformed CSV row width")
                rows.append(dict(zip(expected_header, row)))
            return rows
    except OSError as exc:
        raise BinaryOrbitLedgerError(f"cannot read {path}: {exc}") from exc
    except csv.Error as exc:
        raise BinaryOrbitLedgerError(f"{path}: malformed CSV: {exc}") from exc


def _number(row: dict[str, str], field: str, context: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise BinaryOrbitLedgerError(f"{context}: invalid {field}") from exc
    if not math.isfinite(value):
        raise BinaryOrbitLedgerError(f"{context}: non-finite {field}")
    return value


def _integer(row: dict[str, str], field: str, context: str) -> int:
    try:
        value = int(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise BinaryOrbitLedgerError(f"{context}: invalid {field}") from exc
    return value


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=5e-12, abs_tol=5e-14)


def _relative_close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=5e-12, abs_tol=0.0)


def _read_config(path: Path) -> dict[str, float | int]:
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except OSError as exc:
        raise BinaryOrbitLedgerError(f"cannot read {path}: {exc}") from exc
    config_lines = [line.split() for line in lines if line.startswith("# CONFIG")]
    if len(config_lines) != 1 or len(config_lines[0]) != 17:
        raise BinaryOrbitLedgerError(f"{path}: expected one complete CONFIG record")
    (_, _, stop, interval, cfl, tout_min, tout_max, aout_min, aout_max,
     boxlen, boundary_margin, level_min, level_max, step_cap, foutput,
     fbackup, noutput) = config_lines[0]
    try:
        config: dict[str, float | int] = {
            "stop": float(stop),
            "interval": float(interval),
            "cfl": float(cfl),
            "tout_min": float(tout_min),
            "tout_max": float(tout_max),
            "aout_min": float(aout_min),
            "aout_max": float(aout_max),
            "boxlen": float(boxlen),
            "boundary_margin": float(boundary_margin),
            "level_min": int(level_min),
            "level_max": int(level_max),
            "step_cap": int(step_cap),
            "foutput": int(foutput),
            "fbackup": int(fbackup),
            "noutput": int(noutput),
        }
    except ValueError as exc:
        raise BinaryOrbitLedgerError(f"{path}: malformed CONFIG values") from exc
    if not all(math.isfinite(float(config[key])) for key in (
        "stop", "interval", "cfl", "tout_min", "tout_max", "aout_min", "aout_max",
        "boxlen", "boundary_margin"
    )):
        raise BinaryOrbitLedgerError(f"{path}: non-finite CONFIG value")
    if config["stop"] <= 0 or config["interval"] <= 0 or not 0 < config["cfl"] <= 1:
        raise BinaryOrbitLedgerError(f"{path}: invalid stop, cadence, or collisionless CFL")
    if config["step_cap"] <= 0 or config["foutput"] <= config["step_cap"] or \
            config["fbackup"] <= config["step_cap"]:
        raise BinaryOrbitLedgerError(f"{path}: output/backup clocks are not beyond the safety cap")
    if config["noutput"] < 1 or config["tout_min"] <= config["stop"]:
        raise BinaryOrbitLedgerError(f"{path}: scheduled output is not safely beyond the stop")
    if config["tout_max"] < config["tout_min"] or config["aout_max"] < config["aout_min"]:
        raise BinaryOrbitLedgerError(f"{path}: invalid scheduled output range")
    if config["level_min"] <= 0 or config["level_max"] < config["level_min"]:
        raise BinaryOrbitLedgerError(f"{path}: invalid level range")
    if config["boxlen"] <= 0 or not 0 < config["boundary_margin"] < config["boxlen"] / 2:
        raise BinaryOrbitLedgerError(f"{path}: invalid box or force-stencil margin")
    return config


def validate_binary_orbit_ledger(
    state_path: str | Path,
    terminal_path: str | Path | None = None,
    *,
    expected_stop_time: float | None = None,
) -> BinaryOrbitLedgerSummary:
    """Validate a successful run's state ledger and terminal commit marker."""
    state_path = Path(state_path)
    terminal_path = Path(terminal_path) if terminal_path else state_path.with_name(
        state_path.stem + ".terminal.csv"
    )
    config = _read_config(state_path)
    if expected_stop_time is not None and not _close(float(config["stop"]), expected_stop_time):
        raise BinaryOrbitLedgerError("CONFIG stop time differs from requested stop")

    rows = _read_csv(state_path, STATE_HEADER)
    if not rows:
        raise BinaryOrbitLedgerError("state ledger contains no synchronized sink rows")
    if float(config["aout_min"]) <= _number(rows[0], "aexp", "initial row"):
        raise BinaryOrbitLedgerError("scheduled scale-factor output is not beyond initialization")

    checkpoints: list[tuple[str, int, float, list[dict[str, str]]]] = []
    for row in rows:
        context = f"{row.get('record_type', '?')} row"
        if row["record_type"] not in {"initial", "sample", "terminal"}:
            raise BinaryOrbitLedgerError(f"{context}: unknown record type")
        for field in _FLOAT_FIELDS:
            _number(row, field, context)
        step = _integer(row, "step", context)
        sink_id = _integer(row, "sink_id", context)
        level = _integer(row, "level", context)
        if step < 0 or sink_id <= 0 or level != config["level_max"]:
            raise BinaryOrbitLedgerError(f"{context}: invalid step, sink ID, or level")
        if step > config["step_cap"]:
            raise BinaryOrbitLedgerError(f"{context}: safety step cap was exceeded")
        if row["phase"] != "centered_sink_level_sync":
            raise BinaryOrbitLedgerError(f"{context}: checkpoint is not synchronized")
        if _number(row, "aexp", context) >= float(config["aout_min"]):
            raise BinaryOrbitLedgerError(f"{context}: scheduled scale-factor output was reached")
        if _number(row, "mass", context) <= 0:
            raise BinaryOrbitLedgerError(f"{context}: non-positive sink mass")
        for field in ("x", "y", "z"):
            coordinate = _number(row, field, context)
            if not float(config["boundary_margin"]) <= coordinate < (
                float(config["boxlen"]) - float(config["boundary_margin"])
            ):
                raise BinaryOrbitLedgerError(f"{context}: sink lies outside the interior margin")
        if not _close(_number(row, "requested_stop_time", context), float(config["stop"])):
            raise BinaryOrbitLedgerError(f"{context}: requested stop time changed")

    index = 0
    while index < len(rows):
        first = rows[index]
        kind = first["record_type"]
        step = int(first["step"])
        time = float(first["time_code"])
        pair = rows[index:index + 2]
        if len(pair) != 2 or any(
            row["record_type"] != kind or int(row["step"]) != step or
            not _close(float(row["time_code"]), time)
            for row in pair
        ):
            raise BinaryOrbitLedgerError("each checkpoint must contain exactly two sink rows")
        if int(pair[0]["sink_id"]) == int(pair[1]["sink_id"]):
            raise BinaryOrbitLedgerError("duplicate sink IDs in checkpoint")
        ordered = sorted(pair, key=lambda row: int(row["sink_id"]))
        for field in _COMMON_FIELDS:
            left = ordered[0][field]
            right = ordered[1][field]
            if field in {"level", "phase"}:
                equal = left == right
            else:
                equal = _close(float(left), float(right))
            if not equal:
                raise BinaryOrbitLedgerError(f"sink rows disagree on checkpoint field {field}")
        checkpoints.append((kind, step, time, ordered))
        index += 2

    if checkpoints[0][0] != "initial":
        raise BinaryOrbitLedgerError("first checkpoint must be initial")
    if any(kind == "initial" for kind, *_ in checkpoints[1:]):
        raise BinaryOrbitLedgerError("initial checkpoint may appear only once")
    if checkpoints[-1][0] != "terminal" or sum(kind == "terminal" for kind, *_ in checkpoints) != 1:
        raise BinaryOrbitLedgerError("last checkpoint must be the unique terminal record")

    sink_ids = tuple(int(row["sink_id"]) for row in checkpoints[0][3])
    masses = {int(row["sink_id"]): float(row["mass"]) for row in checkpoints[0][3]}
    interval = float(config["interval"])
    next_sample_time = checkpoints[0][2] + interval
    previous_step = -1
    previous_time = -math.inf
    for checkpoint_index, (kind, step, time, pair) in enumerate(checkpoints):
        if step <= previous_step or time <= previous_time:
            raise BinaryOrbitLedgerError("checkpoint steps/times must increase strictly")
        ids = tuple(int(row["sink_id"]) for row in pair)
        if ids != sink_ids:
            raise BinaryOrbitLedgerError("sink IDs changed across checkpoints")
        for row in pair:
            sink_id = int(row["sink_id"])
            if not _relative_close(float(row["mass"]), masses[sink_id]):
                raise BinaryOrbitLedgerError("sink masses changed across checkpoints")
        since_previous = float(pair[0]["dt_since_previous_checkpoint"])
        completed_dt = float(pair[0]["dt_completed_coarse_step"])
        if checkpoint_index == 0:
            if not _close(since_previous, 0.0) or not _close(completed_dt, 0.0):
                raise BinaryOrbitLedgerError("initial checkpoint timestep fields must be zero")
        else:
            if since_previous <= 0 or not _close(since_previous, time - previous_time):
                raise BinaryOrbitLedgerError("checkpoint interval disagrees with recorded times")
            if completed_dt <= 0:
                raise BinaryOrbitLedgerError("checkpoint lacks a positive completed coarse timestep")
            if completed_dt > interval and not _close(completed_dt, interval):
                raise BinaryOrbitLedgerError("completed timestep exceeds sample interval")
            if since_previous > 2 * interval and not _close(since_previous, 2 * interval):
                raise BinaryOrbitLedgerError("checkpoint gap exceeds the cadence bound")
            if kind == "sample":
                tol = max(5e-14, 5e-12 * max(abs(time), abs(next_sample_time), interval))
                if time + tol < next_sample_time:
                    raise BinaryOrbitLedgerError("sample precedes the origin-anchored cadence")
                if time - next_sample_time > completed_dt + tol:
                    raise BinaryOrbitLedgerError("sample missed its origin-anchored cadence")
                while next_sample_time <= time + tol:
                    next_sample_time += interval
        previous_step, previous_time = step, time

    _, terminal_step, terminal_time, terminal_rows = checkpoints[-1]
    if len(checkpoints) > 1 and checkpoints[-2][2] >= float(config["stop"]):
        raise BinaryOrbitLedgerError("terminal record was not the first checkpoint at/after stop")
    terminal = _read_csv(terminal_path, TERMINAL_HEADER)
    if len(terminal) != 1:
        raise BinaryOrbitLedgerError("terminal file must contain exactly one commit record")
    marker = terminal[0]
    if marker["terminal_state"] != "complete" or marker["reason"] != "physical_time":
        raise BinaryOrbitLedgerError("terminal commit does not confirm physical-time completion")
    if marker["phase"] != "centered_sink_level_sync":
        raise BinaryOrbitLedgerError("terminal commit is not synchronized")
    if _integer(marker, "step_at_record", "terminal") != terminal_step:
        raise BinaryOrbitLedgerError("terminal commit step disagrees with state ledger")
    marker_time = _number(marker, "time_code", "terminal")
    marker_stop = _number(marker, "requested_stop_time", "terminal")
    marker_overshoot = _number(marker, "overshoot", "terminal")
    stop = float(config["stop"])
    overshoot = terminal_time - stop
    if terminal_time < stop or not _close(marker_time, terminal_time):
        raise BinaryOrbitLedgerError("terminal record is before stop or disagrees with commit time")
    if not _close(marker_stop, stop) or not _close(marker_overshoot, overshoot):
        raise BinaryOrbitLedgerError("terminal commit overshoot/stop disagrees with state ledger")
    if not _close(float(terminal_rows[0]["overshoot"]), overshoot):
        raise BinaryOrbitLedgerError("terminal state rows report the wrong overshoot")
    terminal_dt = float(terminal_rows[0]["dt_completed_coarse_step"])
    overshoot_tol = max(5e-14, 5e-12 * max(abs(terminal_time), abs(stop), interval))
    if overshoot > terminal_dt + overshoot_tol or overshoot > interval + overshoot_tol:
        raise BinaryOrbitLedgerError("terminal overshoot exceeds the completed-step bound")
    return BinaryOrbitLedgerSummary(
        checkpoints=len(checkpoints),
        terminal_step=terminal_step,
        terminal_time=terminal_time,
        requested_stop_time=stop,
        overshoot=overshoot,
        sink_ids=sink_ids,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("state_csv", type=Path)
    parser.add_argument("--terminal-csv", type=Path)
    parser.add_argument("--expected-stop-time", type=float)
    args = parser.parse_args()
    try:
        result = validate_binary_orbit_ledger(
            args.state_csv,
            args.terminal_csv,
            expected_stop_time=args.expected_stop_time,
        )
    except BinaryOrbitLedgerError as exc:
        parser.error(str(exc))
    print(
        f"VALID checkpoints={result.checkpoints} terminal_step={result.terminal_step} "
        f"terminal_time={result.terminal_time:.16g} overshoot={result.overshoot:.6g} "
        f"sink_ids={result.sink_ids}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
