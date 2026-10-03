#!/usr/bin/env python3
"""Validate and summarize lagRamses SMBH pre-compaction JSONL ledgers."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 1


@dataclass
class EventBlock:
    begin: dict[str, Any]
    begin_line: int
    rows: list[dict[str, Any]] = field(default_factory=list)
    members: list[dict[str, Any]] = field(default_factory=list)
    pairs: list[dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.rows.append(self.begin)

    @property
    def uid(self) -> str:
        return str(self.begin.get("event_uid", ""))


@dataclass
class BatchBlock:
    begin: dict[str, Any]
    events: list[tuple[dict[str, Any], str, set[int]]] = field(default_factory=list)

    @property
    def uid(self) -> str:
        return str(self.begin.get("batch_uid", ""))


@dataclass
class LedgerReport:
    unique_events: int = 0
    binary_events: int = 0
    multiple_events: int = 0
    duplicate_events: int = 0
    incomplete_events: int = 0
    censored_events: int = 0
    censored_batches: int = 0
    incomplete_batches: int = 0
    invalid_json_lines: int = 0
    committed_batches: int = 0
    run_attempts: int = 0
    superseded_batches: int = 0
    superseded_events: int = 0
    errors: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors and self.invalid_json_lines == 0

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "status": "valid" if self.valid else "invalid",
            "unique_events": self.unique_events,
            "binary_events": self.binary_events,
            "multiple_events": self.multiple_events,
            "duplicate_events": self.duplicate_events,
            "incomplete_events": self.incomplete_events,
            "censored_events": self.censored_events,
            "censored_batches": self.censored_batches,
            "incomplete_batches": self.incomplete_batches,
            "invalid_json_lines": self.invalid_json_lines,
            "committed_batches": self.committed_batches,
            "run_attempts": self.run_attempts,
            "superseded_batches": self.superseded_batches,
            "superseded_events": self.superseded_events,
            "errors": self.errors,
        }


def _close(a: float, b: float, *, rtol: float = 2.0e-12) -> bool:
    return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=1.0e-14)


def _norm(vector: Iterable[float]) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in vector))


def _cross(a: list[float], b: list[float]) -> list[float]:
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def _require_number(
    record: dict[str, Any], key: str, uid: str, errors: list[str]
) -> float | None:
    value = record.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append(f"{uid}: {key} must be a finite number")
        return None
    value = float(value)
    if not math.isfinite(value):
        errors.append(f"{uid}: {key} must be a finite number")
        return None
    return value


def _require_integer(
    record: dict[str, Any], key: str, uid: str, errors: list[str]
) -> int | None:
    value = record.get(key)
    if isinstance(value, bool) or not isinstance(value, int):
        errors.append(f"{uid}: {key} must be an integer")
        return None
    return value


def _require_logical(
    record: dict[str, Any], key: str, uid: str, errors: list[str]
) -> bool | None:
    value = record.get(key)
    if not isinstance(value, bool):
        errors.append(f"{uid}: {key} must be a boolean")
        return None
    return value


def _minimum_image_delta(
    position1: list[float], position2: list[float], box_size: list[float]
) -> list[float]:
    delta = [b - a for a, b in zip(position1, position2)]
    for index, value in enumerate(delta):
        if value > 0.5 * box_size[index]:
            delta[index] -= box_size[index]
        if value < -0.5 * box_size[index]:
            delta[index] += box_size[index]
    return delta


def _event_digest(rows: list[dict[str, Any]]) -> str:
    payload = "\n".join(
        json.dumps(row, sort_keys=True, separators=(",", ":")) for row in rows
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _require_vector(
    record: dict[str, Any], key: str, uid: str, errors: list[str]
) -> list[float] | None:
    value = record.get(key)
    if not isinstance(value, list) or len(value) != 3:
        errors.append(f"{uid}: {key} must be a three-component array")
        return None
    if not all(
        not isinstance(item, bool)
        and isinstance(item, (int, float))
        and math.isfinite(item)
        for item in value
    ):
        errors.append(f"{uid}: {key} contains a non-finite value")
        return None
    return [float(item) for item in value]


def _validate_pair_invariants(
    uid: str,
    pair: dict[str, Any],
    members: dict[int, dict[str, Any]],
    box_size: list[float] | None,
    fact_g: float | None,
    merge_radius: float | None,
    errors: list[str],
) -> None:
    id1 = _require_integer(pair, "sink_id_1", uid, errors)
    id2 = _require_integer(pair, "sink_id_2", uid, errors)
    if id1 is None or id2 is None or id1 not in members or id2 not in members:
        return
    delta_position = _require_vector(pair, "delta_position_code", uid, errors)
    delta_velocity = _require_vector(pair, "delta_velocity_code", uid, errors)
    specific_h = _require_vector(pair, "specific_angular_momentum_code", uid, errors)
    relative_l = _require_vector(pair, "relative_angular_momentum_code", uid, errors)
    if None in (delta_position, delta_velocity, specific_h, relative_l):
        return

    separation = _require_number(pair, "separation_code", uid, errors)
    relative_speed = _require_number(pair, "relative_speed_code", uid, errors)
    reduced_mass = _require_number(pair, "reduced_mass_code", uid, errors)
    relative_kinetic = _require_number(pair, "relative_kinetic_code", uid, errors)
    expected_r = _norm(delta_position)
    expected_v = _norm(delta_velocity)
    if separation is not None and not _close(separation, expected_r):
        errors.append(f"{uid}: pair {id1}-{id2} separation invariant failed")
    if relative_speed is not None and not _close(relative_speed, expected_v):
        errors.append(f"{uid}: pair {id1}-{id2} relative-speed invariant failed")

    mass1 = _require_number(members[id1], "mass_code", uid, errors)
    mass2 = _require_number(members[id2], "mass_code", uid, errors)
    if mass1 is None or mass2 is None or mass1 + mass2 <= 0.0:
        errors.append(f"{uid}: pair {id1}-{id2} has an invalid mass sum")
        return
    expected_mu = mass1 * mass2 / (mass1 + mass2)
    expected_kinetic = 0.5 * expected_mu * expected_v**2
    if reduced_mass is not None and not _close(reduced_mass, expected_mu):
        errors.append(f"{uid}: pair {id1}-{id2} reduced-mass invariant failed")
    if relative_kinetic is not None and not _close(relative_kinetic, expected_kinetic):
        errors.append(f"{uid}: pair {id1}-{id2} kinetic-energy invariant failed")

    position1 = _require_vector(members[id1], "position_code", uid, errors)
    position2 = _require_vector(members[id2], "position_code", uid, errors)
    velocity1 = _require_vector(members[id1], "velocity_code", uid, errors)
    velocity2 = _require_vector(members[id2], "velocity_code", uid, errors)
    if box_size is not None and position1 is not None and position2 is not None:
        expected_delta_position = _minimum_image_delta(position1, position2, box_size)
        if any(
            not _close(got, expected)
            for got, expected in zip(delta_position, expected_delta_position)
        ):
            errors.append(
                f"{uid}: pair {id1}-{id2} minimum-image position invariant failed"
            )
    if velocity1 is not None and velocity2 is not None:
        expected_delta_velocity = [b - a for a, b in zip(velocity1, velocity2)]
        if any(
            not _close(got, expected)
            for got, expected in zip(delta_velocity, expected_delta_velocity)
        ):
            errors.append(
                f"{uid}: pair {id1}-{id2} member-velocity invariant failed"
            )

    expected_h = _cross(delta_position, delta_velocity)
    for got, expected in zip(specific_h, expected_h):
        if not _close(got, expected):
            errors.append(
                f"{uid}: pair {id1}-{id2} specific-angular-momentum invariant failed"
            )
            break
    for got, expected in zip(relative_l, (expected_mu * x for x in expected_h)):
        if not _close(got, expected):
            errors.append(
                f"{uid}: pair {id1}-{id2} relative-angular-momentum invariant failed"
            )
            break

    within_rmerge = _require_logical(pair, "within_rmerge", uid, errors)
    two_body_bound = _require_logical(pair, "two_body_bound", uid, errors)
    legacy_pair_bound = _require_logical(pair, "legacy_pair_bound", uid, errors)
    if merge_radius is not None and within_rmerge is not None:
        expected_within = expected_r <= merge_radius
        if within_rmerge is not expected_within:
            errors.append(f"{uid}: pair {id1}-{id2} within-rmerge invariant failed")

    finite_pair = expected_r > sys.float_info.min
    potential = pair.get("newtonian_potential_1overr_code")
    specific_energy = pair.get("two_body_specific_energy_code")
    legacy_proxy = pair.get("legacy_binding_proxy_1overr2_code")
    if not finite_pair:
        if potential is not None or specific_energy is not None or legacy_proxy is not None:
            errors.append(f"{uid}: pair {id1}-{id2} singular energies must be null")
        if two_body_bound is not False or legacy_pair_bound is not False:
            errors.append(f"{uid}: pair {id1}-{id2} singular binding flags must be false")
        return
    if fact_g is None:
        return

    potential_value = _require_number(
        pair, "newtonian_potential_1overr_code", uid, errors
    )
    specific_energy_value = _require_number(
        pair, "two_body_specific_energy_code", uid, errors
    )
    legacy_proxy_value = _require_number(
        pair, "legacy_binding_proxy_1overr2_code", uid, errors
    )
    expected_potential = -fact_g * mass1 * mass2 / expected_r
    expected_specific_energy = 0.5 * expected_v**2 - fact_g * (mass1 + mass2) / expected_r
    expected_legacy_proxy = fact_g * mass1 * mass2 / expected_r**2
    if potential_value is not None and not _close(potential_value, expected_potential):
        errors.append(f"{uid}: pair {id1}-{id2} potential-energy invariant failed")
    if specific_energy_value is not None and not _close(
        specific_energy_value, expected_specific_energy
    ):
        errors.append(f"{uid}: pair {id1}-{id2} specific-energy invariant failed")
    if legacy_proxy_value is not None and not _close(
        legacy_proxy_value, expected_legacy_proxy
    ):
        errors.append(f"{uid}: pair {id1}-{id2} legacy-binding-proxy invariant failed")
    if two_body_bound is not None and two_body_bound is not (expected_specific_energy < 0.0):
        errors.append(f"{uid}: pair {id1}-{id2} two-body-bound invariant failed")
    if legacy_pair_bound is not None and legacy_pair_bound is not (
        expected_kinetic < expected_legacy_proxy
    ):
        errors.append(f"{uid}: pair {id1}-{id2} legacy-pair-bound invariant failed")


def _validate_event_invariants(
    uid: str,
    begin: dict[str, Any],
    members: list[dict[str, Any]],
    errors: list[str],
) -> tuple[list[float] | None, float | None, float | None]:
    boxlen = _require_number(begin, "boxlen", uid, errors)
    box_size = (
        _require_vector(begin, "periodic_box_size_code", uid, errors)
        if "periodic_box_size_code" in begin
        else None
    )
    fact_g = _require_number(begin, "factG_code", uid, errors)
    merge_radius = _require_number(begin, "merge_radius_code", uid, errors)
    reported_mass = _require_number(begin, "total_mass_code", uid, errors)
    reported_com_position = _require_vector(begin, "com_position_code", uid, errors)
    reported_com_velocity = _require_vector(begin, "com_velocity_code", uid, errors)
    reported_max_separation = _require_number(
        begin, "max_pair_separation_code", uid, errors
    )
    if boxlen is not None and boxlen <= 0.0:
        errors.append(f"{uid}: boxlen must be positive")
        boxlen = None
    if box_size is None and "periodic_box_size_code" not in begin and boxlen is not None:
        box_size = [boxlen] * 3
    if box_size is not None and any(extent <= 0.0 for extent in box_size):
        errors.append(f"{uid}: periodic box extents must be positive")
        box_size = None
    if merge_radius is not None and merge_radius < 0.0:
        errors.append(f"{uid}: merge_radius_code must be non-negative")
        merge_radius = None

    masses: list[float] = []
    positions: list[list[float]] = []
    velocities: list[list[float]] = []
    for member in members:
        mass = _require_number(member, "mass_code", uid, errors)
        position = _require_vector(member, "position_code", uid, errors)
        velocity = _require_vector(member, "velocity_code", uid, errors)
        if mass is None or position is None or velocity is None:
            continue
        masses.append(mass)
        positions.append(position)
        velocities.append(velocity)
    if len(masses) != len(members):
        return box_size, fact_g, merge_radius

    total_mass = sum(masses)
    if total_mass <= 0.0:
        errors.append(f"{uid}: total member mass must be positive")
        return box_size, fact_g, merge_radius
    if reported_mass is not None and not _close(reported_mass, total_mass):
        errors.append(f"{uid}: total-mass invariant failed")

    if box_size is None or not positions:
        return box_size, fact_g, merge_radius
    anchor = positions[0]
    com_position = [0.0, 0.0, 0.0]
    com_velocity = [0.0, 0.0, 0.0]
    for mass, position, velocity in zip(masses, positions, velocities):
        delta = _minimum_image_delta(anchor, position, box_size)
        for index in range(3):
            com_position[index] += mass * delta[index]
            com_velocity[index] += mass * velocity[index]
    for index in range(3):
        com_position[index] = (
            anchor[index] + com_position[index] / total_mass
        ) % box_size[index]
        com_velocity[index] /= total_mass
    if reported_com_position is not None and any(
        not _close(got, expected)
        for got, expected in zip(reported_com_position, com_position)
    ):
        errors.append(f"{uid}: centre-of-mass position invariant failed")
    if reported_com_velocity is not None and any(
        not _close(got, expected)
        for got, expected in zip(reported_com_velocity, com_velocity)
    ):
        errors.append(f"{uid}: centre-of-mass velocity invariant failed")

    max_separation = 0.0
    for position1, position2 in combinations(positions, 2):
        max_separation = max(
            max_separation, _norm(_minimum_image_delta(position1, position2, box_size))
        )
    if reported_max_separation is not None and not _close(
        reported_max_separation, max_separation
    ):
        errors.append(f"{uid}: maximum-separation invariant failed")
    return box_size, fact_g, merge_radius


def _validate_complete_block(
    block: EventBlock, end: dict[str, Any], report: LedgerReport
) -> str | None:
    uid = block.uid
    errors = report.errors
    begin = block.begin
    rows = block.rows + [end]

    if not uid:
        errors.append(f"line {block.begin_line}: event_uid is empty")
        return None
    if begin.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"{uid}: unsupported schema version")
    if any(row.get("schema_version") != SCHEMA_VERSION for row in rows):
        errors.append(f"{uid}: unsupported row schema version")
    if begin.get("complete") is not False or end.get("complete") is not True:
        errors.append(f"{uid}: begin/end completion markers are invalid")
    if end.get("event_uid") != uid:
        errors.append(f"{uid}: event_end UID does not match")

    nmember_value = _require_integer(begin, "nmember", uid, errors)
    nmember = nmember_value if nmember_value is not None else -1
    expected_pairs = nmember * (nmember - 1) // 2
    if nmember < 2:
        errors.append(f"{uid}: an event must contain at least two members")
    if begin.get("expected_pairs") != expected_pairs:
        errors.append(f"{uid}: expected_pairs is inconsistent with nmember")
    if len(block.members) != nmember or end.get("nmember") != nmember:
        errors.append(f"{uid}: member count does not match begin/end records")
    if len(block.pairs) != expected_pairs or end.get("npair") != expected_pairs:
        errors.append(f"{uid}: pair count does not equal nmember choose two")

    member_indices = [row.get("member_index") for row in block.members]
    if member_indices != list(range(1, nmember + 1)):
        errors.append(f"{uid}: member_index sequence is not contiguous")

    member_by_id: dict[int, dict[str, Any]] = {}
    for member in block.members:
        sink_id = _require_integer(member, "sink_id", uid, errors)
        if sink_id is None:
            continue
        if sink_id in member_by_id:
            errors.append(f"{uid}: duplicate sink member ID {sink_id}")
        member_by_id[sink_id] = member
        _require_vector(member, "position_code", uid, errors)
        _require_vector(member, "velocity_code", uid, errors)

    # Schema-v1 ledgers written before the primary-ID extension remain valid.
    # When the extension is present, validate the exact merge_sink survivor
    # rule and every requested (sink_id, primary_sink_id) relation.
    if "primary_sink_id" in begin:
        primary_id = _require_integer(begin, "primary_sink_id", uid, errors)
        expected_primary: int | None = None
        expected_mass = -math.inf
        for member in block.members:
            member_id = member.get("sink_id")
            member_mass = member.get("mass_code")
            if (
                isinstance(member_id, int)
                and not isinstance(member_id, bool)
                and isinstance(member_mass, (int, float))
                and not isinstance(member_mass, bool)
                and math.isfinite(float(member_mass))
                and float(member_mass) > expected_mass
            ):
                expected_primary = member_id
                expected_mass = float(member_mass)
        if primary_id not in member_by_id:
            errors.append(f"{uid}: primary_sink_id does not name a member")
        if primary_id is not None and primary_id != expected_primary:
            errors.append(f"{uid}: primary_sink_id violates the survivor rule")
        for member in block.members:
            member_primary = _require_integer(
                member, "primary_sink_id", uid, errors
            )
            is_primary = _require_logical(member, "is_primary", uid, errors)
            if member_primary != primary_id:
                errors.append(f"{uid}: member primary_sink_id is inconsistent")
            if is_primary is not None and is_primary != (
                member.get("sink_id") == primary_id
            ):
                errors.append(f"{uid}: member is_primary flag is inconsistent")

    box_size, fact_g, merge_radius = _validate_event_invariants(
        uid, begin, block.members, errors
    )

    expected_id_pairs = {
        tuple(sorted(pair)) for pair in combinations(member_by_id.keys(), 2)
    }
    actual_id_pairs: set[tuple[int, int]] = set()
    for pair in block.pairs:
        id1 = _require_integer(pair, "sink_id_1", uid, errors)
        id2 = _require_integer(pair, "sink_id_2", uid, errors)
        if id1 is None or id2 is None:
            continue
        id_pair = tuple(sorted((id1, id2)))
        if id_pair in actual_id_pairs:
            errors.append(f"{uid}: duplicate pair record {id_pair}")
        actual_id_pairs.add(id_pair)
        if id_pair[0] not in member_by_id or id_pair[1] not in member_by_id:
            errors.append(f"{uid}: pair {id_pair} references a non-member sink")
            continue
        _validate_pair_invariants(
            uid, pair, member_by_id, box_size, fact_g, merge_radius, errors
        )
    if actual_id_pairs != expected_id_pairs:
        errors.append(f"{uid}: pair records do not cover every member pair exactly once")

    expected_class = "BINARY" if nmember == 2 else "MULTIPLE"
    if begin.get("classification") != expected_class:
        errors.append(f"{uid}: classification must be {expected_class}")
    return _event_digest(rows)


def validate_ledger(path: Path, *, allow_incomplete_tail: bool = False) -> LedgerReport:
    report = LedgerReport()
    current: EventBlock | None = None
    batch: BatchBlock | None = None
    seen_protocol = False
    legacy: list[tuple[dict[str, Any], str]] = []
    committed: list[tuple[int, BatchBlock, str]] = []
    # Each attempt stores (resume_step, parent_attempt, parent_batch_cutoff).
    attempts: list[tuple[int, int | None, int | None]] = []
    # A checkpoint is tied to the number of committed batches BEFORE it.
    checkpoints: dict[int, tuple[int, int, int]] = {}
    digests: dict[str, str] = {}
    batch_digests: dict[str, str] = {}
    corrupt_restart_tail = False

    def integer(record: dict[str, Any], key: str, minimum: int = 0) -> int | None:
        value = record.get(key)
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            report.errors.append(f"{key} must be an integer >= {minimum}")
            return None
        return value

    def accept(begin: dict[str, Any], digest: str) -> None:
        uid = str(begin.get("event_uid", ""))
        previous = digests.get(uid)
        if previous is None:
            digests[uid] = digest
            report.unique_events += 1
            if begin.get("classification") == "BINARY":
                report.binary_events += 1
            elif begin.get("classification") == "MULTIPLE":
                report.multiple_events += 1
        elif previous == digest:
            report.duplicate_events += 1
        else:
            report.errors.append(f"{uid}: deterministic UID has conflicting event data")

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                continue
            try:
                record = json.loads(raw_line)
            except json.JSONDecodeError as exc:
                if seen_protocol and batch is not None and not corrupt_restart_tail:
                    corrupt_restart_tail = True
                    continue
                report.invalid_json_lines += 1
                report.errors.append(f"line {line_number}: invalid JSON: {exc.msg}")
                if current is not None:
                    report.incomplete_events += 1
                    current = None
                continue

            if not isinstance(record, dict):
                report.errors.append(
                    f"line {line_number}: each ledger record must be a JSON object"
                )
                continue

            record_type = record.get("record_type")
            if corrupt_restart_tail:
                if record_type != "attempt_begin":
                    report.errors.append(
                        f"line {line_number}: corrupt batch tail lacks restart attempt"
                    )
                corrupt_restart_tail = False
            if record_type == "attempt_begin":
                if current is not None:
                    report.incomplete_events += 1
                    if batch is not None:
                        report.censored_events += 1
                    if batch is None:
                        report.errors.append(f"{current.uid}: missing event_end")
                    current = None
                if batch is not None:
                    report.censored_batches += 1
                    batch = None
                if record.get("schema_version") != SCHEMA_VERSION:
                    report.errors.append("unsupported attempt schema")
                restart_output = integer(record, "restart_output")
                resume_step = integer(record, "resume_step")
                if restart_output is None or resume_step is None:
                    continue
                if restart_output == 0:
                    if attempts or committed or legacy:
                        report.errors.append("fresh attempt appended to existing ledger")
                    parent, cutoff = None, None
                else:
                    if legacy:
                        report.errors.append("legacy events cannot prove restart lineage")
                    checkpoint = checkpoints.get(restart_output)
                    if checkpoint is None or checkpoint[1] != resume_step:
                        report.errors.append("restart has no matching ledger checkpoint")
                        parent, cutoff = None, None
                    else:
                        parent, _, cutoff = checkpoint
                attempts.append((resume_step, parent, cutoff))
                seen_protocol = True
            elif record_type == "checkpoint":
                if current is not None or batch is not None or not attempts:
                    report.errors.append("checkpoint outside a completed batch/attempt")
                    continue
                if record.get("schema_version") != SCHEMA_VERSION:
                    report.errors.append("unsupported checkpoint schema")
                output = integer(record, "output_number", 1)
                step = integer(record, "nstep_coarse")
                if output is not None and step is not None:
                    if step < attempts[-1][0]:
                        report.errors.append("checkpoint predates active attempt")
                    checkpoints[output] = (len(attempts) - 1, step, len(committed))
            elif record_type == "batch_begin":
                if current is not None or batch is not None:
                    report.errors.append("missing event_end or batch_commit before next batch")
                    continue
                if not attempts:
                    report.errors.append("batch has no run attempt marker")
                    continue
                if record.get("schema_version") != SCHEMA_VERSION:
                    report.errors.append("unsupported batch schema")
                batch = BatchBlock(record)
                seen_protocol = True
            elif record_type == "event_begin":
                if current is not None:
                    report.incomplete_events += 1
                    report.errors.append(
                        f"{current.uid}: missing event_end before line {line_number}"
                    )
                if seen_protocol and batch is None:
                    report.errors.append("bare event after batch protocol began")
                current = EventBlock(record, line_number)
            elif record_type in {"member", "pair"}:
                if current is None:
                    report.errors.append(
                        f"line {line_number}: {record_type} appears outside an event"
                    )
                    continue
                if record.get("event_uid") != current.uid:
                    report.errors.append(
                        f"line {line_number}: record UID does not match active event"
                    )
                current.rows.append(record)
                if record_type == "member":
                    current.members.append(record)
                else:
                    current.pairs.append(record)
            elif record_type == "event_end":
                if current is None:
                    report.errors.append(
                        f"line {line_number}: event_end appears outside an event"
                    )
                    continue
                if batch is not None:
                    if "periodic_box_size_code" not in current.begin:
                        report.errors.append(f"{current.uid}: batched event lacks periodic box extents")
                    if "primary_sink_id" not in current.begin:
                        report.errors.append(f"{current.uid}: batched event lacks primary sink ID")
                digest = _validate_complete_block(current, record, report)
                if digest is not None:
                    if batch is None:
                        legacy.append((current.begin, digest))
                    else:
                        member_ids = {row.get("sink_id") for row in current.members}
                        if not all(isinstance(item, int) and not isinstance(item, bool)
                                   for item in member_ids):
                            report.errors.append(f"{current.uid}: invalid batched member IDs")
                        else:
                            batch.events.append((current.begin, digest, member_ids))
                current = None
            elif record_type == "batch_commit":
                if batch is None or current is not None:
                    report.errors.append("batch_commit without complete open batch")
                    continue
                begin = batch.begin
                step = integer(begin, "nstep_coarse")
                level = integer(begin, "ilevel", 1)
                before = integer(begin, "nsink_before", 1)
                after = integer(begin, "nsink_after", 1)
                expected = integer(begin, "expected_events", 1)
                committed_after = integer(record, "nsink_after", 1)
                if None not in (step, level, before, after, expected, committed_after):
                    uid = f"{step}-{level}-{before}-{after}"
                    if (begin.get("schema_version") != SCHEMA_VERSION
                        or record.get("schema_version") != SCHEMA_VERSION
                        or batch.uid != uid or record.get("batch_uid") != uid
                        or committed_after != after or before <= after
                        or step < attempts[-1][0] or len(batch.events) != expected):
                        report.errors.append(f"{uid}: invalid batch commit or metadata")
                    reduction = sum(len(ids) - 1 for _, _, ids in batch.events)
                    if reduction != before - after:
                        report.errors.append(f"{uid}: sink-count conservation failed")
                    seen_members: set[int] = set()
                    seen_uids: set[str] = set()
                    for event_begin, _, member_ids in batch.events:
                        event_uid = str(event_begin.get("event_uid", ""))
                        expected_uid = (
                            f"{step}-{level}-{min(member_ids)}-{max(member_ids)}-"
                            f"{len(member_ids)}" if member_ids else ""
                        )
                        if (event_uid != expected_uid
                            or event_begin.get("nstep_coarse") != step
                            or event_begin.get("ilevel") != level
                            or event_uid in seen_uids or seen_members & member_ids):
                            report.errors.append(f"{uid}: inconsistent batched event")
                        seen_members.update(member_ids)
                        seen_uids.add(event_uid)
                    payload = [begin, [(event.get("event_uid"), digest)
                                       for event, digest, _ in batch.events], record]
                    batch_digest = hashlib.sha256(
                        json.dumps(payload, sort_keys=True).encode("utf-8")
                    ).hexdigest()
                    committed.append((len(attempts) - 1, batch, batch_digest))
                batch = None
            else:
                report.errors.append(
                    f"line {line_number}: unknown record_type {record_type!r}"
                )

    if corrupt_restart_tail:
        report.errors.append("corrupt batch tail has no restart attempt")
    if current is not None:
        report.incomplete_events += 1
        if batch is not None:
            report.censored_events += 1
        if not allow_incomplete_tail:
            report.errors.append(f"{current.uid}: incomplete event at end of ledger")
    if batch is not None:
        report.censored_batches += 1
        report.incomplete_batches += 1
        if not allow_incomplete_tail:
            report.errors.append(f"{batch.uid}: incomplete batch at end of ledger")

    active_cutoffs: dict[int, int | None] = {}
    if attempts:
        index: int | None = len(attempts) - 1
        cutoff: int | None = None
        while index is not None:
            if index in active_cutoffs:
                report.errors.append("restart lineage cycle")
                break
            active_cutoffs[index] = cutoff
            _, index, cutoff = attempts[index]

    report.run_attempts = len(attempts)
    for begin, digest in legacy:
        accept(begin, digest)
    for batch_index, (owner, committed_batch, digest) in enumerate(committed):
        if owner not in active_cutoffs:
            report.superseded_batches += 1
            report.superseded_events += len(committed_batch.events)
            continue
        cutoff = active_cutoffs[owner]
        if cutoff is not None and batch_index >= cutoff:
            report.superseded_batches += 1
            report.superseded_events += len(committed_batch.events)
            continue
        report.committed_batches += 1
        previous = batch_digests.get(committed_batch.uid)
        if previous is not None and previous != digest:
            report.errors.append(f"{committed_batch.uid}: conflicting deterministic batch UID")
        batch_digests[committed_batch.uid] = digest
        for begin, event_digest, _ in committed_batch.events:
            accept(begin, event_digest)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument(
        "--allow-incomplete-tail",
        action="store_true",
        help="report but do not fail solely for a final begin block without event_end",
    )
    args = parser.parse_args()
    report = validate_ledger(
        args.ledger, allow_incomplete_tail=args.allow_incomplete_tail
    )
    print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    return 0 if report.valid else 1


if __name__ == "__main__":
    sys.exit(main())
