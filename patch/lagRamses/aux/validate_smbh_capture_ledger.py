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


SCHEMA_VERSION = 2


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
class LedgerReport:
    unique_events: int = 0
    binary_events: int = 0
    multiple_events: int = 0
    duplicate_events: int = 0
    incomplete_events: int = 0
    committed_batches: int = 0
    duplicate_batches: int = 0
    incomplete_batches: int = 0
    invalid_json_lines: int = 0
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
            "committed_batches": self.committed_batches,
            "duplicate_batches": self.duplicate_batches,
            "incomplete_batches": self.incomplete_batches,
            "invalid_json_lines": self.invalid_json_lines,
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
    position1: list[float], position2: list[float], boxlen: float
) -> list[float]:
    delta = [b - a for a, b in zip(position1, position2)]
    half_box = 0.5 * boxlen
    for index, value in enumerate(delta):
        if value > half_box:
            delta[index] -= boxlen
        if value < -half_box:
            delta[index] += boxlen
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
    boxlen: float | None,
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
    if boxlen is not None and position1 is not None and position2 is not None:
        expected_delta_position = _minimum_image_delta(position1, position2, boxlen)
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
) -> tuple[float | None, float | None, float | None]:
    boxlen = _require_number(begin, "boxlen", uid, errors)
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
        return boxlen, fact_g, merge_radius

    total_mass = sum(masses)
    if total_mass <= 0.0:
        errors.append(f"{uid}: total member mass must be positive")
        return boxlen, fact_g, merge_radius
    if reported_mass is not None and not _close(reported_mass, total_mass):
        errors.append(f"{uid}: total-mass invariant failed")

    if boxlen is None or not positions:
        return boxlen, fact_g, merge_radius
    anchor = positions[0]
    com_position = [0.0, 0.0, 0.0]
    com_velocity = [0.0, 0.0, 0.0]
    for mass, position, velocity in zip(masses, positions, velocities):
        delta = _minimum_image_delta(anchor, position, boxlen)
        for index in range(3):
            com_position[index] += mass * delta[index]
            com_velocity[index] += mass * velocity[index]
    for index in range(3):
        com_position[index] = (anchor[index] + com_position[index] / total_mass) % boxlen
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
            max_separation, _norm(_minimum_image_delta(position1, position2, boxlen))
        )
    if reported_max_separation is not None and not _close(
        reported_max_separation, max_separation
    ):
        errors.append(f"{uid}: maximum-separation invariant failed")
    return boxlen, fact_g, merge_radius


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
    if begin.get("schema_version") not in {1, SCHEMA_VERSION}:
        errors.append(f"{uid}: unsupported schema version")
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

    boxlen, fact_g, merge_radius = _validate_event_invariants(
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
            uid, pair, member_by_id, boxlen, fact_g, merge_radius, errors
        )
    if actual_id_pairs != expected_id_pairs:
        errors.append(f"{uid}: pair records do not cover every member pair exactly once")

    expected_class = "BINARY" if nmember == 2 else "MULTIPLE"
    if begin.get("classification") != expected_class:
        errors.append(f"{uid}: classification must be {expected_class}")
    return _event_digest(rows)


def validate_ledger(
    path: Path,
    *,
    allow_incomplete_tail: bool = False,
    output_root: Path | None = None,
    checkpoint_uid: str | None = None,
) -> LedgerReport:
    report = LedgerReport()
    if checkpoint_uid is not None and output_root is None:
        report.errors.append("checkpoint_uid requires output_root")
    current: EventBlock | None = None
    digests: dict[str, str] = {}
    batches: dict[str, str] = {}
    batch: dict[str, Any] | None = None
    attempts: dict[str, Any] = {}
    active_attempt: str | None = None
    committed_payloads: list[tuple[str, int, int, list[tuple[str, str, Any]]]] = []
    batch_sequences: dict[str, set[int]] = {}

    def register_event(uid: str, digest: str, classification: Any) -> None:
        previous = digests.get(uid)
        if previous is None:
            digests[uid] = digest
            report.unique_events += 1
            if classification == "BINARY":
                report.binary_events += 1
            elif classification == "MULTIPLE":
                report.multiple_events += 1
        elif previous == digest:
            report.duplicate_events += 1
        else:
            report.errors.append(
                f"{uid}: deterministic UID has conflicting event data"
            )

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                continue
            try:
                record = json.loads(raw_line)
            except json.JSONDecodeError as exc:
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
            if record_type == "attempt_begin":
                if current is not None or batch is not None:
                    if current is not None:
                        report.incomplete_events += 1
                    if batch is not None:
                        report.incomplete_batches += 1
                    current = None
                    batch = None
                uid = record.get("attempt_uid")
                parent = record.get("parent_checkpoint_uid")
                run_uuid = record.get("run_uuid")
                ledger_file = record.get("ledger_file")
                valid = (
                    record.get("schema_version") == 2
                    and isinstance(uid, str)
                    and bool(uid)
                    and (parent is None or (isinstance(parent, str) and bool(parent)))
                    and isinstance(run_uuid, str)
                    and bool(run_uuid)
                    and isinstance(ledger_file, str)
                    and bool(ledger_file)
                    and isinstance(record.get("restart_output_index"), int)
                    and not isinstance(record.get("restart_output_index"), bool)
                    and isinstance(record.get("restart_nstep_coarse"), int)
                    and not isinstance(record.get("restart_nstep_coarse"), bool)
                    and ((parent is None) == (record.get("restart_output_index") == 0))
                )
                signature = (
                    run_uuid,
                    ledger_file,
                    parent,
                    record.get("restart_output_index"),
                    record.get("restart_nstep_coarse"),
                )
                if not valid:
                    report.errors.append(f"line {line_number}: invalid attempt_begin")
                elif uid in attempts and attempts[uid] != signature:
                    report.errors.append(f"{uid}: conflicting attempt_begin")
                else:
                    attempts[uid] = signature
                    active_attempt = uid
            elif record_type == "batch_begin":
                if current is not None:
                    report.incomplete_events += 1
                    if not (
                        current.begin.get("schema_version") == 2
                        and batch is not None
                        and current.begin.get("batch_uid") == batch["uid"]
                    ):
                        report.errors.append(
                            f"{current.uid}: missing event_end before line {line_number}"
                        )
                    current = None
                if batch is not None:
                    report.incomplete_batches += 1
                    # A restart may append a replay after a crash left an
                    # earlier prepare uncommitted. Its events stay censored;
                    # the new batch marker is an unambiguous recovery boundary.
                uid = record.get("batch_uid")
                expected = record.get("expected_events")
                attempt_uid = record.get("attempt_uid")
                sequence = record.get("committed_batch_seq")
                batch_begin_valid = not (
                    record.get("schema_version") != 2
                    or not isinstance(uid, str)
                    or not uid
                    or isinstance(expected, bool)
                    or not isinstance(expected, int)
                    or expected < 1
                    or attempt_uid != active_attempt
                    or not isinstance(sequence, int)
                    or isinstance(sequence, bool)
                    or sequence < 1
                    or not isinstance(record.get("nstep_coarse"), int)
                    or isinstance(record.get("nstep_coarse"), bool)
                    or record.get("committed") is not False
                )
                if not batch_begin_valid:
                    report.errors.append(f"line {line_number}: invalid batch_begin")
                batch = {
                    "uid": str(uid or ""),
                    "attempt_uid": attempt_uid,
                    "sequence": sequence,
                    "nstep_coarse": record.get("nstep_coarse"),
                    "begin_record": record,
                    "prepared_record": None,
                    "expected": expected,
                    "events": [],
                    "prepared": False,
                    "valid": batch_begin_valid,
                    # Any subsequent structural, JSON, or invariant error
                    # while this batch is active makes it non-promotable.
                    "error_start": len(report.errors),
                }
            elif record_type == "event_begin":
                if current is not None:
                    report.incomplete_events += 1
                    report.errors.append(
                        f"{current.uid}: missing event_end before line {line_number}"
                    )
                if record.get("schema_version") == 2:
                    if batch is None or record.get("batch_uid") != batch["uid"]:
                        report.errors.append(
                            f"line {line_number}: v2 event appears outside its batch"
                        )
                        if batch is not None:
                            batch["valid"] = False
                elif batch is not None:
                    report.errors.append(
                        f"line {line_number}: legacy event appears inside a v2 batch"
                    )
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
                    if batch is not None:
                        batch["valid"] = False
                if record.get("schema_version") != current.begin.get("schema_version"):
                    report.errors.append(
                        f"line {line_number}: record schema does not match active event"
                    )
                    if batch is not None:
                        batch["valid"] = False
                if current.begin.get("schema_version") == 2 and (
                    batch is None
                    or record.get("batch_uid") != batch["uid"]
                    or record.get("batch_uid") != current.begin.get("batch_uid")
                ):
                    report.errors.append(
                        f"line {line_number}: record batch UID does not match active batch"
                    )
                    if batch is not None:
                        batch["valid"] = False
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
                if record.get("schema_version") != current.begin.get("schema_version"):
                    report.errors.append(
                        f"line {line_number}: event_end schema does not match event_begin"
                    )
                errors_before = len(report.errors)
                digest = _validate_complete_block(current, record, report)
                if batch is not None and len(report.errors) != errors_before:
                    batch["valid"] = False
                if digest is not None:
                    if current.begin.get("schema_version") == 2:
                        if batch is None or record.get("batch_uid") != batch["uid"]:
                            report.errors.append(
                                f"{current.uid}: v2 event_end is outside its batch"
                            )
                        else:
                            batch["events"].append(
                                (current.uid, digest, current.begin.get("classification"))
                            )
                    else:
                        register_event(
                            current.uid, digest, current.begin.get("classification")
                        )
                current = None
            elif record_type == "batch_prepared":
                if current is not None or batch is None:
                    report.errors.append(
                        f"line {line_number}: batch_prepared appears outside a complete batch"
                    )
                    continue
                count = record.get("event_count")
                if (
                    record.get("schema_version") != 2
                    or record.get("batch_uid") != batch["uid"]
                    or record.get("attempt_uid") != batch["attempt_uid"]
                    or record.get("committed_batch_seq") != batch["sequence"]
                    or record.get("nstep_coarse") != batch["nstep_coarse"]
                    or record.get("prepared") is not True
                    or record.get("committed") is not False
                    or count != batch["expected"]
                    or count != len(batch["events"])
                ):
                    report.errors.append(f"{batch['uid']}: invalid batch_prepared")
                    batch["valid"] = False
                else:
                    batch["prepared"] = True
                    batch["prepared_record"] = record
            elif record_type == "batch_commit":
                if current is not None or batch is None:
                    report.errors.append(
                        f"line {line_number}: batch_commit appears outside a prepared batch"
                    )
                    continue
                uid = batch["uid"]
                if (
                    record.get("schema_version") != 2
                    or record.get("batch_uid") != uid
                    or record.get("attempt_uid") != batch["attempt_uid"]
                    or record.get("committed_batch_seq") != batch["sequence"]
                    or record.get("nstep_coarse") != batch["nstep_coarse"]
                    or record.get("committed") is not True
                    or record.get("event_count") != batch["expected"]
                    or not batch["prepared"]
                    or not batch["valid"]
                    or len(report.errors) != batch["error_start"]
                ):
                    report.errors.append(f"{uid}: invalid batch_commit")
                else:
                    payload = json.dumps(
                        {
                            "begin": batch["begin_record"],
                            "events": batch["events"],
                            "prepared": batch["prepared_record"],
                            "commit": record,
                        },
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
                    batch_key = f"{batch['attempt_uid']}:{batch['sequence']}"
                    identity_digest = f"{uid}:{digest}"
                    previous = batches.get(batch_key)
                    if previous is None:
                        batches[batch_key] = identity_digest
                        committed_payloads.append(
                            (
                                batch["attempt_uid"],
                                batch["sequence"],
                                batch["nstep_coarse"],
                                batch["events"],
                            )
                        )
                        batch_sequences.setdefault(batch["attempt_uid"], set()).add(
                            batch["sequence"]
                        )
                        report.committed_batches += 1
                        for event_uid, event_digest, classification in batch["events"]:
                            register_event(
                                f"{batch['attempt_uid']}:{event_uid}",
                                event_digest,
                                classification,
                            )
                    elif previous == identity_digest:
                        report.duplicate_batches += 1
                        report.duplicate_events += len(batch["events"])
                    else:
                        report.errors.append(
                            f"{uid}: committed batch sequence has conflicting data"
                        )
                batch = None
            else:
                report.errors.append(
                    f"line {line_number}: unknown record_type {record_type!r}"
                )

    if current is not None:
        report.incomplete_events += 1
        if not allow_incomplete_tail:
            report.errors.append(f"{current.uid}: incomplete event at end of ledger")
    if batch is not None:
        report.incomplete_batches += 1
        if not allow_incomplete_tail:
            report.errors.append(f"{batch['uid']}: uncommitted batch at end of ledger")
    for attempt_uid, sequences in batch_sequences.items():
        if sequences and (min(sequences) != 1 or len(sequences) != max(sequences)):
            report.errors.append(
                f"{attempt_uid}: committed batch sequence is not contiguous from one"
            )

    if output_root is not None:
        markers: dict[str, dict[str, Any]] = {}
        for marker_path in sorted(output_root.glob("output_*/SMBH_CAPTURE_LINEAGE")):
            complete_path = marker_path.parent / "COMPLETE"
            if not complete_path.is_file():
                continue
            try:
                complete_value = complete_path.read_text(encoding="utf-8").strip()
                lines = marker_path.read_text(encoding="utf-8").splitlines()
                if not lines or lines[0] != "LAGRAMSES_SMBH_CAPTURE_LINEAGE_V1":
                    raise ValueError("bad magic")
                pairs = [line.split("=", 1) for line in lines[1:]]
                values = dict(pairs)
                if len(values) != len(pairs):
                    raise ValueError("duplicate marker field")
                marker = {
                    **values,
                    "committed_batch_seq": int(values["committed_batch_seq"]),
                    "nstep_coarse": int(values["nstep_coarse"]),
                    "output_index": int(values["output_index"]),
                }
                if marker["committed_batch_seq"] < 0:
                    raise ValueError("negative committed batch sequence")
                if marker["checkpoint_uid"] in markers:
                    raise ValueError("duplicate checkpoint UID")
                if marker_path.parent.name != f"output_{marker['output_index']:05d}":
                    raise ValueError("output index mismatch")
                if complete_value != f"{marker['output_index']:05d}":
                    raise ValueError("COMPLETE payload mismatch")
                expected_uid = (
                    f"{marker['attempt_uid']}-output-{marker['output_index']:05d}"
                )
                if marker["checkpoint_uid"] != expected_uid:
                    raise ValueError("checkpoint UID binding mismatch")
                markers[marker["checkpoint_uid"]] = marker
            except (OSError, KeyError, ValueError) as exc:
                report.errors.append(f"{marker_path}: invalid lineage marker: {exc}")

        parent_attempts: set[str] = set()
        parent_for: dict[str, tuple[str, int]] = {}
        for uid, signature in attempts.items():
            run_uuid, ledger_file, parent_checkpoint, restart_output, restart_step = signature
            if parent_checkpoint is None:
                continue
            marker = markers.get(parent_checkpoint)
            if marker is None:
                report.errors.append(f"{uid}: parent checkpoint marker is unavailable")
                continue
            parent_uid = marker["attempt_uid"]
            parent_signature = attempts.get(parent_uid)
            if (
                parent_signature is None
                or marker["run_uuid"] != run_uuid
                or marker["ledger_file"] != ledger_file
                or parent_signature[0] != run_uuid
                or parent_signature[1] != ledger_file
                or marker["output_index"] != restart_output
                or marker["nstep_coarse"] != restart_step
            ):
                report.errors.append(f"{uid}: parent checkpoint identity mismatch")
                continue
            parent_for[uid] = (parent_uid, marker["committed_batch_seq"])
            parent_attempts.add(parent_uid)

        leaves = set(attempts) - parent_attempts
        selected_marker = markers.get(checkpoint_uid) if checkpoint_uid else None
        if checkpoint_uid is not None:
            if selected_marker is None or selected_marker["attempt_uid"] not in leaves:
                report.errors.append(
                    "selected checkpoint is unavailable or does not belong to a leaf attempt"
                )
                selected_leaf = None
            else:
                selected_leaf = selected_marker["attempt_uid"]
        elif len(leaves) == 1:
            selected_leaf = next(iter(leaves))
        else:
            report.errors.append(
                f"lineage has {len(leaves)} leaves; select a COMPLETE checkpoint explicitly"
            )
            selected_leaf = None

        if selected_leaf is not None:
            cursor = selected_leaf
            leaf_markers = [
                marker for marker in markers.values() if marker["attempt_uid"] == cursor
            ]
            if checkpoint_uid is not None:
                selected = selected_marker
            elif len(leaf_markers) == 1:
                selected = leaf_markers[0]
            else:
                report.errors.append(
                    f"leaf attempt has {len(leaf_markers)} COMPLETE checkpoints; select one explicitly"
                )
                selected = None

            cutoffs: dict[str, tuple[int, int]] = {}
            if selected is not None:
                cutoffs[cursor] = (
                    selected["committed_batch_seq"], selected["nstep_coarse"]
                )
            while selected is not None and cursor in parent_for:
                parent_uid, cutoff = parent_for[cursor]
                if parent_uid in cutoffs:
                    report.errors.append("capture lineage contains a cycle")
                    break
                parent_marker = markers[attempts[cursor][2]]
                cutoffs[parent_uid] = (cutoff, parent_marker["nstep_coarse"])
                cursor = parent_uid

            for attempt_uid, (cutoff, _) in cutoffs.items():
                available = batch_sequences.get(attempt_uid, set())
                if cutoff > len(available) or (cutoff > 0 and max(available, default=0) < cutoff):
                    report.errors.append(
                        f"{attempt_uid}: checkpoint high-water exceeds committed ledger sequence"
                    )

            report.committed_batches = 0
            report.duplicate_batches = 0
            report.unique_events = 0
            report.binary_events = 0
            report.multiple_events = 0
            report.duplicate_events = 0
            digests.clear()
            for attempt_uid, sequence, batch_step, events in committed_payloads:
                if attempt_uid not in cutoffs:
                    continue
                cutoff, checkpoint_step = cutoffs[attempt_uid]
                if sequence > cutoff:
                    continue
                if batch_step > checkpoint_step:
                    report.errors.append(
                        f"{attempt_uid}: retained batch step exceeds checkpoint step"
                    )
                    continue
                report.committed_batches += 1
                for event_uid, event_digest, classification in events:
                    register_event(
                        f"{attempt_uid}:{event_uid}", event_digest, classification
                    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument(
        "--allow-incomplete-tail",
        action="store_true",
        help="report but do not fail solely for a final begin block without event_end",
    )
    parser.add_argument(
        "--checkpoint-uid",
        help="terminal COMPLETE checkpoint UID; required when the leaf has multiple outputs",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="resolve restart supersession from COMPLETE output lineage markers",
    )
    args = parser.parse_args()
    if args.checkpoint_uid and args.output_root is None:
        parser.error("--checkpoint-uid requires --output-root")
    report = validate_ledger(
        args.ledger,
        allow_incomplete_tail=args.allow_incomplete_tail,
        output_root=args.output_root,
        checkpoint_uid=args.checkpoint_uid,
    )
    print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    return 0 if report.valid else 1


if __name__ == "__main__":
    sys.exit(main())
