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
    begin_line: int
    error_count: int
    invalid_json_count: int
    events: list[tuple[EventBlock, str]] = field(default_factory=list)

    @property
    def uid(self) -> str:
        return str(self.begin.get("batch_uid", ""))


@dataclass
class AttemptBlock:
    resume_step: int | None
    restart_output: int | None
    parent_index: int | None
    parent_batch_cutoff: int | None


@dataclass
class LedgerReport:
    unique_events: int = 0
    binary_events: int = 0
    multiple_events: int = 0
    duplicate_events: int = 0
    incomplete_events: int = 0
    invalid_json_lines: int = 0
    committed_batches: int = 0
    incomplete_batches: int = 0
    censored_events: int = 0
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
            "invalid_json_lines": self.invalid_json_lines,
            "committed_batches": self.committed_batches,
            "incomplete_batches": self.incomplete_batches,
            "censored_events": self.censored_events,
            "run_attempts": self.run_attempts,
            "superseded_batches": self.superseded_batches,
            "superseded_events": self.superseded_events,
            "errors": self.errors,
        }


def _close(
    a: float, b: float, *, rtol: float = 2.0e-12, scale: float = 0.0
) -> bool:
    # Code-unit masses and pair energies can be far below 1e-14.  A fixed
    # absolute floor would accept even sign-flipped values at those scales.
    # For cancellation-dominated expressions, callers pass the size of the
    # uncancelled terms as scale instead of the small residual.
    return abs(float(a) - float(b)) <= rtol * max(
        abs(float(a)), abs(float(b)), abs(float(scale))
    )


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
    for index, (value, extent) in enumerate(zip(delta, box_size)):
        half_box = 0.5 * extent
        if value > half_box:
            delta[index] -= extent
        if value < -half_box:
            delta[index] += extent
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

    position1 = _require_vector(members[id1], "position_code", uid, errors)
    position2 = _require_vector(members[id2], "position_code", uid, errors)
    velocity1 = _require_vector(members[id1], "velocity_code", uid, errors)
    velocity2 = _require_vector(members[id2], "velocity_code", uid, errors)
    separation_scale = (
        _norm(max(abs(a), abs(b)) for a, b in zip(position1, position2))
        if position1 is not None and position2 is not None else 0.0
    )
    speed_scale = (
        _norm(max(abs(a), abs(b)) for a, b in zip(velocity1, velocity2))
        if velocity1 is not None and velocity2 is not None else 0.0
    )

    separation = _require_number(pair, "separation_code", uid, errors)
    relative_speed = _require_number(pair, "relative_speed_code", uid, errors)
    reduced_mass = _require_number(pair, "reduced_mass_code", uid, errors)
    relative_kinetic = _require_number(pair, "relative_kinetic_code", uid, errors)
    expected_r2 = sum(value * value for value in delta_position)
    expected_v2 = sum(value * value for value in delta_velocity)
    expected_r = math.sqrt(expected_r2)
    expected_v = math.sqrt(expected_v2)
    if separation is not None and not _close(separation, expected_r, scale=separation_scale):
        errors.append(f"{uid}: pair {id1}-{id2} separation invariant failed")
    if relative_speed is not None and not _close(relative_speed, expected_v, scale=speed_scale):
        errors.append(f"{uid}: pair {id1}-{id2} relative-speed invariant failed")

    mass1 = _require_number(members[id1], "mass_code", uid, errors)
    mass2 = _require_number(members[id2], "mass_code", uid, errors)
    if mass1 is None or mass2 is None or mass1 <= 0.0 or mass2 <= 0.0:
        errors.append(f"{uid}: pair {id1}-{id2} has a non-positive member mass")
        return
    expected_mu = mass1 * mass2 / (mass1 + mass2)
    expected_kinetic = 0.5 * expected_mu * expected_v2
    if reduced_mass is not None and not _close(reduced_mass, expected_mu):
        errors.append(f"{uid}: pair {id1}-{id2} reduced-mass invariant failed")
    if relative_kinetic is not None and not _close(
        relative_kinetic, expected_kinetic, scale=expected_mu * expected_v * speed_scale
    ):
        errors.append(f"{uid}: pair {id1}-{id2} kinetic-energy invariant failed")

    if box_size is not None and position1 is not None and position2 is not None:
        expected_delta_position = _minimum_image_delta(position1, position2, box_size)
        if any(
            not _close(
                got,
                expected,
                scale=max(abs(position1[i]), abs(position2[i])),
            )
            for i, (got, expected) in enumerate(
                zip(delta_position, expected_delta_position)
            )
        ):
            errors.append(
                f"{uid}: pair {id1}-{id2} minimum-image position invariant failed"
            )
    if velocity1 is not None and velocity2 is not None:
        expected_delta_velocity = [b - a for a, b in zip(velocity1, velocity2)]
        if any(
            not _close(
                got,
                expected,
                scale=max(abs(velocity1[i]), abs(velocity2[i])),
            )
            for i, (got, expected) in enumerate(zip(delta_velocity, expected_delta_velocity))
        ):
            errors.append(
                f"{uid}: pair {id1}-{id2} member-velocity invariant failed"
            )

    expected_h = _cross(delta_position, delta_velocity)
    h_scale = [
        abs(delta_position[1] * delta_velocity[2])
        + abs(delta_position[2] * delta_velocity[1]),
        abs(delta_position[2] * delta_velocity[0])
        + abs(delta_position[0] * delta_velocity[2]),
        abs(delta_position[0] * delta_velocity[1])
        + abs(delta_position[1] * delta_velocity[0]),
    ]
    for got, expected, scale in zip(specific_h, expected_h, h_scale):
        if not _close(got, expected, scale=scale):
            errors.append(
                f"{uid}: pair {id1}-{id2} specific-angular-momentum invariant failed"
            )
            break
    for got, expected, scale in zip(
        relative_l, (expected_mu * x for x in expected_h), h_scale
    ):
        if not _close(got, expected, scale=expected_mu * scale):
            errors.append(
                f"{uid}: pair {id1}-{id2} relative-angular-momentum invariant failed"
            )
            break

    within_rmerge = _require_logical(pair, "within_rmerge", uid, errors)
    two_body_bound = _require_logical(pair, "two_body_bound", uid, errors)
    legacy_pair_bound = _require_logical(pair, "legacy_pair_bound", uid, errors)
    if merge_radius is not None and within_rmerge is not None:
        expected_within = expected_r2 <= merge_radius**2
        if not _close(expected_r2, merge_radius**2,
                      scale=2.0 * expected_r * separation_scale) and within_rmerge is not expected_within:
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
    specific_kinetic = 0.5 * expected_v2
    specific_gravity = fact_g * (mass1 + mass2) / expected_r
    expected_specific_energy = specific_kinetic - specific_gravity
    expected_legacy_proxy = fact_g * mass1 * mass2 / expected_r2
    radius_error_factor = separation_scale / expected_r
    if potential_value is not None and not _close(
        potential_value, expected_potential,
        scale=abs(expected_potential) * radius_error_factor,
    ):
        errors.append(f"{uid}: pair {id1}-{id2} potential-energy invariant failed")
    if specific_energy_value is not None and not _close(
        specific_energy_value,
        expected_specific_energy,
        scale=abs(specific_kinetic) + abs(specific_gravity)
        + expected_v * speed_scale + abs(specific_gravity) * radius_error_factor,
    ):
        errors.append(f"{uid}: pair {id1}-{id2} specific-energy invariant failed")
    if legacy_proxy_value is not None and not _close(
        legacy_proxy_value, expected_legacy_proxy,
        scale=2.0 * abs(expected_legacy_proxy) * radius_error_factor,
    ):
        errors.append(f"{uid}: pair {id1}-{id2} legacy-binding-proxy invariant failed")
    if (
        two_body_bound is not None
        and not _close(
            expected_specific_energy,
            0.0,
            scale=abs(specific_kinetic) + abs(specific_gravity)
            + expected_v * speed_scale + abs(specific_gravity) * radius_error_factor,
        )
        and two_body_bound is not (expected_specific_energy < 0.0)
    ):
        errors.append(f"{uid}: pair {id1}-{id2} two-body-bound invariant failed")
    if (
        legacy_pair_bound is not None
        and not _close(
            expected_kinetic,
            expected_legacy_proxy,
            scale=abs(expected_kinetic) + abs(expected_legacy_proxy)
            + expected_mu * expected_v * speed_scale
            + 2.0 * abs(expected_legacy_proxy) * radius_error_factor,
        )
        and legacy_pair_bound is not (expected_kinetic < expected_legacy_proxy)
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
        else ([boxlen] * 3 if boxlen is not None else None)
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
        if mass <= 0.0:
            errors.append(f"{uid}: member mass must be positive")
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
    velocity_scale = [
        sum(abs(mass * velocity[index]) for mass, velocity in zip(masses, velocities))
        / total_mass
        for index in range(3)
    ]
    if reported_com_position is not None and any(
        not _close(
            _minimum_image_delta([expected], [got], [box_size[index]])[0],
            0.0, scale=box_size[index],
        )
        for index, (got, expected) in enumerate(
            zip(reported_com_position, com_position)
        )
    ):
        errors.append(f"{uid}: centre-of-mass position invariant failed")
    if reported_com_velocity is not None and any(
        not _close(got, expected, scale=velocity_scale[index])
        for index, (got, expected) in enumerate(
            zip(reported_com_velocity, com_velocity)
        )
    ):
        errors.append(f"{uid}: centre-of-mass velocity invariant failed")

    max_separation = 0.0
    max_separation_scale = 0.0
    for position1, position2 in combinations(positions, 2):
        max_separation = max(
            max_separation, _norm(_minimum_image_delta(position1, position2, box_size))
        )
        max_separation_scale = max(
            max_separation_scale,
            _norm(max(abs(a), abs(b)) for a, b in zip(position1, position2)),
        )
    if reported_max_separation is not None and not _close(
        reported_max_separation, max_separation, scale=max_separation_scale
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
    if any(row.get("schema_version") != SCHEMA_VERSION for row in rows):
        errors.append(f"{uid}: unsupported schema version in transaction")
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
    pair_indices = [row.get("pair_index") for row in block.pairs]
    if pair_indices != list(range(1, expected_pairs + 1)):
        errors.append(f"{uid}: pair_index sequence is not contiguous")

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


def validate_ledger(
    path: Path,
    *,
    allow_incomplete_tail: bool = False,
    allow_incomplete_batches: bool = False,
) -> LedgerReport:
    report = LedgerReport()
    current: EventBlock | None = None
    batch: BatchBlock | None = None
    seen_batch_protocol = False
    active_resume_step: int | None = None
    legacy_events: list[tuple[EventBlock, str]] = []
    committed: list[tuple[BatchBlock, str, int | None]] = []
    attempts: list[AttemptBlock] = []
    # A checkpoint marker records the exact committed-batch prefix present in
    # its snapshot.  Coarse-step equality cannot distinguish captures before
    # and after an output written in that same step.
    checkpoint_owner: dict[int, tuple[int, int, int]] = {}
    digests: dict[str, str] = {}
    batch_digests: dict[str, str] = {}

    def accept_event(block: EventBlock, digest: str) -> None:
        previous = digests.get(block.uid)
        if previous is None:
            digests[block.uid] = digest
            report.unique_events += 1
            if block.begin.get("classification") == "BINARY":
                report.binary_events += 1
            elif block.begin.get("classification") == "MULTIPLE":
                report.multiple_events += 1
        elif previous == digest:
            report.duplicate_events += 1
        else:
            report.errors.append(
                f"{block.uid}: deterministic UID has conflicting event data"
            )

    def censor_batch(reason: str, *, allowed: bool) -> None:
        nonlocal batch
        assert batch is not None
        report.incomplete_batches += 1
        report.censored_events += len(batch.events)
        if not allowed:
            report.errors.append(f"{batch.uid}: {reason}")
        batch = None

    def commit_batch(record: dict[str, Any]) -> None:
        nonlocal batch
        assert batch is not None
        uid = batch.uid
        begin = batch.begin
        if begin.get("schema_version") != SCHEMA_VERSION or record.get(
            "schema_version"
        ) != SCHEMA_VERSION:
            report.errors.append(f"{uid}: unsupported batch schema version")
        if not uid or record.get("batch_uid") != uid:
            report.errors.append(f"{uid}: batch_commit UID does not match batch_begin")
        step = _require_integer(begin, "nstep_coarse", uid, report.errors)
        level = _require_integer(begin, "ilevel", uid, report.errors)
        before = _require_integer(begin, "nsink_before", uid, report.errors)
        after = _require_integer(begin, "nsink_after", uid, report.errors)
        expected = _require_integer(begin, "expected_events", uid, report.errors)
        commit_after = _require_integer(record, "nsink_after", uid, report.errors)
        if None not in (step, level, before, after):
            if uid != f"{step}-{level}-{before}-{after}":
                report.errors.append(f"{uid}: batch UID is inconsistent with metadata")
            if before <= after or after < 1:
                report.errors.append(f"{uid}: invalid sink-count transition")
        if step is not None and active_resume_step is not None and step < active_resume_step:
            report.errors.append(f"{uid}: batch predates active attempt resume_step")
        if commit_after != after:
            report.errors.append(f"{uid}: committed sink count does not match batch")
        if expected is None or expected < 1 or len(batch.events) != expected:
            report.errors.append(f"{uid}: batch event count does not match")
        if before is not None and after is not None:
            sink_reduction = sum(
                block.begin["nmember"] - 1
                for block, _ in batch.events
                if isinstance(block.begin.get("nmember"), int)
                and not isinstance(block.begin["nmember"], bool)
            )
            if before - after != sink_reduction:
                report.errors.append(f"{uid}: batch sink-count conservation failed")
        event_uids = [block.uid for block, _ in batch.events]
        if len(set(event_uids)) != len(event_uids):
            report.errors.append(f"{uid}: duplicate event UID inside batch")
        seen_members: set[int] = set()
        for block, _ in batch.events:
            if block.begin.get("nstep_coarse") != step or block.begin.get("ilevel") != level:
                report.errors.append(f"{uid}: event step/level differs from batch")
            for required in ("primary_sink_id", "periodic_box_size_code"):
                if required not in block.begin:
                    report.errors.append(f"{uid}: batched event lacks {required}")
            member_ids = {
                member["sink_id"]
                for member in block.members
                if isinstance(member.get("sink_id"), int)
                and not isinstance(member["sink_id"], bool)
            }
            if step is not None and level is not None and member_ids:
                expected_uid = (
                    f"{step}-{level}-{min(member_ids)}-{max(member_ids)}-"
                    f"{len(block.members)}"
                )
                if block.uid != expected_uid:
                    report.errors.append(f"{uid}: event UID does not match member IDs")
            if member_ids & seen_members:
                report.errors.append(f"{uid}: sink ID occurs in multiple groups")
            seen_members.update(member_ids)
        if (
            len(report.errors) == batch.error_count
            and report.invalid_json_lines == batch.invalid_json_count
        ):
            payload = json.dumps(
                [begin, [(block.uid, digest) for block, digest in batch.events], record],
                sort_keys=True,
                separators=(",", ":"),
            )
            digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
            owner = len(attempts) - 1 if attempts else None
            committed.append((batch, digest, owner))
        batch = None

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
                if current is not None:
                    report.incomplete_events += 1
                    if batch is None or not allow_incomplete_batches:
                        report.errors.append(
                            f"{current.uid}: missing event_end before line {line_number}"
                        )
                    current = None
                if batch is not None:
                    censor_batch(
                        f"missing batch_commit before line {line_number}",
                        allowed=allow_incomplete_batches,
                    )
                if record.get("schema_version") != SCHEMA_VERSION:
                    report.errors.append(f"line {line_number}: unsupported attempt schema")
                resume_step = _require_integer(
                    record, "resume_step", f"line {line_number}", report.errors
                )
                restart_output = _require_integer(
                    record, "restart_output", f"line {line_number}", report.errors
                )
                if resume_step is not None and resume_step < 0:
                    report.errors.append(f"line {line_number}: negative resume_step")
                if restart_output is not None and restart_output < 0:
                    report.errors.append(f"line {line_number}: negative restart_output")
                if restart_output == 0 and (
                    report.run_attempts > 0 or committed or legacy_events
                ):
                    report.errors.append(
                        f"line {line_number}: fresh run appended to an existing ledger"
                    )
                if legacy_events and restart_output is not None and restart_output > 0:
                    report.errors.append(
                        f"line {line_number}: cannot prove restart censoring for legacy bare events"
                    )
                parent_index = None
                parent_batch_cutoff = None
                if restart_output is not None and restart_output > 0:
                    checkpoint = checkpoint_owner.get(restart_output)
                    if checkpoint is None:
                        report.errors.append(
                            f"line {line_number}: restart output {restart_output} has no ledger checkpoint"
                        )
                    else:
                        parent_index, checkpoint_step, parent_batch_cutoff = checkpoint
                        if resume_step != checkpoint_step:
                            report.errors.append(
                                f"line {line_number}: restart step differs from checkpoint"
                            )
                attempts.append(AttemptBlock(
                    resume_step, restart_output, parent_index, parent_batch_cutoff
                ))
                active_resume_step = resume_step
                report.run_attempts += 1
                seen_batch_protocol = True
            elif record_type == "checkpoint":
                if current is not None or batch is not None:
                    report.errors.append(
                        f"line {line_number}: checkpoint appears inside an event or batch"
                    )
                    continue
                if record.get("schema_version") != SCHEMA_VERSION or not attempts:
                    report.errors.append(
                        f"line {line_number}: checkpoint has no valid run attempt"
                    )
                    continue
                output_number = _require_integer(
                    record, "output_number", f"line {line_number}", report.errors
                )
                checkpoint_step = _require_integer(
                    record, "nstep_coarse", f"line {line_number}", report.errors
                )
                if output_number is not None and output_number <= 0:
                    report.errors.append(f"line {line_number}: invalid output_number")
                if (
                    checkpoint_step is not None
                    and active_resume_step is not None
                    and checkpoint_step < active_resume_step
                ):
                    report.errors.append(
                        f"line {line_number}: checkpoint predates active attempt"
                    )
                if (
                    output_number is not None
                    and output_number > 0
                    and checkpoint_step is not None
                ):
                    checkpoint_owner[output_number] = (
                        len(attempts) - 1,
                        checkpoint_step,
                        len(committed),
                    )
            elif record_type == "batch_begin":
                seen_batch_protocol = True
                if current is not None:
                    report.incomplete_events += 1
                    if batch is None or not allow_incomplete_batches:
                        report.errors.append(
                            f"{current.uid}: missing event_end before line {line_number}"
                        )
                    current = None
                if batch is not None:
                    censor_batch(
                        f"missing batch_commit before line {line_number}",
                        allowed=allow_incomplete_batches,
                    )
                batch = BatchBlock(
                    record, line_number, len(report.errors), report.invalid_json_lines
                )
            elif record_type == "event_begin":
                if batch is None and seen_batch_protocol:
                    report.errors.append(
                        f"line {line_number}: bare event after batch protocol began"
                    )
                if current is not None:
                    report.incomplete_events += 1
                    report.errors.append(
                        f"{current.uid}: missing event_end before line {line_number}"
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
                digest = _validate_complete_block(current, record, report)
                if digest is not None:
                    if batch is not None:
                        batch.events.append((current, digest))
                    elif not seen_batch_protocol:
                        legacy_events.append((current, digest))
                current = None
            elif record_type == "batch_commit":
                if batch is None:
                    report.errors.append(
                        f"line {line_number}: batch_commit appears outside a batch"
                    )
                elif current is not None:
                    report.incomplete_events += 1
                    report.errors.append(f"{current.uid}: missing event_end before batch_commit")
                    current = None
                    censor_batch("batch_commit followed an incomplete event", allowed=False)
                else:
                    commit_batch(record)
            else:
                report.errors.append(
                    f"line {line_number}: unknown record_type {record_type!r}"
                )

    if current is not None:
        report.incomplete_events += 1
        if batch is None and not allow_incomplete_tail:
            report.errors.append(f"{current.uid}: incomplete event at end of ledger")
    if batch is not None:
        censor_batch(
            "incomplete batch at end of ledger",
            allowed=allow_incomplete_tail or allow_incomplete_batches,
        )
    active_cutoffs: dict[int, int | None] = {}
    if attempts:
        attempt_index: int | None = len(attempts) - 1
        cutoff: int | None = None
        while attempt_index is not None:
            if attempt_index in active_cutoffs:
                report.errors.append("restart checkpoint lineage contains a cycle")
                break
            active_cutoffs[attempt_index] = cutoff
            attempt = attempts[attempt_index]
            cutoff = attempt.parent_batch_cutoff
            attempt_index = attempt.parent_index
    for block, digest in legacy_events:
        accept_event(block, digest)
    for batch_index, (committed_batch, digest, owner) in enumerate(committed):
        if attempts and owner is None:
            report.errors.append(
                f"{committed_batch.uid}: batch predates first run attempt marker"
            )
            continue
        if attempts and (
            owner not in active_cutoffs
            or (
                active_cutoffs[owner] is not None
                and batch_index >= active_cutoffs[owner]
            )
        ):
            report.superseded_batches += 1
            report.superseded_events += len(committed_batch.events)
            continue
        uid = committed_batch.uid
        previous = batch_digests.get(uid)
        if previous is None:
            batch_digests[uid] = digest
            report.committed_batches += 1
        elif previous != digest:
            report.errors.append(f"{uid}: deterministic batch UID has conflicting data")
            continue
        for block, event_digest in committed_batch.events:
            accept_event(block, event_digest)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument(
        "--allow-incomplete-tail",
        action="store_true",
        help="censor a final incomplete event or batch without treating it as fatal",
    )
    parser.add_argument(
        "--allow-incomplete-batches",
        action="store_true",
        help="censor uncommitted batches superseded by a restart attempt",
    )
    args = parser.parse_args()
    report = validate_ledger(
        args.ledger,
        allow_incomplete_tail=args.allow_incomplete_tail,
        allow_incomplete_batches=args.allow_incomplete_batches,
    )
    print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    return 0 if report.valid else 1


if __name__ == "__main__":
    sys.exit(main())
