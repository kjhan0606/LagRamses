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
class LedgerReport:
    unique_events: int = 0
    binary_events: int = 0
    multiple_events: int = 0
    duplicate_events: int = 0
    incomplete_events: int = 0
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
    if not all(isinstance(item, (int, float)) and math.isfinite(item) for item in value):
        errors.append(f"{uid}: {key} contains a non-finite value")
        return None
    return [float(item) for item in value]


def _validate_pair_invariants(
    uid: str,
    pair: dict[str, Any],
    members: dict[int, dict[str, Any]],
    errors: list[str],
) -> None:
    id1 = int(pair["sink_id_1"])
    id2 = int(pair["sink_id_2"])
    delta_position = _require_vector(pair, "delta_position_code", uid, errors)
    delta_velocity = _require_vector(pair, "delta_velocity_code", uid, errors)
    specific_h = _require_vector(pair, "specific_angular_momentum_code", uid, errors)
    relative_l = _require_vector(pair, "relative_angular_momentum_code", uid, errors)
    if None in (delta_position, delta_velocity, specific_h, relative_l):
        return

    expected_r = _norm(delta_position)
    expected_v = _norm(delta_velocity)
    if not _close(pair["separation_code"], expected_r):
        errors.append(f"{uid}: pair {id1}-{id2} separation invariant failed")
    if not _close(pair["relative_speed_code"], expected_v):
        errors.append(f"{uid}: pair {id1}-{id2} relative-speed invariant failed")

    mass1 = float(members[id1]["mass_code"])
    mass2 = float(members[id2]["mass_code"])
    expected_mu = mass1 * mass2 / (mass1 + mass2)
    expected_kinetic = 0.5 * expected_mu * expected_v**2
    if not _close(pair["reduced_mass_code"], expected_mu):
        errors.append(f"{uid}: pair {id1}-{id2} reduced-mass invariant failed")
    if not _close(pair["relative_kinetic_code"], expected_kinetic):
        errors.append(f"{uid}: pair {id1}-{id2} kinetic-energy invariant failed")

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
    if begin.get("complete") is not False or end.get("complete") is not True:
        errors.append(f"{uid}: begin/end completion markers are invalid")
    if end.get("event_uid") != uid:
        errors.append(f"{uid}: event_end UID does not match")

    nmember = int(begin.get("nmember", -1))
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
        sink_id = int(member["sink_id"])
        if sink_id in member_by_id:
            errors.append(f"{uid}: duplicate sink member ID {sink_id}")
        member_by_id[sink_id] = member
        _require_vector(member, "position_code", uid, errors)
        _require_vector(member, "velocity_code", uid, errors)

    expected_id_pairs = {
        tuple(sorted(pair)) for pair in combinations(member_by_id.keys(), 2)
    }
    actual_id_pairs: set[tuple[int, int]] = set()
    for pair in block.pairs:
        id_pair = tuple(sorted((int(pair["sink_id_1"]), int(pair["sink_id_2"]))))
        if id_pair in actual_id_pairs:
            errors.append(f"{uid}: duplicate pair record {id_pair}")
        actual_id_pairs.add(id_pair)
        if id_pair[0] not in member_by_id or id_pair[1] not in member_by_id:
            errors.append(f"{uid}: pair {id_pair} references a non-member sink")
            continue
        _validate_pair_invariants(uid, pair, member_by_id, errors)
    if actual_id_pairs != expected_id_pairs:
        errors.append(f"{uid}: pair records do not cover every member pair exactly once")

    expected_class = "BINARY" if nmember == 2 else "MULTIPLE"
    if begin.get("classification") != expected_class:
        errors.append(f"{uid}: classification must be {expected_class}")
    return _event_digest(rows)


def validate_ledger(path: Path, *, allow_incomplete_tail: bool = False) -> LedgerReport:
    report = LedgerReport()
    current: EventBlock | None = None
    digests: dict[str, str] = {}

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

            record_type = record.get("record_type")
            if record_type == "event_begin":
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
                    previous = digests.get(current.uid)
                    if previous is None:
                        digests[current.uid] = digest
                        report.unique_events += 1
                        if current.begin.get("classification") == "BINARY":
                            report.binary_events += 1
                        elif current.begin.get("classification") == "MULTIPLE":
                            report.multiple_events += 1
                    elif previous == digest:
                        report.duplicate_events += 1
                    else:
                        report.errors.append(
                            f"{current.uid}: deterministic UID has conflicting event data"
                        )
                current = None
            else:
                report.errors.append(
                    f"line {line_number}: unknown record_type {record_type!r}"
                )

    if current is not None:
        report.incomplete_events += 1
        if not allow_incomplete_tail:
            report.errors.append(f"{current.uid}: incomplete event at end of ledger")
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
