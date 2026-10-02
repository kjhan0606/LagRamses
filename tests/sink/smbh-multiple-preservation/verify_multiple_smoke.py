#!/usr/bin/env python3
"""Check live SMBH IDs and numerical sink state in the tiny MPI fixture."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import struct


EXPECTED_IDS = {1, 2, 3, 5}


def read_record(stream) -> bytes:
    prefix = stream.read(4)
    if len(prefix) != 4:
        raise ValueError("truncated Fortran record marker")
    size = struct.unpack("<i", prefix)[0]
    if size < 0 or size > 100_000_000:
        raise ValueError("invalid Fortran record size")
    data = stream.read(size)
    suffix = stream.read(4)
    if len(data) != size or len(suffix) != 4 or struct.unpack("<i", suffix)[0] != size:
        raise ValueError("truncated or inconsistent Fortran record")
    return data


def read_checkpoint(path: Path) -> dict[int, tuple[float, ...]]:
    with path.open("rb") as stream:
        nsink = struct.unpack("<i", read_record(stream))[0]
        nindsink = struct.unpack("<i", read_record(stream))[0]
        if nsink != 4 or nindsink != 5:
            raise ValueError(f"expected four live sinks and five assigned IDs, got {nsink}/{nindsink}")
        ids = struct.unpack(f"<{nsink}i", read_record(stream))
        if set(ids) != EXPECTED_IDS or len(set(ids)) != nsink:
            raise ValueError(f"live SMBH identities were not preserved: {ids}")
        # The production backup_sink format writes the physical state as one
        # Fortran record per field, after the ID record. Read the first 20
        # state fields through eps_sink; sink_stat follows by level/channel.
        fields = []
        for _ in range(20):
            field = struct.unpack(f"<{nsink}d", read_record(stream))
            if not all(math.isfinite(value) for value in field):
                raise ValueError("non-finite live sink checkpoint field")
            fields.append(field)
    by_id = {sink_id: tuple(field[index] for field in fields)
             for index, sink_id in enumerate(ids)}
    if sum(by_id[sink_id][0] for sink_id in ids) <= 0.0:
        raise ValueError("non-positive retained SMBH mass")
    return by_id


def check_ledger(path: Path) -> tuple[int, int]:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    binary = 0
    multiple = 0
    for index, row in enumerate(rows):
        if row.get("record_type") != "event_begin":
            continue
        members = []
        for following in rows[index + 1:]:
            if following.get("record_type") == "event_end":
                break
            if following.get("record_type") == "member":
                members.append(following["sink_id"])
        if row.get("classification") == "BINARY" and set(members) == {1, 4}:
            binary += 1
        if row.get("classification") == "MULTIPLE" and set(members) == {2, 3, 5}:
            if row.get("multiple_members_preserved") is not True:
                raise ValueError("MULTIPLE ledger omitted live-preservation declaration")
            multiple += 1
    if binary == 0 or multiple == 0:
        raise ValueError("original binary and MULTIPLE groups were not both recorded")
    return binary, multiple


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    args = parser.parse_args()
    root = args.run_directory.resolve()
    ids_by_output = {}
    for output in (1, 2):
        path = root / f"output_{output:05d}" / f"sink_{output:05d}.out"
        ids_by_output[output] = read_checkpoint(path)
    binary, multiple = check_ledger(root / "smbh_capture_ledger_v1.jsonl")
    print(json.dumps({
        "status": "live_multiple_smoke_passed",
        "retained_sink_ids": sorted(ids_by_output[2]),
        "binary_events": binary,
        "multiple_events": multiple,
        "checkpoint_outputs": [1, 2],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
