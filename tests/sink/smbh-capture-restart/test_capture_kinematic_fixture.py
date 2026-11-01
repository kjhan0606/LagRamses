"""Unit tests for the moving-sink capture/restart fixture preflight."""

from __future__ import annotations

import csv
from pathlib import Path
from runpy import run_path

import pytest


HERE = Path(__file__).resolve().parent
IC_SINK = HERE / "ic_sink"
VALIDATOR = run_path(str(HERE / "prepare_and_validate.py"))[
    "validate_fixture_kinematics"
]
SINK_STATE_VALIDATOR = run_path(str(HERE / "prepare_and_validate.py"))[
    "validate_sink_state"
]


def _fixture_rows() -> list[list[float]]:
    with IC_SINK.open(newline="", encoding="utf-8") as stream:
        return [[float(value) for value in row] for row in csv.reader(stream) if row]


def test_fixture_has_nonzero_internal_motion_and_zero_group_com_momentum() -> None:
    rows = _fixture_rows()
    VALIDATOR(rows)


def test_fixture_rejects_stationary_members() -> None:
    rows = _fixture_rows()
    for row in rows:
        row[4:7] = [0.0, 0.0, 0.0]
    with pytest.raises(ValueError, match="nonzero internal angular momentum"):
        VALIDATOR(rows)


def test_fixture_rejects_nonzero_group_com_momentum() -> None:
    rows = _fixture_rows()
    rows[1][5] += 1.0e-4
    with pytest.raises(ValueError, match="zero initial COM momentum"):
        VALIDATOR(rows)


def _write_synthetic_postmerge_outputs(run_dir: Path) -> list[dict]:
    seeds = _fixture_rows()
    event_specs = ((2, (0, 1)), (5, (2, 3, 4)))
    events = []
    for primary, indices in event_specs:
        members = [seeds[index] for index in indices]
        mass = sum(row[0] for row in members)
        position = [
            sum(row[0] * (row[axis + 1] + 0.5) for row in members) / mass
            for axis in range(3)
        ]
        velocity = [sum(row[0] * row[axis + 4] for row in members) / mass
                    for axis in range(3)]
        events.append({
            "primary_sink_id": primary,
            "total_mass_code": mass,
            "com_position_code": position,
            "com_velocity_code": velocity,
            "periodic_box_size_code": [1.0, 1.0, 1.0],
        })

    for output_index in (1, 2):
        output = run_dir / f"output_{output_index:05d}"
        output.mkdir(parents=True)
        with (output / f"sink_{output_index:05d}.csv").open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.writer(stream)
            for event in events:
                writer.writerow([
                    event["primary_sink_id"], event["total_mass_code"],
                    *event["com_position_code"], *event["com_velocity_code"], 0.0, 0.0,
                ])
    (run_dir / "ic_sink").write_text(IC_SINK.read_text(encoding="utf-8"), encoding="utf-8")
    with (run_dir / "sink_log.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        for sink_id in (2, 2, 5, 5):
            writer.writerow([0.0, 0.0, 0.0, sink_id, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return events


def test_postmerge_outputs_preserve_precompaction_com_and_restart_state(
    tmp_path: Path,
) -> None:
    events = _write_synthetic_postmerge_outputs(tmp_path)
    SINK_STATE_VALIDATOR(tmp_path, events)


def test_postmerge_output_rejects_com_velocity_mismatch(tmp_path: Path) -> None:
    events = _write_synthetic_postmerge_outputs(tmp_path)
    snapshot = tmp_path / "output_00002" / "sink_00002.csv"
    with snapshot.open(encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    rows[0][5] = "0.001"
    with snapshot.open("w", newline="", encoding="utf-8") as stream:
        csv.writer(stream).writerows(rows)
    with pytest.raises(SystemExit, match="COM velocity"):
        SINK_STATE_VALIDATOR(tmp_path, events)


def test_postmerge_com_accepts_ramses_sink_csv_precision(tmp_path: Path) -> None:
    events = _write_synthetic_postmerge_outputs(tmp_path)
    for output_index in (1, 2):
        snapshot = tmp_path / f"output_{output_index:05d}" / f"sink_{output_index:05d}.csv"
        with snapshot.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.reader(stream))
        for row in rows:
            row[2:5] = [f"{float(value):.10E}" for value in row[2:5]]
        with snapshot.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerows(rows)
    SINK_STATE_VALIDATOR(tmp_path, events)
