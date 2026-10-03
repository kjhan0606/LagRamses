#!/usr/bin/env python3
"""Evaluate a bounded capture/restart structural smoke, not its physics."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re


EXPECTED_UNITS = {
    "unit_l": 3.085677581282e21,
    "unit_d": 1.66e-24,
    "unit_t": 3.1556926e13,
}


def _event_blocks(rows: list[dict]) -> list[list[dict]]:
    blocks: list[list[dict]] = []
    active: list[dict] | None = None
    for row in rows:
        kind = row["record_type"]
        if kind == "event_begin":
            if active is not None:
                raise ValueError("nested event in capture ledger")
            active = [row]
        elif active is not None:
            active.append(row)
            if kind == "event_end":
                blocks.append(active)
                active = None
    if active is not None:
        raise ValueError("incomplete event in capture ledger")
    return blocks


def evaluate(run: Path) -> dict:
    run = run.resolve()
    fresh = json.loads((run / "fresh_validation.json").read_text())
    restart = json.loads((run / "restart_validation.json").read_text())
    assert fresh["status"] == restart["status"] == "valid"
    assert fresh["multiple_events"] == restart["multiple_events"] == 1
    assert fresh["committed_batches"] == restart["committed_batches"] == 1
    assert restart["superseded_batches"] == restart["superseded_events"] == 1
    assert restart["run_attempts"] == 2 and restart["unique_events"] == 1

    outputs = sorted(p.name for p in run.iterdir() if p.is_dir() and p.name.startswith("output_"))
    assert outputs == ["output_00001"], outputs
    assert (run / outputs[0] / "COMPLETE").is_file()
    fresh_nml = (run / "fresh.nml").read_text()
    restart_nml = (run / "restart.nml").read_text()
    hdf5 = 'outformat="hdf5"' in fresh_nml
    assert hdf5 == ('informat="hdf5"' in restart_nml)
    assert (run / outputs[0] / "data_00001.h5").is_file() == hdf5
    info = (run / outputs[0] / "info_00001.txt").read_text()
    for key, expected in EXPECTED_UNITS.items():
        match = re.search(rf"^{key}\s*=\s*(\S+)", info, re.MULTILINE)
        assert match and math.isclose(float(match.group(1)), expected, rel_tol=1e-12)
    assert re.search(r"^nstep_coarse=\s*0\s*$", info, re.MULTILINE)

    rows = [json.loads(line) for line in (run / "capture.jsonl").read_text().splitlines() if line.strip()]
    kinds = [row["record_type"] for row in rows]
    assert kinds.count("attempt_begin") == 2
    assert kinds.count("batch_begin") == kinds.count("batch_commit") == 2
    assert kinds.count("checkpoint") == 1
    assert kinds.index("checkpoint") < kinds.index("batch_begin")
    blocks = _event_blocks(rows)
    assert len(blocks) == 2 and blocks[0] == blocks[1]
    assert blocks[0][0]["classification"] == "MULTIPLE"
    assert blocks[0][0]["nmember"] == 3
    assert [row["sink_id"] for row in blocks[0] if row["record_type"] == "member"] == [1, 2, 3]
    assert len([row for row in blocks[0] if row["record_type"] == "pair"]) == 3

    diagnostics = {}
    rank_counts = []
    for name in ("fresh", "restart"):
        log = (run / f"{name}.log").read_text()
        assert "Run completed" in log
        assert "FFT direct solve DONE" in log
        assert "Fine multigrid Poisson failed to converge" not in log
        if hdf5 and name == "restart":
            assert "HDF5 restore AMR" in log
        rank_match = re.search(r"Working with nproc\s*=\s*(\d+)", log)
        assert rank_match, name
        rank_counts.append(int(rank_match.group(1)))
        diagnostics[name] = {
            "negative_internal_energy_cells": [
                int(value) for value in re.findall(r"neg_cells=\s*(\d+)", log)
            ],
            "coarse_conservation": [
                {"step": int(step), "mass_error": float(mass), "energy_error": float(energy)}
                for step, mass, energy in re.findall(
                    r"Main step=\s*(\d+) mcons=\s*([+\-\d.Ee]+) econs=\s*([+\-\d.Ee]+)", log
                )
            ],
        }
        assert len(diagnostics[name]["negative_internal_energy_cells"]) == 2
        assert len(diagnostics[name]["coarse_conservation"]) == 2
        assert [item["step"] for item in diagnostics[name]["coarse_conservation"]] == [1, 2]
        assert all(
            math.isfinite(item[key])
            for item in diagnostics[name]["coarse_conservation"]
            for key in ("mass_error", "energy_error")
        )
    assert (diagnostics["fresh"]["negative_internal_energy_cells"] ==
            diagnostics["restart"]["negative_internal_energy_cells"])
    hydro_replay_identical = diagnostics["fresh"] == diagnostics["restart"]
    # The legacy binary checkpoint reproduces these coarse diagnostics in
    # this fixture. A HDF5 restore may follow the cold-Poisson restart path;
    # retain the measured discrepancy and never use this structural smoke as
    # a hydro-conservation or trajectory-equivalence claim.
    if not hdf5:
        assert hydro_replay_identical
    assert rank_counts[0] == rank_counts[1]
    assert blocks[0][0]["ncpu"] == rank_counts[0]
    info_ncpu = re.search(r"^ncpu\s*=\s*(\d+)", info, re.MULTILINE)
    assert info_ncpu and int(info_ncpu.group(1)) == rank_counts[0]
    return {
        "status": "structural_pass",
        "physics_admission": "not_admitted_synthetic_initial_conditions",
        "run_directory": str(run),
        "checkpoint": outputs[0],
        "checkpoint_format": "hdf5" if hdf5 else "original",
        "mpi_ranks": rank_counts[0],
        "active_multiple_events": 1,
        "superseded_batches": 1,
        "replayed_event_records_identical": True,
        "hydro_replay_identical": hydro_replay_identical,
        "diagnostics": diagnostics,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    print(json.dumps(evaluate(args.run), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
