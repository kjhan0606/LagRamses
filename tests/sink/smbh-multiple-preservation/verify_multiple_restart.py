#!/usr/bin/env python3
"""Compare a one-step MPI restart with the uninterrupted sink checkpoint."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from verify_multiple_smoke import read_checkpoint, read_sink_statistics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("continuous_directory", type=Path)
    parser.add_argument("restart_directory", type=Path)
    args = parser.parse_args()
    continuous = read_checkpoint(
        args.continuous_directory / "output_00002/sink_00002.out"
    )
    restarted = read_checkpoint(
        args.restart_directory / "output_00002/sink_00002.out"
    )
    continuous_statistics = read_sink_statistics(
        args.continuous_directory / "output_00002/sink_00002.out"
    )
    restarted_statistics = read_sink_statistics(
        args.restart_directory / "output_00002/sink_00002.out"
    )
    largest_absolute_difference = 0.0
    for sink_id, reference in continuous.items():
        for field_index, (expected, actual) in enumerate(
            zip(reference, restarted[sink_id], strict=True)
        ):
            difference = abs(expected - actual)
            largest_absolute_difference = max(largest_absolute_difference, difference)
            if not math.isclose(expected, actual, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(
                    f"restart state differs for sink {sink_id}, field {field_index}: "
                    f"continuous={expected}, restarted={actual}"
                )
        for field_index, (expected, actual) in enumerate(
            zip(continuous_statistics[sink_id], restarted_statistics[sink_id], strict=True)
        ):
            difference = abs(expected - actual)
            largest_absolute_difference = max(largest_absolute_difference, difference)
            if not math.isclose(expected, actual, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(
                    f"restart sink_stat differs for sink {sink_id}, channel {field_index}: "
                    f"continuous={expected}, restarted={actual}"
                )
    print(json.dumps({
        "status": "multiple_restart_passed",
        "retained_sink_ids": sorted(restarted),
        "compared_checkpoint": "output_00002/sink_00002.out",
        "compared_fields_per_sink": len(next(iter(continuous.values()))),
        "compared_sink_stat_channels_per_sink": len(next(iter(continuous_statistics.values()))),
        "largest_absolute_difference": largest_absolute_difference,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
