#!/usr/bin/env python3
"""Summarize opt-in QUMOND source traces and compare the saved force field."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

from analyze_qumond_spherical_reference import (
    A0_CGS,
    UNITS_LENGTH_CM,
    UNITS_TIME_S,
    discrete_reference,
    error_statistics,
    read_uniform_snapshot,
    source_density,
)


NUM = r"([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)"


def parse_number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_four_values(line: str) -> dict[str, float]:
    values = re.findall(NUM, line.rsplit("=", 1)[1])
    if len(values) != 4:
        raise ValueError(f"expected four diagnostic values: {line}")
    return dict(zip(("mean", "rms", "minimum", "maximum"), map(parse_number, values)))


def parse_log(log: str) -> tuple[list[dict], list[dict], list[dict]]:
    source_records: list[dict] = []
    post_p2_records: list[dict] = []
    fft_records: list[dict] = []
    for index, line in enumerate(log.splitlines()):
        if "[QUMOND_P2_SOURCE_DIAG] stage=post_second_poisson step=" in line:
            match = re.search(
                r"step=\s*(\d+) level=\s*(\d+) active_cells=\s*(\d+)", line
            )
            if not match:
                raise ValueError(f"malformed post-Poisson source header: {line}")
            step, level, cells = match.groups()
            post_p2_records.append({
                "line": index,
                "step": int(step),
                "level": int(level),
                "active_cells": int(cells),
                "fields": {},
            })
        elif "[QUMOND_P2_SOURCE_DIAG]" in line and " mean rms min max=" in line:
            if not post_p2_records:
                raise ValueError(f"post-Poisson values appear before header: {line}")
            name = line.split("[QUMOND_P2_SOURCE_DIAG]", 1)[1].split(" mean rms", 1)[0].strip()
            post_p2_records[-1]["fields"][name] = parse_four_values(line)
        elif "[QUMOND_SOURCE_DIAG] step=" in line:
            match = re.search(
                rf"step=\s*(\d+) level=\s*(\d+) a0_code=\s*{NUM} fourpi=\s*{NUM}",
                line,
            )
            if not match:
                raise ValueError(f"malformed source diagnostic: {line}")
            step, level, a0, fourpi = match.groups()
            source_records.append({
                "line": index,
                "step": int(step),
                "level": int(level),
                "a0_code": parse_number(a0),
                "fourpi": parse_number(fourpi),
                "fields": {},
                "samples": [],
            })
        elif "[QUMOND_SOURCE_DIAG] cells updated skipped missing_faces=" in line:
            values = [int(value) for value in re.findall(r"\d+", line.split("=", 1)[1])]
            if len(values) != 4 or not source_records:
                raise ValueError(f"malformed source cell counts: {line}")
            source_records[-1].update(dict(zip(
                ("cells", "updated_cells", "skipped_cells", "missing_face_links"), values
            )))
        elif "[QUMOND_SOURCE_DIAG]" in line and " mean rms min max=" in line:
            if not source_records:
                raise ValueError(f"source values appear before header: {line}")
            name = line.split("[QUMOND_SOURCE_DIAG]", 1)[1].split(" mean rms", 1)[0].strip()
            source_records[-1]["fields"][name] = parse_four_values(line)
        elif "[QUMOND_SOURCE_DIAG] sample=" in line:
            source_records[-1]["samples"].append(line.strip())

        if "[QUMOND_FFT_DIAG] call=" in line and " cells=" in line:
            match = re.search(
                rf"call=\s*(\d+) step=\s*(\d+) level=\s*(\d+) cells=\s*(\d+) rho_tot=\s*{NUM}",
                line,
            )
            if not match:
                raise ValueError(f"malformed FFT source header: {line}")
            call, step, level, cells, rho_tot = match.groups()
            fft_records.append({
                "line": index,
                "call": int(call),
                "step": int(step),
                "level": int(level),
                "cells": int(cells),
                "rho_tot": parse_number(rho_tot),
            })
        elif "[QUMOND_FFT_DIAG] input_source mean rms min max=" in line:
            if not fft_records:
                raise ValueError(f"FFT source values appear before header: {line}")
            fft_records[-1]["input_source"] = parse_four_values(line)
        elif "[QUMOND_FFT_DIAG] call=" in line and " packed_rhs mean rms min max=" in line:
            if not fft_records:
                raise ValueError(f"FFT RHS values appear before header: {line}")
            fft_records[-1]["packed_rhs"] = parse_four_values(line)

    for source in source_records:
        post_p2 = next((
            item for item in post_p2_records
            if item["line"] > source["line"]
            and item["step"] == source["step"]
            and item["level"] == source["level"]
            and "matched_source" not in item
        ), None)
        if post_p2 is None:
            raise ValueError(
                f"source step {source['step']} level {source['level']} "
                "has no following post-Poisson diagnostic"
            )
        expected_post_fields = {
            "newtonian_pre_source", "phantom_residual_source", "combined_post_source"
        }
        if set(post_p2["fields"]) != expected_post_fields:
            raise ValueError(f"incomplete post-Poisson source statistics: {post_p2}")
        post_p2["matched_source"] = True
        source["post_second_poisson"] = {
            "active_cells": post_p2["active_cells"],
            "fields": post_p2["fields"],
        }
        post_to_pre = {
            "newtonian_pre_source": "newtonian_source",
            "phantom_residual_source": "phantom_source",
            "combined_post_source": "combined_source",
        }
        source["post_minus_pre"] = {
            post_name: {
                key: post_p2["fields"][post_name][key] - source["fields"][pre_name][key]
                for key in ("mean", "rms", "minimum", "maximum")
            }
            for post_name, pre_name in post_to_pre.items()
        }
        if not fft_records:
            continue
        next_fft = next((item for item in fft_records if item["line"] > source["line"]), None)
        if next_fft is None:
            raise ValueError("source diagnostic has no following FFT solve")
        source["following_fft_call"] = next_fft["call"]
        source["combined_to_fft_input"] = {
            key: next_fft["input_source"][key] - source["fields"]["combined_source"][key]
            for key in ("mean", "rms", "minimum", "maximum")
        }
    return source_records, post_p2_records, fft_records


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: analyze_qumond_source_diagnostic.py RUN_DIRECTORY ORIGINAL_QUMOND_SNAPSHOT")
    root = Path(sys.argv[1]).resolve()
    original = Path(sys.argv[2]).resolve()
    snapshots = sorted(root.glob("output_[0-9][0-9][0-9][0-9][0-9]/data_*.h5"))
    if len(snapshots) != 1:
        raise ValueError(f"expected one force snapshot, found {len(snapshots)}")
    actual, positions, boxlen, nmesh = read_uniform_snapshot(snapshots[0])
    baseline, base_positions, base_boxlen, base_nmesh = read_uniform_snapshot(original)
    if nmesh != base_nmesh or boxlen != base_boxlen or not np.array_equal(positions, base_positions):
        raise ValueError("diagnostic and original QUMOND output meshes differ")
    force_mag = np.sqrt(np.sum(baseline**2, axis=0))
    displacement = (positions - 0.5 * boxlen + 0.5 * boxlen) % boxlen - 0.5 * boxlen
    radius = np.sqrt(np.sum(displacement**2, axis=-1))
    mask = (radius > 0.06) & (radius < 0.35) & (force_mag > 1.0e-4 * np.max(force_mag))
    relative = np.sqrt(np.sum((actual - baseline) ** 2, axis=0))[mask] / force_mag[mask]
    density = source_density(radius)
    a0_code = A0_CGS / (UNITS_LENGTH_CM / UNITS_TIME_S**2)
    _, independent_qumond = discrete_reference(density, boxlen / nmesh, a0_code)
    independent_accuracy = error_statistics(actual, independent_qumond, radius)
    log = (root / "run.log").read_text()
    source_records, post_p2_records, fft_records = parse_log(log)
    if not source_records:
        raise ValueError("the opt-in QUMOND source diagnostic is absent")
    if not post_p2_records:
        raise ValueError("the post-second-Poisson source diagnostic is absent")
    result = {
        "status": "ANALYZED",
        "scope": (
            "instrumented pre/post second-Poisson QUMOND source trace on a fixed "
            "spherical density; post-solve phantom source is combined active rho "
            "minus the saved pre-solve Newtonian source; no galaxy RAR result"
        ),
        "snapshot": str(snapshots[0]),
        "baseline_snapshot": str(original),
        "source_records": source_records,
        "post_second_poisson_records": post_p2_records,
        "fft_records": fft_records,
        "comparison_shell_code_radius": [0.06, 0.35],
        "force_reproduction": {
            "cells": int(relative.size),
            "median_relative_vector_difference": float(np.median(relative)),
            "p90_relative_vector_difference": float(np.quantile(relative, 0.90)),
            "max_relative_vector_difference": float(np.max(relative)),
        },
        "independent_discrete_qumond_reference": independent_accuracy,
    }
    (root / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
