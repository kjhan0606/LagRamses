#!/usr/bin/env python3
"""Summarize opt-in SIDM local and AMR-parent timestep-selection records.

These decisions precede child-level recursion and do not determine the final
coarse timestep retroactively. Counts are diagnostic calls, not independent
physical steps or a walltime speedup forecast. Non-SIDM constraints are grouped.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from analyze_mg_convergence_log import analyze_lines as analyze_mg_lines


def record(line: str, parent: bool) -> dict:
    fields = dict(re.findall(r"(\w+)=\s*(\S+)", line))
    integers = ("coarse", "fine", "level")
    booleans = ("bound_active", "parent_limited" if parent else "sidm_limited")
    reals = (("a", "dt_local", "dt_parent", "dt_selected") if parent else
             ("a", "dt_other", "dt_old", "pmax", "courant", "dt_sidm", "dt_local"))
    if set(fields) != set(integers + booleans + reals):
        raise ValueError("missing or unexpected timestep diagnostic fields")
    result = {key: int(fields[key]) for key in integers}
    if result["level"] <= 0 or min(result["coarse"], result["fine"]) < 0:
        raise ValueError("invalid timestep diagnostic indices")
    for key in booleans:
        if fields[key] not in ("T", "F"):
            raise ValueError("invalid logical diagnostic field")
        result[key] = fields[key] == "T"
    for key in reals:
        value = float(fields[key].replace("D", "E").replace("d", "e"))
        if not math.isfinite(value) or value < 0:
            raise ValueError("nonfinite or negative timestep diagnostic value")
        result[key] = value
    if result["a"] <= 0 or result["dt_local"] <= 0:
        raise ValueError("nonpositive scale factor or local timestep")
    return result


def near(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=2e-12, abs_tol=0)


def check_local(item: dict) -> None:
    before = item["dt_other"]
    if before <= 0 or item["courant"] <= 0:
        raise ValueError("nonpositive other-constraint timestep or SIDM courant")
    if item["bound_active"] != (item["pmax"] > 0):
        raise ValueError("SIDM bound activation disagrees with Pmax")
    candidate = item["dt_sidm"]
    if item["bound_active"]:
        if candidate <= 0 or not near(candidate, item["courant"] / item["pmax"] * item["dt_old"]):
            raise ValueError("SIDM candidate disagrees with courant/Pmax * previous dt")
        expected = min(before, candidate)
    else:
        if candidate != 0:
            raise ValueError("inactive SIDM candidate must use the zero sentinel")
        expected = before
    if not near(expected, item["dt_local"]):
        raise ValueError("local timestep disagrees with SIDM minimum")
    if not near(expected, before) and item["sidm_limited"] != (expected < before):
        raise ValueError("SIDM limiting flag disagrees with timestep")


def check_parent(local: dict, parent: dict) -> None:
    for key in ("coarse", "fine", "level"):
        if local[key] != parent[key]:
            raise ValueError("local/parent diagnostic indices do not match")
    if not near(local["a"], parent["a"]) or not near(local["dt_local"], parent["dt_local"]):
        raise ValueError("local/parent epoch or local timestep does not match")
    candidate = parent["dt_parent"]
    if parent["bound_active"]:
        if candidate <= 0:
            raise ValueError("active parent bound is nonpositive")
        expected = min(local["dt_local"], candidate)
    else:
        if candidate != 0:
            raise ValueError("inactive parent bound must use the zero sentinel")
        expected = local["dt_local"]
    if parent["dt_selected"] <= 0 or not near(expected, parent["dt_selected"]):
        raise ValueError("selected timestep disagrees with parent minimum")
    if not near(expected, local["dt_local"]) and parent["parent_limited"] != (expected < local["dt_local"]):
        raise ValueError("parent limiting flag disagrees with timestep")


def paired_records(lines) -> list:
    pending = None
    pairs = []
    for line in lines:
        if "SIDM_DT_PARENT:" in line:
            if pending is None:
                raise ValueError("parent diagnostic without a preceding local record")
            parent = record(line, True)
            check_parent(pending, parent)
            pairs.append((pending, parent))
            pending = None
        elif "SIDM_DT:" in line:
            if pending is not None:
                raise ValueError("local diagnostic without matching parent record")
            pending = record(line, False)
            check_local(pending)
    if pending is not None:
        raise ValueError("final local diagnostic has no matching parent record")
    if not pairs:
        raise ValueError("no paired SIDM timestep diagnostic evidence")
    return pairs


def summarize_pairs(pairs) -> dict:
    grouped = defaultdict(list)
    for local, parent in pairs:
        grouped[local["level"]].append((local, parent))
    levels = []
    for level, pairs in sorted(grouped.items()):
        levels.append({
            "level": level, "diagnostic_calls": len(pairs),
            "sidm_local_limited_calls": sum(local["sidm_limited"] for local, _ in pairs),
            "parent_further_limited_calls": sum(parent["parent_limited"] for _, parent in pairs),
            "minimum_local_over_other_ratio": min(local["dt_local"] / local["dt_other"] for local, _ in pairs),
            "minimum_selected_over_other_ratio": min(parent["dt_selected"] / local["dt_other"] for local, parent in pairs),
            "maximum_previous_Pmax": max(local["pmax"] for local, _ in pairs),
            "aexp_min": min(local["a"] for local, _ in pairs),
            "aexp_max": max(local["a"] for local, _ in pairs),
        })
    return {
        "paired_diagnostic_calls": sum(item["diagnostic_calls"] for item in levels),
        "levels": levels,
        "scope": "per-level selection before child recursion; not final coarse dt or independent samples",
        "walltime_speedup_inferred": False,
    }


def analyze_lines(lines) -> dict:
    return summarize_pairs(paired_records(lines))


def analyze_pilot_lines(lines, seed_coarse: int, seed_aexp: float,
                        base_level: int, expected_steps: int = 20,
                        warmup_steps: int = 2) -> dict:
    if not (seed_coarse >= 0 and 0 < seed_aexp < 1 and base_level > 0 and
            0 <= warmup_steps < expected_steps - 2):
        raise ValueError("invalid pilot seed, step count, or warmup")
    text = "\n".join(lines)
    number_pattern = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    restarts = re.findall(r"Restarting at t=.*?nstep_coarse=\s*(\d+)", text)
    mains = [int(x) for x in re.findall(r"Main step=\s*(\d+)", text)]
    timings = [float(x.replace("D", "E")) for x in re.findall(
        rf"Time elapsed since last coarse step:\s*({number_pattern})\s+s", text)]
    elapsed = re.findall(rf"Total elapsed time:\s*({number_pattern})", text)
    if restarts != [str(seed_coarse)] or not mains or mains[-1] != seed_coarse + expected_steps:
        raise ValueError("pilot restart or stopping step does not match requested interval")
    if text.count("Run completed") != 1 or len(elapsed) != 1 or len(timings) != expected_steps:
        raise ValueError("pilot completion/timing evidence is missing or incomplete; ncontrol=1 required")
    if any(not math.isfinite(x) or x <= 0 for x in timings):
        raise ValueError("nonpositive or nonfinite coarse walltime")
    solver = analyze_mg_lines(lines, 1e-4, 80)
    if not solver["convergence_gate_pass"] or not solver["fftw_solver_verified"]:
        raise ValueError("pilot has not passed strict MG/actual FFT solve gates")
    pairs = paired_records(lines)
    base = {}
    for local, _ in pairs:
        if local["level"] == base_level:
            step = local["coarse"]
            if step in base:
                raise ValueError("duplicate base-level timestep decision for a coarse step")
            base[step] = local["a"]
    if seed_coarse not in base or not near(base[seed_coarse], seed_aexp):
        raise ValueError("base-level diagnostic epoch does not match the seed")
    start = seed_coarse + warmup_steps
    end = max(base)
    # The final stop can occur inside update_time. Use only exact epochs from
    # base-level diagnostics and their matching completed coarse-step times,
    # never the rounded four-significant-digit Fine-step epoch.
    if not (start in base and start + 2 <= end <= mains[-1]) or set(base) != set(range(seed_coarse, end + 1)):
        raise ValueError("insufficient or discontinuous precise base-level epoch evidence")
    a_start, a_end = base[start], base[end]
    if not 0 < a_start < a_end < 1:
        raise ValueError("pilot interval must precede z=0")
    if any(not base[step] < base[step + 1] for step in range(seed_coarse, end)):
        raise ValueError("base-level scale factors are not strictly increasing")
    seconds = sum(timings[step - seed_coarse - 1] for step in range(start + 1, end + 1))
    total_elapsed = float(elapsed[0].replace("D", "E"))
    if not math.isfinite(total_elapsed) or total_elapsed <= 0 or total_elapsed < seconds - 0.01 * expected_steps:
        raise ValueError("total elapsed time is inconsistent with measured coarse-step times")
    seconds_per_a = seconds / (a_end - a_start)
    coarse_records = [pair for pair in pairs if start <= pair[0]["coarse"] < end]
    remaining = seconds_per_a * (1 - a_end) / 86400
    return {
        "completed_coarse_steps": expected_steps, "warmup_coarse_steps_excluded": warmup_steps,
        "measured_start_coarse_step": start, "measured_end_coarse_step": end,
        "measured_coarse_steps": end - start, "a_start": a_start, "a_end": a_end,
        "measured_coarse_wall_seconds": seconds,
        "mean_seconds_per_measured_coarse_step": seconds / (end - start),
        "total_logged_elapsed_seconds": total_elapsed,
        "seconds_per_delta_a": seconds_per_a,
        "z0_compute_days_at_fixed_local_rate": remaining,
        "speedup_for_90_days_at_full_allocation": remaining / 90,
        "solver_gate": solver, "steady_interval_diagnostics": summarize_pairs(coarse_records),
        "accepted_as_science_trajectory": False,
        "limitations": [
            "A short pilot is an operational local-rate measurement, not a stationary-rate guarantee.",
            "Legacy seed trajectory remains uncertified; passing reconstruction does not repair its history.",
            "Diagnostic logging and 0.01-second printed timing precision contribute measurement overhead/error.",
            "Final control work, initialization, queue delay, future AMR growth and further restarts are not in the measured interval.",
            "No CDM matched-epoch speedup or scientific convergence is inferred from this single-model pilot.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--pilot-seed-info", type=Path)
    parser.add_argument("--pilot-steps", type=int, default=20)
    parser.add_argument("--warmup-steps", type=int, default=2)
    args = parser.parse_args()
    lines = args.log.read_text().splitlines()
    if args.pilot_seed_info:
        seed = args.pilot_seed_info.read_text()
        def info_value(key):
            found = re.findall(rf"^{key}\s*=\s*(\S+)", seed, re.M)
            if len(found) != 1:
                raise ValueError(f"missing or ambiguous seed {key}")
            return found[0]
        result = analyze_pilot_lines(lines, int(info_value("nstep_coarse")),
                                    float(info_value("aexp").replace("D", "E")),
                                    int(info_value("levelmin")), args.pilot_steps, args.warmup_steps)
        result["seed_info"] = str(args.pilot_seed_info)
    else:
        result = analyze_lines(lines)
    result["log"] = str(args.log)
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    print(text, end="")
    if args.json:
        if args.json.exists():
            parser.error("refusing to overwrite a previous diagnostic analysis")
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text)


if __name__ == "__main__":
    main()
