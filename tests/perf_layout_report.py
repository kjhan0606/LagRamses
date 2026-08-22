#!/usr/bin/env python3
"""Summarize the fixed ABBA layout performance characterization."""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
from pathlib import Path


TOTAL_RE = re.compile(r"Total elapsed time:\s*([0-9]+(?:\.[0-9]*)?(?:[Ee][+-]?\d+)?)")
TIMER_RE = re.compile(
    r"^\s*([0-9]+\.\d+)\s+([0-9]+\.\d+)\s+([0-9]+\.\d+)\s+"
    r"([0-9]+\.\d+)\s+([0-9]+\.\d+)\s+([0-9]+\.\d+)\s+"
    r"(\d+)\s+(\d+)\s+(.+?)\s*$"
)
TIMER_TOTAL_RE = re.compile(r"^\s*([0-9]+\.\d+)\s+100\.0\s+TOTAL\s*$")


def parse_timer_report(text: str, source: object) -> tuple[float, dict[str, float]]:
    marker = "PERF_TIMER_REPORT step=8 interval=8"
    if text.count(marker) != 1:
        raise ValueError(f"{source}: expected exactly one step-8 timer report")
    timer_rows: dict[str, float] = {}
    timer_total = None
    in_step8_report = False
    for line in text.splitlines():
        if marker in line:
            in_step8_report = True
            continue
        if not in_step8_report:
            continue
        match = TIMER_RE.match(line)
        if match:
            # The third field is the maximum time over ranks and therefore
            # determines parallel wall time for this phase.
            timer_rows[match.group(9).strip()] = float(match.group(3))
        total_match = TIMER_TOTAL_RE.match(line)
        if total_match:
            timer_total = float(total_match.group(1))
            in_step8_report = False
            break
    if timer_total is None or timer_total <= 0 or not timer_rows:
        raise ValueError(f"{source}: missing complete step-8 timer report")
    for required in ("particles", "poisson - mg base", "hydro - godunov"):
        if timer_rows.get(required, 0.0) <= 0:
            raise ValueError(f"{source}: missing positive timer phase {required!r}")
    return timer_total, timer_rows


def parse_run(run_dir: Path) -> dict[str, object]:
    text = (run_dir / "run.log").read_text(errors="replace")
    totals = [float(value) for value in TOTAL_RE.findall(text)]
    if len(totals) != 1 or not math.isfinite(totals[0]) or totals[0] <= 0:
        raise ValueError(f"{run_dir}: expected one positive Total elapsed time")
    timer_total, timer_rows = parse_timer_report(text, run_dir)

    external_ns = int((run_dir / "external_ns.txt").read_text().strip())
    if external_ns <= 0:
        raise ValueError(f"{run_dir}: invalid external wall time")
    return {
        "program_seconds": totals[0],
        "external_seconds": external_ns / 1.0e9,
        "timer_total_seconds": timer_total,
        "timer_max_seconds": timer_rows,
    }


def summarize(values: list[float]) -> dict[str, float]:
    mean = statistics.fmean(values)
    spread = abs(values[0] - values[1]) / mean * 100.0
    return {"samples": values, "mean": mean, "spread_percent": spread}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("job_root", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()

    names = ("block_1", "legacy_1", "legacy_2", "block_2")
    runs = {name: parse_run(args.job_root / "runs" / name) for name in names}
    grouped = {
        "block": [runs["block_1"], runs["block_2"]],
        "legacy": [runs["legacy_1"], runs["legacy_2"]],
    }
    metrics: dict[str, object] = {}
    for metric in ("program_seconds", "external_seconds", "timer_total_seconds"):
        metrics[metric] = {
            variant: summarize([float(run[metric]) for run in variant_runs])
            for variant, variant_runs in grouped.items()
        }

    program = metrics["program_seconds"]
    block_mean = program["block"]["mean"]
    legacy_mean = program["legacy"]["mean"]
    ratio = block_mean / legacy_mean
    co_tenant_snapshots = [
        (args.job_root / "runs" / name / "co_tenants.txt").read_text().strip()
        for name in names
    ]
    co_tenant_sets = set(co_tenant_snapshots)
    co_tenant_changed = any(
        "CO_TENANTS_CHANGED" in snapshot for snapshot in co_tenant_snapshots
    )
    noisy = (
        program["block"]["spread_percent"] > 5.0
        or program["legacy"]["spread_percent"] > 5.0
        or len(co_tenant_sets) != 1
        or co_tenant_changed
    )
    if noisy:
        assessment = "NOISY_INCONCLUSIVE"
    elif ratio > 1.05:
        assessment = "BLOCK_REGRESSION_OVER_5_PERCENT"
    elif ratio < 0.95:
        assessment = "BLOCK_FASTER_OVER_5_PERCENT"
    else:
        assessment = "WITHIN_5_PERCENT_INCONCLUSIVE"

    common_labels = sorted(
        set.intersection(
            *(set(run["timer_max_seconds"]) for run in runs.values())
        )
    )
    phases: dict[str, object] = {}
    for label in common_labels:
        block_values = [float(run["timer_max_seconds"][label]) for run in grouped["block"]]
        legacy_values = [float(run["timer_max_seconds"][label]) for run in grouped["legacy"]]
        block_summary = summarize(block_values)
        legacy_summary = summarize(legacy_values)
        phases[label] = {
            "block": block_summary,
            "legacy": legacy_summary,
            "block_over_legacy": block_summary["mean"] / legacy_summary["mean"]
            if legacy_summary["mean"] > 0
            else None,
        }

    report = {
        "schema": "lagramses-layout-perf-v1",
        "order": list(names),
        "metrics": metrics,
        "program_block_over_legacy": ratio,
        "observed_endpoint_co_tenant_sets_unchanged": (
            len(co_tenant_sets) == 1 and not co_tenant_changed
        ),
        "assessment": assessment,
        "phase_timers": phases,
    }
    args.json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    lines = [
        "# Block grid-major performance characterization",
        "",
        f"Assessment: `{assessment}`",
        "",
        "| metric | legacy samples (s) | block samples (s) | block/legacy |",
        "|---|---:|---:|---:|",
    ]
    for metric in ("program_seconds", "timer_total_seconds", "external_seconds"):
        left = metrics[metric]["legacy"]
        right = metrics[metric]["block"]
        lines.append(
            f"| {metric} | {left['samples'][0]:.6f}, {left['samples'][1]:.6f} "
            f"| {right['samples'][0]:.6f}, {right['samples'][1]:.6f} "
            f"| {right['mean'] / left['mean']:.6f} |"
        )
    lines.extend(
        [
            "",
            "The primary metric is the program-reported elapsed time. External wall time",
            "also includes initial output I/O. A per-variant spread over 5% or a co-tenant",
            "set change makes this shared-node result noisy/inconclusive.",
        ]
    )
    args.markdown.write_text("\n".join(lines) + "\n")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
