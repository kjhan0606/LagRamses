#!/usr/bin/env python3
"""Attribute hybrid FDM mass drift to instrumented solver transitions."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path


NUMBER = r"[+\-0-9.EDed]+"
MCHK_RE = re.compile(
    rf"MCHK:\s*(?P<label>\S+)\s+parent=\s*(?P<level>\d+)"
    rf"\s+total=\s*(?P<total>{NUMBER})"
)
CN_RE = re.compile(
    rf"FDM_CN_BUDGET\s+level=\s*(?P<level>\d+)"
    rf"\s+wave=(?P<wave>{NUMBER})"
    rf"\s+requested=(?P<requested>{NUMBER})"
    rf"\s+applied=(?P<applied>{NUMBER})"
    rf"\s+closure=(?P<closure>{NUMBER})"
    rf"\s+unapplied=(?P<unapplied>{NUMBER})"
)


def ffloat(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def analyze(case_dir: Path) -> dict[str, object]:
    log_path = case_dir / "run.log"
    result: dict[str, object] = {"case": case_dir.name, "log": str(log_path)}
    if not log_path.is_file():
        result.update({"completed": False, "error": "missing run.log"})
        return result

    text = log_path.read_text(errors="replace")
    mchk = [
        (match.group("label"), int(match.group("level")), ffloat(match.group("total")))
        for match in MCHK_RE.finditer(text)
    ]
    transitions: dict[str, dict[str, float | int]] = defaultdict(
        lambda: {"count": 0, "sum": 0.0, "sum_abs": 0.0, "max_abs": 0.0}
    )
    for previous, current in zip(mchk, mchk[1:]):
        key = f"{previous[0]}(L{previous[1]})->{current[0]}(L{current[1]})"
        delta = current[2] - previous[2]
        item = transitions[key]
        item["count"] = int(item["count"]) + 1
        item["sum"] = float(item["sum"]) + delta
        item["sum_abs"] = float(item["sum_abs"]) + abs(delta)
        item["max_abs"] = max(float(item["max_abs"]), abs(delta))

    cn_rows = [
        {
            "level": int(match.group("level")),
            **{
                field: ffloat(match.group(field))
                for field in ("wave", "requested", "applied", "closure", "unapplied")
            },
        }
        for match in CN_RE.finditer(text)
    ]
    cn_summary = {
        field: sum(float(row[field]) for row in cn_rows)
        for field in ("wave", "requested", "applied", "closure", "unapplied")
    }
    initial = mchk[0][2] if mchk else math.nan
    final = mchk[-1][2] if mchk else math.nan
    result.update(
        {
            "completed": "Run completed" in text,
            "mchk_count": len(mchk),
            "mass_initial": initial,
            "mass_final": final,
            "mass_drift": final - initial,
            "mass_drift_rel": (final - initial) / initial if initial else math.nan,
            "transitions": dict(
                sorted(
                    transitions.items(),
                    key=lambda item: abs(float(item[1]["sum"])),
                    reverse=True,
                )
            ),
            "cn_budget_count": len(cn_rows),
            "cn_budget": cn_summary,
            "fatal": bool(re.search(r"FATAL|No more free|NaN_CHK.*=[ ]*[1-9]", text)),
        }
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    cases = [
        analyze(path)
        for path in sorted(args.run_root.iterdir())
        if path.is_dir() and (path / "cosmo.nml").is_file()
    ]
    report = {"run_root": str(args.run_root), "cases": cases}

    for case in cases:
        print(
            f"{case['case']}: completed={case.get('completed')} "
            f"drift_rel={float(case.get('mass_drift_rel', math.nan)):+.6e} "
            f"CN_closure={float(dict(case.get('cn_budget', {})).get('closure', math.nan)):+.6e}"
        )
        transitions = dict(case.get("transitions", {}))
        for name, values in list(transitions.items())[:8]:
            print(
                f"  {name:48s} sum={float(values['sum']):+.6e} "
                f"max={float(values['max_abs']):.3e} n={int(values['count'])}"
            )

    if args.output:
        args.output.write_text(json.dumps(report, indent=2, allow_nan=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
