#!/usr/bin/env python3
"""Compare the first coarse-entry owned-cloud and sink-stat census after restart."""

import json
import re
import sys


PATTERN = re.compile(
    r"CAPTURE_RESTART_PROBE pre_output step=2 rank=1 "
    r"local_owned=(\d+) global_owned=(\d+) "
    r"global_sumx=\s*([+-]?[\d.Ee+-]+) "
    r"global_sink_stat=\s*([+-]?[\d.Ee+-]+)"
)


def read_census(path):
    with open(path, encoding="utf-8") as handle:
        matches = PATTERN.findall(handle.read())
    if len(matches) != 1:
        raise ValueError(f"expected exactly one step-2 rank-1 census in {path}, got {len(matches)}")
    local, total, sumx, stat = matches[0]
    return {
        "local_owned": int(local),
        "global_owned": int(total),
        "global_sumx": float(sumx),
        "global_sink_stat": float(stat),
    }


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: evaluate_restart_ownership.py FRESH.LOG RESTART.LOG")
    fresh = read_census(sys.argv[1])
    restarted = read_census(sys.argv[2])
    report = {
        "fresh_step2_pre_output": fresh,
        "restarted_step2_pre_output": restarted,
        "owned_cloud_count_equal": fresh["global_owned"] == restarted["global_owned"],
        "sink_stat_equal": fresh["global_sink_stat"] == restarted["global_sink_stat"],
        "interpretation": "Diagnostic ownership gate only; matching counts would not prove trajectory equivalence.",
    }
    print(json.dumps(report, indent=2, sort_keys=True))
