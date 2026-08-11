#!/usr/bin/env python3
"""Evaluate the low-cost FDM fluid--wave interface regression matrix."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path


MCHK_RE = re.compile(
    r"MCHK:\s*(?P<label>\S+)\s+parent=\s*(?P<level>\d+)"
    r"\s+total=\s*(?P<total>[+\-0-9.EDed]+)"
)
LEVEL_RE = re.compile(r"Level\s+(?P<level>\d+)\s+has\s+(?P<count>\d+)\s+grids")
MTOT_RE = re.compile(r"FDM:\s+M_tot=\s*(?P<total>[+\-0-9.EDed]+)")
ELAPSED_RE = re.compile(r"Total elapsed time:\s*(?P<seconds>[+\-0-9.EDed]+)")
STEP_DT_RE = re.compile(
    r"Fine step=\s*\d+.*?\bdt=\s*(?P<dt>[+\-0-9.EDed]+)"
)
FATAL_RE = re.compile(
    r"No more free|FATAL|NaN_CHK.*=\s*[1-9]|mismatch=\s*[1-9]",
    re.IGNORECASE,
)


def fortran_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def nml_scalar(text: str, name: str, default: str) -> str:
    match = re.search(
        rf"(?im)^\s*{re.escape(name)}\s*=\s*([^!\n/]+)",
        text,
    )
    return match.group(1).strip() if match else default


def json_safe(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def analyze_case(case_dir: Path, abs_tol: float) -> dict[str, object]:
    log_path = case_dir / "run.log"
    nml_path = case_dir / "cosmo.nml"
    result: dict[str, object] = {
        "case": case_dir.name,
        "log": str(log_path),
        "pass": False,
        "reasons": [],
    }
    reasons: list[str] = result["reasons"]  # type: ignore[assignment]

    if not log_path.is_file():
        reasons.append("missing run.log")
        return result
    if not nml_path.is_file():
        reasons.append("missing cosmo.nml")
        return result

    log = log_path.read_text(errors="replace")
    nml = nml_path.read_text(errors="replace")
    levelmax = int(nml_scalar(nml, "levelmax", "0"))
    first_wave = int(nml_scalar(nml, "fdm_first_wave_level", "0"))
    use_hjm = nml_scalar(nml, "fdm_use_hjm", ".false.").lower() == ".true."
    courant = fortran_float(nml_scalar(nml, "fdm_courant", "nan"))

    result.update(
        {
            "levelmax": levelmax,
            "first_wave_level": first_wave,
            "use_hjm": use_hjm,
            "fdm_courant": courant,
            "completed": "Run completed" in log,
        }
    )

    if not result["completed"]:
        reasons.append("run did not complete")
    fatal_matches = sorted(set(match.group(0) for match in FATAL_RE.finditer(log)))
    result["fatal_diagnostics"] = fatal_matches
    if fatal_matches:
        reasons.append("fatal diagnostic in log")

    active: dict[int, int] = {}
    for match in LEVEL_RE.finditer(log):
        level = int(match.group("level"))
        active[level] = max(active.get(level, 0), int(match.group("count")))
    result["max_grids_by_level"] = {str(key): active[key] for key in sorted(active)}
    if active.get(levelmax, 0) <= 0:
        reasons.append(f"level {levelmax} was not active")

    entries = [
        {
            "label": match.group("label"),
            "level": int(match.group("level")),
            "total": fortran_float(match.group("total")),
        }
        for match in MCHK_RE.finditer(log)
    ]
    result["mchk_count"] = len(entries)
    if not entries:
        reasons.append("no MCHK diagnostics")
        result["mchk_span_abs"] = math.nan
        result["mchk_span_rel"] = math.nan
    else:
        totals = [float(entry["total"]) for entry in entries]
        span_abs = max(totals) - min(totals)
        scale = abs(totals[0])
        span_rel = span_abs / scale if scale else math.inf
        result["mchk_first"] = totals[0]
        result["mchk_last"] = totals[-1]
        result["mchk_span_abs"] = span_abs
        result["mchk_span_rel"] = span_rel
        if span_abs > abs_tol:
            reasons.append(f"MCHK span {span_abs:.3e} exceeds {abs_tol:.3e}")

    pending_pre: dict[int, float] = {}
    refine_deltas: list[float] = []
    for entry in entries:
        label = str(entry["label"])
        level = int(entry["level"])
        total = float(entry["total"])
        if label == "pre-refine":
            pending_pre[level] = total
        elif label == "post-refine" and level in pending_pre:
            refine_deltas.append(abs(total - pending_pre.pop(level)))
    max_refine_delta = max(refine_deltas, default=math.nan)
    result["refine_pair_count"] = len(refine_deltas)
    result["max_pre_post_refine_abs"] = max_refine_delta
    if not refine_deltas:
        reasons.append("no paired pre/post-refine diagnostics")
    elif max_refine_delta > abs_tol:
        reasons.append(
            f"pre/post-refine delta {max_refine_delta:.3e} exceeds {abs_tol:.3e}"
        )

    mtot = [fortran_float(match.group("total")) for match in MTOT_RE.finditer(log)]
    result["mtot_count"] = len(mtot)
    result["mtot_first"] = mtot[0] if mtot else math.nan
    result["mtot_last"] = mtot[-1] if mtot else math.nan
    result["mtot_drift_abs"] = abs(mtot[-1] - mtot[0]) if mtot else math.nan

    elapsed = [fortran_float(match.group("seconds")) for match in ELAPSED_RE.finditer(log)]
    result["elapsed_seconds"] = elapsed[-1] if elapsed else math.nan
    step_dt = [fortran_float(match.group("dt")) for match in STEP_DT_RE.finditer(log)]
    result["fine_step_count"] = len(step_dt)
    result["last_fine_dt"] = step_dt[-1] if step_dt else math.nan
    result["pass"] = not reasons
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--abs-tol", type=float, default=2.0e-10)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    case_dirs = sorted(
        path
        for path in args.run_root.iterdir()
        if path.is_dir() and (path / "cosmo.nml").is_file()
    )
    results = [analyze_case(path, args.abs_tol) for path in case_dirs]
    by_name = {str(result["case"]): result for result in results}
    timestep_pairs = []
    for family in ("wave_l6", "hybrid_l7", "hybrid_l8"):
        coarse = by_name.get(f"{family}_dt1")
        fine = by_name.get(f"{family}_dt2")
        if coarse is None or fine is None:
            continue
        dt1 = float(coarse.get("last_fine_dt", math.nan))
        dt2 = float(fine.get("last_fine_dt", math.nan))
        ratio = dt2 / dt1 if dt1 > 0.0 and math.isfinite(dt2) else math.nan
        pair_pass = math.isfinite(ratio) and ratio <= 0.75
        timestep_pairs.append(
            {
                "family": family,
                "dt1": dt1,
                "dt2": dt2,
                "ratio": ratio,
                "pass": pair_pass,
            }
        )
        if not pair_pass:
            fine["reasons"].append(
                f"Courant-factor pair did not reduce the last fine dt: ratio={ratio:.3g}"
            )
            fine["pass"] = False

    report = {
        "run_root": str(args.run_root),
        "absolute_mass_tolerance": args.abs_tol,
        "pass": bool(results) and all(result["pass"] for result in results),
        "cases": results,
        "timestep_pairs": timestep_pairs,
    }

    print(
        f"{'case':20s} {'HJM':>5s} {'Lwave':>5s} {'CFL':>6s} "
        f"{'MCHK span':>12s} {'refine delta':>13s} {'status':>7s}"
    )
    for result in results:
        span = float(result.get("mchk_span_abs", math.nan))
        refine = float(result.get("max_pre_post_refine_abs", math.nan))
        print(
            f"{str(result['case']):20s} "
            f"{str(result.get('use_hjm', '?')):>5s} "
            f"{int(result.get('first_wave_level', 0)):5d} "
            f"{float(result.get('fdm_courant', math.nan)):6.3f} "
            f"{span:12.3e} {refine:13.3e} "
            f"{'PASS' if result['pass'] else 'FAIL':>7s}"
        )
        for reason in result["reasons"]:
            print(f"  {result['case']}: {reason}")

    encoded = json.dumps(json_safe(report), indent=2, sort_keys=True, allow_nan=False)
    if args.output:
        args.output.write_text(encoded + "\n")
    else:
        print(encoded)
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
