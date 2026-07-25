#!/usr/bin/env python3
"""Plot matched-phase DMO power-spectrum ratios and write a validation report."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_MODELS = ("q1", "cde10", "sym_a")
COLORS = {
    "q1": "#2878b5",
    "cde10": "#d95319",
    "f5": "#d62728",
    "f6": "#ff9896",
    "n1": "#1f77b4",
    "n5": "#9ecae1",
    "sym_a": "#6f4c9b",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--models", nargs="+", default=list(DEFAULT_MODELS))
    parser.add_argument("--redshifts", nargs="+", type=float, default=[2.0, 1.0, 0.5, 0.0])
    parser.add_argument("--kmax", type=float, default=0.8, help="maximum plotted k [h/Mpc]")
    parser.add_argument("--output-stem", default="pk_ratio_stageA")
    parser.add_argument(
        "--nearest",
        action="store_true",
        help="use the nearest dump instead of interpolating P(k) to each target redshift",
    )
    return parser.parse_args()


def read_pk(path: Path) -> dict:
    header = "\n".join(path.read_text().splitlines()[:7])
    match = re.search(r"a_exp\s*=\s*([0-9.Ee+-]+)", header)
    if not match:
        raise ValueError(f"missing a_exp header in {path}")
    aexp = float(match.group(1))
    data = np.loadtxt(path)
    return {
        "path": str(path),
        "aexp": aexp,
        "z": 1.0 / aexp - 1.0,
        "k": data[:, 0],
        "pk": data[:, 1],
        "nmodes": data[:, 4],
    }


def snapshots(model_dir: Path) -> list[dict]:
    return [
        read_pk(path)
        for path in sorted(model_dir.glob("output_*/pk_*.dat"))
    ]


def nearest(items: list[dict], redshift: float) -> dict:
    if not items:
        raise RuntimeError("no power-spectrum snapshots found")
    item = min(items, key=lambda value: abs(value["z"] - redshift)).copy()
    item["source_aexp"] = [item["aexp"]]
    item["interpolated"] = False
    return item


def at_redshift(items: list[dict], redshift: float, use_nearest: bool) -> dict:
    """Return P(k) at a common epoch, using log-P/log-a interpolation."""
    if use_nearest:
        return nearest(items, redshift)
    if not items:
        raise RuntimeError("no power-spectrum snapshots found")

    target_a = 1.0 / (1.0 + redshift)
    ordered = sorted(items, key=lambda item: item["aexp"])
    if target_a <= ordered[0]["aexp"] or target_a >= ordered[-1]["aexp"]:
        return nearest(ordered, redshift)

    lower, upper = next(
        (left, right)
        for left, right in zip(ordered[:-1], ordered[1:])
        if left["aexp"] <= target_a <= right["aexp"]
    )
    if not np.allclose(lower["k"], upper["k"], rtol=1.0e-7, atol=0.0):
        raise ValueError(
            f"incompatible k grids in {lower['path']} and {upper['path']}"
        )

    weight = (
        (np.log(target_a) - np.log(lower["aexp"]))
        / (np.log(upper["aexp"]) - np.log(lower["aexp"]))
    )
    positive = (lower["pk"] > 0.0) & (upper["pk"] > 0.0)
    pk = np.zeros_like(lower["pk"])
    pk[positive] = np.exp(
        (1.0 - weight) * np.log(lower["pk"][positive])
        + weight * np.log(upper["pk"][positive])
    )
    return {
        "path": f"{lower['path']}::{upper['path']}",
        "aexp": target_a,
        "z": redshift,
        "k": lower["k"],
        "pk": pk,
        "nmodes": np.minimum(lower["nmodes"], upper["nmodes"]),
        "source_aexp": [lower["aexp"], upper["aexp"]],
        "interpolated": True,
    }


def log_warnings(model_dir: Path) -> dict:
    logs = sorted(model_dir.glob("run-*.out"))
    text = "\n".join(path.read_text(errors="replace") for path in logs)
    return {
        "logs": [str(path) for path in logs],
        "nonconvergence_warnings": len(re.findall(r"NOT converged", text)),
        "fatal_markers": len(
            re.findall(r"\b(?:FATAL|SIGSEGV|segmentation fault|forrtl: severe)\b", text, re.I)
        ),
    }


def main() -> int:
    args = arguments()
    campaign = args.campaign.resolve()
    figure_dir = campaign / "figures"
    figure_dir.mkdir(exist_ok=True)

    model_names = ["lcdm", *args.models]
    all_snapshots = {name: snapshots(campaign / name) for name in model_names}
    missing = [name for name, values in all_snapshots.items() if not values]
    if missing:
        raise RuntimeError(f"missing P(k) outputs for: {', '.join(missing)}")

    ncols = 2
    nrows = (len(args.redshifts) + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(7.2, 2.65 * nrows), sharex=True, sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).ravel()
    rows = []

    for ax, target_z in zip(axes, args.redshifts):
        reference = at_redshift(all_snapshots["lcdm"], target_z, args.nearest)
        for name in args.models:
            sample = at_redshift(all_snapshots[name], target_z, args.nearest)
            good_ref = (reference["pk"] > 0.0) & (reference["nmodes"] > 0.0)
            good_mod = (sample["pk"] > 0.0) & (sample["nmodes"] > 0.0)
            good = good_ref & good_mod & (reference["k"] <= args.kmax)
            if not np.any(good):
                continue
            k = reference["k"][good]
            if np.array_equal(reference["k"], sample["k"]):
                ratio = sample["pk"][good] / reference["pk"][good]
            else:
                interp = np.interp(
                    np.log(k),
                    np.log(sample["k"][good_mod]),
                    np.log(sample["pk"][good_mod]),
                )
                ratio = np.exp(interp) / reference["pk"][good]
            ax.plot(k, ratio, lw=1.55, color=COLORS.get(name), label=name)
            for kval, value in zip(k, ratio):
                rows.append(
                    (
                        target_z,
                        reference["z"],
                        name,
                        sample["z"],
                        kval,
                        value,
                        ";".join(f"{value:.9g}" for value in reference["source_aexp"]),
                        ";".join(f"{value:.9g}" for value in sample["source_aexp"]),
                    )
                )

        ax.axhline(1.0, color="0.35", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_title(f"$z={target_z:.2f}$", fontsize=10)
        ax.grid(alpha=0.18)

    for ax in axes[len(args.redshifts):]:
        ax.set_visible(False)
    for ax in axes[-ncols:]:
        if ax.get_visible():
            ax.set_xlabel(r"$k\ [h\,{\rm Mpc}^{-1}]$")
    for index, ax in enumerate(axes):
        if ax.get_visible() and index % ncols == 0:
            ax.set_ylabel(r"$P(k)/P_{\Lambda{\rm CDM}}(k)$")
    axes[0].legend(frameon=False, fontsize=9, ncol=min(3, len(args.models)))

    png = figure_dir / f"{args.output_stem}.png"
    pdf = figure_dir / f"{args.output_stem}.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)

    csv_path = figure_dir / f"{args.output_stem}.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "target_z",
                "lcdm_z",
                "model",
                "model_z",
                "k_h_mpc",
                "pk_ratio",
                "lcdm_source_aexp",
                "model_source_aexp",
            )
        )
        writer.writerows(rows)

    report = {
        "campaign": str(campaign),
        "models": {
            name: {
                "snapshot_count": len(all_snapshots[name]),
                "redshifts": [item["z"] for item in all_snapshots[name]],
                **log_warnings(campaign / name),
            }
            for name in model_names
        },
        "kmax_h_mpc": args.kmax,
        "epoch_sampling": "nearest dump" if args.nearest else "log-P/log-a interpolation",
        "figure_png": str(png),
        "figure_pdf": str(pdf),
        "ratio_csv": str(csv_path),
    }
    report_path = figure_dir / "validation_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
