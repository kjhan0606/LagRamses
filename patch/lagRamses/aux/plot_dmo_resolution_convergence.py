#!/usr/bin/env python3
"""Compare z=0 DMO power ratios across a uniform-resolution ladder."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_dmo_pk_ratios import (
    linear_theory,
    snapshots,
    spectrum_ratio,
)


MODELS = ("f5", "f6", "n1", "n5", "sym_a")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--levels", nargs="+", type=int, default=[6, 7, 8])
    parser.add_argument("--redshift", type=float, default=0.0)
    parser.add_argument("--kmax", type=float, default=0.5)
    parser.add_argument(
        "--large-scale-kmax",
        type=float,
        default=0.2,
        help="upper k used for the separate large-scale residual metric",
    )
    parser.add_argument(
        "--aexp-tolerance",
        type=float,
        default=1.0e-6,
        help="absolute scale-factor tolerance for accepting a completed epoch",
    )
    parser.add_argument("--output-stem", default="z0_resolution_convergence")
    return parser.parse_args()


def exact_snapshot(model_dir: Path, redshift: float, tolerance: float) -> dict:
    """Return a measured spectrum only when it matches the requested epoch."""
    target_aexp = 1.0 / (1.0 + redshift)
    candidates = snapshots(model_dir)
    if not candidates:
        raise RuntimeError(f"no P(k) snapshots in {model_dir}")
    result = min(candidates, key=lambda item: abs(item["aexp"] - target_aexp))
    delta = abs(result["aexp"] - target_aexp)
    if delta > tolerance:
        raise RuntimeError(
            f"{model_dir} has no spectrum at a={target_aexp:.9g}; "
            f"nearest is a={result['aexp']:.9g}"
        )
    return result


def main() -> int:
    args = arguments()
    root = args.root.resolve()
    campaigns = {
        level: root / f"L{level}_{2**level:03d}" for level in args.levels
    }
    missing = [str(path) for path in campaigns.values() if not path.is_dir()]
    if missing:
        raise FileNotFoundError(f"missing campaigns: {', '.join(missing)}")

    available = {}
    snapshots_by_level = {}
    incomplete = {}
    for level, campaign in campaigns.items():
        try:
            selected = {
                model: exact_snapshot(
                    campaign / model, args.redshift, args.aexp_tolerance
                )
                for model in ("lcdm", *MODELS)
            }
        except RuntimeError as error:
            incomplete[level] = str(error)
            continue
        available[level] = campaign
        snapshots_by_level[level] = selected
    if not available:
        details = "; ".join(f"L{level}: {reason}" for level, reason in incomplete.items())
        raise RuntimeError(f"no resolution has complete exact-epoch P(k): {details}")

    finest_level = max(available)
    campaign_metadata = json.loads(
        (available[finest_level] / "campaign.json").read_text()
    )
    boxlen = float(campaign_metadata["boxlen_mpc_h"])
    theory = linear_theory(
        available[finest_level], ["lcdm", *MODELS], [args.redshift], args.kmax
    )
    fig, axes = plt.subplots(
        2, 3, figsize=(10.0, 5.8), sharex=True, sharey=True,
        constrained_layout=True,
    )
    axes = axes.ravel()
    rows = []
    metrics = []
    level_colors = plt.cm.viridis(np.linspace(0.18, 0.88, len(available)))
    level_color = dict(zip(sorted(available), level_colors))

    for ax, model in zip(axes, MODELS):
        theory_reference = theory["lcdm"][args.redshift]
        theory_sample = theory[model][args.redshift]
        theory_k = np.geomspace(
            max(0.98 * 2.0 * np.pi / boxlen, theory_reference["k"].min()),
            args.kmax,
            200,
        )
        theory_ratio = spectrum_ratio(theory_sample, theory_reference, theory_k)
        ax.plot(theory_k, theory_ratio, "k--", lw=1.4, label="lagCAMB linear")

        level_ratios = {}
        for level in sorted(available):
            reference = snapshots_by_level[level]["lcdm"]
            sample = snapshots_by_level[level][model]
            good = (
                (reference["pk"] > 0.0)
                & (sample["pk"] > 0.0)
                & (reference["nmodes"] > 0.0)
                & (sample["nmodes"] > 0.0)
                & (reference["k"] <= args.kmax)
            )
            k = reference["k"][good]
            ratio = sample["pk"][good] / reference["pk"][good]
            level_ratios[level] = (k, ratio)
            ax.plot(
                k,
                ratio,
                color=level_color[level],
                lw=1.3,
                label=rf"$L_{{\rm min}}={level}$ ({2**level}$^3$)",
            )
            predicted = spectrum_ratio(theory_sample, theory_reference, k)
            residual = ratio / predicted - 1.0
            large_scale = k <= min(args.large_scale_kmax, args.kmax)
            large_scale_residual = residual[large_scale]
            metrics.append(
                {
                    "level": level,
                    "particle_load": 2**level,
                    "model": model,
                    "redshift": args.redshift,
                    "kmax_h_mpc": args.kmax,
                    "n_bins": int(k.size),
                    "median_abs_theory_residual": float(
                        np.median(np.abs(residual))
                    ),
                    "rms_theory_residual": float(np.sqrt(np.mean(residual**2))),
                    "large_scale_kmax_h_mpc": min(
                        args.large_scale_kmax, args.kmax
                    ),
                    "large_scale_n_bins": int(np.count_nonzero(large_scale)),
                    "large_scale_median_abs_theory_residual": (
                        float(np.median(np.abs(large_scale_residual)))
                        if large_scale_residual.size
                        else None
                    ),
                    "large_scale_rms_theory_residual": (
                        float(np.sqrt(np.mean(large_scale_residual**2)))
                        if large_scale_residual.size
                        else None
                    ),
                }
            )
            rows.extend(
                (level, 2**level, model, args.redshift, kval, value, pred)
                for kval, value, pred in zip(k, ratio, predicted)
            )

        if len(level_ratios) > 1:
            fine_k, fine_ratio = level_ratios[finest_level]
            for level, (k, ratio) in level_ratios.items():
                if level == finest_level:
                    continue
                fine_on_k = np.interp(np.log(k), np.log(fine_k), fine_ratio)
                delta = ratio / fine_on_k - 1.0
                metric = next(
                    item
                    for item in metrics
                    if item["level"] == level and item["model"] == model
                )
                metric["rms_residual_to_finest"] = float(
                    np.sqrt(np.mean(delta**2))
                )
                low = k <= min(args.large_scale_kmax, args.kmax)
                metric["large_scale_rms_residual_to_finest"] = (
                    float(np.sqrt(np.mean(delta[low] ** 2)))
                    if np.any(low)
                    else None
                )

        ax.axhline(1.0, color="0.5", lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_title(model)
        ax.grid(alpha=0.18)

    axes[-1].set_visible(False)
    for ax in axes[3:5]:
        ax.set_xlabel(r"$k\ [h\,{\rm Mpc}^{-1}]$")
    for ax in (axes[0], axes[3]):
        ax.set_ylabel(r"$P(k)/P_{\Lambda{\rm CDM}}(k)$")
    axes[0].legend(frameon=False, fontsize=8)

    figure_dir = root / "figures"
    figure_dir.mkdir(exist_ok=True)
    png = figure_dir / f"{args.output_stem}.png"
    pdf = figure_dir / f"{args.output_stem}.pdf"
    csv_path = figure_dir / f"{args.output_stem}.csv"
    report_path = figure_dir / f"{args.output_stem}.json"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)

    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "level",
                "particle_load",
                "model",
                "redshift",
                "k_h_mpc",
                "simulation_ratio",
                "lagcamb_linear_ratio",
            )
        )
        writer.writerows(rows)
    report = {
        "root": str(root),
        "available_levels": sorted(available),
        "redshift": args.redshift,
        "boxlen_mpc_h": boxlen,
        "kmax_h_mpc": args.kmax,
        "large_scale_kmax_h_mpc": min(args.large_scale_kmax, args.kmax),
        "incomplete_levels": incomplete,
        "metrics": metrics,
        "figure_png": str(png),
        "figure_pdf": str(pdf),
        "csv": str(csv_path),
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
