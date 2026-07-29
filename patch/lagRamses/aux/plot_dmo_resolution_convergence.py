#!/usr/bin/env python3
"""Compare z=0 DMO power ratios across a uniform-resolution ladder."""

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

from plot_dmo_pk_ratios import (
    linear_theory,
    snapshots,
    spectrum_ratio,
)


MODELS = ("f5", "f6", "n1", "n5", "sym_a")
MODEL_CHOICES = (*MODELS, "phicdm_a01", "cde10")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--levels", nargs="+", type=int, default=[6, 7, 8])
    parser.add_argument(
        "--models",
        nargs="+",
        choices=MODEL_CHOICES,
        default=list(MODELS),
        help="model subset to validate (default: all benchmark models)",
    )
    parser.add_argument("--redshift", type=float, default=0.0)
    parser.add_argument("--kmax", type=float, default=0.5)
    parser.add_argument(
        "--pk-estimator",
        choices=("cic", "runtime-ngp"),
        default="cic",
        help="spectrum source; CIC is required for precision validation",
    )
    parser.add_argument(
        "--shot-noise",
        choices=("none", "poisson"),
        default="none",
        help=(
            "CIC spectrum convention; matched perturbed-lattice DMO "
            "validation should normally use 'none'"
        ),
    )
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
    parser.add_argument(
        "--residual-target",
        type=float,
        default=1.0e-3,
        help="fractional PASS threshold (default: 1e-3 = 0.1 percent)",
    )
    parser.add_argument(
        "--require-resolution-pass",
        action="store_true",
        help="return a nonzero status unless every finest adjacent pair passes",
    )
    parser.add_argument(
        "--require-theory-pass",
        action="store_true",
        help="return nonzero unless every finest-level large-scale theory test passes",
    )
    parser.add_argument("--output-stem", default="z0_resolution_convergence")
    return parser.parse_args()


def cic_snapshots(model_dir: Path, shot_noise: str) -> list[dict]:
    """Read common-mesh, CIC-deconvolved spectra."""
    result = []
    for path in sorted(model_dir.glob("output_*/pk_cic_*.dat")):
        header = "\n".join(path.read_text().splitlines()[:9])
        match = re.search(r"a_exp\s*=\s*([0-9.Ee+-]+)", header)
        if not match:
            raise ValueError(f"missing a_exp header in {path}")
        shot_match = re.search(
            r"shot_noise.*=\s*([0-9.Ee+-]+)", header
        )
        data = np.loadtxt(path)
        result.append(
            {
                "path": str(path),
                "aexp": float(match.group(1)),
                "z": 1.0 / float(match.group(1)) - 1.0,
                "k": data[:, 0],
                "pk_raw": data[:, 1],
                "pk": data[:, 2] if shot_noise == "poisson" else data[:, 1],
                "nmodes": data[:, 3],
                "shot_noise": (
                    float(shot_match.group(1)) if shot_match else None
                ),
            }
        )
    return result


def exact_snapshot(
    model_dir: Path,
    redshift: float,
    tolerance: float,
    estimator: str,
    shot_noise: str,
) -> dict:
    """Return a measured spectrum only when it matches the requested epoch."""
    target_aexp = 1.0 / (1.0 + redshift)
    candidates = (
        cic_snapshots(model_dir, shot_noise)
        if estimator == "cic"
        else snapshots(model_dir)
    )
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
                    campaign / model,
                    args.redshift,
                    args.aexp_tolerance,
                    args.pk_estimator,
                    args.shot_noise,
                )
                for model in ("lcdm", *args.models)
            }
        except RuntimeError as error:
            incomplete[level] = str(error)
            continue
        available[level] = campaign
        snapshots_by_level[level] = selected
    if not available:
        details = "; ".join(f"L{level}: {reason}" for level, reason in incomplete.items())
        raise RuntimeError(f"no resolution has complete exact-epoch P(k): {details}")

    phase_records = {}
    for level, campaign in available.items():
        metadata = json.loads((campaign / "campaign.json").read_text())
        phase_records[level] = {
            "seed": metadata.get("seed"),
            # Legacy campaigns placed the seed independently at levelmin.
            "phase_anchor_level": metadata.get("phase_anchor_level", level),
            "explicit_phase_anchor": "phase_anchor_level" in metadata,
        }
    phase_seeds = {record["seed"] for record in phase_records.values()}
    phase_anchors = {
        record["phase_anchor_level"] for record in phase_records.values()
    }
    phase_matched = (
        len(phase_seeds) == 1
        and None not in phase_seeds
        and len(phase_anchors) == 1
        and next(iter(phase_anchors)) >= max(available)
        and all(
            record["explicit_phase_anchor"]
            for record in phase_records.values()
        )
    )

    finest_level = max(available)
    campaign_metadata = json.loads(
        (available[finest_level] / "campaign.json").read_text()
    )
    boxlen = float(campaign_metadata["boxlen_mpc_h"])
    theory = linear_theory(
        available[finest_level],
        ["lcdm", *args.models],
        [args.redshift],
        args.kmax,
    )
    fig, axes = plt.subplots(
        2, 3, figsize=(10.0, 5.8), sharex=True, sharey=True,
        constrained_layout=True,
    )
    axes = axes.ravel()
    rows = []
    metrics = []
    resolution_pairs = []
    level_colors = plt.cm.viridis(np.linspace(0.18, 0.88, len(available)))
    level_color = dict(zip(sorted(available), level_colors))

    for ax, model in zip(axes, args.models):
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
            if not np.allclose(
                reference["k"], sample["k"], rtol=1.0e-10, atol=0.0
            ):
                raise ValueError(
                    f"incompatible k grids: {reference['path']} and "
                    f"{sample['path']}"
                )
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
            worst_index = int(np.argmax(np.abs(residual)))
            max_abs_residual = float(abs(residual[worst_index]))
            large_scale_indices = np.flatnonzero(large_scale)
            large_scale_worst_index = (
                int(
                    large_scale_indices[
                        np.argmax(np.abs(residual[large_scale_indices]))
                    ]
                )
                if large_scale_indices.size
                else None
            )
            large_scale_max_abs_residual = (
                float(abs(residual[large_scale_worst_index]))
                if large_scale_worst_index is not None
                else None
            )
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
                    "max_abs_theory_residual": max_abs_residual,
                    "worst_k_h_mpc": float(k[worst_index]),
                    "worst_signed_theory_residual": float(
                        residual[worst_index]
                    ),
                    "theory_full_range_pass": (
                        max_abs_residual <= args.residual_target
                    ),
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
                    "large_scale_max_abs_theory_residual": (
                        large_scale_max_abs_residual
                    ),
                    "large_scale_worst_k_h_mpc": (
                        float(k[large_scale_worst_index])
                        if large_scale_worst_index is not None
                        else None
                    ),
                    "large_scale_worst_signed_theory_residual": (
                        float(residual[large_scale_worst_index])
                        if large_scale_worst_index is not None
                        else None
                    ),
                    "theory_large_scale_pass": (
                        large_scale_max_abs_residual <= args.residual_target
                        if large_scale_max_abs_residual is not None
                        else False
                    ),
                }
            )
            rows.extend(
                (level, 2**level, model, args.redshift, kval, value, pred)
                for kval, value, pred in zip(k, ratio, predicted)
            )

        ordered_levels = sorted(level_ratios)
        for coarse_level, fine_level in zip(
            ordered_levels[:-1], ordered_levels[1:]
        ):
            k, coarse_ratio = level_ratios[coarse_level]
            fine_k, fine_ratio = level_ratios[fine_level]
            overlap = (k >= fine_k.min()) & (k <= fine_k.max())
            k = k[overlap]
            coarse_ratio = coarse_ratio[overlap]
            fine_on_k = np.interp(np.log(k), np.log(fine_k), fine_ratio)
            delta = coarse_ratio / fine_on_k - 1.0
            low = k <= min(args.large_scale_kmax, args.kmax)
            absolute_delta = np.abs(delta)
            prefix_pass = (
                np.maximum.accumulate(absolute_delta)
                <= args.residual_target
            )
            worst_index = int(np.argmax(np.abs(delta)))
            max_abs = float(abs(delta[worst_index]))
            low_indices = np.flatnonzero(low)
            low_worst_index = (
                int(low_indices[np.argmax(np.abs(delta[low_indices]))])
                if low_indices.size
                else None
            )
            low_max_abs = (
                float(abs(delta[low_worst_index]))
                if low_worst_index is not None
                else None
            )
            resolution_pairs.append(
                {
                    "model": model,
                    "coarse_level": coarse_level,
                    "fine_level": fine_level,
                    "kmax_h_mpc": args.kmax,
                    "n_bins": int(k.size),
                    "n_failing_bins": int(
                        np.count_nonzero(
                            absolute_delta > args.residual_target
                        )
                    ),
                    "largest_contiguous_kmax_pass_h_mpc": (
                        float(k[np.flatnonzero(prefix_pass)[-1]])
                        if np.any(prefix_pass)
                        else None
                    ),
                    "rms_fractional_residual": float(
                        np.sqrt(np.mean(delta**2))
                    ),
                    "max_abs_fractional_residual": max_abs,
                    "worst_k_h_mpc": float(k[worst_index]),
                    "worst_signed_fractional_residual": float(
                        delta[worst_index]
                    ),
                    "full_range_pass": max_abs <= args.residual_target,
                    "large_scale_kmax_h_mpc": min(
                        args.large_scale_kmax, args.kmax
                    ),
                    "large_scale_n_bins": int(np.count_nonzero(low)),
                    "large_scale_rms_fractional_residual": (
                        float(np.sqrt(np.mean(delta[low] ** 2)))
                        if np.any(low)
                        else None
                    ),
                    "large_scale_max_abs_fractional_residual": low_max_abs,
                    "large_scale_worst_k_h_mpc": (
                        float(k[low_worst_index])
                        if low_worst_index is not None
                        else None
                    ),
                    "large_scale_worst_signed_fractional_residual": (
                        float(delta[low_worst_index])
                        if low_worst_index is not None
                        else None
                    ),
                    "large_scale_pass": (
                        low_max_abs <= args.residual_target
                        if low_max_abs is not None
                        else False
                    ),
                }
            )

        ax.axhline(1.0, color="0.5", lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_title(model)
        ax.grid(alpha=0.18)

    for ax in axes[len(args.models):]:
        ax.set_visible(False)
    if not phase_matched:
        fig.suptitle(
            "Diagnostic only: resolution ICs do not share one explicit "
            "white-noise phase anchor",
            color="firebrick",
            fontsize=10,
        )
    for index, ax in enumerate(axes[:len(args.models)]):
        if index >= 3 or index + 3 >= len(args.models):
            ax.set_xlabel(r"$k\ [h\,{\rm Mpc}^{-1}]$")
        if index % 3 == 0:
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
    finest_metrics = [item for item in metrics if item["level"] == finest_level]
    finest_pairs = [
        item for item in resolution_pairs if item["fine_level"] == finest_level
    ]
    report = {
        "root": str(root),
        "available_levels": sorted(available),
        "models": list(args.models),
        "redshift": args.redshift,
        "pk_estimator": args.pk_estimator,
        "shot_noise": args.shot_noise,
        "boxlen_mpc_h": boxlen,
        "kmax_h_mpc": args.kmax,
        "large_scale_kmax_h_mpc": min(args.large_scale_kmax, args.kmax),
        "residual_target_fraction": args.residual_target,
        "residual_target_percent": 100.0 * args.residual_target,
        "incomplete_levels": incomplete,
        "phase_matching": {
            "certified": phase_matched,
            "levels": phase_records,
            "requirement": (
                "one explicit common phase_anchor_level at or above the "
                "finest compared resolution, with one common seed"
            ),
        },
        "metrics": metrics,
        "resolution_pairs": resolution_pairs,
        "acceptance": {
            "theory_large_scale_all_levels_pass": all(
                item["theory_large_scale_pass"] for item in metrics
            ),
            "finest_theory_large_scale_all_pass": all(
                item["theory_large_scale_pass"] for item in finest_metrics
            ),
            "adjacent_resolution_full_range_all_pass": (
                bool(resolution_pairs)
                and all(item["full_range_pass"] for item in resolution_pairs)
            ),
            "adjacent_resolution_large_scale_all_pass": (
                bool(resolution_pairs)
                and all(item["large_scale_pass"] for item in resolution_pairs)
            ),
            "finest_pair_full_range_all_pass": (
                bool(finest_pairs)
                and all(item["full_range_pass"] for item in finest_pairs)
            ),
            "finest_pair_large_scale_all_pass": (
                bool(finest_pairs)
                and all(item["large_scale_pass"] for item in finest_pairs)
            ),
            "phase_matched_across_resolutions": phase_matched,
            "resolution_convergence_certified": (
                phase_matched
                and bool(finest_pairs)
                and all(item["full_range_pass"] for item in finest_pairs)
            ),
            "finest_resolution_independently_certified": False,
            "note": (
                "The finest level is not independently converged unless an "
                "additional finer level or a validated extrapolation is used."
            ),
        },
        "figure_png": str(png),
        "figure_pdf": str(pdf),
        "csv": str(csv_path),
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if (
        args.require_resolution_pass
        and not report["acceptance"]["resolution_convergence_certified"]
    ):
        return 2
    if (
        args.require_theory_pass
        and not report["acceptance"]["finest_theory_large_scale_all_pass"]
    ):
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
