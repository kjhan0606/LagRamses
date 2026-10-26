#!/usr/bin/env python3
"""Direct statistical regression test for the production SIDM angular sampler."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "patch" / "cuRamses" / "sidm_angular.f90"
DRIVER = Path(__file__).with_suffix(".f90")
EPS2 = np.array((2.0e-3, 2.0e-2, 2.0e-1, 2.0, 20.0))
ARCHIVE_PRODUCTION_REVISION = "0a9d76850aca265e53cebc36843e00c184592318"
ISOTROPIC_MOMENTS = {"momentum": 1.0, "viscosity": 2.0 / 3.0,
                     "identical_particle": 0.5}


def compact_fortran(source: str) -> str:
    """Remove insignificant spacing and comments for a source-formula check."""
    return "".join(
        "".join(line.split("!", maxsplit=1)[0].split()).replace("&", "").lower()
        for line in source.splitlines()
    )


def verify_archived_production_formula() -> None:
    """Tie the refactored helper to the formula in the archived run binary."""
    archive = subprocess.run(
        [
            "git",
            "-C",
            str(ROOT),
            "show",
            f"{ARCHIVE_PRODUCTION_REVISION}:patch/cuRamses/sidm_scatter.f90",
        ],
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    formula = compact_fortran(archive)
    required = (
        "eps2=2.0d0*sidm_epsilon",
        "a_ruth=1.0d0/eps2",
        "b_ruth=1.0d0/(2.0d0+eps2)",
        "cos_theta=1.0d0+eps2-1.0d0/(b_ruth+r1*(a_ruth-b_ruth))",
    )
    missing = [term for term in required if term not in formula]
    if missing:
        raise AssertionError(
            "archived production Rutherford formula differs from the "
            f"tested helper: {missing}"
        )


def rutherford_cdf(cos_theta: np.ndarray, eps2: float) -> np.ndarray:
    lower = 1.0 / (2.0 + eps2)
    upper = 1.0 / eps2
    return (1.0 / (1.0 - cos_theta + eps2) - lower) / (upper - lower)


def kolmogorov_smirnov_uniform(values: np.ndarray) -> float:
    ordered = np.sort(values)
    sample_size = ordered.size
    upper = np.arange(1, sample_size + 1, dtype=float) / sample_size - ordered
    lower = ordered - np.arange(sample_size, dtype=float) / sample_size
    return float(max(np.max(upper), np.max(lower)))


def mean_momentum_transfer(eps2: float) -> float:
    """Analytic <1-cos(theta)> for the sampled Rutherford kernel."""
    require_positive_finite(eps2, "eps2")
    if eps2 >= 20.0:
        return weak_anisotropy_moment(eps2, viscosity=False)
    numerator = (
        math.log1p(2.0 / eps2)
        + eps2 / (2.0 + eps2)
        - 1.0
    )
    denominator = 1.0 / eps2 - 1.0 / (2.0 + eps2)
    return numerator / denominator


def mean_viscosity_transfer(eps2: float) -> float:
    """Integral of 1-mu^2 against the normalized sampled angular density."""
    require_positive_finite(eps2, "eps2")
    if eps2 >= 20.0:
        return weak_anisotropy_moment(eps2, viscosity=True)
    return eps2 * (eps2 + 2.0) * ((eps2 + 1.0) * math.log1p(2.0 / eps2) - 2.0)


def mean_identical_particle_transfer(eps2: float) -> float:
    """Integral of 1-|mu|, without a conventional factor-of-two rescaling."""
    require_positive_finite(eps2, "eps2")
    inverse_product = (1.0 / eps2) / (eps2 + 2.0)
    if inverse_product == 0.0:
        return 0.5
    return 0.5 * math.log1p(inverse_product) / inverse_product


def require_positive_finite(value: float, name: str) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")


def weak_anisotropy_moment(eps2: float, *, viscosity: bool) -> float:
    """Evaluate the exact moment series without large-b cancellation.

    With t=(1-mu)/2 and r=2/b, p(t)=(1+r)/(1+r*t)^2.
    For b>=20, 24 integrated series terms have a truncation error <1e-22.
    """
    r = 2.0 / eps2
    if viscosity:
        terms = ((k + 1.0) * (-r)**k / ((k + 2.0) * (k + 3.0))
                 for k in range(24))
        return 4.0 * (1.0 + r) * math.fsum(terms)
    terms = ((k + 1.0) * (-r)**k / (k + 2.0) for k in range(24))
    return 2.0 * (1.0 + r) * math.fsum(terms)


def matched_isotropic_normalizations(eps2: float) -> dict[str, dict[str, float]]:
    """Match weighted cross sections, not inverse weighted fractions alone.

    The ratio applies at fixed b and with the same velocity law. Equal angular
    moments do not establish equal halo evolution or measured wall time.
    """
    moments = {"momentum": mean_momentum_transfer(eps2),
               "viscosity": mean_viscosity_transfer(eps2),
               "identical_particle": mean_identical_particle_transfer(eps2)}
    return {
        name: {"rutherford_weight_mean": moment,
               "isotropic_weight_mean": ISOTROPIC_MOMENTS[name],
               "isotropically_normalized_weight_mean": moment / ISOTROPIC_MOMENTS[name],
               "total_over_weighted_cross_section": 1.0 / moment,
               "rutherford_total_over_isotropic_total": ISOTROPIC_MOMENTS[name] / moment}
        for name, moment in moments.items()
    }


def yukawa_model_contract() -> dict[str, object]:
    """Audit implemented formulas against a Born reference without changing physics."""
    source = SOURCE.with_name("sidm_scatter.f90")
    formula = compact_fortran(source.read_text())
    required = (
        "sigma_over_m=sidm_cross_section/(1.0d0+(v_rel_kms/sidm_v0)**2)**2",
        "eps2=2.0d0*sidm_v0**2/max(v_rel_kms**2,1.0d0)",
        "callsidm_pair_probability(sigma_over_m,m1,m2,v_rel_mag,",
    )
    missing = [term for term in required if term not in formula]
    if missing:
        raise AssertionError(f"production Yukawa formula contract changed: {missing}")
    v0 = 100.0
    cases = []
    for velocity in (0.1, 1.0, 10.0, 100.0, 1000.0):
        ratio_squared = (velocity / v0)**2
        eps2 = 2.0 * v0**2 / max(velocity**2, 1.0)
        moments = matched_isotropic_normalizations(eps2)
        total_fraction = 1.0 / (1.0 + ratio_squared)**2
        cases.append({
            "velocity_kms": velocity, "eps2_with_code_floor": eps2,
            "angular_floor_active": velocity < 1.0,
            "implemented_total_over_sigma0": total_fraction,
            "born_total_over_sigma0_without_floor": 1.0 / (1.0 + ratio_squared),
            "implemented_over_born_total": 1.0 / (1.0 + ratio_squared),
            "transport_over_sigma0": {
                name: total_fraction * item["rutherford_weight_mean"]
                for name, item in moments.items()},
            "matched_isotropic_normalizations": moments,
        })
    return {"scope": "source-formula audit and analytic integrals, not a runtime test",
            "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "v0_kms": v0, "sigma0_convention": "low-velocity total cross section",
            "born_reference": "https://arxiv.org/html/1612.03906#S2.SS4",
            "born_equations": [12, 13],
            "matching_scope": "pointwise velocity only; no single fixed-b factor for all velocities",
            "cases": cases}


def save_figure(destination: Path, results: list[dict[str, object]]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.7), constrained_layout=True)
    eps2 = np.array([item["eps2"] for item in results])
    definitions = (("momentum", r"$\langle1-\mu\rangle$"),
                   ("viscosity", r"$\langle1-\mu^2\rangle$"),
                   ("identical_particle", r"$\langle1-|\mu|\rangle$"))
    for name, label in definitions:
        analytic = np.array([item["transport_moments"][name]["analytic"] for item in results])
        measured = np.array([item["transport_moments"][name]["sample_mean"] for item in results])
        error = np.array([item["transport_moments"][name]["standard_error"] for item in results])
        line, = axes[0].plot(eps2, analytic, label=label)
        axes[0].errorbar(eps2, measured, yerr=5.0 * error, fmt="o", color=line.get_color(),
                        markersize=4, capsize=2)
        matching = [item["isotropic_matching"][name]["rutherford_total_over_isotropic_total"]
                    for item in results]
        axes[2].plot(eps2, matching, "o-", label=label, color=line.get_color())
    axes[0].set(xscale="log", yscale="log", xlabel=r"$2\epsilon$",
                ylabel="Normalized transport moment")
    axes[0].legend(frameon=False, fontsize=9)
    axes[1].plot(eps2, [item["KS"] for item in results], "o-", label="Measured KS")
    axes[1].plot(eps2, [item["KS_limit"] for item in results], "--", label=r"$1.95/\sqrt{N}$")
    axes[1].set(xscale="log", xlabel=r"$2\epsilon$", ylabel="CDF deviation")
    axes[1].legend(frameon=False, fontsize=9)
    axes[2].set(xscale="log", yscale="log", xlabel=r"$2\epsilon$",
                ylabel=r"$\sigma_{\mathrm{tot,R}}/\sigma_{\mathrm{tot,iso}}$ at equal transport")
    axes[2].legend(frameon=False, fontsize=9)
    fig.savefig(destination / "angular_validation.pdf")
    fig.savefig(destination / "angular_validation.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        help="save measured JSON and PDF/PNG verification figures")
    args = parser.parse_args()
    verify_archived_production_formula()
    velocity_contract = yukawa_model_contract()
    compiler = shutil.which("gfortran")
    if compiler is None:
        raise SystemExit("gfortran is required for the SIDM angular regression test")
    with tempfile.TemporaryDirectory(prefix="sidm-angular-") as temporary:
        executable = Path(temporary) / "test_sidm_angular"
        subprocess.run(
            [compiler, "-O2", "-std=f2008", "-Wall", "-Wextra", str(SOURCE),
             str(DRIVER), "-o", str(executable)],
            check=True,
        )
        result = subprocess.run(
            [str(executable)], check=True, text=True, capture_output=True
        )
    samples = np.fromstring(result.stdout, sep=" ").reshape(-1, 2)
    if samples.shape != (50000 * EPS2.size, 2) or not np.all(np.isfinite(samples)):
        raise AssertionError("angular driver must emit 50000 finite samples per regulator")
    results: list[dict[str, object]] = []
    for case_index, eps2 in enumerate(EPS2, start=1):
        cos_theta = samples[samples[:, 0] == case_index, 1]
        if cos_theta.size != 50000:
            raise AssertionError(f"missing angular samples for case {case_index}")
        if np.any(cos_theta < -1.0) or np.any(cos_theta > 1.0):
            raise AssertionError(f"cos(theta) escaped [-1, 1] for eps2={eps2}")
        cdf = rutherford_cdf(cos_theta, eps2)
        statistic = kolmogorov_smirnov_uniform(cdf)
        limit = 1.95 / math.sqrt(cdf.size)
        if statistic > limit:
            raise AssertionError(
                f"KS={statistic:.5f} exceeds {limit:.5f} for eps2={eps2}"
            )
        sampled_momentum_transfer = float(np.mean(1.0 - cos_theta))
        expected_momentum_transfer = mean_momentum_transfer(eps2)
        standard_error = float(np.std(1.0 - cos_theta, ddof=1) / math.sqrt(cdf.size))
        if abs(sampled_momentum_transfer - expected_momentum_transfer) > 5.0 * standard_error:
            raise AssertionError(
                "momentum-transfer mean disagrees with the Rutherford kernel "
                f"for eps2={eps2}"
            )
        moments: dict[str, dict[str, float]] = {}
        for name, values, analytic in (
            ("momentum", 1.0 - cos_theta, expected_momentum_transfer),
            ("viscosity", 1.0 - cos_theta**2, mean_viscosity_transfer(eps2)),
            ("identical_particle", 1.0 - np.abs(cos_theta),
             mean_identical_particle_transfer(eps2)),
        ):
            mean = float(np.mean(values))
            error = float(np.std(values, ddof=1) / math.sqrt(values.size))
            if abs(mean - analytic) > 5.0 * error:
                raise AssertionError(f"{name} transport moment failed for eps2={eps2}")
            moments[name] = {"sample_mean": mean, "analytic": analytic,
                             "standard_error": error, "sigma_deviation": (mean-analytic)/error}
        results.append({"eps2": float(eps2), "epsilon": float(eps2/2.0),
                        "sample_size": int(cdf.size), "KS": statistic, "KS_limit": limit,
                        "transport_moments": moments,
                        "isotropic_matching": matched_isotropic_normalizations(eps2)})
        print(
            f"eps2={eps2:g} epsilon={eps2 / 2.0:g} n={cdf.size} "
            f"KS={statistic:.5f} limit={limit:.5f} "
            f"<1-cos>={sampled_momentum_transfer:.5f} "
            f"analytic={expected_momentum_transfer:.5f} "
            f"sigma_total/sigma_T={1.0 / expected_momentum_transfer:.5f}"
        )
    print(
        "archived-production-formula="
        f"{ARCHIVE_PRODUCTION_REVISION[:12]} matches tested Rutherford CDF"
    )
    print("rotation geometry: 200 cases across both basis branches passed")
    print("three normalized transport moments: all five regulators within five standard errors")
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        payload = {"status": "passed", "scope": "Rutherford helper and lab-frame rotation",
                   "archive_revision": ARCHIVE_PRODUCTION_REVISION,
                   "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                   "driver_sha256": hashlib.sha256(DRIVER.read_bytes()).hexdigest(),
                   "test_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "compiler": subprocess.check_output([compiler, "--version"], text=True).splitlines()[0],
                   "rotation_geometry_cases": 200,
                   "moment_error_bars_standard_errors": 5,
                   "identical_particle_weight": "1-abs(mu), no factor-of-two rescaling",
                   "matching_scope": "equal weighted cross sections at fixed epsilon, not halo-response equivalence",
                   "matching_cost_scope": "cross-section/event-rate ratio at fixed particle configuration; no wall-time forecast",
                   "yukawa_model_contract": velocity_contract,
                   "cases": results}
        (args.output_dir / "angular_validation.json").write_text(
            json.dumps(payload, indent=2, allow_nan=False) + "\n")
        save_figure(args.output_dir, results)


if __name__ == "__main__":
    main()
