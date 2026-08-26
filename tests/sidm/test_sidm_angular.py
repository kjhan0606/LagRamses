#!/usr/bin/env python3
"""Direct statistical regression test for the production SIDM angular sampler."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "patch" / "cuRamses" / "sidm_angular.f90"
DRIVER = Path(__file__).with_suffix(".f90")
EPS2 = np.array((2.0e-3, 2.0e-2, 2.0e-1))
ARCHIVE_PRODUCTION_REVISION = "0a9d76850aca265e53cebc36843e00c184592318"


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
    numerator = (
        math.log((2.0 + eps2) / eps2)
        + eps2 / (2.0 + eps2)
        - 1.0
    )
    denominator = 1.0 / eps2 - 1.0 / (2.0 + eps2)
    return numerator / denominator


def main() -> None:
    verify_archived_production_formula()
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
    for case_index, eps2 in enumerate(EPS2, start=1):
        cos_theta = samples[samples[:, 0] == case_index, 1]
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


if __name__ == "__main__":
    main()
