#!/usr/bin/env python3
"""Validate the AGN photon-ledger converter on the P0 and legacy group tables."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.primordial import group_spectral_closure_from_metadata
from snrt_core.source_ledger import read_photon_source_ledger_csv


def _write_candidates(path: Path) -> None:
    fields = (
        "source_id",
        "source_kind",
        "x_code",
        "y_code",
        "z_code",
        "aexp",
        "mass_msun",
        "inflow_mdot_msun_per_year",
        "bolometric_luminosity_erg_s",
        "radiative_efficiency",
        "accretion_rate_convention",
    )
    rows = (
        (10, "agn", 0.2, 0.3, 0.4, 0.2085, 2.0e7, 0.15, 2.0e45, 0.32, "synthetic"),
        (11, "agn", 0.6, 0.7, 0.8, 0.2085, 1.0e7, 0.05, 5.0e44, 0.30, "synthetic"),
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)


def _run_converter(candidates: Path, output: Path, metadata: Path, *extra: str) -> str:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT)
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "tools" / "p4_build_agn_photon_ledger.py"),
            "--candidates",
            str(candidates),
            "--output",
            str(output),
            "--metadata-output",
            str(metadata),
            *extra,
        ],
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    return result.stdout


def main() -> int:
    with TemporaryDirectory(prefix="agn-ledger-test-") as directory:
        work = Path(directory)
        candidates = work / "candidates.csv"
        _write_candidates(candidates)

        p0_output = work / "p0.csv"
        p0_metadata_path = work / "p0.json"
        stdout = _run_converter(
            candidates,
            p0_output,
            p0_metadata_path,
        )
        assert "P4_AGN_PHOTON_LEDGER_OK" in stdout
        p0_ledger = read_photon_source_ledger_csv(p0_output)
        p0_metadata = json.loads(p0_metadata_path.read_text(encoding="utf-8"))
        p0_closure = group_spectral_closure_from_metadata(p0_metadata)
        assert p0_ledger.photon_luminosity_s.shape == (2, 9)
        assert np.all(p0_ledger.photon_luminosity_s[:, :2] == 0.0)
        assert np.all(p0_ledger.photon_luminosity_s[:, 2:] > 0.0)
        assert p0_metadata["group_table_mode"] == "p0_default"
        assert p0_metadata["groups"][0]["closure_status"] == "agn_sed_below_support_zero_photons"
        assert p0_closure.cross_sections.hydrogen_i.shape == (9,)
        assert np.all(np.isfinite(p0_closure.photoelectron_excess_energy_ev))

        legacy_output = work / "legacy.csv"
        legacy_metadata_path = work / "legacy.json"
        _run_converter(candidates, legacy_output, legacy_metadata_path, "--legacy-five-groups")
        legacy_ledger = read_photon_source_ledger_csv(legacy_output)
        legacy_metadata = json.loads(legacy_metadata_path.read_text(encoding="utf-8"))
        assert legacy_ledger.photon_luminosity_s.shape == (2, 5)
        assert np.all(legacy_ledger.photon_luminosity_s > 0.0)
        assert legacy_metadata["group_table_mode"] == "legacy_five_group_control"

    print("AGN_PHOTON_LEDGER_TEST_OK p0_groups=9 legacy_groups=5")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
