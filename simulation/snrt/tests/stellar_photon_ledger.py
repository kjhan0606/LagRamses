#!/usr/bin/env python3
"""Synthetic contract test for the stellar SED-to-ledger converter."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.primordial import group_spectral_closure_from_metadata
from snrt_core.source_ledger import read_photon_source_ledger_csv


TOOL = ROOT / "tools" / "p4_build_stellar_photon_ledger.py"
SPEC = importlib.util.spec_from_file_location("p4_build_stellar_photon_ledger", TOOL)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load stellar ledger builder")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write_catalogue(path: Path) -> None:
    fields = [
        "source_id",
        "position_x_code",
        "position_y_code",
        "position_z_code",
        "age_myr",
        "birth_metallicity_mass_fraction",
        "mass_msun",
        "initial_mass_msun",
    ]
    rows = [
        (101, 0.1, 0.2, 0.3, 1.0, 0.01, 8.0, 10.0),
        (102, 0.4, 0.5, 0.6, 10.0, 0.02, 16.0, 20.0),
        (103, 0.7, 0.8, 0.9, 100.0, 0.04, 24.0, 30.0),
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)


def _write_sed(path: Path) -> None:
    ages = (1.0, 10.0, 100.0)
    metallicities = (0.5, 1.0, 2.0)
    energies = np.unique(
        np.concatenate((np.geomspace(0.01, 10000.0, 64), MODULE._read_group_edges(MODULE.DEFAULT_GROUP_EDGES)))
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(MODULE.REQUIRED_SED_COLUMNS)
        for age in ages:
            for metallicity in metallicities:
                for energy in energies:
                    rate = 0.0
                    if energy < 500.0:
                        rate = 1.0e38 * (age / 10.0) ** -0.2 * metallicity**0.1 * (
                            energy / 1.0
                        ) ** -1.1
                    writer.writerow((age, metallicity, energy, rate))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="stellar-ledger-test-") as directory:
        work = Path(directory)
        catalogue = work / "catalogue.csv"
        sed = work / "sed.csv"
        output = work / "ledger.csv"
        metadata_path = work / "ledger.json"
        _write_catalogue(catalogue)
        _write_sed(sed)
        metadata = MODULE.build_ledger(
            catalogue_path=catalogue,
            sed_table_path=sed,
            group_edges_path=MODULE.DEFAULT_GROUP_EDGES,
            output_csv=output,
            output_metadata=metadata_path,
            solar_mass_fraction=0.02,
            scale_factor=0.15,
            mass_field="initial_mass_msun",
            escape_fraction=0.8,
        )
        with output.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 3
        q_fields = [f"q_group_{index}_s" for index in range(9)]
        q = np.asarray([[float(row[field]) for field in q_fields] for row in rows])
        assert np.all(np.isfinite(q)) and np.all(q[:, :7] > 0.0)
        assert np.all(q[:, 7:] == 0.0)
        assert metadata["source_count"] == 3
        assert len(metadata["group_edges_ev"]) == 10
        assert metadata["closure_complete"] is False
        assert metadata["group_spectral_closure"]["group_status"][7:] == [
            "empty_source_group_zero_photons",
            "empty_source_group_zero_photons",
        ]
        ledger = read_photon_source_ledger_csv(output)
        assert ledger.photon_luminosity_s.shape == (3, 9)
        assert np.allclose(ledger.photon_luminosity_s.sum(axis=0), metadata["group_photon_rate_total_s"])
        closure = group_spectral_closure_from_metadata(metadata)
        assert closure.cross_sections.hydrogen_i.shape == (9,)
        saved = json.loads(metadata_path.read_text(encoding="utf-8"))
        assert saved["status"] == "complete_stellar_photon_ledger"
        assert saved["normalization"]["mass_field_used"] == "initial_mass_msun"
        assert saved["normalization"]["escape_fraction"] == 0.8
    print("STELLAR_PHOTON_LEDGER_TEST_OK sources=3 groups=9")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
