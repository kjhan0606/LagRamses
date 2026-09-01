#!/usr/bin/env python3
"""Validate strict STAR+AGN ledger merging and aggregate spectral closure."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.primordial import group_spectral_closure_from_metadata
from snrt_core.source_ledger import read_photon_source_ledger_csv
from tools.merge_photon_source_ledgers import merge_ledgers


EDGES = np.loadtxt(ROOT / "config" / "p0_photon_group_edges_ev.txt", comments="#")


def _write_ledger(path: Path, source_id: int, source_kind: str, luminosity: np.ndarray) -> None:
    fields = [
        "source_id",
        "source_kind",
        "x_code",
        "y_code",
        "z_code",
        *[f"q_group_{index}_s" for index in range(luminosity.size)],
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerow((source_id, source_kind, 0.2, 0.3, 0.4, *luminosity))


def _write_metadata(
    path: Path,
    source_kind: str,
    luminosity: np.ndarray,
    hydrogen_sigma: np.ndarray,
    excess: np.ndarray,
    scale_factor: float | None = None,
) -> None:
    mean_energy = np.sqrt(EDGES[:-1] * EDGES[1:])
    metadata = {
        "schema": f"{source_kind}_photon_source_ledger_v1",
        "status": "complete_test_ledger",
        "source_kind": source_kind,
        "source_count": 1,
        "group_edges_ev": EDGES.tolist(),
        "groups": [
            {
                "index": index,
                "energy_interval_ev": [float(EDGES[index]), float(EDGES[index + 1])],
                "photon_weighted_mean_energy_ev": float(mean_energy[index]),
                "total_photon_rate_s": float(luminosity[index]),
            }
            for index in range(luminosity.size)
        ],
        "group_photon_rate_total_s": luminosity.tolist(),
        "group_spectral_closure": {
            "cross_sections_cm2": {
                "hydrogen_i": hydrogen_sigma.tolist(),
                "helium_i": (0.5 * hydrogen_sigma).tolist(),
                "helium_ii": (0.25 * hydrogen_sigma).tolist(),
            },
            "photoelectron_excess_energy_ev": {
                "hydrogen_i": excess.tolist(),
                "helium_i": (0.5 * excess).tolist(),
                "helium_ii": (0.25 * excess).tolist(),
            },
        },
    }
    if scale_factor is not None:
        metadata["source_scale_factor"] = scale_factor
    path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    with TemporaryDirectory(prefix="merge-ledger-test-") as directory:
        work = Path(directory)
        star_luminosity = np.arange(1.0, 10.0) * 1.0e50
        agn_luminosity = np.arange(9.0, 0.0, -1.0) * 1.0e50
        star_sigma = np.arange(1.0, 10.0) * 1.0e-18
        agn_sigma = 2.0 * star_sigma
        star_excess = np.arange(1.0, 10.0)
        agn_excess = 3.0 * star_excess
        star_ledger = work / "star.csv"
        star_metadata = work / "star.json"
        agn_ledger = work / "agn.csv"
        agn_metadata = work / "agn.json"
        _write_ledger(star_ledger, 1001, "star", star_luminosity)
        _write_metadata(star_metadata, "star", star_luminosity, star_sigma, star_excess)
        _write_ledger(agn_ledger, 2001, "agn", agn_luminosity)
        _write_metadata(agn_metadata, "agn", agn_luminosity, agn_sigma, agn_excess)

        output = work / "star_plus_agn.csv"
        output_metadata = work / "star_plus_agn.json"
        merged_metadata = merge_ledgers(
            [star_ledger, agn_ledger],
            [star_metadata, agn_metadata],
            output,
            output_metadata,
        )
        merged = read_photon_source_ledger_csv(output)
        closure = group_spectral_closure_from_metadata(merged_metadata)
        totals = star_luminosity + agn_luminosity
        expected_sigma_group_3 = (
            star_luminosity[3] * star_sigma[3] + agn_luminosity[3] * agn_sigma[3]
        ) / totals[3]
        expected_excess_group_3 = (
            star_luminosity[3] * star_sigma[3] * star_excess[3]
            + agn_luminosity[3] * agn_sigma[3] * agn_excess[3]
        ) / (star_luminosity[3] * star_sigma[3] + agn_luminosity[3] * agn_sigma[3])
        assert merged.photon_luminosity_s.shape == (2, 9)
        assert np.allclose(merged.photon_luminosity_s.sum(axis=0), totals)
        assert np.allclose(closure.cross_sections.hydrogen_i[3], expected_sigma_group_3)
        assert np.allclose(closure.photoelectron_excess_energy_ev[0, 3], expected_excess_group_3)
        assert merged_metadata["source_kind"] == "agn+star"
        assert merged_metadata["group_count"] == 9
        assert merged_metadata["source_id_policy"]["input_offsets"] == [0, 0]

        star_epoch_metadata = work / "star_epoch.json"
        agn_epoch_metadata = work / "agn_epoch.json"
        _write_metadata(star_epoch_metadata, "star", star_luminosity, star_sigma, star_excess, 0.15)
        _write_metadata(agn_epoch_metadata, "agn", agn_luminosity, agn_sigma, agn_excess, 0.20)
        try:
            merge_ledgers(
                [star_ledger, agn_ledger],
                [star_epoch_metadata, agn_epoch_metadata],
                work / "mixed.csv",
                work / "mixed.json",
            )
        except ValueError as error:
            assert "different source scale factors" in str(error)
        else:
            raise AssertionError("mixed-epoch ledger merge was not rejected")

        collision_ledger = work / "collision.csv"
        collision_metadata = work / "collision.json"
        _write_ledger(collision_ledger, 1001, "agn", agn_luminosity)
        _write_metadata(collision_metadata, "agn", agn_luminosity, agn_sigma, agn_excess)
        try:
            merge_ledgers(
                [star_ledger, collision_ledger],
                [star_metadata, collision_metadata],
                work / "collision_fail.csv",
                work / "collision_fail.json",
            )
        except ValueError as error:
            assert "source_id values collide" in str(error)
        else:
            raise AssertionError("source-ID collision was not rejected")

    print("MERGE_PHOTON_LEDGER_TEST_OK sources=2 groups=9 closure=aggregate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
