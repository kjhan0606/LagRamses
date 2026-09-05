#!/usr/bin/env python3
"""Test deterministic STAR+AGN photon-ledger merging and provenance."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools import merge_photon_source_ledgers as merger  # noqa: E402


def _write_component(
    directory: Path,
    name: str,
    source_id: int,
    source_kind: str,
    values: list[list[float]],
    identity: str | None,
) -> tuple[Path, Path]:
    ledger_path = directory / f"{name}.csv"
    metadata_path = directory / f"{name}.json"
    with ledger_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("source_id", "source_kind", "x_code", "y_code", "z_code", "q_group_0_s", "q_group_1_s"))
        for index, row in enumerate(values):
            writer.writerow((source_id + index, source_kind, 0.25 + 0.1 * index, 0.5, 0.75, *row))
    totals = np.asarray(values, dtype=np.float64).sum(axis=0)
    metadata = {
        "schema": "component_photon_source_ledger_v1",
        "group_edges_ev": [1.0, 10.0, 100.0],
        "group_photon_rate_total_s": totals.tolist(),
        "source_sed_identity": identity,
        "source_sed_sha256": None if identity is None else "a" * 64,
        "source_sed_contract": {
            "status": "candidate_component_sed" if identity else "reference_control",
            "identity": identity,
            "input_sha256": None if identity is None else "a" * 64,
        },
        "groups": [
            {
                "index": 0,
                "energy_interval_ev": [1.0, 10.0],
                "photon_weighted_mean_energy_ev": 4.0,
            },
            {
                "index": 1,
                "energy_interval_ev": [10.0, 100.0],
                "photon_weighted_mean_energy_ev": 40.0,
            },
        ],
        "group_spectral_closure": {
            "cross_sections_cm2": {
                "hydrogen_i": [1.0e-18, 1.0e-19],
                "helium_i": [2.0e-18, 2.0e-19],
                "helium_ii": [3.0e-18, 3.0e-19],
            },
            "photoelectron_excess_energy_ev": {
                "hydrogen_i": [1.0, 10.0],
                "helium_i": [2.0, 20.0],
                "helium_ii": [3.0, 30.0],
            },
        },
    }
    metadata_path.write_text(json.dumps(metadata) + "\n", encoding="utf-8")
    return ledger_path, metadata_path


def main() -> int:
    with TemporaryDirectory(prefix="merge-source-ledger-test-") as directory:
        work = Path(directory)
        star_ledger, star_metadata = _write_component(
            work, "star", 1, "star", [[2.0, 3.0], [1.0, 2.0]], "1" * 64
        )
        agn_ledger, agn_metadata = _write_component(
            work, "agn", 100, "agn", [[4.0, 5.0]], "2" * 64
        )
        output = work / "merged.csv"
        output_metadata = work / "merged.json"
        result = merger.merge_ledgers(
            [star_ledger, agn_ledger],
            [star_metadata, agn_metadata],
            output,
            output_metadata,
        )
        assert result["source_count"] == 3
        assert result["source_sed_identity"]
        assert result["dust_binding"]["component_only_sidecars_allowed"] is False
        assert np.allclose(result["group_photon_rate_total_s"], [7.0, 10.0])
        assert json.loads(output_metadata.read_text(encoding="utf-8"))["source_sed_identity"] == result["source_sed_identity"]

        duplicate_ledger, duplicate_metadata = _write_component(
            work, "duplicate", 1, "agn", [[1.0, 1.0]], "3" * 64
        )
        duplicate_output = work / "duplicate-merged.csv"
        duplicate_output_metadata = work / "duplicate-merged.json"
        try:
            merger.merge_ledgers(
                [star_ledger, duplicate_ledger],
                [star_metadata, duplicate_metadata],
                duplicate_output,
                duplicate_output_metadata,
            )
        except ValueError as error:
            assert "globally unique" in str(error)
        else:
            raise AssertionError("duplicate source IDs were accepted")
        assert not duplicate_output.exists()
        assert not duplicate_output_metadata.exists()

        wrong_edges = json.loads(agn_metadata.read_text(encoding="utf-8"))
        wrong_edges["group_edges_ev"] = [1.0, 11.0, 100.0]
        wrong_edges["groups"][0]["energy_interval_ev"] = [1.0, 11.0]
        wrong_edges["groups"][1]["energy_interval_ev"] = [11.0, 100.0]
        wrong_metadata = work / "wrong.json"
        wrong_metadata.write_text(json.dumps(wrong_edges) + "\n", encoding="utf-8")
        wrong_output = work / "wrong-merged.csv"
        wrong_output_metadata = work / "wrong-merged.json"
        try:
            merger.merge_ledgers(
                [star_ledger, agn_ledger],
                [star_metadata, wrong_metadata],
                wrong_output,
                wrong_output_metadata,
            )
        except ValueError as error:
            assert "identical group edges" in str(error)
        else:
            raise AssertionError("mismatched group edges were accepted")
        assert not wrong_output.exists()
        assert not wrong_output_metadata.exists()

    print("MERGE_PHOTON_SOURCE_LEDGERS_TEST_OK components=2 sources=3 mixed_dust_gate=1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
