#!/usr/bin/env python3
"""Validate the staged BPASS candidate ledger and its provenance sidecar."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.source_ledger import read_photon_source_ledger_csv


def main() -> int:
    ledger_path = ROOT / "data" / "feedback_transition_phase0_output_00011_bpass_stellar_photon_ledger.csv"
    metadata_path = ledger_path.with_suffix(".json")
    if not ledger_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError(ledger_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    ledger = read_photon_source_ledger_csv(ledger_path)
    assert metadata["status"] == "candidate_bpass_stellar_photon_ledger"
    assert metadata["source_count"] == 42342
    assert ledger.source_id.shape == (42342,)
    assert ledger.photon_luminosity_s.shape == (42342, 9)
    assert np.isfinite(ledger.photon_luminosity_s).all()
    assert np.all(ledger.photon_luminosity_s >= 0.0)
    assert metadata["group_edges_ev"] == [0.01, 1.0, 5.6, 11.2, 13.6, 24.59, 54.42, 500.0, 2000.0, 10000.0]
    assert metadata["normalization"]["metallicity_floor_count"] == 338
    assert metadata["interpolation_clamped_sources"]["age"] == 178
    assert metadata["interpolation_clamped_sources"]["metallicity"] == 42004
    with ledger_path.open("r", encoding="utf-8", newline="") as handle:
        row_count = sum(1 for _ in csv.DictReader(handle))
    assert row_count == 42342
    print("BPASS_STELLAR_LEDGER_TEST_OK sources=42342 groups=9 candidate=true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
