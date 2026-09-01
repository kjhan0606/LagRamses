#!/usr/bin/env python3
"""Tests for the deterministic normalized-row yield converter."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from convert_yield_rows_to_canonical import ConversionError, convert  # noqa: E402


def _row(
    age: float, *, semantics: str = "cumulative", returned: float = 0.0,
    tracked: float | None = None
) -> dict:
    del semantics
    if tracked is None:
        tracked = returned
    return {
        "channel": 1,
        "initial_mass_msun_per_star": 1.0,
        "birth_metallicity_mass_fraction": 0.001,
        "age_yr": age,
        "returned_mass_msun_per_star": returned,
        "remnant_mass_msun_per_star": 0.0,
        "energy_erg_per_star": 0.0,
        "momentum_g_cm_s_per_star": [0.0, 0.0, 0.0],
        "ejecta_msun_per_star": [tracked] + [0.0] * 10,
        "net_yield_msun_per_star": [0.0] * 11,
    }


def _source(**overrides: object) -> dict:
    source = {
        "citation": "test source",
        "source_version": "test-v1",
        "source_sha256": "a" * 64,
        "release_history_semantics": "cumulative",
        "approval_id": "TEST-G2-001",
        "license_status": "approved",
        "provenance_status": "approved",
        "units": "canonical",
        "IMF": "Kroupa",
        "population_model": "single_star_ssp",
        "channel_boundaries": {"1": [0.8, 120.0]},
        "metallicity_definition": "mass fraction",
        "solar_abundance_set": "test",
        "remnant_model": "test",
        "untracked_ejecta_policy": (
            "returned_mass_minus_sum_tracked_ejecta_deposited_as_generic_metals"
        ),
    }
    source.update(overrides)
    return source


def main() -> int:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source_json = root / "source.json"
        output = root / "yield.dat"
        sidecar = root / "yield.dat.json"
        source_json.write_text(
            json.dumps(
                {
                    "source": _source(),
                    "rows": [_row(1.0, returned=0.1, tracked=0.075), _row(0.0)],
                }
            ),
            encoding="utf-8",
        )
        metadata = convert(source_json, output, sidecar)
        lines = [line for line in output.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
        assert len(lines) == 2
        assert float(lines[0].split()[3]) == 0.0
        assert metadata["asset_sha256"] == hashlib.sha256(output.read_bytes()).hexdigest()
        sidecar_data = json.loads(sidecar.read_text(encoding="utf-8"))
        assert sidecar_data["conversion_code_sha256"]
        assert sidecar_data["rows_with_untracked_ejecta"] == 1
        assert math.isclose(
            sidecar_data["maximum_untracked_ejecta_msun_per_star"], 0.025
        )

        overfull_json = root / "overfull.json"
        overfull_json.write_text(
            json.dumps(
                {"source": _source(), "rows": [_row(0.0, returned=0.1, tracked=0.2)]}
            ),
            encoding="utf-8",
        )
        try:
            convert(overfull_json, root / "overfull.dat", root / "overfull.dat.json")
        except ConversionError as exc:
            assert "tracked ejecta exceeding" in str(exc)
        else:
            raise AssertionError("overfull tracked ejecta were not rejected")

        rate_json = root / "rate.json"
        rate_json.write_text(
            json.dumps({"source": _source(release_history_semantics="rate"), "rows": [_row(0.0)]}),
            encoding="utf-8",
        )
        try:
            convert(rate_json, root / "rate.dat", root / "rate.dat.json")
        except ConversionError as exc:
            assert "rate tables" in str(exc)
        else:
            raise AssertionError("rate input was not rejected")
    print("YIELD_CONVERTER_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
