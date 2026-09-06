#!/usr/bin/env python3
"""Tests for the manifest-scoped G2 candidate package fingerprints."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_source_package_fingerprints import audit_fingerprints  # noqa: E402


EXPECTED = {
    "limongi_chieffi_2018_cds": "01873ea4ca2746b5f4a44941237a711777e6e35b3f8e1ced65b524e27bfdb228",
    "nugrid_set1ext_mesaonly_fryer12_delay": "fd0cc5d2dced86807756fc4c39b94ed1af221aaf8434c4778fe6a6ecab0f936d",
    "huscher2025_agb": "39cfcf0477039954c6a2d20ba40c49f9d97163a1a747c1e236316c6678224b46",
    "boccioli_roberti2026_neutrino_ccsn": "3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b",
    "doherty2014_sagb": "ff4dda0d45f61e0c8342379a65afef857c6d66a4c5671fdaa860eb007912aae6",
    "stockinger2020_low_mass_ccsn": "84960add39cf27ae97ce27655a5f162c6441deedf2bfdb8598881c9539f27731",
    "sukhbold2016_ccsn": "11e2988969f69da9fb3f6939af9bbd3f0f26afd4420cd676d27a05b30486aa60",
    "limongi2024_transition_fates": "86442e0fa7e65482fa0476b9f0da5c895b98672c95ca35c9944340f2a9e46d85",
    "roberti2024_ultralowz_ccsn": "3085e6ad9044e6f35a37cb2eb95a06da0cfb9f79361eb175d3d7ffbd3f3699e2",
    "heger_woosley2010_popiii": "c7a352d9f162a5378b5d171f1d57423b6c0f1c388e04526c996767ef87ba3d5e",
    "nubase2020_decay_projection_data": "71184c8651296b4703bc49f2bc62c0b849490d3a5541c150ede7c2660dbee27c",
}


def _assert_clean_report() -> None:
    root = ROOT.parents[1] / "external" / "g2_candidates"
    report = audit_fingerprints(root)
    assert report["status"] == "candidate_fingerprint_review_only"
    assert report["production_ready"] is False
    assert report["input_integrity_passed"] is True
    assert report["candidate_count"] == 11
    assert report["file_count"] == 65
    observed = {
        candidate["candidate_id"]: candidate["composite_sha256"]
        for candidate in report["candidates"]
    }
    assert observed == EXPECTED
    assert all(candidate["input_integrity_passed"] for candidate in report["candidates"])


def _assert_manifest_mutation_is_fatal() -> None:
    project = ROOT.parents[1]
    source_root = project / "external" / "g2_candidates"
    manifest = source_root / "acquisition_manifest_v1.json"
    with tempfile.TemporaryDirectory(prefix="snrt-g2-fingerprint-") as directory:
        mutated = Path(directory) / manifest.name
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        payload["candidates"][0]["files"][0]["sha256"] = "0" * 64
        mutated.write_text(json.dumps(payload), encoding="utf-8")
        report = audit_fingerprints(source_root, mutated)
        assert report["status"] == "candidate_fingerprint_blocked_input_integrity"
        assert report["input_integrity_passed"] is False
        assert report["audit_failures"]
        result = subprocess.run(
            [
                sys.executable,
                str(TOOLS / "audit_g2_source_package_fingerprints.py"),
                "--root",
                str(source_root),
                "--manifest",
                str(mutated),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 2
        cli_report = json.loads(result.stdout)
        assert cli_report["status"] == "candidate_fingerprint_blocked_input_integrity"


def _assert_unsafe_path_is_fatal() -> None:
    project = ROOT.parents[1]
    source_root = project / "external" / "g2_candidates"
    manifest = source_root / "acquisition_manifest_v1.json"
    with tempfile.TemporaryDirectory(prefix="snrt-g2-fingerprint-path-") as directory:
        mutated = Path(directory) / manifest.name
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        payload["candidates"][0]["files"][0]["path"] = "../outside"
        mutated.write_text(json.dumps(payload), encoding="utf-8")
        report = audit_fingerprints(source_root, mutated)
        assert report["status"] == "candidate_fingerprint_blocked_input_integrity"
        assert any(failure["reason"] == "unsafe_file_path" for failure in report["audit_failures"])


def main() -> int:
    _assert_clean_report()
    _assert_manifest_mutation_is_fatal()
    _assert_unsafe_path_is_fatal()
    print("G2_SOURCE_PACKAGE_FINGERPRINT_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
