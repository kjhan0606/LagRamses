#!/usr/bin/env python3
"""Tests for the G0 production manifest audit."""

from __future__ import annotations

import json
import hashlib
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

SNRT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SNRT_ROOT / "tools"))

from audit_production_manifest import audit_manifest  # noqa: E402


def write(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


def base_payloads() -> tuple[dict, dict, dict]:
    registry = {
        "schema": "lrd-jwst-external-assets",
        "assets": [
            {
                "id": "approved_binary",
                "path": "/tmp/approved_binary",
                "status": "available_local_approved",
                "sha256": "a" * 64,
                "simulation": {},
            }
        ],
    }
    production = {
        "schema": "lrd-jwst-production-readiness-manifest",
        "project_root": "/gpfs/kjhan/LRD_JWST",
        "repository": {
            "working_tree_at_recording": "clean",
            "production_clean_tree_required": True,
        },
        "asset_metadata_policy": {
            "require_per_asset_record": True,
            "required_fields": ["license_status", "provenance_status", "owner"],
            "accepted_license_statuses": ["verified"],
            "accepted_provenance_statuses": ["approved"],
        },
        "asset_metadata": {
            "approved_binary": {
                "license_status": "verified",
                "provenance_status": "approved",
                "owner": "test",
            }
        },
        "required_production_assets": [{"id": "approved_binary"}],
        "production_rules": {
            "reject_embedded_yield_fallback": True,
            "reject_legacy_and_transitional_assets": True,
        },
    }
    environment = {
        "schema": "snrt_environment_v1",
        "project_root": "/gpfs/kjhan/LRD_JWST",
        "backend": {"jax": "0.11.1", "jaxlib": "0.11.1"},
        "dependency_files": {"requirements": {"sha256": "b" * 64}},
        "reproducibility_status": "locked",
    }
    return registry, production, environment


def main() -> None:
    with TemporaryDirectory(prefix="g0-manifest-") as directory:
        root = Path(directory)
        registry, production, environment = base_payloads()
        paths = [root / "registry.json", root / "production.json", root / "environment.json"]
        asset_path = root / "approved_binary"
        asset_path.write_bytes(b"approved test asset\n")
        registry["assets"][0]["path"] = str(asset_path)
        lock_path = root / "requirements.lock.txt"
        lock_path.write_text("jax==0.11.1\njaxlib==0.11.1\n", encoding="utf-8")
        environment["dependency_files"]["resolved_lock"] = {
            "path": str(lock_path),
            "sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
        }
        for path, payload in zip(paths, (registry, production, environment)):
            write(path, payload)
        report = audit_manifest(*paths)
        assert report["production_gate_pass"] is True

        production["repository"]["working_tree_at_recording"] = "dirty"
        write(paths[1], production)
        report = audit_manifest(*paths)
        assert report["status"] == "blocked"
        assert "production_repository_dirty" in report["blocking_reasons"]

        registry["assets"][0]["status"] = "not_migrated"
        production["repository"]["working_tree_at_recording"] = "clean"
        write(paths[0], registry)
        write(paths[1], production)
        report = audit_manifest(*paths)
        assert "asset_status_blocked:approved_binary:not_migrated" in report["blocking_reasons"]

    print("PRODUCTION_MANIFEST_TEST_OK fail_closed=true")


if __name__ == "__main__":
    main()
