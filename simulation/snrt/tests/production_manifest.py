#!/usr/bin/env python3
"""Tests for the G0 production manifest audit."""

from __future__ import annotations

import json
import hashlib
import sys
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

SNRT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SNRT_ROOT / "tools"))

from audit_production_manifest import audit_manifest  # noqa: E402


def write(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


def init_repository(path: Path) -> str:
    path.mkdir()
    subprocess.run(["git", "-C", str(path), "init", "--quiet"], check=True)
    (path / "source.txt").write_text("clean source\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(path), "add", "source.txt"], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(path),
            "-c",
            "user.name=Manifest Test",
            "-c",
            "user.email=manifest-test@example.invalid",
            "commit",
            "--quiet",
            "-m",
            "initial source",
        ],
        check=True,
    )
    return subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def base_payloads(repository_path: Path, head: str) -> tuple[dict, dict, dict]:
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
        "project_root": str(repository_path),
        "repository": {
            "path": str(repository_path),
            "head_at_recording": head,
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
        "project_root": str(repository_path),
        "backend": {"jax": "0.11.1", "jaxlib": "0.11.1"},
        "dependency_files": {"requirements": {"sha256": "b" * 64}},
        "reproducibility_status": "locked",
    }
    return registry, production, environment


def main() -> None:
    with TemporaryDirectory(prefix="g0-manifest-") as directory:
        root = Path(directory)
        repository_path = root / "source"
        source_head = init_repository(repository_path)
        registry, production, environment = base_payloads(repository_path, source_head)
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
        registry["assets"].append(
            {
                "id": "optional_legacy_comparison",
                "path": str(root / "legacy_comparison"),
                "status": "not_migrated",
                "sha256": "c" * 64,
                "simulation": {"comparison_role": "legacy_comparison_only"},
            }
        )
        for path, payload in zip(paths, (registry, production, environment)):
            write(path, payload)
        report = audit_manifest(*paths)
        # The external registry is a project inventory; unselected controls do
        # not need production approval metadata or block the selected runtime.
        assert report["production_gate_pass"] is True

        registry["assets"][0]["simulation"] = {"comparison_role": "legacy_comparison_only"}
        write(paths[0], registry)
        report = audit_manifest(*paths)
        assert "comparison_asset_selected_for_production:approved_binary" in report["blocking_reasons"]
        registry["assets"][0]["simulation"] = {}

        production["repository"]["working_tree_at_recording"] = "dirty"
        write(paths[0], registry)
        write(paths[1], production)
        report = audit_manifest(*paths)
        assert report["status"] == "blocked"
        assert "production_repository_dirty" in report["blocking_reasons"]

        registry["assets"][0]["status"] = "not_migrated"
        production["repository"]["working_tree_at_recording"] = "clean"
        write(paths[0], registry)
        write(paths[1], production)
        report = audit_manifest(*paths)
        assert "required_asset_status_blocked:approved_binary:not_migrated" in report["blocking_reasons"]

        registry["assets"][0]["status"] = "available_local_approved"
        production["asset_metadata"] = {}
        write(paths[0], registry)
        write(paths[1], production)
        report = audit_manifest(*paths)
        assert "asset_metadata_missing:approved_binary" in report["blocking_reasons"]

        production["asset_metadata"] = {
            "approved_binary": {
                "license_status": "verified",
                "provenance_status": "approved",
                "owner": "test",
            }
        }
        production["required_production_assets"].append({"id": "absent_asset"})
        write(paths[1], production)
        report = audit_manifest(*paths)
        assert "required_asset_not_registered:absent_asset" in report["blocking_reasons"]
        assert "asset_metadata_missing:absent_asset" not in report["blocking_reasons"]

        (repository_path / "source.txt").write_text("dirty source\n", encoding="utf-8")
        production["required_production_assets"] = [{"id": "approved_binary"}]
        write(paths[1], production)
        report = audit_manifest(*paths)
        assert "production_repository_dirty" in report["blocking_reasons"]

        (repository_path / "source.txt").write_text("clean source\n", encoding="utf-8")
        production["repository"]["head_at_recording"] = "0" * 40
        write(paths[1], production)
        report = audit_manifest(*paths)
        assert "production_source_commit_mismatch" in report["blocking_reasons"]

    print("PRODUCTION_MANIFEST_TEST_OK fail_closed=true")


if __name__ == "__main__":
    main()
