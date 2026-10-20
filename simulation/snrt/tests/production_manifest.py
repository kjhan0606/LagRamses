#!/usr/bin/env python3
"""Tests for the G0 production manifest audit."""

from __future__ import annotations

import json
import hashlib
import sys
import subprocess
import os
from unittest.mock import patch
from pathlib import Path
from tempfile import TemporaryDirectory

SNRT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SNRT_ROOT / "tools"))

from audit_production_manifest import audit_manifest, audit_native_release_candidate  # noqa: E402


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

        # Launch admission, not scientific validation: all approvals here are
        # synthetic test records. No real physical asset receives approval.
        production["repository"]["head_at_recording"] = source_head
        production["production_rules"]["production_status"] = "ready"
        registry["assets"][0]["sha256"] = hashlib.sha256(asset_path.read_bytes()).hexdigest()
        asset_path.chmod(0o755)
        namelist = root / "run.nml"
        namelist.write_text("&RUN_PARAMS\n hydro=.true.\n/\n", encoding="utf-8")
        native_lock = root / "native-build.lock"
        native_lock.write_text("compiler=TEST\nmpi=TEST\nSNRT=1\n", encoding="utf-8")
        environment["backend"] = {"runtime": "native_ramses", "compiler": "TEST", "mpi": "TEST"}
        environment["dependency_files"]["resolved_lock"] = {
            "path": str(native_lock), "sha256": hashlib.sha256(native_lock.read_bytes()).hexdigest(),
        }
        production["release_record"] = {
            "status": "approved", "approval_id": "TEST_ONLY",
            "executable_asset_id": "approved_binary",
            "namelist_sha256": hashlib.sha256(namelist.read_bytes()).hexdigest(),
            "environment_variables": {},
        }
        production["source_registry"] = str(paths[0])
        production["environment"] = str(paths[2])
        for path, payload in zip(paths, (registry, production, environment)):
            write(path, payload)

        def launch_report(**overrides):
            arguments = {"launch_root": repository_path, "binary": asset_path, "namelist": namelist}
            arguments.update(overrides)
            return audit_manifest(*paths, **arguments)

        with patch.dict(os.environ, {}, clear=True):
            report = launch_report()
            assert report["production_gate_pass"] is True
            assert report["launch_identity_checked"] is True
            cli = subprocess.run(
                [sys.executable, "-B", str(SNRT_ROOT / "tools/audit_production_manifest.py"),
                 "--production-manifest", str(paths[1]), "--launch-root", str(repository_path),
                 "--binary", str(asset_path), "--namelist", str(namelist)],
                capture_output=True, text=True, check=False,
            )
            assert cli.returncode == 0, cli.stderr + cli.stdout

            # An otherwise approved different checkout is not this launcher.
            assert "production_launch_root_mismatch" in launch_report(launch_root=root)["blocking_reasons"]
            production["production_rules"]["production_status"] = "blocked"
            write(paths[1], production)
            assert "physical_production_manifest_blocked" in launch_report()["blocking_reasons"]
            production["production_rules"]["production_status"] = "ready"
            production["release_record"]["approval_id"] = ""
            write(paths[1], production)
            assert "production_release_unapproved" in launch_report()["blocking_reasons"]
            production["release_record"]["approval_id"] = "TEST_ONLY"
            production["release_record"]["executable_asset_id"] = "optional_legacy_comparison"
            write(paths[1], production)
            assert "production_launch_binary_not_selected" in launch_report()["blocking_reasons"]
            production["release_record"]["executable_asset_id"] = "approved_binary"
            write(paths[1], production)
            namelist.write_text("changed namelist\n", encoding="utf-8")
            assert "production_launch_namelist_sha256_mismatch" in launch_report()["blocking_reasons"]
            namelist.write_text("&RUN_PARAMS\n hydro=.true.\n/\n", encoding="utf-8")
            asset_path.write_bytes(b"changed binary\n")
            assert "production_launch_binary_sha256_mismatch" in launch_report()["blocking_reasons"]
            asset_path.write_bytes(b"approved test asset\n")
            with patch.dict(os.environ, {"SNRT_STELLAR_SED": "changed-source"}):
                assert "production_launch_environment_undeclared:SNRT_STELLAR_SED" in launch_report()["blocking_reasons"]
            production["release_record"]["environment_variables"] = {"SNRT_STELLAR_SED": "approved-source"}
            write(paths[1], production)
            assert "production_launch_environment_mismatch:SNRT_STELLAR_SED" in launch_report()["blocking_reasons"]
            with patch.dict(os.environ, {"SNRT_STELLAR_SED": "approved-source"}):
                assert launch_report()["production_gate_pass"] is True
            production["release_record"]["environment_variables"] = {}
            production["repository"]["production_clean_tree_required"] = False
            write(paths[1], production)
            assert "production_launch_clean_source_policy_missing" in launch_report()["blocking_reasons"]
            production["repository"]["production_clean_tree_required"] = True
            production["asset_metadata_policy"]["require_per_asset_record"] = False
            write(paths[1], production)
            assert "production_launch_asset_metadata_policy_missing" in launch_report()["blocking_reasons"]
            production["asset_metadata_policy"]["require_per_asset_record"] = True
            environment["backend"].pop("mpi")
            write(paths[1], production)
            write(paths[2], environment)
            assert "native_environment_field_missing:mpi" in launch_report()["blocking_reasons"]

            # A mismatched launch never needs the other workspace's inventory.
            report = audit_manifest(root / "absent-registry.json", paths[1], root / "absent-environment.json",
                                    launch_root=root, binary=asset_path, namelist=namelist)
            assert report["blocking_reasons"] == ["production_launch_root_mismatch"]
            assert report["asset_count"] is None

            # Exercise the actual shell frontend in a tiny clean synthetic
            # checkout, without launching MPI, a simulation, or a Slurm job.
            environment["backend"]["mpi"] = "TEST"
            script = repository_path / "simulation/snrt/runs/simulation_ready_preflight.sh"
            tool = repository_path / "simulation/snrt/tools/audit_production_manifest.py"
            for destination, original in (
                (script, SNRT_ROOT / "runs/simulation_ready_preflight.sh"),
                (tool, SNRT_ROOT / "tools/audit_production_manifest.py"),
            ):
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(original.read_bytes())
            subprocess.run(["git", "-C", str(repository_path), "add", "simulation"], check=True)
            subprocess.run(["git", "-C", str(repository_path), "-c", "user.name=Manifest Test",
                            "-c", "user.email=manifest-test@example.invalid", "commit", "--quiet",
                            "-m", "synthetic admission source"], check=True)
            production["repository"]["head_at_recording"] = subprocess.run(
                ["git", "-C", str(repository_path), "rev-parse", "HEAD"],
                capture_output=True, text=True, check=True,
            ).stdout.strip()
            write(paths[1], production)
            write(paths[2], environment)
            command = ["bash", str(script), "--mode", "physical-production",
                       "--production-manifest", str(paths[1]), "--binary", str(asset_path),
                       "--namelist", str(namelist)]
            admitted = subprocess.run(command, capture_output=True, text=True, check=False)
            assert admitted.returncode == 0, admitted.stdout + admitted.stderr
            assert "SIMULATION_READY_PASS mode=physical-production" in admitted.stdout
            production["release_record"]["status"] = "pending"
            write(paths[1], production)
            rejected = subprocess.run(command, capture_output=True, text=True, check=False)
            assert rejected.returncode == 2
            assert "production_release_unapproved" in rejected.stdout

    # The current M5 candidate schema must admit only a committed, clean,
    # explicitly licensed release whose exact binary/input/evidence is used.
    with TemporaryDirectory(prefix="native-release-audit-") as directory:
        root = Path(directory)
        repository = root / "repository"
        source_head = init_repository(repository)
        subprocess.run(["git", "-C", str(repository), "remote", "add", "origin",
                        "git@github.com:kjhan0606/LagRamses.git"], check=True)
        binary = root / "native-binary"
        binary.write_bytes(b"native test executable\n")
        binary.chmod(0o755)
        namelist = root / "run.nml"
        namelist.write_text("&RUN_PARAMS\n hydro=.true.\n/\n", encoding="utf-8")
        inventory = root / "inputs.sha256"
        inventory.write_text("test inventory\n", encoding="utf-8")
        restart_artifact = root / "restart.json"
        restart_artifact.write_text("restart pass\n", encoding="utf-8")
        amr_artifact = root / "amr.json"
        amr_artifact.write_text("AMR pass\n", encoding="utf-8")
        digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
        candidate = {
            "schema": "lagramses-native-release-candidate-v1",
            "project_root": str(repository),
            "status": "approved_for_physical_production",
            "launch_admission_manifest": True,
            "repository": {
                "origin": "git@github.com:kjhan0606/LagRamses.git",
                "implementation_commit": source_head,
                "release_commit": source_head,
            },
            "implementation_snapshot": {
                "source.txt": hashlib.sha256(b"clean source\n").hexdigest(),
            },
            "executable": {"path": str(binary), "sha256": digest(binary)},
            "inputs": {
                "namelist": str(namelist), "namelist_sha256": digest(namelist),
                "complete_tested_input_inventory": str(inventory),
                "inventory_sha256": digest(inventory),
            },
            "environment_bindings": {},
            "source_contracts": [{
                "role": "test_source",
                "license_record_for_this_release": {
                    "license_status": "verified", "provenance_status": "recorded",
                },
            }],
            "selected_data_use_review": {"explicit_data_license": "verified"},
            "release_conditions_remaining": [],
            "verification": {
                "status": "M5_COUPLED_RESTART_PASS", "artifact": str(restart_artifact),
                "sha256": digest(restart_artifact),
            },
            "mixed_amr_verification": {
                "status": "M5_COUPLED_AMR_INTERFACE_PASS", "artifact": str(amr_artifact),
                "sha256": digest(amr_artifact),
            },
        }
        candidate_path = root / "candidate.json"
        write(candidate_path, candidate)
        with patch.dict(os.environ, {}, clear=True):
            report = audit_native_release_candidate(
                candidate_path, launch_root=repository, binary=binary, namelist=namelist,
            )
            assert report["production_gate_pass"] is True, report["blocking_reasons"]

            candidate["launch_admission_manifest"] = False
            candidate["source_contracts"][0]["license_record_for_this_release"]["license_status"] = "pending"
            candidate["release_conditions_remaining"] = ["upstream data terms unresolved"]
            write(candidate_path, candidate)
            report = audit_native_release_candidate(
                candidate_path, launch_root=repository, binary=binary, namelist=namelist,
            )
            assert "native_candidate_launch_admission_disabled" in report["blocking_reasons"]
            assert "native_source_license_unapproved:test_source" in report["blocking_reasons"]
            assert any(reason.startswith("native_release_condition_unresolved:")
                       for reason in report["blocking_reasons"])

    print("PRODUCTION_MANIFEST_TEST_OK fail_closed=true")


if __name__ == "__main__":
    main()
