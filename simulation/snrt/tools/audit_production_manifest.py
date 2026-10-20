#!/usr/bin/env python3
"""Fail-closed audit for the G0 production-readiness manifest.

The migration registry is intentionally left untouched.  This tool audits the
production overlay against that registry and an environment fingerprint.  It
does not copy, hash large files, edit manifests, or make missing assets appear
available.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


BLOCKING_STATUSES = {
    "missing",
    "not_migrated",
    "hash_pending",
    "available_external_dirty",
}


def _load(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"JSON object required: {path}")
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _audit_repository_identity(
    repository: dict[str, Any],
    production_root_value: Any,
    environment_root_value: Any,
    blockers: list[str],
) -> None:
    repository_path_value = repository.get("path")
    if not isinstance(repository_path_value, str) or not Path(repository_path_value).is_absolute():
        blockers.append("production_repository_path_invalid")
        return
    repository_path = Path(repository_path_value)
    try:
        repository_path = repository_path.resolve(strict=True)
    except OSError:
        blockers.append("production_repository_path_missing")
        return
    if not repository_path.is_dir():
        blockers.append("production_repository_path_missing")
        return

    for value, reason in (
        (production_root_value, "production_root_invalid"),
        (environment_root_value, "environment_root_invalid"),
    ):
        if not isinstance(value, str) or not Path(value).is_absolute():
            blockers.append(reason)
            continue
        try:
            if Path(value).resolve(strict=True) != repository_path:
                blockers.append(reason)
        except OSError:
            blockers.append(reason)

    def git_output(*args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", "-C", str(repository_path), *args],
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )

    try:
        top = git_output("rev-parse", "--show-toplevel")
        head = git_output("rev-parse", "HEAD")
        status = git_output("status", "--porcelain", "--untracked-files=normal")
    except (OSError, subprocess.TimeoutExpired):
        blockers.append("production_repository_probe_failed")
        return
    if top.returncode != 0 or head.returncode != 0 or status.returncode != 0:
        blockers.append("production_repository_probe_failed")
        return
    try:
        actual_top = Path(top.stdout.strip()).resolve(strict=True)
    except OSError:
        blockers.append("production_repository_root_mismatch")
        return
    if actual_top != repository_path:
        blockers.append("production_repository_root_mismatch")

    expected_head = repository.get("head_at_recording")
    actual_head = head.stdout.strip().lower()
    if not isinstance(expected_head, str) or re.fullmatch(r"[0-9a-fA-F]{40}", expected_head) is None:
        blockers.append("production_source_commit_invalid")
    elif actual_head != expected_head.lower():
        blockers.append("production_source_commit_mismatch")

    if repository.get("production_clean_tree_required") and (
        repository.get("working_tree_at_recording") != "clean" or status.stdout.strip()
    ):
        blockers.append("production_repository_dirty")


def audit_manifest(
    registry_path: Path,
    production_path: Path,
    environment_path: Path,
    *,
    launch_root: Path | None = None,
    binary: Path | None = None,
    namelist: Path | None = None,
) -> dict[str, Any]:
    production = _load(production_path)
    blockers: list[str] = []
    warnings: list[str] = []

    if production.get("schema") != "lrd-jwst-production-readiness-manifest":
        blockers.append("production_schema_invalid")
    repository = production.get("repository", {})
    if not isinstance(repository, dict):
        blockers.append("production_repository_record_invalid")
        repository = {}
    launch_root_matches = True
    if launch_root is not None:
        try:
            launch_root_matches = launch_root.resolve(strict=True) == Path(repository.get("path", "")).resolve(strict=True)
        except (OSError, TypeError, ValueError):
            launch_root_matches = False
    if not launch_root_matches:
        # Do not probe Git in another worktree when a launcher supplied its
        # actual checkout, or recast its historical assets as current blockers.
        blockers.append("production_launch_root_mismatch")
        declared = production.get("required_production_assets", [])
        return {
            "schema": "snrt_production_manifest_audit_v1",
            "production_manifest": str(production_path),
            "registry": str(registry_path), "environment": str(environment_path),
            "asset_count": None,
            "required_asset_ids": [r["id"] for r in declared if isinstance(r, dict) and isinstance(r.get("id"), str)] if isinstance(declared, list) else [],
            "status": "blocked", "production_gate_pass": False,
            "blocking_reasons": sorted(set(blockers)),
            "warnings": ["historical_assets_not_audited_for_current_launch"],
            "launch_identity_checked": True,
        }
    registry = _load(registry_path)
    environment = _load(environment_path)
    if registry.get("schema") != "lrd-jwst-external-assets":
        blockers.append("registry_schema_invalid")
    if environment.get("schema") != "snrt_environment_v1":
        blockers.append("environment_schema_invalid")
    _audit_repository_identity(
        repository, production.get("project_root"), environment.get("project_root"), blockers,
    )

    assets = registry.get("assets", [])
    if not isinstance(assets, list) or not assets:
        blockers.append("registry_assets_missing")
        assets = []

    required_assets = production.get("required_production_assets", [])
    required_ids: list[str] = []
    if not isinstance(required_assets, list):
        blockers.append("required_production_assets_invalid")
        required_assets = []
    for requirement in required_assets:
        if not isinstance(requirement, dict) or not isinstance(requirement.get("id"), str):
            blockers.append("required_asset_record_invalid")
            continue
        required_ids.append(requirement["id"])
    required_id_set = set(required_ids)

    asset_by_id: dict[str, dict[str, Any]] = {}
    duplicate_ids: list[str] = []
    for asset in assets:
        if not isinstance(asset, dict) or not isinstance(asset.get("id"), str):
            blockers.append("asset_record_invalid")
            continue
        asset_id = asset["id"]
        if asset_id in asset_by_id:
            duplicate_ids.append(asset_id)
        asset_by_id[asset_id] = asset
        path_value = asset.get("path")
        if not isinstance(path_value, str) or not Path(path_value).is_absolute():
            blockers.append(f"asset_path_not_absolute:{asset_id}")
        if asset_id in required_id_set:
            status = asset.get("status")
            if status in BLOCKING_STATUSES:
                blockers.append(f"required_asset_status_blocked:{asset_id}:{status}")
            simulation = asset.get("simulation")
            comparison_role = simulation.get("comparison_role") if isinstance(simulation, dict) else None
            if comparison_role in {
                "legacy_comparison_only",
                "phase0_preselector_transitional",
            }:
                blockers.append(f"comparison_asset_selected_for_production:{asset_id}")
    if duplicate_ids:
        blockers.extend(f"duplicate_asset_id:{asset_id}" for asset_id in duplicate_ids)

    policy = production.get("asset_metadata_policy", {})
    metadata = production.get("asset_metadata", {})
    required_fields = policy.get("required_fields", [])
    if policy.get("require_per_asset_record"):
        if not isinstance(metadata, dict):
            blockers.append("asset_metadata_not_object")
            metadata = {}
        for asset_id in required_ids:
            if asset_id not in asset_by_id:
                # The missing registry record is reported by the required-
                # asset gate below; metadata cannot meaningfully be audited.
                continue
            record = metadata.get(asset_id)
            if not isinstance(record, dict):
                blockers.append(f"asset_metadata_missing:{asset_id}")
                continue
            for field in required_fields:
                if not record.get(field):
                    blockers.append(f"asset_metadata_field_missing:{asset_id}:{field}")
            if record.get("license_status") not in policy.get("accepted_license_statuses", []):
                blockers.append(f"asset_license_unapproved:{asset_id}")
            if record.get("provenance_status") not in policy.get("accepted_provenance_statuses", []):
                blockers.append(f"asset_provenance_unapproved:{asset_id}")

    for requirement in required_assets:
        if not isinstance(requirement, dict) or not isinstance(requirement.get("id"), str):
            continue
        asset_id = requirement["id"]
        asset = asset_by_id.get(asset_id)
        if asset is None:
            blockers.append(f"required_asset_not_registered:{asset_id}")
            continue
        asset_path = asset.get("path")
        if not isinstance(asset_path, str) or not Path(asset_path).exists():
            blockers.append(f"required_asset_path_missing:{asset_id}")
        asset_sha256 = asset.get("sha256")
        if not isinstance(asset_sha256, str) or re.fullmatch(r"[0-9a-fA-F]{64}", asset_sha256) is None:
            blockers.append(f"required_asset_sha256_missing:{asset_id}")

    environment_backend = environment.get("backend", {})
    if not isinstance(environment_backend, dict):
        environment_backend = {}
        blockers.append("environment_backend_invalid")
    runtime = environment_backend.get("runtime", "jax")
    if runtime == "native_ramses":
        # Native Fortran/C/CUDA execution does not depend on the P0 JAX stack.
        # Its exact build/runtime recipe still requires the existing hash lock.
        for field in ("compiler", "mpi"):
            if not isinstance(environment_backend.get(field), str) or not environment_backend[field].strip():
                blockers.append(f"native_environment_field_missing:{field}")
    elif runtime == "jax":
        if environment_backend.get("jax") != "0.11.1" or environment_backend.get("jaxlib") != "0.11.1":
            blockers.append("jax_cpu_environment_version_mismatch")
    else:
        blockers.append("environment_runtime_unsupported")
    dependency_files = environment.get("dependency_files")
    if not isinstance(dependency_files, dict) or not dependency_files:
        blockers.append("environment_dependency_files_missing")
        dependency_files = {}
    resolved_lock = dependency_files.get("resolved_lock")
    if not isinstance(resolved_lock, dict):
        blockers.append("environment_resolved_lock_missing")
    else:
        lock_path_value = resolved_lock.get("path")
        lock_sha256 = resolved_lock.get("sha256")
        if not isinstance(lock_path_value, str) or not Path(lock_path_value).is_file():
            blockers.append("environment_resolved_lock_path_missing")
        elif not isinstance(lock_sha256, str) or len(lock_sha256) != 64:
            blockers.append("environment_resolved_lock_sha256_missing")
        elif _sha256(Path(lock_path_value)) != lock_sha256:
            blockers.append("environment_resolved_lock_sha256_mismatch")
    if environment.get("reproducibility_status") != "locked":
        blockers.append("environment_not_locked")

    if not production.get("production_rules", {}).get("reject_embedded_yield_fallback"):
        blockers.append("embedded_fallback_rejection_not_enabled")
    if not production.get("production_rules", {}).get("reject_legacy_and_transitional_assets"):
        blockers.append("legacy_transitional_rejection_not_enabled")

    launch_requested = any(value is not None for value in (launch_root, binary, namelist))
    if launch_requested:
        if launch_root is None or binary is None or namelist is None:
            blockers.append("production_launch_arguments_incomplete")
        else:
            if repository.get("production_clean_tree_required") is not True:
                blockers.append("production_launch_clean_source_policy_missing")
            if policy.get("require_per_asset_record") is not True:
                blockers.append("production_launch_asset_metadata_policy_missing")
            for asset_id in required_ids:
                if asset_id not in asset_by_id:
                    continue
                record = metadata.get(asset_id) if isinstance(metadata, dict) else None
                if not isinstance(record, dict):
                    continue  # Existing per-asset admission reports its absence.
                if not record.get("owner"):
                    blockers.append(f"asset_metadata_field_missing:{asset_id}:owner")
                if record.get("license_status") not in {"verified", "approved"}:
                    blockers.append(f"asset_license_unapproved:{asset_id}")
                if record.get("provenance_status") not in {"recorded", "approved"}:
                    blockers.append(f"asset_provenance_unapproved:{asset_id}")
            if production.get("production_rules", {}).get("production_status") != "ready":
                blockers.append("physical_production_manifest_blocked")
            release = production.get("release_record")
            if not isinstance(release, dict):
                blockers.append("production_release_record_missing")
            else:
                if release.get("status") != "approved" or not isinstance(release.get("approval_id"), str) or not release["approval_id"].strip():
                    blockers.append("production_release_unapproved")
                asset_id = release.get("executable_asset_id")
                asset = asset_by_id.get(asset_id) if isinstance(asset_id, str) else None
                if asset is None or asset_id not in required_id_set:
                    blockers.append("production_launch_binary_not_selected")
                else:
                    try:
                        if binary.resolve(strict=True) != Path(asset["path"]).resolve(strict=True) or not binary.is_file() or not os.access(binary, os.X_OK):
                            blockers.append("production_launch_binary_path_mismatch")
                        elif _sha256(binary) != asset.get("sha256"):
                            blockers.append("production_launch_binary_sha256_mismatch")
                    except (OSError, KeyError, TypeError, ValueError):
                        blockers.append("production_launch_binary_path_mismatch")
                try:
                    if not namelist.is_file() or _sha256(namelist) != release.get("namelist_sha256"):
                        blockers.append("production_launch_namelist_sha256_mismatch")
                except OSError:
                    blockers.append("production_launch_namelist_sha256_mismatch")
                variables = release.get("environment_variables")
                if not isinstance(variables, dict) or any(not isinstance(k, str) or not isinstance(v, str) for k, v in variables.items()):
                    blockers.append("production_launch_environment_record_invalid")
                else:
                    for key, expected in variables.items():
                        if os.environ.get(key, "") != expected:
                            blockers.append(f"production_launch_environment_mismatch:{key}")
                    # Source contracts and reference-control switches may be
                    # selected through the environment, not only namelist text.
                    for key, value in os.environ.items():
                        if value and key.startswith(("SNRT_", "PHASE0_", "CHIMES_")) and key not in variables:
                            blockers.append(f"production_launch_environment_undeclared:{key}")

    if not blockers:
        status = "pass"
    else:
        status = "blocked"
    return {
        "schema": "snrt_production_manifest_audit_v1",
        "production_manifest": str(production_path),
        "registry": str(registry_path),
        "environment": str(environment_path),
        "asset_count": len(asset_by_id),
        "required_asset_ids": required_ids,
        "status": status,
        "production_gate_pass": status == "pass",
        "blocking_reasons": sorted(set(blockers)),
        "warnings": sorted(set(warnings)),
        "launch_identity_checked": launch_requested,
    }


def audit_native_release_candidate(
    manifest_path: Path,
    *,
    launch_root: Path | None = None,
    binary: Path | None = None,
    namelist: Path | None = None,
) -> dict[str, Any]:
    """Fail closed on the current LagRamses native M5 candidate schema."""
    candidate = _load(manifest_path)
    blockers: list[str] = []
    if candidate.get("schema") != "lagramses-native-release-candidate-v1":
        blockers.append("native_candidate_schema_invalid")

    root_value = candidate.get("project_root")
    try:
        root = Path(root_value).resolve(strict=True) if isinstance(root_value, str) else None
    except OSError:
        root = None
    try:
        actual_root = launch_root.resolve(strict=True) if launch_root is not None else None
    except OSError:
        actual_root = None
    if root is None or actual_root is None or root != actual_root:
        blockers.append("native_candidate_launch_root_mismatch")
        return {
            "schema": "snrt_native_release_candidate_audit_v1",
            "manifest": str(manifest_path),
            "status": "blocked",
            "production_gate_pass": False,
            "blocking_reasons": blockers,
            "launch_identity_checked": any(value is not None for value in (launch_root, binary, namelist)),
        }

    repository = candidate.get("repository", {})
    if not isinstance(repository, dict):
        repository = {}
        blockers.append("native_candidate_repository_record_invalid")
    try:
        top = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
            capture_output=True, text=True, check=False, timeout=30,
        )
        head = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=False, timeout=30,
        )
        status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=normal"],
            capture_output=True, text=True, check=False, timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        blockers.append("native_candidate_repository_probe_failed")
        top = head = status = None
    if top is not None and head is not None and status is not None:
        if top.returncode or head.returncode or status.returncode:
            blockers.append("native_candidate_repository_probe_failed")
        else:
            if Path(top.stdout.strip()).resolve() != root:
                blockers.append("native_candidate_repository_root_mismatch")
            expected_origin = repository.get("origin")
            if isinstance(expected_origin, str):
                origin = subprocess.run(
                    ["git", "-C", str(root), "remote", "get-url", "origin"],
                    capture_output=True, text=True, check=False, timeout=30,
                )
                if origin.returncode or origin.stdout.strip() != expected_origin:
                    blockers.append("native_candidate_repository_origin_mismatch")
            else:
                blockers.append("native_candidate_repository_origin_missing")
            actual_head = head.stdout.strip().lower()
            for field, reason in (
                ("implementation_commit", "native_candidate_implementation_commit_invalid"),
                ("release_commit", "native_candidate_release_commit_missing"),
            ):
                expected = repository.get(field)
                if not isinstance(expected, str) or re.fullmatch(r"[0-9a-fA-F]{40}", expected) is None:
                    blockers.append(reason)
                    continue
                ancestor = subprocess.run(
                    ["git", "-C", str(root), "merge-base", "--is-ancestor", expected, actual_head],
                    capture_output=True, text=True, check=False, timeout=30,
                )
                if ancestor.returncode:
                    blockers.append(f"{field}_not_in_current_history")
            if status.stdout.strip():
                blockers.append("native_candidate_repository_dirty")

    snapshot = candidate.get("implementation_snapshot")
    if not isinstance(snapshot, dict) or not snapshot:
        blockers.append("native_candidate_source_snapshot_missing")
    else:
        for relative_path, expected_hash in snapshot.items():
            source = root / relative_path if isinstance(relative_path, str) else root
            if not isinstance(expected_hash, str) or re.fullmatch(r"[0-9a-fA-F]{64}", expected_hash) is None:
                blockers.append(f"native_candidate_source_hash_invalid:{relative_path}")
            elif not source.is_file() or _sha256(source) != expected_hash:
                blockers.append(f"native_candidate_source_snapshot_mismatch:{relative_path}")

    if candidate.get("launch_admission_manifest") is not True:
        blockers.append("native_candidate_launch_admission_disabled")
    if candidate.get("status") != "approved_for_physical_production":
        blockers.append("native_candidate_physical_release_not_approved")
    conditions = candidate.get("release_conditions_remaining")
    if not isinstance(conditions, list):
        blockers.append("native_candidate_release_conditions_invalid")
    elif conditions:
        blockers.extend(f"native_release_condition_unresolved:{i + 1}" for i in range(len(conditions)))

    # Candidate-use approval does not certify upstream data terms. Require a
    # per-source release license/provenance record and reject unresolved
    # license/terms declarations anywhere in the selected-data review.
    contracts = candidate.get("source_contracts")
    if not isinstance(contracts, list) or not contracts:
        blockers.append("native_candidate_source_contracts_missing")
    else:
        for contract in contracts:
            if not isinstance(contract, dict):
                blockers.append("native_candidate_source_contract_invalid")
                continue
            role = contract.get("role", "unknown")
            record = contract.get("license_record_for_this_release")
            if not isinstance(record, dict):
                blockers.append(f"native_source_license_record_missing:{role}")
                continue
            if record.get("license_status") not in {"verified", "approved"}:
                blockers.append(f"native_source_license_unapproved:{role}")
            if record.get("provenance_status") not in {"recorded", "approved"}:
                blockers.append(f"native_source_provenance_unapproved:{role}")

    def check_license_text(value: Any, path: str = "selected_data_use_review") -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                child_path = f"{path}.{key}"
                if ("license" in key.lower() or "terms" in key.lower()) and isinstance(child, str):
                    lowered = child.lower()
                    if any(marker in lowered for marker in ("unresolved", "pending", "not located", "not established")):
                        blockers.append(f"native_data_terms_unresolved:{child_path}")
                else:
                    check_license_text(child, child_path)
        elif isinstance(value, list):
            for index, child in enumerate(value):
                check_license_text(child, f"{path}[{index}]")

    check_license_text(candidate.get("selected_data_use_review", {}))

    executable = candidate.get("executable", {})
    if not isinstance(executable, dict):
        executable = {}
        blockers.append("native_candidate_executable_record_invalid")
    expected_binary = executable.get("path")
    expected_binary_hash = executable.get("sha256")
    if binary is None:
        blockers.append("native_candidate_launch_binary_missing")
    else:
        try:
            if not binary.is_file() or not os.access(binary, os.X_OK):
                blockers.append("native_candidate_launch_binary_unavailable")
            elif not isinstance(expected_binary, str) or binary.resolve(strict=True) != Path(expected_binary).resolve(strict=True):
                blockers.append("native_candidate_launch_binary_path_mismatch")
            elif not isinstance(expected_binary_hash, str) or _sha256(binary) != expected_binary_hash:
                blockers.append("native_candidate_launch_binary_sha256_mismatch")
        except (OSError, TypeError, ValueError):
            blockers.append("native_candidate_launch_binary_unavailable")

    inputs = candidate.get("inputs", {})
    if not isinstance(inputs, dict):
        inputs = {}
        blockers.append("native_candidate_inputs_invalid")
    declared_namelists = []
    for path_key, hash_key in (("namelist", "namelist_sha256"), ("restart_namelist", "restart_namelist_sha256")):
        path_value = inputs.get(path_key)
        digest = inputs.get(hash_key)
        if isinstance(path_value, str) and isinstance(digest, str):
            declared_namelists.append((Path(path_value), digest))
    if namelist is None or not namelist.is_file():
        blockers.append("native_candidate_launch_namelist_unavailable")
    else:
        try:
            if not any(namelist.resolve(strict=True) == path.resolve(strict=True) and _sha256(namelist) == digest
                       for path, digest in declared_namelists):
                blockers.append("native_candidate_launch_namelist_identity_mismatch")
        except OSError:
            blockers.append("native_candidate_launch_namelist_unavailable")

    inventory_path = inputs.get("complete_tested_input_inventory")
    inventory_hash = inputs.get("inventory_sha256")
    if not isinstance(inventory_path, str) or not isinstance(inventory_hash, str):
        blockers.append("native_candidate_input_inventory_missing")
    else:
        inventory = Path(inventory_path)
        if not inventory.is_file():
            blockers.append("native_candidate_input_inventory_unavailable")
        elif _sha256(inventory) != inventory_hash:
            blockers.append("native_candidate_input_inventory_sha256_mismatch")

    environment = candidate.get("environment_bindings", {})
    if not isinstance(environment, dict) or any(not isinstance(k, str) or not isinstance(v, str) for k, v in environment.items()):
        blockers.append("native_candidate_environment_record_invalid")
    else:
        for key, expected in environment.items():
            if os.environ.get(key, "") != expected:
                blockers.append(f"native_candidate_environment_mismatch:{key}")
        for key, value in os.environ.items():
            if value and key.startswith(("SNRT_", "PHASE0_", "CHIMES_")) and key not in environment:
                blockers.append(f"native_candidate_environment_undeclared:{key}")

    for section, expected_status in (
        ("verification", "M5_COUPLED_RESTART_PASS"),
        ("mixed_amr_verification", "M5_COUPLED_AMR_INTERFACE_PASS"),
    ):
        record = candidate.get(section, {})
        if not isinstance(record, dict) or record.get("status") != expected_status:
            blockers.append(f"native_candidate_verification_missing:{section}")
            continue
        artifact_value = record.get("artifact")
        artifact_hash = record.get("sha256", record.get("artifact_sha256"))
        if not isinstance(artifact_value, str) or not isinstance(artifact_hash, str):
            blockers.append(f"native_candidate_verification_artifact_missing:{section}")
            continue
        artifact = Path(artifact_value)
        if not artifact.is_file() or _sha256(artifact) != artifact_hash:
            blockers.append(f"native_candidate_verification_artifact_mismatch:{section}")

    return {
        "schema": "snrt_native_release_candidate_audit_v1",
        "manifest": str(manifest_path),
        "status": "pass" if not blockers else "blocked",
        "production_gate_pass": not blockers,
        "blocking_reasons": sorted(set(blockers)),
        "launch_identity_checked": all(value is not None for value in (launch_root, binary, namelist)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path)
    parser.add_argument("--production-manifest", type=Path, required=True)
    parser.add_argument("--environment", type=Path)
    parser.add_argument("--launch-root", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--namelist", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    try:
        production = _load(args.production_manifest)
        if production.get("schema") == "lagramses-native-release-candidate-v1":
            report = audit_native_release_candidate(
                args.production_manifest, launch_root=args.launch_root,
                binary=args.binary, namelist=args.namelist,
            )
        else:
            def input_path(explicit: Path | None, key: str) -> Path:
                if explicit is not None:
                    return explicit
                value = production.get(key)
                if not isinstance(value, str) or not Path(value).is_absolute():
                    raise ValueError(f"absolute {key} path required in production manifest")
                return Path(value)
            report = audit_manifest(
                input_path(args.registry, "source_registry"), args.production_manifest,
                input_path(args.environment, "environment"), launch_root=args.launch_root,
                binary=args.binary, namelist=args.namelist,
            )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"manifest audit error: {exc}", file=sys.stderr)
        return 1
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        if args.report.exists():
            print(f"refusing to overwrite existing report: {args.report}", file=sys.stderr)
            return 1
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if report["production_gate_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
