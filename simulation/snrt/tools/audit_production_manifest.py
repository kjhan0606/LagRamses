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


def audit_manifest(
    registry_path: Path,
    production_path: Path,
    environment_path: Path,
) -> dict[str, Any]:
    registry = _load(registry_path)
    production = _load(production_path)
    environment = _load(environment_path)
    blockers: list[str] = []
    warnings: list[str] = []

    if registry.get("schema") != "lrd-jwst-external-assets":
        blockers.append("registry_schema_invalid")
    if production.get("schema") != "lrd-jwst-production-readiness-manifest":
        blockers.append("production_schema_invalid")
    if environment.get("schema") != "snrt_environment_v1":
        blockers.append("environment_schema_invalid")
    if production.get("project_root") != "/gpfs/kjhan/LRD_JWST":
        blockers.append("production_root_invalid")
    if environment.get("project_root") != "/gpfs/kjhan/LRD_JWST":
        blockers.append("environment_root_invalid")

    assets = registry.get("assets", [])
    if not isinstance(assets, list) or not assets:
        blockers.append("registry_assets_missing")
        assets = []
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
        status = asset.get("status")
        if status in BLOCKING_STATUSES:
            blockers.append(f"asset_status_blocked:{asset_id}:{status}")
        simulation = asset.get("simulation")
        comparison_role = simulation.get("comparison_role") if isinstance(simulation, dict) else None
        if comparison_role in {
            "legacy_comparison_only",
            "phase0_preselector_transitional",
        }:
            blockers.append(f"comparison_asset_in_production_registry:{asset_id}")
    if duplicate_ids:
        blockers.extend(f"duplicate_asset_id:{asset_id}" for asset_id in duplicate_ids)

    policy = production.get("asset_metadata_policy", {})
    metadata = production.get("asset_metadata", {})
    required_fields = policy.get("required_fields", [])
    if policy.get("require_per_asset_record"):
        if not isinstance(metadata, dict):
            blockers.append("asset_metadata_not_object")
            metadata = {}
        for asset_id in asset_by_id:
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

    required_assets = production.get("required_production_assets", [])
    required_ids: list[str] = []
    if not isinstance(required_assets, list):
        blockers.append("required_production_assets_invalid")
        required_assets = []
    for requirement in required_assets:
        if not isinstance(requirement, dict) or not isinstance(requirement.get("id"), str):
            blockers.append("required_asset_record_invalid")
            continue
        asset_id = requirement["id"]
        required_ids.append(asset_id)
        asset = asset_by_id.get(asset_id)
        if asset is None:
            blockers.append(f"required_asset_not_registered:{asset_id}")
            continue
        asset_path = asset.get("path")
        if not isinstance(asset_path, str) or not Path(asset_path).exists():
            blockers.append(f"required_asset_path_missing:{asset_id}")
        if asset.get("status") in BLOCKING_STATUSES:
            blockers.append(f"required_asset_status_blocked:{asset_id}")
        asset_sha256 = asset.get("sha256")
        if not isinstance(asset_sha256, str) or len(asset_sha256) != 64:
            blockers.append(f"required_asset_sha256_missing:{asset_id}")

    environment_backend = environment.get("backend", {})
    if environment_backend.get("jax") != "0.11.1" or environment_backend.get("jaxlib") != "0.11.1":
        blockers.append("jax_cpu_environment_version_mismatch")
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

    repository = production.get("repository", {})
    if repository.get("production_clean_tree_required") and repository.get("working_tree_at_recording") != "clean":
        blockers.append("production_repository_dirty")

    if not production.get("production_rules", {}).get("reject_embedded_yield_fallback"):
        blockers.append("embedded_fallback_rejection_not_enabled")
    if not production.get("production_rules", {}).get("reject_legacy_and_transitional_assets"):
        blockers.append("legacy_transitional_rejection_not_enabled")

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
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--production-manifest", type=Path, required=True)
    parser.add_argument("--environment", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    try:
        report = audit_manifest(args.registry, args.production_manifest, args.environment)
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
