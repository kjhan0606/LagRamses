#!/usr/bin/env python3
"""Audit the manifest-scoped, review-only F-P2 SNIa event-yield asset."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from fp2_provenance import project_relative


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_ROOT = PROJECT_ROOT / "assets" / "review_only" / "fp2_snia" / "keegans2023"
DEFAULT_MANIFEST = PROJECT_ROOT / "manifests" / "fp2_snia_keegans2023_review_v1.json"
HASH_POLICY = "sha256(path + NUL + bytes + NUL + file_sha256 + LF) over sorted relative paths"


def _sha256(path: Path) -> tuple[int, str] | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read manifest {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"manifest must be a JSON object: {path}")
    return value


def _composite(records: list[dict[str, Any]]) -> str:
    payload = bytearray()
    for record in sorted(records, key=lambda item: item["path"].encode("utf-8")):
        payload.extend(record["path"].encode("utf-8"))
        payload.extend(b"\0")
        payload.extend(str(record["bytes"]).encode("ascii"))
        payload.extend(b"\0")
        payload.extend(record["sha256"].encode("ascii"))
        payload.extend(b"\n")
    return hashlib.sha256(payload).hexdigest()


def audit_asset(root: Path = DEFAULT_ROOT, manifest_path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    root = Path(root).resolve()
    manifest_path = Path(manifest_path).resolve()
    failures: list[dict[str, Any]] = []
    if not manifest_path.is_file():
        failures.append({"reason": "manifest_missing"})
        manifest: dict[str, Any] = {}
    else:
        try:
            manifest = _read_json(manifest_path)
        except ValueError as exc:
            failures.append({"reason": "manifest_invalid", "detail": str(exc)})
            manifest = {}
    if manifest:
        if manifest.get("schema") != "snrt-fp2-snia-event-yield-asset-manifest":
            failures.append({"reason": "schema_mismatch"})
        if manifest.get("gate") != "F-P2":
            failures.append({"reason": "gate_mismatch"})
        if manifest.get("status") != "review_only_asset":
            failures.append({"reason": "asset_must_remain_review_only"})
        if manifest.get("asset_root") != "assets/review_only/fp2_snia/keegans2023":
            failures.append({"reason": "asset_root_identity_mismatch"})
        source = manifest.get("source", {})
        if source.get("license") != "CC-BY-4.0":
            failures.append({"reason": "license_metadata_mismatch"})
        if source.get("approval_id") is not None:
            failures.append({"reason": "source_approval_must_remain_null"})
        approval = manifest.get("approval", {})
        if approval.get("approval_id") is not None:
            failures.append({"reason": "approval_id_must_remain_null"})
        if approval.get("canonical_conversion_allowed") is not False:
            failures.append({"reason": "canonical_conversion_must_remain_disabled"})
        if approval.get("runtime_activation_allowed") is not False:
            failures.append({"reason": "runtime_activation_must_remain_disabled"})

    entries = manifest.get("files") if manifest else None
    if not isinstance(entries, list) or not entries:
        failures.append({"reason": "manifest_file_list_missing"})
        entries = []
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    file_reports: list[dict[str, Any]] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str) or not entry["path"]:
            failures.append({"index": index, "reason": "file_record_invalid"})
            continue
        relative = Path(entry["path"])
        if relative.is_absolute() or ".." in relative.parts or entry["path"] in seen:
            failures.append({"path": entry["path"], "reason": "unsafe_or_duplicate_path"})
            continue
        seen.add(entry["path"])
        path = (root / relative).resolve()
        try:
            path.relative_to(root)
        except ValueError:
            failures.append({"path": entry["path"], "reason": "path_outside_asset_root"})
            continue
        observed = _sha256(path)
        report: dict[str, Any] = {"path": entry["path"], "exists": observed is not None}
        if observed is None:
            failures.append({"path": entry["path"], "reason": "file_missing"})
            file_reports.append(report)
            continue
        observed_bytes, observed_hash = observed
        passed = entry.get("bytes") == observed_bytes and str(entry.get("sha256", "")).lower() == observed_hash
        report.update({"bytes": observed_bytes, "sha256": observed_hash, "input_integrity_passed": passed})
        file_reports.append(report)
        if not passed:
            failures.append({
                "path": entry["path"],
                "reason": "fingerprint_mismatch",
                "recorded_bytes": entry.get("bytes"),
                "observed_bytes": observed_bytes,
                "recorded_sha256": entry.get("sha256"),
                "observed_sha256": observed_hash,
            })
            continue
        records.append({"path": entry["path"], "bytes": observed_bytes, "sha256": observed_hash})

    observed_package_hash = _composite(records) if len(records) == len(entries) else None
    if observed_package_hash != manifest.get("package_sha256"):
        failures.append({
            "reason": "package_fingerprint_mismatch",
            "recorded": manifest.get("package_sha256"),
            "observed": observed_package_hash,
        })
    return {
        "schema": "snrt-fp2-snia-event-yield-asset-audit",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "review_only_asset_integrity_passed" if not failures else "blocked_asset_integrity",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "runtime_activation_allowed": False,
        "root": project_relative(root),
        "manifest": project_relative(manifest_path),
        "hash_policy": HASH_POLICY,
        "candidate_id": manifest.get("candidate_id"),
        "file_count": len(entries),
        "files": file_reports,
        "package_sha256": observed_package_hash,
        "manifest_sha256": _sha256(manifest_path)[1] if _sha256(manifest_path) else None,
        "audit_code_sha256": _sha256(TOOL_PATH)[1] if _sha256(TOOL_PATH) else None,
        "audit_failures": failures,
        "interpretation": "Review-only source bytes are intact; this report does not approve isotope conversion, event physics, or runtime activation.",
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    report = audit_asset(args.root, args.manifest)
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 2 if report["audit_failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
