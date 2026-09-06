#!/usr/bin/env python3
"""Build reproducible, manifest-scoped fingerprints for G2 candidates.

The acquisition manifest remains the source of truth for staged files. This
tool verifies every recorded byte count and SHA256, then hashes a canonical
UTF-8 serialization of each candidate's manifest entries. The result is a
provenance aid for candidate review only: it is not a physical source
approval and is never accepted as a canonical-yield input by itself.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"
DEFAULT_MANIFEST_NAME = "acquisition_manifest_v1.json"
SCHEME = (
    "sha256(candidate_id + NUL + sorted UTF-8 records "
    "path + NUL + bytes + NUL + file_sha256 + LF)"
)


class FingerprintAuditError(ValueError):
    """The acquisition manifest cannot be read."""


def _sha256(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise FingerprintAuditError(f"cannot read staged file {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_manifest(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FingerprintAuditError(f"cannot read acquisition manifest {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise FingerprintAuditError("acquisition manifest must be a JSON object")
    return value


def _resolve_relative_file(root: Path, raw_path: Any) -> tuple[str, Path]:
    if not isinstance(raw_path, str) or not raw_path:
        raise FingerprintAuditError("manifest file path must be a non-empty string")
    relative = Path(raw_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise FingerprintAuditError(
            f"manifest file path must be relative and confined to the candidate root: {raw_path!r}"
        )
    resolved = (root / relative).resolve(strict=False)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise FingerprintAuditError(
            f"manifest file path resolves outside the candidate root: {raw_path!r}"
        ) from exc
    return raw_path, resolved


def _candidate_fingerprint(candidate_id: str, records: list[dict[str, Any]]) -> str:
    ordered = sorted(records, key=lambda record: record["path"].encode("utf-8"))
    payload = bytearray(candidate_id.encode("utf-8"))
    payload.extend(b"\0")
    for record in ordered:
        payload.extend(record["path"].encode("utf-8"))
        payload.extend(b"\0")
        payload.extend(str(record["bytes"]).encode("ascii"))
        payload.extend(b"\0")
        payload.extend(record["sha256"].encode("ascii"))
        payload.extend(b"\n")
    return hashlib.sha256(payload).hexdigest()


def _blocked_report(root: Path, manifest_path: Path, failures: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "schema": "snrt-g2-source-package-fingerprint-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "candidate_fingerprint_blocked_input_integrity",
        "production_ready": False,
        "root": str(root),
        "manifest": str(manifest_path),
        "scheme": SCHEME,
        "input_integrity_passed": False,
        "audit_failures": failures,
        "candidate_count": 0,
        "file_count": 0,
        "candidates": [],
        "approval_id": None,
    }


def audit_fingerprints(root: Path, manifest_path: Path | None = None) -> dict[str, Any]:
    root = Path(root).resolve()
    manifest_path = (
        Path(manifest_path).resolve()
        if manifest_path is not None
        else root / DEFAULT_MANIFEST_NAME
    )
    if not manifest_path.is_file():
        return _blocked_report(root, manifest_path, [{"reason": "acquisition_manifest_missing"}])
    try:
        manifest = _load_manifest(manifest_path)
    except FingerprintAuditError as exc:
        return _blocked_report(root, manifest_path, [{"reason": "manifest_read_error", "detail": str(exc)}])

    failures: list[dict[str, Any]] = []
    candidate_reports: list[dict[str, Any]] = []
    candidates = manifest.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        failures.append({"reason": "candidate_manifest_empty_or_invalid"})
        candidates = []

    seen_candidate_ids: set[str] = set()
    for candidate_index, candidate in enumerate(candidates):
        if not isinstance(candidate, dict):
            failures.append({
                "candidate_index": candidate_index,
                "reason": "candidate_record_not_object",
            })
            continue
        candidate_id = candidate.get("candidate_id")
        if not isinstance(candidate_id, str) or not candidate_id:
            failures.append({"candidate_index": candidate_index, "reason": "candidate_id_missing"})
            candidate_id = f"<candidate-{candidate_index}>"
        elif candidate_id in seen_candidate_ids:
            failures.append({"candidate_id": candidate_id, "reason": "duplicate_candidate_id"})
        else:
            seen_candidate_ids.add(candidate_id)

        entries = candidate.get("files")
        if not isinstance(entries, list) or not entries:
            failures.append({"candidate_id": candidate_id, "reason": "candidate_file_coverage_empty"})
            candidate_reports.append({
                "candidate_id": candidate_id,
                "file_count": 0,
                "files": [],
                "composite_sha256": None,
                "input_integrity_passed": False,
            })
            continue

        records: list[dict[str, Any]] = []
        file_reports: list[dict[str, Any]] = []
        seen_paths: set[str] = set()
        for file_index, entry in enumerate(entries):
            if not isinstance(entry, dict):
                failures.append({
                    "candidate_id": candidate_id,
                    "file_index": file_index,
                    "reason": "file_record_not_object",
                })
                continue
            try:
                relative, path = _resolve_relative_file(root, entry.get("path"))
            except FingerprintAuditError as exc:
                failures.append({
                    "candidate_id": candidate_id,
                    "file_index": file_index,
                    "reason": "unsafe_file_path",
                    "detail": str(exc),
                })
                continue
            if relative in seen_paths:
                failures.append({"candidate_id": candidate_id, "path": relative, "reason": "duplicate_file_path"})
                continue
            seen_paths.add(relative)
            file_report: dict[str, Any] = {"path": relative, "exists": path.is_file()}
            if not path.is_file():
                failures.append({"candidate_id": candidate_id, "path": relative, "reason": "missing_file"})
                file_reports.append(file_report)
                continue
            observed_bytes, observed_sha256 = _sha256(path)
            passed = (
                entry.get("bytes") == observed_bytes
                and str(entry.get("sha256", "")).lower() == observed_sha256
            )
            file_report.update({
                "bytes": observed_bytes,
                "sha256": observed_sha256,
                "input_integrity_passed": passed,
            })
            file_reports.append(file_report)
            if not passed:
                failures.append({
                    "candidate_id": candidate_id,
                    "path": relative,
                    "reason": "fingerprint_mismatch",
                    "recorded_bytes": entry.get("bytes"),
                    "observed_bytes": observed_bytes,
                    "recorded_sha256": entry.get("sha256"),
                    "observed_sha256": observed_sha256,
                })
                continue
            records.append({
                "path": relative,
                "bytes": observed_bytes,
                "sha256": observed_sha256,
            })

        candidate_ok = len(records) == len(entries)
        candidate_reports.append({
            "candidate_id": candidate_id,
            "file_count": len(entries),
            "files": sorted(file_reports, key=lambda record: record["path"].encode("utf-8")),
            "composite_sha256": _candidate_fingerprint(candidate_id, records) if candidate_ok else None,
            "input_integrity_passed": candidate_ok,
        })

    report = {
        "schema": "snrt-g2-source-package-fingerprint-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "candidate_fingerprint_review_only" if not failures else "candidate_fingerprint_blocked_input_integrity",
        "production_ready": False,
        "root": str(root),
        "manifest": str(manifest_path),
        "scheme": SCHEME,
        "manifest_hash_policy": manifest.get("policy", {}).get("hash_policy"),
        "input_integrity_passed": not failures,
        "audit_failures": failures,
        "candidate_count": len(candidate_reports),
        "file_count": sum(item["file_count"] for item in candidate_reports),
        "candidates": candidate_reports,
        "approval_id": None,
        "interpretation": (
            "Composite hashes identify the exact manifest-scoped staged file set. "
            "They do not assert scientific validity, licensing approval, or canonical conversion."
        ),
        "audit_code_sha256": _sha256(TOOL_PATH)[1],
    }
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help=f"candidate root (default: {DEFAULT_ROOT})")
    parser.add_argument("--manifest", type=Path, help="acquisition manifest (default: ROOT/acquisition_manifest_v1.json)")
    parser.add_argument("--json-out", type=Path, help="write the JSON audit report")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    report = audit_fingerprints(args.root, args.manifest)
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 2 if report["audit_failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
