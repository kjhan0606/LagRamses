#!/usr/bin/env python3
"""Execute the F-P1 source-identity and redistribution-rights gate."""

from __future__ import annotations

from datetime import date
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Any


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
REPOSITORY_ROOT = SNRT_ROOT.parents[1]
DEFAULT_CANDIDATE_ROOT = REPOSITORY_ROOT / "external" / "g2_candidates"
DEFAULT_MANIFEST = DEFAULT_CANDIDATE_ROOT / "acquisition_manifest_v1.json"
DEFAULT_TERMS = SNRT_ROOT / "config" / "g2_source_use_terms_evidence_v1.json"
DEFAULT_SOURCE_CONTRACT = (
    SNRT_ROOT / "config" / "g2_boccioli_roberti2026_candidate_contract_v1.json"
)
VALIDATOR_ID = "fp1.source_identity_and_rights.v1"
GATE_ID = "source_identity_and_rights"
REQUIREMENTS = {
    "citation_and_data_version",
    "per_file_and_composite_sha256",
    "machine_readable_license",
    "redistribution_permission",
    "hash_locked_local_source_mirror",
}


# This profile is the independent trust anchor. External manifests and source
# contracts are checked against it, rather than being allowed to authenticate
# mutually rewritten bytes. Updating it is a reviewed code change.
LOCKED_CANDIDATE_PROFILES: dict[str, dict[str, Any]] = {
    "boccioli_roberti2026_lc18": {
        "source_candidate_id": "boccioli_roberti2026_neutrino_ccsn",
        "release_root_relative_path": "boccioli_roberti2026_ccsn",
        "article_citation": (
            "Boccioli & Roberti 2026, Astronomy & Astrophysics 709, A201"
        ),
        "manifest_citation": (
            "Boccioli & Roberti 2026, Astronomy & Astrophysics 709, A201, "
            "DOI 10.1051/0004-6361/202557714"
        ),
        "attribution_citation": (
            "Boccioli & Roberti 2026, A&A 709, A201; article DOI "
            "10.1051/0004-6361/202557714; data DOI "
            "10.5281/zenodo.19503168"
        ),
        "article_doi": "10.1051/0004-6361/202557714",
        "data_doi": "10.5281/zenodo.19503168",
        "source_url": "https://zenodo.org/records/19503168",
        "zenodo_record_id": 19503168,
        "zenodo_record_filename": "zenodo_record_19503168.json",
        "zenodo_record_sha256": (
            "bc04d8a528e5a852524be3e0cd2d973a12b39090ad0bf2aaf71305016d17eb67"
        ),
        "zenodo_title": "Nucleosynthesis Yields from Boccioli & Roberti (2026)",
        "zenodo_creators": ["Boccioli, Luca", "Roberti, Lorenzo"],
        "license_id": "cc-by-4.0",
        "license_reference_url": "https://creativecommons.org/licenses/by/4.0/",
        "retrieved_date": "2026-09-01",
        "expected_composite_sha256": (
            "3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b"
        ),
        "files": {
            "README": {
                "relative_path": "boccioli_roberti2026_ccsn/README",
                "bytes": 3987,
                "sha256": (
                    "66dae04f90bf7b96460199a7ebbedf0126c1e70ab840e1693a8a326cd7ae2316"
                ),
                "md5": "4a42dda7587574978b7730cc45d02c8e",
            },
            "LC18.zip": {
                "relative_path": "boccioli_roberti2026_ccsn/LC18.zip",
                "bytes": 816270,
                "sha256": (
                    "249aea46713ab41cad7e8d7406835c205e4f02d36958e113a8f2231f81ebef5e"
                ),
                "md5": "c1d21fcbdf7ed200344881f8a47a211b",
            },
            "WH07.zip": {
                "relative_path": "boccioli_roberti2026_ccsn/WH07.zip",
                "bytes": 377318,
                "sha256": (
                    "1576e96d9366dd7b36a75abbbc6e9ced55c1157891417062a5cbbf439dbed84c"
                ),
                "md5": "9f6f436fbb1ffdd5ccc535f2a0976a79",
            },
            "F23.zip": {
                "relative_path": "boccioli_roberti2026_ccsn/F23.zip",
                "bytes": 444104,
                "sha256": (
                    "41b72c45e7743d9855637ae56cf1ef264c000b4639ed4805b177ab7fe5079a10"
                ),
                "md5": "43f568e8ba7a94248a7627b24f04a7c4",
            },
            "zenodo_record_19503168.json": {
                "relative_path": (
                    "boccioli_roberti2026_ccsn/zenodo_record_19503168.json"
                ),
                "bytes": 4618,
                "sha256": (
                    "bc04d8a528e5a852524be3e0cd2d973a12b39090ad0bf2aaf71305016d17eb67"
                ),
                "md5": None,
            },
        },
    }
}


class SourceIdentityRightsError(ValueError):
    """The source-identity evidence is malformed or cannot be inspected."""


def _absolute(path: Path) -> Path:
    """Return an absolute lexical path without resolving symlinks."""

    return Path(os.path.abspath(os.fspath(path)))


def _digest(path: Path, algorithm: str) -> tuple[int, str]:
    if algorithm == "md5":
        digest = hashlib.new(algorithm, usedforsecurity=False)
    else:
        digest = hashlib.new(algorithm)
    size = 0
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                size += len(block)
                digest.update(block)
    except OSError as exc:
        raise SourceIdentityRightsError(f"cannot hash {path}: {exc}") from exc
    return size, digest.hexdigest()


def _sha256(path: Path) -> tuple[int, str]:
    return _digest(path, "sha256")


def _valid_digest(value: Any, length: int) -> bool:
    return (
        isinstance(value, str)
        and len(value) == length
        and all(character in "0123456789abcdef" for character in value.lower())
    )


def _regular_file_without_symlink(path: Path, label: str) -> Path:
    path = _absolute(path)
    current = path
    while True:
        try:
            mode = current.lstat().st_mode
        except OSError as exc:
            raise SourceIdentityRightsError(
                f"{label}_missing_or_unreadable:{exc}"
            ) from exc
        if stat.S_ISLNK(mode):
            raise SourceIdentityRightsError(f"{label}_symlink_forbidden")
        if current == path:
            if not stat.S_ISREG(mode):
                raise SourceIdentityRightsError(f"{label}_not_regular_file")
        elif not stat.S_ISDIR(mode):
            raise SourceIdentityRightsError(f"{label}_parent_not_directory")
        if current.parent == current:
            break
        current = current.parent
    return path


def _directory_without_symlink(path: Path, label: str) -> Path:
    path = _absolute(path)
    try:
        mode = path.lstat().st_mode
    except OSError as exc:
        raise SourceIdentityRightsError(
            f"{label}_missing_or_unreadable:{exc}"
        ) from exc
    if stat.S_ISLNK(mode):
        raise SourceIdentityRightsError(f"{label}_symlink_forbidden")
    if not stat.S_ISDIR(mode):
        raise SourceIdentityRightsError(f"{label}_not_directory")
    return path


def _confined_regular_file(root: Path, relative: str, label: str) -> Path:
    root = _absolute(root)
    raw = Path(relative)
    if raw.is_absolute() or not raw.parts or ".." in raw.parts:
        raise SourceIdentityRightsError(f"{label}_path_unsafe")
    try:
        root_mode = root.lstat().st_mode
    except OSError as exc:
        raise SourceIdentityRightsError(
            f"candidate_root_missing_or_unreadable:{exc}"
        ) from exc
    if stat.S_ISLNK(root_mode):
        raise SourceIdentityRightsError("candidate_root_symlink_forbidden")
    if not stat.S_ISDIR(root_mode):
        raise SourceIdentityRightsError("candidate_root_not_directory")

    current = root
    for index, part in enumerate(raw.parts):
        current = current / part
        try:
            mode = current.lstat().st_mode
        except OSError as exc:
            raise SourceIdentityRightsError(
                f"{label}_missing_or_unreadable:{exc}"
            ) from exc
        if stat.S_ISLNK(mode):
            raise SourceIdentityRightsError(f"{label}_symlink_forbidden")
        final = index == len(raw.parts) - 1
        if final and not stat.S_ISREG(mode):
            raise SourceIdentityRightsError(f"{label}_not_regular_file")
        if not final and not stat.S_ISDIR(mode):
            raise SourceIdentityRightsError(f"{label}_parent_not_directory")
    try:
        current.resolve(strict=True).relative_to(root.resolve(strict=True))
    except (OSError, ValueError) as exc:
        raise SourceIdentityRightsError(f"{label}_escapes_candidate_root") from exc
    return current


def _read_json(path: Path, label: str) -> dict[str, Any]:
    path = _regular_file_without_symlink(path, label)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceIdentityRightsError(f"{label}_json_unreadable:{exc}") from exc
    if not isinstance(value, dict):
        raise SourceIdentityRightsError(f"{label}_must_be_json_object")
    return value


def _calendar_date(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    try:
        return date.fromisoformat(value).isoformat() == value
    except ValueError:
        return False


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


def _blocked_report(candidate_id: Any, reason: str) -> dict[str, Any]:
    reported_id = candidate_id if isinstance(candidate_id, str) else repr(candidate_id)
    try:
        validator_code_sha256: str | None = _sha256(TOOL_PATH)[1]
    except Exception:
        validator_code_sha256 = None
    return {
        "schema": "snrt-fp1-executable-gate-validation",
        "schema_version": 1,
        "validator_id": VALIDATOR_ID,
        "gate_id": GATE_ID,
        "candidate_id": reported_id,
        "status": "blocked",
        "passed": False,
        "requirements": {name: False for name in sorted(REQUIREMENTS)},
        "blockers": [reason],
        "package_fingerprint_sha256": None,
        "artifacts": {},
        "validator_code_sha256": validator_code_sha256,
    }


def _manifest_inventory(
    manifest_candidate: dict[str, Any], profile: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    records = manifest_candidate.get("files")
    if not isinstance(records, list) or not records:
        raise SourceIdentityRightsError("manifest_file_inventory_empty_or_malformed")
    by_path: dict[str, dict[str, Any]] = {}
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise SourceIdentityRightsError("manifest_file_record_malformed")
        relative = record["path"]
        raw = Path(relative)
        if raw.is_absolute() or ".." in raw.parts:
            raise SourceIdentityRightsError("manifest_file_path_unsafe")
        if relative in by_path:
            raise SourceIdentityRightsError("manifest_file_path_duplicate")
        by_path[relative] = record
    expected_paths = {
        record["relative_path"] for record in profile["files"].values()
    }
    if set(by_path) != expected_paths:
        raise SourceIdentityRightsError("manifest_file_inventory_not_lock_pinned")
    return by_path


def _contract_inventory(
    source: dict[str, Any], profile: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    records = source.get("files")
    if not isinstance(records, dict) or not records:
        raise SourceIdentityRightsError("source_contract_file_inventory_empty_or_malformed")
    if set(records) != set(profile["files"]):
        raise SourceIdentityRightsError("source_contract_file_inventory_not_lock_pinned")
    if any(not isinstance(record, dict) for record in records.values()):
        raise SourceIdentityRightsError("source_contract_file_record_malformed")
    return records


def _non_authoritative_terms_report(
    terms_path: Path, source_id: str, expected_citation: str
) -> dict[str, Any]:
    """Report local terms alignment without using it to establish rights."""

    try:
        terms = _read_json(terms_path, "source_use_terms")
        source = terms.get("sources", {}).get(source_id)
        if not isinstance(source, dict):
            raise SourceIdentityRightsError("source_use_terms_profile_missing")
        path = _absolute(terms_path)
        return {
            "path": str(path),
            "sha256": _sha256(path)[1],
            "authoritative_for_verdict": False,
            "citation_matches_lock": source.get("citation") == expected_citation,
            "reported_license_status": source.get("production_license_status"),
            "reported_redistribution_status": source.get("redistribution_status"),
        }
    except SourceIdentityRightsError as exc:
        return {
            "path": str(_absolute(terms_path)),
            "sha256": None,
            "authoritative_for_verdict": False,
            "citation_matches_lock": False,
            "error": str(exc),
        }


def _validate_source_identity_and_rights(
    candidate_id: str,
    *,
    candidate_root: Path,
    manifest_path: Path,
    terms_path: Path,
    source_contract_path: Path,
) -> dict[str, Any]:
    if not isinstance(candidate_id, str):
        raise SourceIdentityRightsError("candidate_id_must_be_string")
    profile = LOCKED_CANDIDATE_PROFILES.get(candidate_id)
    if profile is None:
        return _blocked_report(
            candidate_id, "candidate_has_no_code_registered_rights_profile"
        )

    candidate_root = _directory_without_symlink(candidate_root, "candidate_root")
    manifest_path = _absolute(manifest_path)
    source_contract_path = _absolute(source_contract_path)
    manifest = _read_json(manifest_path, "acquisition_manifest")
    source_contract = _read_json(source_contract_path, "candidate_source_contract")
    blockers: list[str] = []
    requirements = {name: False for name in REQUIREMENTS}
    source_id = profile["source_candidate_id"]

    candidates = manifest.get("candidates")
    if not isinstance(candidates, list):
        raise SourceIdentityRightsError("acquisition_manifest_candidate_set_malformed")
    matches = [
        record
        for record in candidates
        if isinstance(record, dict) and record.get("candidate_id") == source_id
    ]
    if len(matches) != 1:
        raise SourceIdentityRightsError(
            "acquisition_manifest_candidate_identity_not_unique"
        )
    manifest_candidate = matches[0]
    source = source_contract.get("source")
    if not isinstance(source, dict):
        raise SourceIdentityRightsError("source_contract_profile_missing")

    identity_checks = {
        "source_contract_candidate_identity_mismatch": (
            source.get("candidate_id") == source_id
        ),
        "release_root_identity_mismatch": (
            source.get("release_root_relative_path")
            == profile["release_root_relative_path"]
        ),
        "article_citation_identity_mismatch": (
            source.get("article_citation") == profile["article_citation"]
        ),
        "article_doi_identity_mismatch": (
            source.get("article_doi") == profile["article_doi"]
        ),
        "data_doi_identity_mismatch": (
            source.get("data_doi") == profile["data_doi"]
            and manifest_candidate.get("data_doi") == profile["data_doi"]
        ),
        "zenodo_record_identity_mismatch": (
            type(source.get("zenodo_record_id")) is int
            and source.get("zenodo_record_id") == profile["zenodo_record_id"]
        ),
        "manifest_citation_identity_mismatch": (
            manifest_candidate.get("citation") == profile["manifest_citation"]
        ),
        "manifest_source_url_identity_mismatch": (
            manifest_candidate.get("source_url") == profile["source_url"]
        ),
        "manifest_retrieval_date_invalid_or_unpinned": (
            _calendar_date(manifest.get("retrieved_date"))
            and manifest.get("retrieved_date") == profile["retrieved_date"]
        ),
    }
    blockers.extend(name for name, okay in identity_checks.items() if not okay)
    identity_ok = all(identity_checks.values())

    manifest_by_path = _manifest_inventory(manifest_candidate, profile)
    contract_by_name = _contract_inventory(source, profile)
    actual_records: list[dict[str, Any]] = []
    file_reports: list[dict[str, Any]] = []
    file_identity_ok = True
    mirror_ok = True
    for name, locked in profile["files"].items():
        relative = locked["relative_path"]
        manifest_record = manifest_by_path[relative]
        contract_record = contract_by_name[name]
        strict_integer_fields = (
            type(manifest_record.get("bytes")) is int
            and type(contract_record.get("bytes")) is int
            and type(locked["bytes"]) is int
        )
        declared_identity_ok = (
            strict_integer_fields
            and manifest_record.get("bytes") == locked["bytes"]
            and contract_record.get("bytes") == locked["bytes"]
            and manifest_record.get("sha256") == locked["sha256"]
            and contract_record.get("sha256") == locked["sha256"]
            and _valid_digest(locked["sha256"], 64)
        )
        expected_md5 = locked["md5"]
        if expected_md5 is not None:
            declared_identity_ok = (
                declared_identity_ok
                and manifest_record.get("md5") == expected_md5
                and contract_record.get("md5") == expected_md5
                and _valid_digest(expected_md5, 32)
            )
        elif "md5" in manifest_record or "md5" in contract_record:
            declared_identity_ok = False

        try:
            path = _confined_regular_file(candidate_root, relative, f"source_file:{name}")
            observed_bytes, observed_sha = _sha256(path)
            observed_md5 = _digest(path, "md5")[1] if expected_md5 else None
        except SourceIdentityRightsError as exc:
            file_identity_ok = False
            mirror_ok = False
            blockers.append(str(exc))
            continue
        observed_ok = (
            observed_bytes == locked["bytes"]
            and observed_sha == locked["sha256"]
            and observed_md5 == expected_md5
        )
        passed = declared_identity_ok and observed_ok
        file_reports.append(
            {
                "path": relative,
                "bytes": observed_bytes,
                "sha256": observed_sha,
                "md5": observed_md5,
                "identity_passed": passed,
            }
        )
        if not passed:
            file_identity_ok = False
            mirror_ok = False
            blockers.append(f"source_file_not_lock_pinned:{name}")
            continue
        actual_records.append(
            {"path": relative, "bytes": observed_bytes, "sha256": observed_sha}
        )

    package_fingerprint = None
    if file_identity_ok and len(actual_records) == len(profile["files"]):
        package_fingerprint = _candidate_fingerprint(source_id, actual_records)
    if package_fingerprint != profile["expected_composite_sha256"]:
        blockers.append("package_composite_fingerprint_not_lock_pinned")
        file_identity_ok = False
    requirements["per_file_and_composite_sha256"] = file_identity_ok
    requirements["hash_locked_local_source_mirror"] = mirror_ok and file_identity_ok

    version_name = profile["zenodo_record_filename"]
    version_path = _confined_regular_file(
        candidate_root,
        profile["files"][version_name]["relative_path"],
        "zenodo_version_record",
    )
    version_record = _read_json(version_path, "zenodo_version_record")
    _, version_sha = _sha256(version_path)
    metadata = version_record.get("metadata")
    creators = metadata.get("creators") if isinstance(metadata, dict) else None
    creator_names = (
        [record.get("name") for record in creators]
        if isinstance(creators, list) and all(isinstance(record, dict) for record in creators)
        else None
    )
    version_identity_ok = (
        version_sha == profile["zenodo_record_sha256"]
        and type(version_record.get("id")) is int
        and version_record.get("id") == profile["zenodo_record_id"]
        and version_record.get("doi") == profile["data_doi"]
        and isinstance(metadata, dict)
        and metadata.get("doi") == profile["data_doi"]
        and metadata.get("title") == profile["zenodo_title"]
        and creator_names == profile["zenodo_creators"]
        and identity_ok
    )
    requirements["citation_and_data_version"] = version_identity_ok
    if not version_identity_ok:
        blockers.append("citation_or_version_specific_data_identity_not_verified")

    license_record = metadata.get("license") if isinstance(metadata, dict) else None
    machine_license_ok = (
        isinstance(license_record, dict)
        and license_record.get("id") == profile["license_id"]
        and source.get("license") == profile["license_id"]
        and manifest_candidate.get("license_status")
        == "verified_cc_by_4.0_from_zenodo_record"
        and version_sha == profile["zenodo_record_sha256"]
    )
    requirements["machine_readable_license"] = machine_license_ok
    if not machine_license_ok:
        blockers.append("machine_readable_license_not_verified")

    zenodo_files = version_record.get("files")
    published_identity_ok = isinstance(zenodo_files, list)
    published_by_name: dict[str, dict[str, Any]] = {}
    if published_identity_ok:
        for record in zenodo_files:
            if not isinstance(record, dict) or not isinstance(record.get("key"), str):
                published_identity_ok = False
                break
            if record["key"] in published_by_name:
                published_identity_ok = False
                break
            published_by_name[record["key"]] = record
    published_names = {
        name for name, locked in profile["files"].items() if locked["md5"] is not None
    }
    if set(published_by_name) != published_names:
        published_identity_ok = False
    if published_identity_ok:
        for name in published_names:
            locked = profile["files"][name]
            record = published_by_name[name]
            if (
                type(record.get("size")) is not int
                or record.get("size") != locked["bytes"]
                or record.get("checksum") != f"md5:{locked['md5']}"
            ):
                published_identity_ok = False
                break
    if not published_identity_ok:
        blockers.append("zenodo_published_file_identity_not_verified")

    redistribution_ok = (
        machine_license_ok
        and published_identity_ok
        and profile["license_reference_url"]
        == "https://creativecommons.org/licenses/by/4.0/"
        and profile["attribution_citation"]
        == (
            "Boccioli & Roberti 2026, A&A 709, A201; article DOI "
            "10.1051/0004-6361/202557714; data DOI 10.5281/zenodo.19503168"
        )
    )
    requirements["redistribution_permission"] = redistribution_ok
    if not redistribution_ok:
        blockers.append("redistribution_permission_not_verified_from_lock")

    blockers = sorted(set(blockers))
    passed = all(requirements.values()) and not blockers
    terms_report = _non_authoritative_terms_report(
        terms_path, source_id, profile["attribution_citation"]
    )
    return {
        "schema": "snrt-fp1-executable-gate-validation",
        "schema_version": 1,
        "validator_id": VALIDATOR_ID,
        "gate_id": GATE_ID,
        "candidate_id": candidate_id,
        "status": "pass" if passed else "blocked",
        "passed": passed,
        "requirements": dict(sorted(requirements.items())),
        "blockers": blockers,
        "package_fingerprint_sha256": package_fingerprint,
        "artifacts": {
            "acquisition_manifest": {
                "path": str(manifest_path),
                "sha256": _sha256(manifest_path)[1],
            },
            "source_contract": {
                "path": str(source_contract_path),
                "sha256": _sha256(source_contract_path)[1],
            },
            "source_use_terms_report": terms_report,
            "source_files": file_reports,
            "retrieved_date": manifest.get("retrieved_date"),
            "version_record_id": version_record.get("id"),
            "data_doi": profile["data_doi"],
            "license_id": license_record.get("id")
            if isinstance(license_record, dict)
            else None,
            "license_reference_url": profile["license_reference_url"],
            "attribution_citation": profile["attribution_citation"],
            "lock_profile_sha256": hashlib.sha256(
                json.dumps(profile, sort_keys=True, separators=(",", ":")).encode(
                    "utf-8"
                )
            ).hexdigest(),
        },
        "validator_code_sha256": _sha256(TOOL_PATH)[1],
    }


def validate_source_identity_and_rights(
    candidate_id: str,
    *,
    candidate_root: Path = DEFAULT_CANDIDATE_ROOT,
    manifest_path: Path = DEFAULT_MANIFEST,
    terms_path: Path = DEFAULT_TERMS,
    source_contract_path: Path = DEFAULT_SOURCE_CONTRACT,
) -> dict[str, Any]:
    """Recompute the gate, converting malformed evidence into a block."""

    try:
        return _validate_source_identity_and_rights(
            candidate_id,
            candidate_root=Path(candidate_root),
            manifest_path=Path(manifest_path),
            terms_path=Path(terms_path),
            source_contract_path=Path(source_contract_path),
        )
    except Exception as exc:
        reason = (
            str(exc)
            if isinstance(exc, SourceIdentityRightsError)
            else f"validator_exception:{type(exc).__name__}"
        )
        return _blocked_report(candidate_id, reason)
