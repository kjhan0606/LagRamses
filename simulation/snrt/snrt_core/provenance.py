"""Small, deterministic provenance primitives for source-closure sidecars."""

from __future__ import annotations

from collections.abc import Mapping
import hashlib
import json
from pathlib import Path
import re


_SHA256 = re.compile(r"^[0-9a-f]{64}$")
PAYLOAD_HASH_SCHEME = "sha256_canonical_json_without_payload_sha256_v1"


def sha256_file(path: str | Path) -> str:
    """Return the raw-byte SHA-256 of a file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_sha256(value: object, label: str, context: str | Path) -> str:
    """Validate a serialized lowercase SHA-256 digest."""

    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"{context}: {label} must be a lowercase SHA-256")
    return value


def canonical_payload_sha256(payload: Mapping[str, object]) -> str:
    """Hash a JSON object after removing its self-hash field.

    JSON formatting and dictionary insertion order are intentionally excluded;
    values, paths, hashes, and closure arrays remain covered.  ``allow_nan``
    is disabled so a non-JSON numeric payload cannot acquire an ambiguous
    digest.
    """

    body = {key: value for key, value in payload.items() if key != "payload_sha256"}
    encoded = json.dumps(
        body,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_code_manifest(expected: Mapping[str, str | Path]) -> list[dict[str, str]]:
    """Build a sorted role/path/hash manifest for a closure's code inputs."""

    manifest: list[dict[str, str]] = []
    for role, path in sorted(expected.items()):
        resolved = Path(path).resolve()
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        manifest.append(
            {
                "role": str(role),
                "path": str(resolved),
                "sha256": sha256_file(resolved),
            }
        )
    return manifest


def validate_code_manifest(
    metadata: Mapping[str, object],
    expected: Mapping[str, str | Path],
    *,
    context: str | Path,
) -> tuple[dict[str, str], ...]:
    """Require exactly the expected closure-code roles and current file hashes."""

    raw_manifest = metadata.get("closure_code_manifest")
    if not isinstance(raw_manifest, list):
        raise ValueError(f"{context}: closure_code_manifest must be a list")
    entries: dict[str, dict[str, str]] = {}
    for raw_entry in raw_manifest:
        if not isinstance(raw_entry, Mapping):
            raise ValueError(f"{context}: closure code manifest entries must be objects")
        role = raw_entry.get("role")
        path = raw_entry.get("path")
        digest = raw_entry.get("sha256")
        if not isinstance(role, str) or not role.strip() or role in entries:
            raise ValueError(f"{context}: closure code manifest has invalid or duplicate role")
        if not isinstance(path, str) or not path.strip():
            raise ValueError(f"{context}: closure code manifest role {role!r} lacks path")
        entries[role] = {
            "role": role,
            "path": path,
            "sha256": require_sha256(digest, f"closure_code_manifest[{role}].sha256", context),
        }
    expected_roles = set(expected)
    if set(entries) != expected_roles:
        raise ValueError(
            f"{context}: closure code manifest roles {sorted(entries)} do not match "
            f"required {sorted(expected_roles)}"
        )
    for role, expected_path in expected.items():
        entry = entries[role]
        resolved = Path(expected_path).resolve()
        if Path(entry["path"]).resolve() != resolved:
            raise ValueError(f"{context}: closure code path for role {role!r} is not the expected file")
        try:
            actual_hash = sha256_file(resolved)
        except (FileNotFoundError, OSError) as error:
            raise ValueError(f"{context}: closure code file for role {role!r} is unavailable") from error
        if actual_hash != entry["sha256"]:
            raise ValueError(f"{context}: closure code hash for role {role!r} does not match its file")
    return tuple(entries[role] for role in sorted(entries))


def validate_payload_hash(metadata: Mapping[str, object], *, context: str | Path) -> str:
    """Validate the canonical self-hash used by source-bound sidecars."""

    if metadata.get("payload_hash_scheme") != PAYLOAD_HASH_SCHEME:
        raise ValueError(f"{context}: unsupported or missing payload hash scheme")
    declared = require_sha256(metadata.get("payload_sha256"), "payload_sha256", context)
    actual = canonical_payload_sha256(metadata)
    if actual != declared:
        raise ValueError(f"{context}: payload hash does not match the serialized payload")
    return declared
