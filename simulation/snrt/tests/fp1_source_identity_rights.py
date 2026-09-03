#!/usr/bin/env python3
"""Adversarial tests for the executable F-P1 source-rights validator."""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from fp1_gate_validator_registry import (  # noqa: E402
    GateValidatorRegistryError,
    REGISTERED_VALIDATORS,
    run_registered_validator,
)
from validate_fp1_source_identity_rights import (  # noqa: E402
    DEFAULT_CANDIDATE_ROOT,
    DEFAULT_MANIFEST,
    DEFAULT_SOURCE_CONTRACT,
    DEFAULT_TERMS,
    GATE_ID,
    VALIDATOR_ID,
    validate_source_identity_and_rights,
)


CANDIDATE_ID = "boccioli_roberti2026_lc18"
SOURCE_ID = "boccioli_roberti2026_neutrino_ccsn"
RELEASE = "boccioli_roberti2026_ccsn"
EXPECTED_FINGERPRINT = (
    "3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b"
)


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value), encoding="utf-8")


def _fixture(directory: Path) -> dict[str, Path]:
    candidate_root = directory / "g2_candidates"
    candidate_root.mkdir()
    manifest_path = candidate_root / DEFAULT_MANIFEST.name
    source_contract_path = directory / "source_contract.json"
    terms_path = directory / "terms.json"
    shutil.copy2(DEFAULT_MANIFEST, manifest_path)
    shutil.copy2(DEFAULT_SOURCE_CONTRACT, source_contract_path)
    shutil.copy2(DEFAULT_TERMS, terms_path)
    shutil.copytree(
        DEFAULT_CANDIDATE_ROOT / RELEASE,
        candidate_root / RELEASE,
    )
    return {
        "candidate_root": candidate_root,
        "manifest_path": manifest_path,
        "source_contract_path": source_contract_path,
        "terms_path": terms_path,
    }


def _validate(paths: dict[str, Path]) -> dict[str, Any]:
    return validate_source_identity_and_rights(CANDIDATE_ID, **paths)


def _manifest_candidate(paths: dict[str, Path]) -> tuple[dict[str, Any], dict[str, Any]]:
    manifest = _load(paths["manifest_path"])
    candidate = next(
        value for value in manifest["candidates"] if value.get("candidate_id") == SOURCE_ID
    )
    return manifest, candidate


def _contract_source(paths: dict[str, Path]) -> tuple[dict[str, Any], dict[str, Any]]:
    contract = _load(paths["source_contract_path"])
    return contract, contract["source"]


def _expect_blocked(
    mutate: Callable[[dict[str, Path]], None], fragment: str
) -> None:
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-rights-") as directory:
        paths = _fixture(Path(directory))
        mutate(paths)
        report = _validate(paths)
        assert report["status"] == "blocked", report
        assert report["passed"] is False
        assert report["blockers"]
        assert any(fragment in value for value in report["blockers"]), report


def _mutate_json(
    path: Path, mutate: Callable[[dict[str, Any]], None]
) -> None:
    value = _load(path)
    mutate(value)
    _write(path, value)


def _candidate_identity(paths: dict[str, Path]) -> None:
    _mutate_json(
        paths["source_contract_path"],
        lambda value: value["source"].__setitem__("candidate_id", "substitute"),
    )


def _release_root(paths: dict[str, Path]) -> None:
    _mutate_json(
        paths["source_contract_path"],
        lambda value: value["source"].__setitem__("release_root_relative_path", RELEASE + "_other"),
    )


def _empty_inventory(paths: dict[str, Path]) -> None:
    manifest, candidate = _manifest_candidate(paths)
    candidate["files"] = []
    _write(paths["manifest_path"], manifest)
    contract, source = _contract_source(paths)
    source["files"] = {}
    _write(paths["source_contract_path"], contract)


def _extra_inventory(paths: dict[str, Path]) -> None:
    manifest, candidate = _manifest_candidate(paths)
    candidate["files"].append(
        {
            "path": f"{RELEASE}/extra",
            "bytes": 0,
            "sha256": "0" * 64,
        }
    )
    _write(paths["manifest_path"], manifest)


def _duplicate_inventory(paths: dict[str, Path]) -> None:
    manifest, candidate = _manifest_candidate(paths)
    candidate["files"].append(copy.deepcopy(candidate["files"][0]))
    _write(paths["manifest_path"], manifest)


def _coherent_rewrite(paths: dict[str, Path]) -> None:
    target = paths["candidate_root"] / RELEASE / "README"
    target.write_bytes(target.read_bytes() + b"\ncoherent rewrite\n")
    payload = target.read_bytes()
    size = len(payload)
    sha256 = hashlib.sha256(payload).hexdigest()
    md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    manifest, candidate = _manifest_candidate(paths)
    record = next(value for value in candidate["files"] if value["path"].endswith("/README"))
    record.update({"bytes": size, "sha256": sha256, "md5": md5})
    _write(paths["manifest_path"], manifest)
    contract, source = _contract_source(paths)
    source["files"]["README"].update(
        {"bytes": size, "sha256": sha256, "md5": md5}
    )
    _write(paths["source_contract_path"], contract)


def _internal_symlink(paths: dict[str, Path]) -> None:
    release = paths["candidate_root"] / RELEASE
    target = release / "README.real"
    (release / "README").rename(target)
    (release / "README").symlink_to(target.name)


def _external_symlink(paths: dict[str, Path]) -> None:
    release = paths["candidate_root"] / RELEASE
    external = paths["candidate_root"].parent / "README.external"
    (release / "README").rename(external)
    (release / "README").symlink_to(external)


def _non_regular_file(paths: dict[str, Path]) -> None:
    target = paths["candidate_root"] / RELEASE / "README"
    target.unlink()
    os.mkfifo(target)


def _traversal(paths: dict[str, Path]) -> None:
    manifest, candidate = _manifest_candidate(paths)
    candidate["files"][0]["path"] = f"{RELEASE}/../README"
    _write(paths["manifest_path"], manifest)


def _missing_file(paths: dict[str, Path]) -> None:
    (paths["candidate_root"] / RELEASE / "LC18.zip").unlink()


def _malformed_json(paths: dict[str, Path]) -> None:
    paths["manifest_path"].write_text("{", encoding="utf-8")


def _float_bytes(paths: dict[str, Path]) -> None:
    manifest, candidate = _manifest_candidate(paths)
    candidate["files"][0]["bytes"] = float(candidate["files"][0]["bytes"])
    _write(paths["manifest_path"], manifest)


def _boolean_bytes(paths: dict[str, Path]) -> None:
    contract, source = _contract_source(paths)
    source["files"]["README"]["bytes"] = True
    _write(paths["source_contract_path"], contract)


def _invalid_date(paths: dict[str, Path]) -> None:
    _mutate_json(
        paths["manifest_path"],
        lambda value: value.__setitem__("retrieved_date", "2026-99-99"),
    )


def _null_doi(paths: dict[str, Path]) -> None:
    _mutate_json(
        paths["source_contract_path"],
        lambda value: value["source"].__setitem__("article_doi", None),
    )


def _substitute_record_and_license(paths: dict[str, Path]) -> None:
    record_path = paths["candidate_root"] / RELEASE / "zenodo_record_19503168.json"
    record = _load(record_path)
    record["id"] = 1
    record["doi"] = "10.5281/zenodo.1"
    record["metadata"]["doi"] = "10.5281/zenodo.1"
    record["metadata"]["license"]["id"] = "cc0-1.0"
    _write(record_path, record)


def _manifest_symlink(paths: dict[str, Path]) -> None:
    target = paths["manifest_path"].with_suffix(".real.json")
    paths["manifest_path"].rename(target)
    paths["manifest_path"].symlink_to(target.name)


def _registry_exception_tests(baseline: dict[str, Any]) -> None:
    record = REGISTERED_VALIDATORS[VALIDATOR_ID]
    original = record["runner"]
    try:
        def raises(_: str) -> dict[str, Any]:
            raise RuntimeError("synthetic runner failure")

        record["runner"] = raises
        try:
            run_registered_validator(
                validator_id=VALIDATOR_ID,
                gate_id=GATE_ID,
                candidate_id=CANDIDATE_ID,
            )
        except GateValidatorRegistryError as exc:
            assert "validator raised RuntimeError" in str(exc)
        else:
            raise AssertionError("runner exception escaped the registry boundary")

        synthetic = copy.deepcopy(baseline)
        synthetic["requirements"] = {
            name: True for name in synthetic["requirements"]
        }
        synthetic["blockers"] = ["synthetic_identity_blocker"]
        synthetic["passed"] = False
        synthetic["status"] = "blocked"
        record["runner"] = lambda _: copy.deepcopy(synthetic)
        blocked = run_registered_validator(
            validator_id=VALIDATOR_ID,
            gate_id=GATE_ID,
            candidate_id=CANDIDATE_ID,
        )
        assert all(blocked["requirements"].values())
        assert blocked["passed"] is False
        assert blocked["blockers"] == ["synthetic_identity_blocker"]
    finally:
        record["runner"] = original

    try:
        run_registered_validator(
            validator_id=[],  # type: ignore[arg-type]
            gate_id=GATE_ID,
            candidate_id=CANDIDATE_ID,
        )
    except GateValidatorRegistryError as exc:
        assert "validator id must be a string" in str(exc)
    else:
        raise AssertionError("unhashable validator id was not controlled")


def _terms_are_non_authoritative() -> None:
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-terms-") as directory:
        paths = _fixture(Path(directory))
        terms = _load(paths["terms_path"])
        terms["sources"][SOURCE_ID]["redistribution_status"] = "substituted"
        terms["sources"][SOURCE_ID]["citation"] = "substituted"
        _write(paths["terms_path"], terms)
        report = _validate(paths)
        assert report["passed"] is True
        assert report["requirements"]["redistribution_permission"] is True
        evidence = report["artifacts"]["source_use_terms_report"]
        assert evidence["authoritative_for_verdict"] is False
        assert evidence["citation_matches_lock"] is False


def main() -> int:
    before = validate_source_identity_and_rights(CANDIDATE_ID)
    assert before["status"] == "pass", before
    assert before["passed"] is True
    assert before["blockers"] == []
    assert all(before["requirements"].values())
    assert before["package_fingerprint_sha256"] == EXPECTED_FINGERPRINT
    assert "hash_locked_local_source_mirror" in before["requirements"]
    assert "immutable_local_source_mirror" not in before["requirements"]

    unsupported = validate_source_identity_and_rights("sukhbold2016_w18_n20")
    assert unsupported["status"] == "blocked"
    assert unsupported["passed"] is False
    assert unsupported["blockers"] == [
        "candidate_has_no_code_registered_rights_profile"
    ]

    cases = (
        (_candidate_identity, "source_contract_candidate_identity_mismatch"),
        (_release_root, "release_root_identity_mismatch"),
        (_empty_inventory, "manifest_file_inventory_empty_or_malformed"),
        (_extra_inventory, "manifest_file_inventory_not_lock_pinned"),
        (_duplicate_inventory, "manifest_file_path_duplicate"),
        (_coherent_rewrite, "source_file_not_lock_pinned:README"),
        (_internal_symlink, "source_file:README_symlink_forbidden"),
        (_external_symlink, "source_file:README_symlink_forbidden"),
        (_non_regular_file, "source_file:README_not_regular_file"),
        (_traversal, "manifest_file_path_unsafe"),
        (_missing_file, "source_file:LC18.zip_missing_or_unreadable"),
        (_malformed_json, "acquisition_manifest_json_unreadable"),
        (_float_bytes, "source_file_not_lock_pinned:README"),
        (_boolean_bytes, "source_file_not_lock_pinned:README"),
        (_invalid_date, "manifest_retrieval_date_invalid_or_unpinned"),
        (_null_doi, "article_doi_identity_mismatch"),
        (_substitute_record_and_license, "source_file_not_lock_pinned"),
        (_manifest_symlink, "acquisition_manifest_symlink_forbidden"),
    )
    for mutate, fragment in cases:
        _expect_blocked(mutate, fragment)

    _registry_exception_tests(before)
    _terms_are_non_authoritative()

    after = validate_source_identity_and_rights(CANDIDATE_ID)
    assert after["passed"] is True
    assert after["package_fingerprint_sha256"] == EXPECTED_FINGERPRINT
    assert after["artifacts"]["source_files"] == before["artifacts"]["source_files"]
    print("FP1_SOURCE_IDENTITY_RIGHTS_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
