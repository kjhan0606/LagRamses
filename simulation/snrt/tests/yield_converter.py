#!/usr/bin/env python3
"""Tests for the deterministic normalized-row yield converter."""

from __future__ import annotations

import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import convert_yield_rows_to_canonical as converter  # noqa: E402
from convert_yield_rows_to_canonical import ConversionError  # noqa: E402
from audit_g2_source_package_fingerprints import audit_fingerprints  # noqa: E402
from fp1_source_node_fixture import (  # noqa: E402
    APPROVAL_ID,
    NODE_ID,
    approved_source_node_contract,
)
from validate_fp1_source_identity_rights import (  # noqa: E402
    LOCKED_CANDIDATE_PROFILES,
)


STAGED_SOURCE_ROOT = ROOT.parents[1] / "external" / "g2_candidates"
STAGED_SOURCE_MANIFEST = STAGED_SOURCE_ROOT / "acquisition_manifest_v1.json"


def _row(
    age: float, *, semantics: str = "cumulative", returned: float = 0.0,
    tracked: float | None = None
) -> dict:
    del semantics
    if tracked is None:
        tracked = returned
    return {
        "source_node_id": NODE_ID,
        "channel": 1,
        "initial_mass_msun_per_star": 60.0,
        "birth_metallicity_mass_fraction": 0.001,
        "age_yr": age,
        "returned_mass_msun_per_star": returned,
        "remnant_mass_msun_per_star": 0.0,
        "energy_erg_per_star": 0.0,
        "momentum_g_cm_s_per_star": [0.0, 0.0, 0.0],
        "ejecta_msun_per_star": [tracked] + [0.0] * 10,
        "net_yield_msun_per_star": [0.0] * 11,
    }


def _source(**overrides: object) -> dict:
    source = {
        "citation": "test source",
        "source_version": "test-v1",
        "source_sha256": "a" * 64,
        "release_history_semantics": "cumulative",
        "approval_id": APPROVAL_ID,
        "license_status": "approved",
        "provenance_status": "approved",
        "units": "canonical",
        "IMF": "Kroupa",
        "population_model": "single_star_ssp",
        "channel_boundaries": {"1": [40.0, 120.0]},
        "metallicity_definition": "mass fraction",
        "solar_abundance_set": "test",
        "remnant_model": "test",
        "untracked_ejecta_policy": (
            "returned_mass_minus_sum_tracked_ejecta_deposited_as_generic_metals"
        ),
        "axis_reduction_policy": {"mode": "none"},
        "energy_semantics": "cumulative_physical_erg_per_initial_star",
        "momentum_deposition_contract": "source_frame_vector_only_no_scalar_radial_deposition",
    }
    source.update(overrides)
    return source


def _repository_artifact_hashes() -> dict[str, dict[str, str]]:
    """Hash every file in the repository config/data inputs used by SNRT."""

    snapshots: dict[str, dict[str, str]] = {}
    for name in ("config", "data"):
        base = ROOT / name
        snapshots[name] = {
            str(path.relative_to(base)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(base.rglob("*"))
            if path.is_file()
        }
    return snapshots


def _staged_source_hashes() -> dict[str, object]:
    """Snapshot only manifest-listed staged files and their existing composites."""

    report = audit_fingerprints(STAGED_SOURCE_ROOT, STAGED_SOURCE_MANIFEST)
    assert report["status"] == "candidate_fingerprint_review_only"
    assert report["input_integrity_passed"] is True
    assert report["candidate_count"] == 11
    assert report["file_count"] == 65
    files: dict[str, dict[str, object]] = {}
    composites: dict[str, str] = {}
    for candidate in report["candidates"]:
        candidate_id = candidate["candidate_id"]
        composites[candidate_id] = candidate["composite_sha256"]
        for record in candidate["files"]:
            assert record["input_integrity_passed"] is True
            files[f"{candidate_id}/{record['path']}"] = {
                "bytes": record["bytes"],
                "sha256": record["sha256"],
            }

    profile = LOCKED_CANDIDATE_PROFILES["boccioli_roberti2026_lc18"]
    source_id = profile["source_candidate_id"]
    assert composites[source_id] == profile["expected_composite_sha256"]
    for locked in profile["files"].values():
        observed = files[f"{source_id}/{locked['relative_path']}"]
        assert observed["bytes"] == locked["bytes"]
        assert observed["sha256"] == locked["sha256"]

    return {
        "manifest_sha256": hashlib.sha256(STAGED_SOURCE_MANIFEST.read_bytes()).hexdigest(),
        "files": files,
        "composites": composites,
        "code_owned_lc18": {
            "source_candidate_id": source_id,
            "files": {
                name: {
                    "relative_path": locked["relative_path"],
                    "bytes": locked["bytes"],
                    "sha256": locked["sha256"],
                }
                for name, locked in profile["files"].items()
            },
            "expected_composite_sha256": profile["expected_composite_sha256"],
        },
    }


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _mapping_from_proposal(
    proposal: dict[str, object], *, node_contract_hash: str, package_hash: str
) -> dict[str, object]:
    """Construct the admitted mapping from the non-writing proposal rows."""

    return {
        "schema": "snrt-fp1-source-node-row-mapping",
        "schema_version": 1,
        "source_node_contract_sha256": node_contract_hash,
        "source_node_contract_approval_id": APPROVAL_ID,
        "physical_package_approval_id": APPROVAL_ID,
        "physical_package_sha256": package_hash,
        "canonical_asset_sha256": proposal["canonical_asset_sha256"],
        "canonical_row_count": proposal["canonical_row_count"],
        "rows": copy.deepcopy(proposal["canonical_mapping_rows"]),
    }


def _synthetic_admitted_converter_path(
    *, source_json: Path, proposal: dict[str, object], root: Path
) -> None:
    """Exercise converter writes without making a synthetic package production state."""

    node_contract_path = root / "synthetic-approved-source-nodes.json"
    package_contract_path = root / "synthetic-admitted-package.json"
    _write_json(node_contract_path, approved_source_node_contract())
    node_contract_hash = hashlib.sha256(node_contract_path.read_bytes()).hexdigest()
    package_hash = "a" * 64
    admitted_mapping = _mapping_from_proposal(
        proposal, node_contract_hash=node_contract_hash, package_hash=package_hash
    )
    package_contract = {
        "evidence_artifacts": {
            "source_node_contract": {
                "path": "config/fp1_source_node_contract_v1.json",
                "sha256": node_contract_hash,
            }
        },
        "selection": {
            "selected_package_id": "synthetic-admitted-package",
            "selected_package_sha256": package_hash,
            "source_node_mapping_sha256": converter.mapping_sha256(admitted_mapping),
            "source_node_mapping": admitted_mapping,
            "approval_id": APPROVAL_ID,
        },
    }
    _write_json(package_contract_path, package_contract)

    # Capture every tracked input before any module seam is patched. The
    # staged-source snapshot is deliberately manifest-scoped and reuses the
    # existing package-fingerprint implementation and LC18 code-owned lock.
    before = {
        "repository": _repository_artifact_hashes(),
        "staged": _staged_source_hashes(),
    }
    original_node_contract = converter.DEFAULT_SOURCE_NODE_CONTRACT
    original_package_contract = converter.DEFAULT_PHYSICAL_PACKAGE_CONTRACT
    original_node_audit = converter.audit_source_node_contract
    original_package_audit = converter.audit_physical_package_admission
    node_audit_calls: list[Path] = []
    package_audit_calls: list[Path] = []

    def synthetic_node_audit(*, node_contract_path: Path) -> dict[str, object]:
        node_audit_calls.append(Path(node_contract_path))
        return {
            "status": "approved_physical_nodes",
            "production_ready": True,
        }

    def synthetic_package_audit(*, contract_path: Path) -> dict[str, object]:
        package_audit_calls.append(Path(contract_path))
        return {
            "status": "admitted_physical_package",
            "production_ready": True,
            "canonical_conversion_allowed": True,
        }

    try:
        # These are test-only module seams. No converter production override or
        # CLI contract switch is introduced, and the real repository contracts
        # remain untouched and blocked.
        converter.DEFAULT_SOURCE_NODE_CONTRACT = node_contract_path
        converter.DEFAULT_PHYSICAL_PACKAGE_CONTRACT = package_contract_path
        converter.audit_source_node_contract = synthetic_node_audit
        converter.audit_physical_package_admission = synthetic_package_audit

        output = root / "positive" / "yield.dat"
        sidecar = root / "positive" / "yield.dat.json"
        mapping_path = root / "positive" / "yield.nodes.json"
        metadata = converter.convert(
            source_json, output, sidecar, mapping_path, node_contract_path
        )
        assert output.exists() and sidecar.exists() and mapping_path.exists()
        assert metadata["asset_sha256"] == hashlib.sha256(output.read_bytes()).hexdigest()
        assert metadata["source_node_mapping_sha256"] == converter.mapping_sha256(
            json.loads(mapping_path.read_text(encoding="utf-8"))
        )
        sidecar_data = json.loads(sidecar.read_text(encoding="utf-8"))
        assert sidecar_data["asset_sha256"] == metadata["asset_sha256"]
        assert sidecar_data["source_node_mapping_sha256"] == metadata[
            "source_node_mapping_sha256"
        ]
        assert node_audit_calls == [node_contract_path]
        assert package_audit_calls == [package_contract_path]

        def expect_mapping_rejection(
            label: str, mutated_package: dict[str, object]
        ) -> None:
            _write_json(package_contract_path, mutated_package)
            case_root = root / label
            case_output = case_root / "yield.dat"
            case_sidecar = case_root / "yield.dat.json"
            case_mapping = case_root / "yield.nodes.json"
            try:
                converter.convert(
                    source_json,
                    case_output,
                    case_sidecar,
                    case_mapping,
                    node_contract_path,
                )
            except ConversionError as exc:
                assert "source-node mapping" in str(exc), str(exc)
            else:
                raise AssertionError(f"{label} mapping mutation was accepted")
            assert not case_output.exists()
            assert not case_sidecar.exists()
            assert not case_mapping.exists()

        changed_mapping = copy.deepcopy(package_contract)
        changed_mapping["selection"]["source_node_mapping"]["rows"][0][
            "source_node_id"
        ] = "mutated-source-node"
        changed_mapping["selection"]["source_node_mapping_sha256"] = converter.mapping_sha256(
            changed_mapping["selection"]["source_node_mapping"]
        )
        expect_mapping_rejection("mapping-content-recomputed", changed_mapping)

        changed_mapping_without_hash = copy.deepcopy(package_contract)
        changed_mapping_without_hash["selection"]["source_node_mapping"]["rows"][0][
            "source_node_id"
        ] = "mutated-source-node"
        expect_mapping_rejection(
            "mapping-content-unrecomputed", changed_mapping_without_hash
        )

        changed_hash = copy.deepcopy(package_contract)
        changed_hash["selection"]["source_node_mapping_sha256"] = "f" * 64
        expect_mapping_rejection("mapping-hash-only", changed_hash)
    finally:
        converter.DEFAULT_SOURCE_NODE_CONTRACT = original_node_contract
        converter.DEFAULT_PHYSICAL_PACKAGE_CONTRACT = original_package_contract
        converter.audit_source_node_contract = original_node_audit
        converter.audit_physical_package_admission = original_package_audit

    assert converter.DEFAULT_SOURCE_NODE_CONTRACT == original_node_contract
    assert converter.DEFAULT_PHYSICAL_PACKAGE_CONTRACT == original_package_contract
    assert converter.audit_source_node_contract is original_node_audit
    assert converter.audit_physical_package_admission is original_package_audit

    # The genuine repository path remains fail-closed after all synthetic seams
    # have been restored. Call the Python audit function, not its writing CLI,
    # so this check cannot rewrite the tracked audit artifact.
    real_package_report = converter.audit_physical_package_admission()
    assert real_package_report["status"] == "blocked_no_qualified_physical_package"
    assert real_package_report["canonical_conversion_allowed"] is False
    assert real_package_report["runtime_deposition_allowed"] is False
    assert real_package_report["production_ready"] is False
    assert real_package_report["publication_ready"] is False
    assert real_package_report["physical_node_count"] == 0
    assert real_package_report["selected_package_id"] is None

    try:
        blocked_output = root / "repository-blocked" / "yield.dat"
        blocked_sidecar = root / "repository-blocked" / "yield.dat.json"
        blocked_mapping = root / "repository-blocked" / "yield.nodes.json"
        converter.convert(
            source_json,
            blocked_output,
            blocked_sidecar,
            blocked_mapping,
        )
    except ConversionError as exc:
        assert (
            "requires the approved repository source-node contract" in str(exc)
            or "physical package" in str(exc)
        )
    else:
        raise AssertionError("real repository conversion was not fail-closed")
    assert not blocked_output.exists()
    assert not blocked_sidecar.exists()
    assert not blocked_mapping.exists()

    # These comparisons intentionally occur after every post-restore check,
    # covering the entire synthetic fixture window.
    assert _repository_artifact_hashes() == before["repository"]
    assert _staged_source_hashes() == before["staged"]


def main() -> int:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source_json = root / "source.json"
        output = root / "yield.dat"
        sidecar = root / "yield.dat.json"
        mapping = root / "yield.nodes.json"
        source_json.write_text(
            json.dumps(
                {
                    "source": _source(),
                    "rows": [_row(1.0, returned=6.0), _row(0.0)],
                }
            ),
            encoding="utf-8",
        )
        proposal = converter.build_conversion_proposal(source_json)
        assert proposal["status"] == "proposal_only_not_admitted"
        assert proposal["writes_performed"] is False
        assert proposal["canonical_row_count"] == 2
        assert len(proposal["canonical_mapping_rows"]) == 2
        assert not output.exists() and not sidecar.exists() and not mapping.exists()
        _synthetic_admitted_converter_path(
            source_json=source_json, proposal=proposal, root=root
        )
        original_argv = list(sys.argv)
        proposal_stdout = io.StringIO()
        sys.argv = [
            "convert_yield_rows_to_canonical.py",
            str(source_json),
            "--proposal",
            "--output",
            str(output),
            "--sidecar",
            str(sidecar),
            "--node-mapping",
            str(mapping),
        ]
        try:
            with contextlib.redirect_stdout(proposal_stdout):
                assert converter.main() == 0
        finally:
            sys.argv = original_argv
        cli_proposal = json.loads(proposal_stdout.getvalue())
        assert cli_proposal["status"] == "proposal_only_not_admitted"
        assert not output.exists() and not sidecar.exists() and not mapping.exists()
        node_contract = root / "approved-source-nodes.json"
        node_contract.write_text(
            json.dumps(approved_source_node_contract()), encoding="utf-8"
        )
        repository_contract = converter.DEFAULT_SOURCE_NODE_CONTRACT
        converter.DEFAULT_SOURCE_NODE_CONTRACT = node_contract
        normalized = converter._normalize_rows(  # noqa: SLF001 - bounded unit seam
            json.loads(source_json.read_text())["rows"]
        )
        lines = [
            line
            for line in converter._format_table(normalized).splitlines()  # noqa: SLF001
            if not line.startswith("#")
        ]
        assert len(lines) == 2
        assert float(lines[0].split()[3]) == 0.0
        try:
            converter.convert(source_json, output, sidecar, mapping, node_contract)
        except ConversionError as exc:
            assert "physical package" in str(exc)
        else:
            raise AssertionError("conversion bypassed blocked F-P1H-E package admission")
        assert not output.exists() and not sidecar.exists() and not mapping.exists()

        inconsistent_projection_json = root / "inconsistent-projection.json"
        inconsistent_projection = _row(5.0e6, returned=50.0)
        inconsistent_projection["channel"] = 3
        inconsistent_projection["remnant_mass_msun_per_star"] = 0.0
        inconsistent_projection["energy_erg_per_star"] = 1.0e51
        inconsistent_projection_json.write_text(
            json.dumps(
                {"source": _source(), "rows": [inconsistent_projection]}
            ),
            encoding="utf-8",
        )
        try:
            converter.convert(
                inconsistent_projection_json,
                root / "inconsistent.dat",
                root / "inconsistent.dat.json",
                root / "inconsistent.nodes.json",
                node_contract,
            )
        except ConversionError as exc:
            assert "violates its source-node projection" in str(exc)
        else:
            raise AssertionError("canonical payload was not bound to source-node physics")

        overfull_json = root / "overfull.json"
        overfull_json.write_text(
            json.dumps(
                {"source": _source(), "rows": [_row(0.0, returned=0.1, tracked=0.2)]}
            ),
            encoding="utf-8",
        )
        try:
            converter.convert(
                overfull_json,
                root / "overfull.dat",
                root / "overfull.dat.json",
                root / "overfull.nodes.json",
            )
        except ConversionError as exc:
            assert "tracked ejecta exceeding" in str(exc)
        else:
            raise AssertionError("overfull tracked ejecta were not rejected")

        rate_json = root / "rate.json"
        rate_json.write_text(
            json.dumps({"source": _source(release_history_semantics="rate"), "rows": [_row(0.0)]}),
            encoding="utf-8",
        )
        try:
            converter.convert(
                rate_json,
                root / "rate.dat",
                root / "rate.dat.json",
                root / "rate.nodes.json",
            )
        except ConversionError as exc:
            assert "rate tables" in str(exc)
        else:
            raise AssertionError("rate input was not rejected")

        missing_node_id = root / "missing-node-id.json"
        bad_row = _row(0.0)
        del bad_row["source_node_id"]
        missing_node_id.write_text(
            json.dumps({"source": _source(), "rows": [bad_row]}), encoding="utf-8"
        )
        try:
            converter.convert(
                missing_node_id,
                root / "missing-node-id.dat",
                root / "missing-node-id.dat.json",
                root / "missing-node-id.nodes.json",
            )
        except ConversionError as exc:
            assert "source_node_id" in str(exc)
        else:
            raise AssertionError("canonical row without source_node_id was accepted")

        unknown_node_json = root / "unknown-node.json"
        unknown_node_row = _row(0.0)
        unknown_node_row["source_node_id"] = "does-not-exist"
        unknown_node_json.write_text(
            json.dumps({"source": _source(), "rows": [unknown_node_row]}),
            encoding="utf-8",
        )
        try:
            converter.convert(
                unknown_node_json,
                root / "unknown-node.dat",
                root / "unknown-node.dat.json",
                root / "unknown-node.nodes.json",
                node_contract,
            )
        except ConversionError as exc:
            assert "absent from the approved contract" in str(exc)
        else:
            raise AssertionError("unknown source_node_id was accepted")

        converter.DEFAULT_SOURCE_NODE_CONTRACT = repository_contract
        blocked_by_review_contract = root / "blocked-review-contract.json"
        blocked_by_review_contract.write_text(
            json.dumps({"source": _source(), "rows": [_row(0.0)]}),
            encoding="utf-8",
        )
        try:
            converter.convert(
                blocked_by_review_contract,
                root / "blocked.dat",
                root / "blocked.dat.json",
                root / "blocked.nodes.json",
            )
        except ConversionError as exc:
            assert "requires the approved repository source-node contract" in str(exc)
        else:
            raise AssertionError("review-only source-node contract allowed conversion")
    print("YIELD_CONVERTER_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
