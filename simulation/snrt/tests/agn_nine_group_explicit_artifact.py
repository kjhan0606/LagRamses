#!/usr/bin/env python3
"""Validate the canonical explicit-SED nine-group engineering control."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parents[1]
ARTIFACT = ROOT / "data" / "agn_nine_group_explicit_validation.json"
LEDGER = ROOT / "data" / "p4_explicit_agn_photon_ledger.csv"
METADATA = ROOT / "data" / "p4_explicit_agn_photon_ledger.json"
GROUP_EDGES = ROOT / "config" / "p0_photon_group_edges_ev.txt"
GENERATOR = ROOT / "tools" / "p4_build_agn_photon_ledger.py"
VALIDATOR = ROOT / "tools" / "validate_agn_nine_group_ledger.py"
SOURCE_REBIND_TOOL = ROOT / "tools" / "p4_attach_pilot_sources.py"
P4_RUNNER = ROOT / "tools" / "p4_run_transport_pilot.py"
STATIC_INPUT = ROOT / "data" / "p4_coeval_static_rt_input_agn9_explicit.h5"
STATIC_METADATA = ROOT / "data" / "p4_coeval_static_rt_input_agn9_explicit.json"
TRANSPORT_CONTROL = ROOT / "data" / "p4_validation" / "p4_agn9_explicit_stage4_0p001myr.h5"
EXPLICIT_SED = ROOT / "data" / "p4_explicit_agn_sed_control.csv"
EXTERNAL_ASSET_MANIFEST = ROOT / "data" / "agn_nine_group_external_assets.json"
ATTESTATION_SCOPE = "simulation/snrt"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_current_provenance(payload: dict[str, object], metadata: dict[str, object]) -> None:
    provenance = payload["provenance"]
    assert isinstance(provenance, dict)
    current_head = subprocess.check_output(
        ("git", "rev-parse", "HEAD"), cwd=REPOSITORY_ROOT, text=True
    ).strip()
    working_tree_status = subprocess.check_output(
        ("git", "status", "--short", "--untracked-files=all", "--", ATTESTATION_SCOPE),
        cwd=REPOSITORY_ROOT,
        text=True,
    )
    assert provenance["git_head"] == current_head
    assert provenance["working_tree_attestation_scope"] == ATTESTATION_SCOPE
    assert provenance["working_tree_clean"] is (working_tree_status == "")
    assert provenance["working_tree_status_sha256"] == hashlib.sha256(
        working_tree_status.encode("utf-8")
    ).hexdigest()
    candidates = Path(metadata["candidates"])
    for key, path in (
        ("validator_sha256", VALIDATOR),
        ("generator_sha256", GENERATOR),
        ("group_edges_sha256", GROUP_EDGES),
        ("ledger_sha256", LEDGER),
        ("metadata_sha256", METADATA),
        ("candidates_sha256", candidates),
        ("primordial_sha256", ROOT / "snrt_core" / "primordial.py"),
        ("source_ledger_sha256", ROOT / "snrt_core" / "source_ledger.py"),
        ("source_rebind_tool_sha256", SOURCE_REBIND_TOOL),
        ("p4_runner_sha256", P4_RUNNER),
        ("static_input_sha256", STATIC_INPUT),
        ("static_metadata_sha256", STATIC_METADATA),
        ("transport_control_sha256", TRANSPORT_CONTROL),
        ("external_asset_manifest_sha256", EXTERNAL_ASSET_MANIFEST),
    ):
        assert provenance[key] == sha256(path)
    assert provenance["source_sed_sha256"] == sha256(EXPLICIT_SED)


def main() -> int:
    validator = ROOT / "tools" / "validate_agn_nine_group_ledger.py"
    canonical = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    metadata = json.loads(METADATA.read_text(encoding="utf-8"))
    assert canonical["passed"] is True
    assert all(canonical["criteria"].values())
    assert canonical["configuration"]["source_mode"] == "explicit"
    assert canonical["configuration"]["working_tree_attestation_scope"] == ATTESTATION_SCOPE
    assert canonical["ledger"]["source_count"] == 10
    assert metadata["source_model_status"] == "synthetic_non_physical_wiring_fixture"
    assert_current_provenance(canonical, metadata)
    with TemporaryDirectory(prefix="agn-explicit-artifact-test-") as directory:
        output = Path(directory) / "validation.json"
        result = subprocess.run(
            [
                sys.executable,
                str(validator),
                "--source-mode",
                "explicit",
                "--output",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stdout + result.stderr)
        payload = json.loads(output.read_text(encoding="utf-8"))
        assert payload["passed"] is True
        assert all(payload["criteria"].values())
        assert payload["configuration"]["source_mode"] == "explicit"
        assert_current_provenance(payload, metadata)
        assert payload["criteria"] == canonical["criteria"]
    print("AGN_NINE_GROUP_EXPLICIT_ARTIFACT_OK mode=explicit criteria=all_true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
