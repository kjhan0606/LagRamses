#!/usr/bin/env python3
"""Record auditable evidence for a forced production-linked RAMSES build."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
TOOLS = ROOT / "simulation/snrt/tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from validate_stellar_source_parity import (  # noqa: E402
    DEFAULT_CONFIG,
    audit,
    binary_linkage_contract,
    git_head,
    git_status_without,
    production_smoke_contract,
    production_build_log_contract,
    sha256,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--build-log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()

    config_path = args.config.resolve()
    payload = audit(config_path)
    production = payload["production"]
    config = json.loads(config_path.read_text(encoding="utf-8"))
    harness = config["production_linked_harness"]
    patch_root = ROOT / config["production"]["patch_root"]
    required_sources = [
        object_name[:-2] + ".f90"
        for object_name in config["production"]["required_objects"]
    ]
    binary = args.binary.resolve()
    build_log = args.build_log.resolve()
    smoke_log = (ROOT / harness["smoke_log"]).resolve()
    expected_binary = (ROOT / harness["binary"]).resolve()
    expected_build_log = (ROOT / harness["build_log"]).resolve()
    output = args.output.resolve()
    expected_output = (ROOT / harness["evidence_file"]).resolve()
    if binary != expected_binary:
        raise SystemExit(f"unexpected production binary path: {binary}")
    if build_log != expected_build_log:
        raise SystemExit(f"unexpected production build-log path: {build_log}")
    if output != expected_output:
        raise SystemExit(f"unexpected production evidence path: {output}")
    if not binary.is_file():
        raise SystemExit(f"production binary is missing: {binary}")
    if not build_log.is_file():
        raise SystemExit(f"build log is missing: {build_log}")
    if not smoke_log.is_file():
        raise SystemExit(f"smoke log is missing: {smoke_log}")
    if payload["production"]["declared_patch_repo_relative"] != config["production"]["patch_root"]:
        raise SystemExit("Makefile PATCH does not resolve to the configured production tree")
    if production["missing_objects"] or production["missing_sources"]:
        raise SystemExit("production source/object contract is incomplete")
    binary_sha256 = sha256(binary)

    build_log_contract = production_build_log_contract(
        build_log.read_text(encoding="utf-8", errors="replace"),
        config["production"]["required_objects"],
        config["production"]["compile_parameters"],
        config["production"]["required_compile_flags"],
        config["production"]["embedded_yield_policy"]["macro"],
        harness["build_command"],
        str(Path("..") / config["production"]["patch_root"]),
        config["production"]["compile_policy"]["required_optimization_flag"],
        tuple(config["production"]["compile_policy"]["forbidden_compile_flags"]),
        harness["link_output"],
        binary_sha256,
    )
    if build_log_contract["status"] != "pass":
        raise SystemExit(f"production build log contract failed: {build_log_contract}")
    linkage_contract = binary_linkage_contract(
        binary, harness["linkage_symbol_patterns"]
    )
    if linkage_contract["status"] != "pass":
        raise SystemExit(f"production binary linkage contract failed: {linkage_contract}")
    smoke_contract = production_smoke_contract(
        smoke_log.read_text(encoding="utf-8", errors="replace"),
        harness["smoke_command"],
        harness["smoke_expected_output"],
        harness["smoke_expected_exit_code"],
        str(Path("..") / config["production"]["patch_root"]),
        git_head(),
        binary_sha256,
    )
    if smoke_contract["status"] != "pass":
        raise SystemExit(f"production smoke contract failed: {smoke_contract}")

    source_hashes = {
        source: sha256(patch_root / source) for source in required_sources
    }
    evidence = {
        "schema": "snrt-p0-production-linked-build-evidence-v1",
        "status": "pass",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository_head": git_head(),
        "worktree_status": git_status_without([output]),
        "source_root": config["production"]["patch_root"],
        "makefile": config["production"]["makefile"],
        "makefile_sha256": sha256(ROOT / config["production"]["makefile"]),
        "config_sha256": sha256(config_path),
        "harness_sha256": sha256(ROOT / harness["script"]),
        "validator_sha256": sha256(
            ROOT / "simulation/snrt/tools/validate_stellar_source_parity.py"
        ),
        "recorder_sha256": sha256(Path(__file__).resolve()),
        "runner_sha256": sha256(ROOT / config["g1_runner"]["script"]),
        "source_manifest_sha256": sha256(
            ROOT / config["production"]["source_manifest"]["path"]
        ),
        "required_objects": config["production"]["required_objects"],
        "production_source_hashes": source_hashes,
        "build_input_source_hashes": production["build_input_manifest"]["source_hashes"],
        "build_input_tree_sha256": production["build_input_manifest"]["tree_sha256"],
        "build_log_contract": build_log_contract,
        "binary_linkage_contract": linkage_contract,
        "smoke_contract": smoke_contract,
        "compile_parameters": build_log_contract["compile_parameters"],
        "required_compile_flags": config["production"]["required_compile_flags"],
        "build_command": harness["build_command"],
        "forced_rebuild": build_log_contract["forced_rebuild"],
        "binary_path": str(binary),
        "binary_sha256": binary_sha256,
        "build_log": str(build_log),
        "build_log_sha256": sha256(build_log),
        "smoke_log": str(smoke_log),
        "smoke_log_sha256": sha256(smoke_log),
        "stale_build_objects_are_not_evidence": True,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(f"P0_PRODUCTION_BUILD_EVIDENCE_RECORDED {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
