#!/usr/bin/env python3
"""Regression tests for F-P1 terminal deposition and ownership semantics."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_fp1_terminal_deposition_contract import (  # noqa: E402
    TerminalDepositionContractError,
    audit_terminal_deposition_contract,
)


def _expect_error(payload: dict, fragment: str) -> None:
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-deposition-") as directory:
        path = Path(directory) / "contract.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        try:
            audit_terminal_deposition_contract(contract_path=path)
        except TerminalDepositionContractError as exc:
            assert fragment in str(exc), str(exc)
        else:
            raise AssertionError(f"expected TerminalDepositionContractError containing {fragment!r}")


def main() -> int:
    path = ROOT / "config" / "fp1_terminal_deposition_contract_v1.json"
    current = json.loads(path.read_text(encoding="utf-8"))
    report = audit_terminal_deposition_contract()
    assert report["status"] == "review_only_contract_complete_physical_policy_unselected"
    assert report["production_ready"] is False
    assert report["runtime_deposition_allowed"] is False
    assert report["candidate_mass_range_msun"] == [8.0, 120.0]
    assert report["fate_filtered"] is True
    assert report["scalar_radial_momentum"] is None
    assert report["selected_deposition_mode"] is None
    assert report["ownership_closed"] is True

    unresolved = copy.deepcopy(current)
    unresolved["channel"]["unresolved_outcome_deposition_allowed"] = True
    _expect_error(unresolved, "must not deposit")

    inferred_energy = copy.deepcopy(current)
    inferred_energy["energy"]["infer_injected_from_asymptotic_or_diagnostic_allowed"] = True
    _expect_error(inferred_energy, "unsafe energy policy")

    inferred_momentum = copy.deepcopy(current)
    inferred_momentum["momentum"]["infer_scalar_from_energy_allowed"] = True
    _expect_error(inferred_momentum, "must not be inferred")

    arbitrary_vector = copy.deepcopy(current)
    arbitrary_vector["momentum"]["isotropic_source_vector"] = [1.0, 0.0, 0.0]
    _expect_error(arbitrary_vector, "exactly zero")

    selected_without_approval = copy.deepcopy(current)
    selected_without_approval["deposition"]["selected_mode"] = "thermal_pure"
    _expect_error(selected_without_approval, "must not select selected_mode")

    reused = copy.deepcopy(current)
    reused["ownership"]["cross_channel_packet_reuse_allowed"] = True
    _expect_error(reused, "reuse is forbidden")

    print("FP1_TERMINAL_DEPOSITION_CONTRACT_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
