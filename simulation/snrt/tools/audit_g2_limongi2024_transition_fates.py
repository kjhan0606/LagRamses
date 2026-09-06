#!/usr/bin/env python3
"""Audit the Limongi et al. (2024) transition-fate evidence fail closed."""

from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import html
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_limongi2024_transition_fate_contract_v1.json"


class LimongiTransitionAuditError(ValueError):
    """The staged fate evidence violates its review-only contract."""


class _VisibleText(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self._hidden_depth = 0
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() in {"script", "style"}:
            self._hidden_depth += 1

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() in {"script", "style"} and self._hidden_depth:
            self._hidden_depth -= 1

    def handle_data(self, data: str) -> None:
        if not self._hidden_depth:
            self.parts.append(data)


def _hash(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise LimongiTransitionAuditError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise LimongiTransitionAuditError(f"cannot read contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-limongi2024-transition-fate-contract"
        or contract.get("schema_version") != 1
    ):
        raise LimongiTransitionAuditError("unsupported Limongi transition-fate contract")
    policy = contract.get("project_transition_policy", {})
    required_false = (
        "runtime_boundary_supported_as_universal_explosion_threshold",
        "continuous_fate_interpolation_allowed",
        "cross_metallicity_extrapolation_allowed",
        "cross_rotation_extrapolation_allowed",
        "stockinger_e8p8_anchor_may_define_population_fate_law",
        "canonical_yields_available",
        "canonical_energy_available",
        "canonical_momentum_available",
        "age_resolved_composition_history_available",
        "runtime_promotion_allowed",
    )
    if any(policy.get(key) is not False for key in required_false):
        raise LimongiTransitionAuditError("transition-fate policy is not fail closed")
    if policy.get("canonical_rows_emitted") != 0:
        raise LimongiTransitionAuditError("review contract unexpectedly emits canonical rows")
    approval = contract.get("approval", {})
    if approval.get("physics_fate_policy_selected") is not False:
        raise LimongiTransitionAuditError("review contract unexpectedly selects a fate policy")
    if approval.get("runtime_channel_assignment_approved") is not False:
        raise LimongiTransitionAuditError("review contract unexpectedly approves a runtime channel")
    return contract


def _article_text(path: Path) -> tuple[str, str]:
    try:
        raw = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise LimongiTransitionAuditError(f"cannot read article HTML: {exc}") from exc
    parser = _VisibleText()
    try:
        parser.feed(raw)
        parser.close()
    except Exception as exc:
        raise LimongiTransitionAuditError(f"cannot parse article HTML: {exc}") from exc
    visible = html.unescape(" ".join(parser.parts)).replace("−", "-")
    normalized = re.sub(r"\s+", " ", visible).strip()
    return raw, normalized


def _require_pattern(text: str, pattern: str, label: str) -> None:
    if re.search(pattern, text, flags=re.IGNORECASE) is None:
        raise LimongiTransitionAuditError(f"article evidence drifted: {label}")


def _audit_article(path: Path) -> dict[str, Any]:
    raw, text = _article_text(path)
    required_patterns = {
        "solar_nonrotating_scope": r"solar metallicity nonrotating stars.*?7.?15",
        "white_dwarf_interval": r"7\.50.*?8\.00.*?end their lives as hybrid CO/ONeMg or ONeMg WDs",
        "potential_ecsn_interval": r"8\.50.*?9\.20.*?may potentially explode as electron-capture supernovae",
        "ordinary_ccsn_boundary": r"initial mass.*?9\.22.*?end their lives as core-collapse supernovae",
        "extrapolated_tp_fate": r"several thousand thermal pulses.*?not feasible.*?extrapolated",
        "reference_ec_threshold": r"M\s*CO-ec\s*=\s*1\.415",
        "reference_ecsn_lower_mass": r"minimum mass that can potentially explode as an ECSN is.*?8\.5.?8\.8",
        "alternate_ec_threshold_and_mass": r"Zha et al.*?M\s*CO-ec\s*=\s*1\.36.*?minimum mass.*?8\.3",
        "high_sagb_fate_uncertainty": r"9\.10.*?9\.20.*?final fate.*?difficult to predict|9\.10.*?9\.20.*?difficult to predict.*?final fate",
    }
    for label, pattern in required_patterns.items():
        _require_pattern(text, pattern, label)
    license_markers = (
        "http://creativecommons.org/licenses/by/4.0/",
        "Creative Commons Attribution 4.0",
    )
    if any(marker not in raw for marker in license_markers):
        raise LimongiTransitionAuditError("CC BY 4.0 license evidence drifted")
    return {
        "scope_evidence_pass": True,
        "fate_statement_evidence_pass": True,
        "thermal_pulse_extrapolation_evidence_pass": True,
        "threshold_sensitivity_evidence_pass": True,
        "license_evidence_pass": True,
        "required_evidence_labels": sorted(required_patterns),
    }


def _finite(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise LimongiTransitionAuditError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise LimongiTransitionAuditError(f"{field} is not finite: {token!r}")
    return value


def _audit_table(path: Path, contract: dict[str, Any]) -> dict[str, Any]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise LimongiTransitionAuditError(f"cannot read Table 4: {exc}") from exc
    header = "\n".join(lines[:20])
    required_header = (
        "Table: Main properties of the TP phase for all models",
        "massModel  Mass of model",
        "TP         Pulse number",
        "dMco       Increase in mass of CO core",
    )
    if any(fragment not in header for fragment in required_header):
        raise LimongiTransitionAuditError("Table 4 header semantics drifted")

    rows_by_mass: dict[float, list[dict[str, float | int]]] = {}
    for line_number, raw in enumerate(lines, start=1):
        if re.match(r"^\d+\.\d{2}\s+\d+\s", raw) is None:
            continue
        fields = raw.split()
        if len(fields) != 10:
            raise LimongiTransitionAuditError(f"Table 4 line {line_number}: expected 10 fields")
        mass = _finite(fields[0], "model mass")
        try:
            pulse = int(fields[1])
        except ValueError as exc:
            raise LimongiTransitionAuditError(f"Table 4 line {line_number}: invalid pulse") from exc
        numeric = [_finite(value, f"Table 4 column {index}") for index, value in enumerate(fields[2:], start=3)]
        rows_by_mass.setdefault(mass, []).append(
            {
                "pulse": pulse,
                "pulse_duration_yr": numeric[0],
                "log_lhe": numeric[1],
                "maximum_he_shell_extent_msun": numeric[2],
                "maximum_he_shell_temperature_1e8k": numeric[3],
                "third_dredge_up_efficiency": numeric[4],
                "interpulse_time_yr": numeric[5],
                "he_core_growth_msun": numeric[6],
                "co_core_growth_msun": numeric[7],
            }
        )

    scope = contract["review_scope"]
    expected_masses = [float(value) for value in scope["table4_tp_model_masses_msun"]]
    if sorted(rows_by_mass) != expected_masses:
        raise LimongiTransitionAuditError("Table 4 model-mass grid drifted")
    expected_counts = {float(key): int(value) for key, value in scope["table4_expected_pulse_count_by_mass"].items()}
    observed_counts = {mass: len(rows) for mass, rows in rows_by_mass.items()}
    if observed_counts != expected_counts:
        raise LimongiTransitionAuditError("Table 4 pulse counts drifted")
    if sum(observed_counts.values()) != scope["table4_expected_data_row_count"]:
        raise LimongiTransitionAuditError("Table 4 data-row count drifted")
    for mass, rows in rows_by_mass.items():
        pulses = [int(row["pulse"]) for row in rows]
        if pulses != list(range(1, expected_counts[mass] + 1)):
            raise LimongiTransitionAuditError(f"Table 4 pulse sequence drifted at {mass} Msun")
    return {
        "title": "Main properties of the TP phase for all models",
        "data_row_count": sum(observed_counts.values()),
        "model_count": len(rows_by_mass),
        "model_masses_msun": sorted(rows_by_mass),
        "pulse_count_by_mass": {f"{mass:.2f}": observed_counts[mass] for mass in sorted(observed_counts)},
        "pulse_sequences_contiguous": True,
        "table_contains_terminal_event_yields": False,
    }


def audit_limongi2024_transition_fates(
    *, root: Path = DEFAULT_ROOT, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    root = Path(root).resolve()
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    source = contract["source"]
    base = root / source["release_root_relative_path"]
    fingerprints: dict[str, dict[str, Any]] = {}
    for name, expected in source["files"].items():
        size, sha256 = _hash(base / name)
        if size != expected["bytes"] or sha256 != expected["sha256"]:
            raise LimongiTransitionAuditError(f"staged transition-fate source fingerprint drifted: {name}")
        fingerprints[name] = {"bytes": size, "sha256": sha256}

    article = _audit_article(base / "ARTICLE.html")
    table = _audit_table(base / "TABLE4_TP_MRT.txt", contract)
    fate = contract["source_reported_fate_statements"]
    policy = contract["project_transition_policy"]
    blockers = [
        "solar_nonrotating_fate_evidence_only",
        "thermal_pulse_table_is_not_a_yield_or_event_table",
        "eight_to_eight_point_eight_msun_terminal_fate_is_not_deterministic",
        "discrete_fate_models_may_not_be_interpolated",
        "stockinger_e8p8_event_anchor_cannot_define_a_population_fate_law",
        "no_event_yields_energy_momentum_or_age_resolved_composition",
        "runtime_channel3_lower_boundary_at_eight_msun_is_not_approved",
        "project_fate_policy_and_runtime_channel_assignment_are_not_approved",
    ]
    return {
        "schema": "snrt-g2-limongi2024-transition-fate-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "review_only_fate_policy_blocked",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "contract_path": str(contract_path),
        "source_identity": {
            "candidate_id": source["candidate_id"],
            "article_doi": source["article_doi"],
            "metallicity": contract["review_scope"]["metallicity"],
            "rotation": contract["review_scope"]["rotation"],
            "license": contract["use_terms"]["article_license"],
            "license_verified": article["license_evidence_pass"],
        },
        "fingerprints": fingerprints,
        "article_evidence": article,
        "machine_readable_tp_table": table,
        "source_reported_fate_statements": fate,
        "project_transition_policy": policy,
        "interpretation": (
            "At source-solar metallicity and zero rotation, Limongi et al. constrain the transition "
            "but do not define a deterministic continuous explosion law. The 8--8.8 Msun runtime "
            "edge remains a fate-policy gap, not a yield interval that may be filled by interpolation."
        ),
        "blockers": blockers,
        "audit_code_sha256": _hash(TOOL_PATH)[1],
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_limongi2024_transition_fates(root=args.root, contract_path=args.contract)
    except LimongiTransitionAuditError as exc:
        report = {"schema": "snrt-g2-limongi2024-transition-fate-audit", "status": "error", "error": str(exc)}
        text = json.dumps(report, indent=2) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(text, encoding="utf-8")
        sys.stdout.write(text)
        return 1
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
