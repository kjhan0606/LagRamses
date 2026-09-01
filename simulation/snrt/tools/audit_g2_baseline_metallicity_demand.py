#!/usr/bin/env python3
"""Audit inherited stellar-particle metallicity demand against staged G2 domains."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_baseline_metallicity_demand_contract_v1.json"


class BaselineMetallicityDemandError(ValueError):
    """The baseline metallicity-demand evidence violates its contract."""


def _hash(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise BaselineMetallicityDemandError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BaselineMetallicityDemandError(f"cannot read contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-baseline-metallicity-demand-contract"
        or contract.get("schema_version") != 1
    ):
        raise BaselineMetallicityDemandError("unsupported baseline metallicity-demand contract")
    policy = contract.get("policy", {})
    required_false = (
        "comparison_population_defines_production_domain",
        "metallicity_floor_or_clamp_allowed",
        "solar_source_extrapolation_to_ultra_low_z_allowed",
        "discrete_zero_metallicity_event_anchor_counts_as_full_grid",
        "unstaged_article_values_may_fill_domain",
        "required_production_birth_metallicity_domain_selected",
        "production_ready",
    )
    if any(policy.get(key) is not False for key in required_false):
        raise BaselineMetallicityDemandError("baseline-demand policy is not fail closed")
    if policy.get("canonical_rows_emitted") != 0:
        raise BaselineMetallicityDemandError("review contract unexpectedly emits canonical rows")
    reference = contract.get("staged_source_domain_reference", {})
    if reference.get("stockinger_zero_metallicity_model_is_discrete_event_anchor") is not True:
        raise BaselineMetallicityDemandError("Stockinger zero-metallicity scope is not preserved")
    if reference.get("jost_primordial_article_candidate_has_staged_yield_asset") is not False:
        raise BaselineMetallicityDemandError("unstaged primordial yield asset is asserted")
    return contract


def _finite(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise BaselineMetallicityDemandError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise BaselineMetallicityDemandError(f"{field} is not finite: {token!r}")
    return value


def _close(observed: float, expected: float, field: str) -> None:
    if not math.isclose(observed, expected, rel_tol=2e-13, abs_tol=1e-20):
        raise BaselineMetallicityDemandError(
            f"metadata/CSV {field} mismatch: observed={observed:.17g}, expected={expected:.17g}"
        )


def _quantile(values: list[float], fraction: float) -> float:
    index = math.floor((len(values) - 1) * fraction)
    return values[index]


def audit_g2_baseline_metallicity_demand(
    *, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    baseline = contract["baseline"]
    metadata_path = SNRT_ROOT / baseline["metadata_path"]
    catalogue_path = SNRT_ROOT / baseline["catalogue_path"]
    fingerprints: dict[str, dict[str, Any]] = {}
    for role, path, bytes_key, hash_key in (
        ("metadata", metadata_path, "metadata_bytes", "metadata_sha256"),
        ("catalogue", catalogue_path, "catalogue_bytes", "catalogue_sha256"),
    ):
        size, sha256 = _hash(path)
        if size != baseline[bytes_key] or sha256 != baseline[hash_key]:
            raise BaselineMetallicityDemandError(f"{role} fingerprint drifted")
        fingerprints[role] = {"path": str(path), "bytes": size, "sha256": sha256}

    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BaselineMetallicityDemandError(f"cannot parse stellar catalogue metadata: {exc}") from exc
    expected = contract["expected_metadata"]
    if metadata.get("status") != expected["status"]:
        raise BaselineMetallicityDemandError("stellar catalogue metadata status drifted")
    if metadata.get("catalogue_csv_sha256") != baseline["catalogue_sha256"]:
        raise BaselineMetallicityDemandError("metadata catalogue hash disagrees with contract")
    if metadata.get("field_provenance", {}).get("birth_metallicity") != expected["metallicity_source_field"]:
        raise BaselineMetallicityDemandError("birth-metallicity field provenance drifted")
    sanitization = metadata.get("sanitization", {})
    if sanitization.get("negative_birth_metallicity_clamped_count") != expected["negative_birth_metallicity_clamped_count"]:
        raise BaselineMetallicityDemandError("negative-metallicity sanitization count drifted")
    _close(
        float(sanitization.get("negative_birth_metallicity_clamp_tolerance")),
        expected["negative_birth_metallicity_clamp_tolerance"],
        "sanitization tolerance",
    )

    metallicities: list[float] = []
    ages: list[float] = []
    source_ids: set[int] = set()
    zero_count = 0
    near_floor_count = 0
    try:
        with catalogue_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            required_fields = {
                "source_id",
                "source_kind",
                baseline["birth_metallicity_field"],
                baseline["age_field"],
            }
            if reader.fieldnames is None or not required_fields.issubset(reader.fieldnames):
                raise BaselineMetallicityDemandError("stellar catalogue header drifted")
            for line_number, row in enumerate(reader, start=2):
                if row["source_kind"] != "star":
                    raise BaselineMetallicityDemandError(f"line {line_number}: non-star source")
                try:
                    source_id = int(row["source_id"])
                except ValueError as exc:
                    raise BaselineMetallicityDemandError(f"line {line_number}: invalid source_id") from exc
                if source_id in source_ids:
                    raise BaselineMetallicityDemandError(f"line {line_number}: duplicate source_id")
                source_ids.add(source_id)
                metallicity = _finite(row[baseline["birth_metallicity_field"]], "birth metallicity")
                age = _finite(row[baseline["age_field"]], "stellar age")
                if metallicity < 0.0 or age < 0.0:
                    raise BaselineMetallicityDemandError(f"line {line_number}: negative physical coordinate")
                metallicities.append(metallicity)
                ages.append(age)
                zero_count += metallicity == 0.0
                near_floor_count += 0.0 < metallicity <= 1.01e-50
    except OSError as exc:
        raise BaselineMetallicityDemandError(f"cannot parse stellar catalogue: {exc}") from exc

    if len(metallicities) != baseline["expected_star_count"] or len(source_ids) != len(metallicities):
        raise BaselineMetallicityDemandError("stellar catalogue count drifted")
    ordered_z = sorted(metallicities)
    z_min = ordered_z[0]
    z_max = ordered_z[-1]
    z_mean = math.fsum(ordered_z) / len(ordered_z)
    age_min = min(ages)
    age_max = max(ages)
    for observed, key in (
        (z_min, "birth_metallicity_minimum"),
        (z_max, "birth_metallicity_maximum"),
        (z_mean, "birth_metallicity_mean"),
        (age_min, "age_myr_minimum"),
        (age_max, "age_myr_maximum"),
    ):
        _close(observed, expected[key], key)

    candidate_minimum = float(
        contract["staged_source_domain_reference"]["lowest_positive_full_grid_candidate_metallicity_mass_fraction"]
    )
    below_count = sum(value < candidate_minimum for value in ordered_z)
    if below_count != len(ordered_z):
        raise BaselineMetallicityDemandError("baseline demand unexpectedly enters a positive-Z staged full grid")
    max_to_candidate_offset_dex = math.log10(candidate_minimum / z_max)
    return {
        "schema": "snrt-g2-baseline-metallicity-demand-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "comparison_population_ultra_low_z_uncovered",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "contract_path": str(contract_path),
        "baseline_identity": {
            "baseline_id": baseline["baseline_id"],
            "role": baseline["role"],
            "output": baseline["output"],
            "defines_production_domain": False,
        },
        "fingerprints": fingerprints,
        "stellar_population": {
            "star_count": len(ordered_z),
            "unique_source_id_count": len(source_ids),
            "birth_metallicity_mass_fraction": {
                "minimum": z_min,
                "maximum": z_max,
                "mean": z_mean,
                "quantiles_nearest_rank": {
                    "p01": _quantile(ordered_z, 0.01),
                    "p50": _quantile(ordered_z, 0.50),
                    "p90": _quantile(ordered_z, 0.90),
                    "p99": _quantile(ordered_z, 0.99),
                },
                "zero_count_after_recorded_sanitization": zero_count,
                "positive_at_or_below_1p01e_minus_50_count": near_floor_count,
            },
            "age_myr": {"minimum": age_min, "maximum": age_max},
            "metadata_csv_statistics_close": True,
        },
        "candidate_domain_comparison": {
            "lowest_positive_full_grid_candidate_metallicity_mass_fraction": candidate_minimum,
            "stars_below_lowest_positive_full_grid_candidate": below_count,
            "fraction_below_lowest_positive_full_grid_candidate": below_count / len(ordered_z),
            "maximum_baseline_z_to_candidate_lower_edge_offset_dex": max_to_candidate_offset_dex,
            "stockinger_zero_metallicity_model_is_discrete_event_anchor": True,
            "jost_primordial_yield_asset_staged": False,
            "positive_z_full_grid_covers_comparison_population": False,
        },
        "policy": contract["policy"],
        "interpretation": (
            "The inherited comparison population probes ultra-low birth metallicity, but it does "
            "not by itself select the future production domain. Every comparison star lies below "
            "the lowest staged positive-Z full-grid candidate node, and neither solar extrapolation "
            "nor a discrete zero-metallicity event anchor is accepted as coverage."
        ),
        "blockers": [
            "required_production_birth_metallicity_domain_is_not_selected",
            "comparison_population_is_below_all_staged_positive_z_full_grid_candidate_nodes",
            "no_staged_primordial_or_ultra_low_z_full_channel_yield_asset",
            "stockinger_zero_metallicity_event_anchor_is_not_a_full_grid",
            "metallicity_floor_clamp_and_solar_extrapolation_are_forbidden",
        ],
        "audit_code_sha256": _hash(TOOL_PATH)[1],
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_g2_baseline_metallicity_demand(contract_path=args.contract)
    except BaselineMetallicityDemandError as exc:
        report = {"schema": "snrt-g2-baseline-metallicity-demand-audit", "status": "error", "error": str(exc)}
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
