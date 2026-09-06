#!/usr/bin/env python3
"""Quantify mass omitted by the canonical eleven-species chemistry scope."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable

from adapt_g2_candidate_sources import (
    DEFAULT_CONTRACT,
    DEFAULT_ROOT,
    DEFAULT_SEMANTICS,
    LIMONGI_ID,
    NUGRID_ID,
    SourceAdapterError,
    adapt_candidate,
)


TOOL_PATH = Path(__file__).resolve()
TRACKED_ELEMENTS = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
_ISOTOPE = re.compile(r"^([A-Z][a-z]?)")


class ChemistryScopeError(ValueError):
    """Source species cannot be classified without an explicit policy."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _limongi_element(isotope: str) -> str:
    match = _ISOTOPE.match(isotope)
    if match is None:
        if isotope == "n":
            return "free_neutron"
        raise ChemistryScopeError(f"cannot identify parent element for Limongi isotope {isotope!r}")
    return match.group(1)


def _aggregate_species(species: dict[str, float], *, isotope_labels: bool) -> dict[str, float]:
    elements: dict[str, float] = defaultdict(float)
    for label, value in species.items():
        if not math.isfinite(value) or value < 0.0:
            raise ChemistryScopeError(f"invalid source yield for {label!r}")
        element = _limongi_element(label) if isotope_labels else label
        elements[element] += value
    return dict(elements)


def _summarize(records: list[dict[str, Any]], *, component: str) -> dict[str, Any]:
    omitted_totals: dict[str, float] = defaultdict(float)
    tracked_fractions: list[float] = []
    omitted_fractions: list[float] = []
    closure_residuals: list[float] = []
    zero_mass_records = 0
    for record in records:
        elements = record["elements"]
        source_total = math.fsum(elements.values())
        tracked = math.fsum(elements.get(element, 0.0) for element in TRACKED_ELEMENTS)
        omitted = math.fsum(value for element, value in elements.items() if element not in TRACKED_ELEMENTS)
        closure_residuals.append(source_total - tracked - omitted)
        if source_total == 0.0:
            zero_mass_records += 1
            continue
        tracked_fractions.append(tracked / source_total)
        omitted_fractions.append(omitted / source_total)
        for element, value in elements.items():
            if element not in TRACKED_ELEMENTS:
                omitted_totals[element] += value
    return {
        "component": component,
        "record_count": len(records),
        "zero_source_mass_record_count": zero_mass_records,
        "minimum_tracked_mass_fraction": min(tracked_fractions),
        "maximum_tracked_mass_fraction": max(tracked_fractions),
        "minimum_omitted_mass_fraction": min(omitted_fractions),
        "maximum_omitted_mass_fraction": max(omitted_fractions),
        "records_with_nonzero_omitted_mass": sum(value > 0.0 for value in omitted_fractions),
        "maximum_absolute_tracked_plus_omitted_closure_residual": max(
            abs(value) for value in closure_residuals
        ),
        "omitted_elements_by_unweighted_source_grid_sum": dict(
            sorted(omitted_totals.items(), key=lambda item: (-item[1], item[0]))
        ),
        "grid_sum_is_not_an_imf_weighted_population_quantity": True,
    }


def _limongi_components(report: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    recommended = {
        (
            record["source_model_coordinate"]["rotation_velocity_km_s"],
            record["source_model_coordinate"]["metallicity_feh"],
            record["source_model_coordinate"]["initial_mass_msun"],
        ): record["source_reported_isotopic_yields"]
        for record in report["source_components"]["recommended_yields"]["records"]
    }
    winds_13_25 = {
        (
            record["source_model_coordinate"]["rotation_velocity_km_s"],
            record["source_model_coordinate"]["metallicity_feh"],
            record["source_model_coordinate"]["initial_mass_msun"],
        ): record["source_reported_isotopic_yields"]
        for record in report["source_components"]["wind_yields"]["records"]
    }
    wind_records: list[dict[str, Any]] = []
    terminal_records: list[dict[str, Any]] = []
    for coordinate, final_yields in sorted(recommended.items()):
        mass = coordinate[2]
        if mass <= 25.0:
            if coordinate not in winds_13_25:
                raise ChemistryScopeError(f"missing Limongi wind coordinate {coordinate}")
            wind_yields = winds_13_25[coordinate]
            terminal_yields = {
                isotope: final_yields[isotope] - wind_yields[isotope]
                for isotope in final_yields
            }
            if any(value < 0.0 for value in terminal_yields.values()):
                raise ChemistryScopeError(f"negative Limongi terminal component at {coordinate}")
            terminal_records.append(
                {
                    "coordinate": coordinate,
                    "elements": _aggregate_species(terminal_yields, isotope_labels=True),
                }
            )
        else:
            wind_yields = final_yields
        wind_records.append(
            {
                "coordinate": coordinate,
                "elements": _aggregate_species(wind_yields, isotope_labels=True),
            }
        )
    return {"source_supported_wind": wind_records, "source_supported_terminal_set_R": terminal_records}


def _nugrid_components(report: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    components = report["source_components"]
    total_records = components["total"]["records"]
    wind_records = components["winds"]["records"]
    result: dict[str, list[dict[str, Any]]] = {
        "integrated_agb_ejecta_candidate": [],
        "massive_star_wind": [],
        "delayed_explosion_terminal_ejecta": [],
    }
    for total, winds in zip(total_records, wind_records, strict=True):
        total_coordinate = total["source_coordinate"]
        wind_coordinate = winds["source_coordinate"]
        if total_coordinate != wind_coordinate:
            raise ChemistryScopeError("NuGrid total/wind coordinates are not aligned")
        mass = total_coordinate["initial_mass_msun"]
        total_yields = total["source_reported_element_yields_msun"]
        wind_yields = winds["source_reported_element_yields_msun"]
        if mass <= 7.0:
            result["integrated_agb_ejecta_candidate"].append(
                {
                    "coordinate": total_coordinate,
                    "elements": _aggregate_species(total_yields, isotope_labels=False),
                }
            )
            continue
        terminal = {element: total_yields[element] - wind_yields[element] for element in total_yields}
        if any(value < 0.0 for value in terminal.values()):
            raise ChemistryScopeError(f"negative NuGrid terminal component at {total_coordinate}")
        result["massive_star_wind"].append(
            {
                "coordinate": wind_coordinate,
                "elements": _aggregate_species(wind_yields, isotope_labels=False),
            }
        )
        result["delayed_explosion_terminal_ejecta"].append(
            {
                "coordinate": total_coordinate,
                "elements": _aggregate_species(terminal, isotope_labels=False),
            }
        )
    return result


def audit_reduced_chemistry_scope(
    *,
    root: Path = DEFAULT_ROOT,
    manifest_path: Path | None = None,
    contract_path: Path = DEFAULT_CONTRACT,
    semantics_path: Path = DEFAULT_SEMANTICS,
) -> dict[str, Any]:
    limongi = adapt_candidate(
        LIMONGI_ID,
        root=root,
        manifest_path=manifest_path,
        contract_path=contract_path,
        semantics_path=semantics_path,
        include_records=True,
    )
    nugrid = adapt_candidate(
        NUGRID_ID,
        root=root,
        manifest_path=manifest_path,
        contract_path=contract_path,
        semantics_path=semantics_path,
        include_records=True,
    )
    candidate_components = {
        LIMONGI_ID: _limongi_components(limongi),
        NUGRID_ID: _nugrid_components(nugrid),
    }
    summaries = {
        candidate_id: {
            component: _summarize(records, component=component)
            for component, records in components.items()
        }
        for candidate_id, components in candidate_components.items()
    }
    maximum_omitted_fraction = max(
        component["maximum_omitted_mass_fraction"]
        for candidate in summaries.values()
        for component in candidate.values()
    )
    return {
        "schema": "snrt-g2-reduced-chemistry-scope-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "blocked_decay_horizon_and_source_approval_required",
        "production_ready": False,
        "tracked_elements": list(TRACKED_ELEMENTS),
        "radioactive_decay_applied": False,
        "candidate_component_summaries": summaries,
        "maximum_observed_omitted_mass_fraction": maximum_omitted_fraction,
        "canonical_ejecta_sum_equals_returned_mass_possible_with_only_tracked_elements": False,
        "untracked_ejecta_residual_contract_implemented": True,
        "runtime_mass_preservation_policy": (
            "returned_mass minus the eleven tracked-element masses is deposited into "
            "generic metallicity and never assigned to an individual element"
        ),
        "required_resolution": [
            "Approve one quantified radioactive-decay horizon before mapping Limongi parent isotopes into the reduced chemistry fields.",
            "Include the omitted-species fraction in publication uncertainty and abundance-scope statements."
        ],
        "interpretation": (
            "The source-supported component decomposition is used only to measure chemistry truncation. "
            "No IMF weighting, age history, canonical row, or runtime channel approval is produced."
        ),
        "adapter_code_sha256": limongi["adapter_code_sha256"],
        "source_semantics_evidence_sha256": limongi["source_semantics_evidence_sha256"],
        "audit_code_sha256": _sha256(TOOL_PATH),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--semantics", type=Path, default=DEFAULT_SEMANTICS)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_reduced_chemistry_scope(
            root=args.root,
            manifest_path=args.manifest,
            contract_path=args.contract,
            semantics_path=args.semantics,
        )
    except (SourceAdapterError, ChemistryScopeError) as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, indent=2), file=sys.stderr)
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
