#!/usr/bin/env python3
"""Merge photon ledgers without losing source-SED provenance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
if str(SNRT_ROOT) not in sys.path:
    sys.path.insert(0, str(SNRT_ROOT))

from snrt_core.primordial import group_spectral_closure_from_metadata  # noqa: E402
from snrt_core.source_ledger import read_photon_source_ledger_csv  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _metadata_edges(metadata: dict[str, object], path: Path) -> np.ndarray:
    try:
        edges = np.asarray(metadata["group_edges_ev"], dtype=np.float64)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: metadata lacks group_edges_ev") from error
    if (
        edges.ndim != 1
        or edges.size < 2
        or not np.isfinite(edges).all()
        or np.any(edges <= 0.0)
        or np.any(np.diff(edges) <= 0.0)
    ):
        raise ValueError(f"{path}: group edges are invalid")
    return edges


def _metadata_group_totals(metadata: dict[str, object], path: Path, count: int) -> np.ndarray:
    try:
        totals = np.asarray(metadata["group_photon_rate_total_s"], dtype=np.float64)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: metadata lacks group_photon_rate_total_s") from error
    if totals.shape != (count,) or not np.isfinite(totals).all() or np.any(totals < 0.0):
        raise ValueError(f"{path}: group photon totals are invalid")
    return totals


def _component_identity(metadata: dict[str, object], ledger_path: Path) -> dict[str, object]:
    contract = metadata.get("source_sed_contract")
    identity = metadata.get("source_sed_identity")
    source_hash = metadata.get("source_sed_sha256")
    if identity is not None and (not isinstance(identity, str) or len(identity) != 64):
        raise ValueError(f"{ledger_path}: source_sed_identity must be a SHA-256 identity or null")
    if source_hash is not None and (not isinstance(source_hash, str) or len(source_hash) != 64):
        raise ValueError(f"{ledger_path}: source_sed_sha256 must be a SHA-256 or null")
    return {
        "ledger_sha256": _sha256(ledger_path),
        "source_sed_identity": identity,
        "source_sed_sha256": source_hash,
        "source_sed_status": (
            contract.get("status") if isinstance(contract, dict) else "unspecified"
        ),
    }


def merge_ledgers(
    ledger_paths: list[Path],
    metadata_paths: list[Path],
    output: Path,
    metadata_output: Path,
) -> dict[str, object]:
    """Merge ledgers and construct an aggregate group spectral closure."""

    if len(ledger_paths) < 2 or len(ledger_paths) != len(metadata_paths):
        raise ValueError("supply at least two ledgers and one metadata file for each ledger")
    if output.exists() or metadata_output.exists():
        raise FileExistsError("refusing to overwrite merged ledger outputs")

    ledgers = []
    metadata_list: list[dict[str, object]] = []
    edges_reference: np.ndarray | None = None
    component_totals: list[np.ndarray] = []
    component_closures = []
    components: list[dict[str, object]] = []
    component_scale_factors: list[float | None] = []
    all_ids: list[int] = []
    for ledger_path, metadata_path in zip(ledger_paths, metadata_paths, strict=True):
        ledger = read_photon_source_ledger_csv(ledger_path)
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if not isinstance(metadata, dict):
            raise ValueError(f"{metadata_path}: metadata root must be an object")
        edges = _metadata_edges(metadata, metadata_path)
        if edges_reference is None:
            edges_reference = edges
        elif not np.array_equal(edges_reference, edges):
            raise ValueError("all component photon ledgers must use exactly identical group edges")
        if ledger.photon_luminosity_s.shape[1] != edges.size - 1:
            raise ValueError(f"{ledger_path}: CSV group count does not match metadata")
        totals = _metadata_group_totals(metadata, metadata_path, edges.size - 1)
        csv_totals = ledger.photon_luminosity_s.sum(axis=0, dtype=np.float64)
        if not np.allclose(csv_totals, totals, rtol=2.0e-12, atol=1.0e-30):
            raise ValueError(f"{ledger_path}: CSV and metadata group totals disagree")
        component_closure = group_spectral_closure_from_metadata(metadata)
        scale_factor_value = metadata.get("source_scale_factor")
        if scale_factor_value is None:
            component_scale_factors.append(None)
        else:
            try:
                scale_factor = float(scale_factor_value)
            except (TypeError, ValueError) as error:
                raise ValueError(f"{metadata_path}: source_scale_factor is not numeric") from error
            if not np.isfinite(scale_factor) or scale_factor <= 0.0:
                raise ValueError(f"{metadata_path}: source_scale_factor is invalid")
            component_scale_factors.append(scale_factor)
        ledgers.append(ledger)
        metadata_list.append(metadata)
        component_totals.append(totals)
        component_closures.append(component_closure)
        components.append(_component_identity(metadata, ledger_path))
        all_ids.extend(int(value) for value in ledger.source_id)

    if len(set(all_ids)) != len(all_ids):
        raise ValueError("component source_id values collide; IDs must be globally unique")
    if any(value is not None for value in component_scale_factors):
        if not all(value is not None for value in component_scale_factors):
            raise ValueError("all component photon ledgers must declare source scale factors together")
        scale_factors = np.asarray(component_scale_factors, dtype=np.float64)
        if not np.allclose(scale_factors, scale_factors[0], rtol=1.0e-12, atol=1.0e-14):
            raise ValueError("component photon ledgers have different source scale factors")
        merged_scale_factor: float | None = float(scale_factors[0])
        epoch_status = "verified_coeval_source_scale_factor"
    else:
        merged_scale_factor = None
        epoch_status = "missing_source_scale_factor_not_production_eligible"
    assert edges_reference is not None
    group_count = edges_reference.size - 1
    total_photons = np.sum(np.asarray(component_totals), axis=0, dtype=np.float64)
    aggregate_energy = np.zeros(group_count, dtype=np.float64)
    aggregate_sigma = np.zeros((3, group_count), dtype=np.float64)
    aggregate_excess = np.zeros((3, group_count), dtype=np.float64)
    for totals, closure in zip(component_totals, component_closures, strict=True):
        aggregate_energy += totals * np.asarray(closure.photon_weighted_energy_ev, dtype=np.float64)
        sigma = np.asarray(
            [closure.cross_sections.hydrogen_i, closure.cross_sections.helium_i, closure.cross_sections.helium_ii],
            dtype=np.float64,
        )
        excess = np.asarray(closure.photoelectron_excess_energy_ev, dtype=np.float64)
        aggregate_sigma += totals[None, :] * sigma
        aggregate_excess += totals[None, :] * sigma * excess
    mean_energy = np.divide(
        aggregate_energy,
        total_photons,
        out=np.sqrt(edges_reference[:-1] * edges_reference[1:]),
        where=total_photons > 0.0,
    )
    excess_energy = np.divide(
        aggregate_excess,
        aggregate_sigma,
        out=np.zeros_like(aggregate_excess),
        where=aggregate_sigma > 0.0,
    )
    aggregate_sigma = np.divide(
        aggregate_sigma,
        total_photons[None, :],
        out=np.zeros_like(aggregate_sigma),
        where=total_photons[None, :] > 0.0,
    )

    rows = []
    for ledger in ledgers:
        rows.extend(
            (
                int(source_id),
                str(source_kind),
                *[float(value) for value in position],
                *[float(value) for value in luminosity],
            )
            for source_id, source_kind, position, luminosity in zip(
                ledger.source_id,
                ledger.source_kind,
                ledger.position_code,
                ledger.photon_luminosity_s,
                strict=True,
            )
        )
    rows.sort(key=lambda row: row[0])
    group_fields = [f"q_group_{index}_s" for index in range(group_count)]
    output_fields = ["source_id", "source_kind", "x_code", "y_code", "z_code", *group_fields]
    source_kinds = sorted({str(row[1]) for row in rows})
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.parent.mkdir(parents=True, exist_ok=True)
    mixed_identity = _canonical_sha256(
        {
            "schema": "snrt_mixed_source_identity_v1",
            "group_edges_ev": edges_reference.tolist(),
            "components": components,
            "group_photon_rate_total_s": total_photons.tolist(),
        }
    )
    metadata: dict[str, object] = {
        "schema": "snrt_mixed_photon_source_ledger_v1",
        "status": "candidate_mixed_source_ledger",
        "source_kind": "+".join(source_kinds),
        "source_kinds": source_kinds,
        "source_count": len(rows),
        "group_count": group_count,
        "group_edges_ev": edges_reference.tolist(),
        "group_interval_convention": "left_closed_right_open_except_final_closed",
        "group_photon_rate_total_s": total_photons.tolist(),
        "source_epoch": {
            "status": epoch_status,
            "scale_factor": merged_scale_factor,
            "component_scale_factors": component_scale_factors,
        },
        "source_id_policy": {
            "mode": "preserve_input_source_ids",
            "input_offsets": [0 for _ in ledgers],
            "collision_policy": "reject_before_write",
            "sort_order": "ascending_source_id",
        },
        "source_sed_identity": mixed_identity,
        "source_sed_sha256": None,
        "source_sed_contract": {
            "schema": "snrt_source_sed_v1",
            "status": "candidate_mixed_group_aggregate",
            "identity": mixed_identity,
            "normalization": "component_ledgers_in_source_photons_per_s",
            "component_identities": components,
            "requires_aggregate_dust_closure": True,
        },
        "component_ledgers": components,
        "groups": [
            {
                "index": int(index),
                "energy_interval_ev": [float(edges_reference[index]), float(edges_reference[index + 1])],
                "photon_weighted_mean_energy_ev": float(mean_energy[index]),
                "total_photon_rate_s": float(total_photons[index]),
                "closure_status": (
                    "mixed_source_aggregate" if total_photons[index] > 0.0 else "empty_source_group_zero_photons"
                ),
            }
            for index in range(group_count)
        ],
        "group_spectral_closure": {
            "method": "photon-rate-weighted aggregate of component absorber closures",
            "species_order": ["hydrogen_i", "helium_i", "helium_ii"],
            "cross_sections_cm2": {
                "hydrogen_i": aggregate_sigma[0].tolist(),
                "helium_i": aggregate_sigma[1].tolist(),
                "helium_ii": aggregate_sigma[2].tolist(),
            },
            "photoelectron_excess_energy_ev": {
                "hydrogen_i": excess_energy[0].tolist(),
                "helium_i": excess_energy[1].tolist(),
                "helium_ii": excess_energy[2].tolist(),
            },
            "group_status": [
                "mixed_source_aggregate" if value > 0.0 else "empty_source_group_zero_photons"
                for value in total_photons
            ],
        },
        "dust_binding": {
            "status": "requires_aggregate_source_sed_weighted_sidecar",
            "component_only_sidecars_allowed": False,
        },
        "limits": [
            "The aggregate group closure cannot reconstruct within-group spectral shape.",
            "A source-matched dust sidecar must be built from the aggregate continuous SED or an explicitly approved equivalent.",
            "This merge does not approve stellar/AGN physics or enable live feedback.",
        ],
    }
    if merged_scale_factor is not None:
        metadata["source_scale_factor"] = merged_scale_factor
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(output_fields)
        for row in rows:
            writer.writerow(row[:5] + row[5:])
    metadata["output_csv_sha256"] = _sha256(output)
    metadata["generator_sha256"] = _sha256(TOOL_PATH)
    metadata_output.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", action="append", type=Path, required=True)
    parser.add_argument("--metadata", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    args = parser.parse_args()
    try:
        metadata = merge_ledgers(args.ledger, args.metadata, args.output, args.metadata_output)
    except (FileExistsError, FileNotFoundError, OSError, ValueError) as error:
        parser.error(str(error))
    print(
        "MERGE_PHOTON_SOURCE_LEDGERS_OK "
        f"sources={metadata['source_count']} groups={len(metadata['group_edges_ev']) - 1} "
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
