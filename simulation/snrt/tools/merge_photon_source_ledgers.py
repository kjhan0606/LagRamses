#!/usr/bin/env python3
"""Merge compatible stellar/AGN photon ledgers into one aggregate source ledger."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.primordial import GroupSpectralClosure, group_spectral_closure_from_metadata
from snrt_core.source_ledger import PhotonSourceLedger, read_photon_source_ledger_csv


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _group_edges(metadata: dict[str, object], path: Path) -> np.ndarray:
    raw_edges = metadata.get("group_edges_ev")
    if raw_edges is not None:
        edges = np.asarray(raw_edges, dtype=np.float64)
    else:
        raw_groups = metadata.get("groups")
        if not isinstance(raw_groups, list) or not raw_groups:
            raise ValueError(f"{path}: metadata has no group_edges_ev or groups")
        try:
            intervals = np.asarray(
                [group["energy_interval_ev"] for group in raw_groups], dtype=np.float64  # type: ignore[index]
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"{path}: invalid group energy intervals") from error
        if intervals.shape != (len(raw_groups), 2):
            raise ValueError(f"{path}: group intervals have inconsistent dimensions")
        if not np.allclose(intervals[1:, 0], intervals[:-1, 1], rtol=0.0, atol=1.0e-12):
            raise ValueError(f"{path}: group intervals are not contiguous")
        edges = np.concatenate((intervals[:1, 0], intervals[:, 1]))
    if edges.ndim != 1 or edges.size < 2 or not np.isfinite(edges).all() or np.any(edges <= 0.0):
        raise ValueError(f"{path}: group edges must be finite and positive")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError(f"{path}: group edges must be strictly increasing")
    return edges


def _declared_group_totals(metadata: dict[str, object], path: Path, number_of_groups: int) -> np.ndarray:
    raw_totals = metadata.get("group_photon_rate_total_s")
    if raw_totals is not None:
        totals = np.asarray(raw_totals, dtype=np.float64)
    else:
        raw_groups = metadata.get("groups")
        if not isinstance(raw_groups, list):
            raise ValueError(f"{path}: metadata has no group photon totals")
        try:
            totals = np.asarray(
                [group["total_photon_rate_s"] for group in raw_groups], dtype=np.float64  # type: ignore[index]
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"{path}: invalid group photon totals") from error
    if totals.shape != (number_of_groups,) or not np.isfinite(totals).all() or np.any(totals < 0.0):
        raise ValueError(f"{path}: group photon totals are invalid")
    return totals


def _read_input(ledger_path: Path, metadata_path: Path) -> tuple[PhotonSourceLedger, np.ndarray, GroupSpectralClosure, dict[str, object]]:
    ledger = read_photon_source_ledger_csv(ledger_path)
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"could not read ledger metadata {metadata_path}") from error
    if not isinstance(metadata, dict):
        raise ValueError(f"{metadata_path}: metadata root must be an object")
    edges = _group_edges(metadata, metadata_path)
    if ledger.photon_luminosity_s.shape[1] != edges.size - 1:
        raise ValueError(
            f"{ledger_path}: CSV has {ledger.photon_luminosity_s.shape[1]} groups but metadata has {edges.size - 1}"
        )
    closure = group_spectral_closure_from_metadata(metadata)
    if closure.cross_sections.hydrogen_i.shape != (edges.size - 1,):
        raise ValueError(f"{metadata_path}: spectral closure does not match its group table")
    actual_totals = ledger.photon_luminosity_s.sum(axis=0, dtype=np.float64)
    declared_totals = _declared_group_totals(metadata, metadata_path, edges.size - 1)
    if not np.allclose(actual_totals, declared_totals, rtol=2.0e-10, atol=0.0):
        raise ValueError(
            f"{ledger_path}: CSV photon totals disagree with metadata; refusing to combine a stale closure"
        )
    return ledger, edges, closure, metadata


def _combine_closures(
    totals: np.ndarray,
    edges: np.ndarray,
    inputs: list[tuple[np.ndarray, GroupSpectralClosure]],
) -> GroupSpectralClosure:
    number_of_groups = edges.size - 1
    mean_energy = np.sqrt(edges[:-1] * edges[1:])
    sigma_numerator = np.zeros((3, number_of_groups), dtype=np.float64)
    excess_numerator = np.zeros_like(sigma_numerator)
    mean_energy_numerator = np.zeros(number_of_groups, dtype=np.float64)
    for source_totals, closure in inputs:
        source_mean_energy = np.asarray(closure.photon_weighted_energy_ev, dtype=np.float64)
        source_sigma = np.asarray(
            [closure.cross_sections.hydrogen_i, closure.cross_sections.helium_i, closure.cross_sections.helium_ii],
            dtype=np.float64,
        )
        source_excess = np.asarray(closure.photoelectron_excess_energy_ev, dtype=np.float64)
        mean_energy_numerator += source_totals * source_mean_energy
        sigma_numerator += source_totals[None, :] * source_sigma
        excess_numerator += source_totals[None, :] * source_sigma * source_excess
    populated = totals > 0.0
    mean_energy[populated] = mean_energy_numerator[populated] / totals[populated]
    sigma = np.zeros_like(sigma_numerator)
    sigma[:, populated] = sigma_numerator[:, populated] / totals[None, populated]
    excess = np.zeros_like(excess_numerator)
    nonzero_sigma = sigma_numerator > 0.0
    excess[nonzero_sigma] = excess_numerator[nonzero_sigma] / sigma_numerator[nonzero_sigma]
    return GroupSpectralClosure(
        cross_sections=type(closure.cross_sections)(
            hydrogen_i=sigma[0],
            helium_i=sigma[1],
            helium_ii=sigma[2],
        ),
        photon_weighted_energy_ev=mean_energy,
        photoelectron_excess_energy_ev=excess,
    )


def _source_kind(ledgers: list[PhotonSourceLedger]) -> str:
    kinds = sorted({str(kind) for ledger in ledgers for kind in ledger.source_kind})
    if not kinds:
        raise ValueError("input ledgers contain no source kinds")
    return "+".join(kinds)


def merge_ledgers(
    ledger_paths: list[Path],
    metadata_paths: list[Path],
    output_path: Path,
    metadata_output_path: Path,
    allow_mixed_epochs: bool = False,
    source_id_offsets: list[int] | None = None,
) -> dict[str, object]:
    if len(ledger_paths) < 2:
        raise ValueError("at least two source ledgers are required")
    if len(ledger_paths) != len(metadata_paths):
        raise ValueError("--ledger and --metadata must contain the same number of paths")
    if source_id_offsets is None:
        source_id_offsets = [0] * len(ledger_paths)
    if len(source_id_offsets) != len(ledger_paths):
        raise ValueError("--source-id-offset must contain one integer per --ledger input")
    if any(not isinstance(offset, (int, np.integer)) for offset in source_id_offsets):
        raise ValueError("source ID offsets must be integers")
    if output_path.exists() or metadata_output_path.exists():
        raise FileExistsError("refusing to overwrite an existing merged ledger or metadata file")
    for path in (*ledger_paths, *metadata_paths):
        if not path.is_file():
            raise FileNotFoundError(path)

    records = [_read_input(ledger_path, metadata_path) for ledger_path, metadata_path in zip(ledger_paths, metadata_paths, strict=True)]
    reference_edges = records[0][1]
    for ledger_path, edges in zip(ledger_paths[1:], (record[1] for record in records[1:]), strict=True):
        if edges.shape != reference_edges.shape or not np.allclose(edges, reference_edges, rtol=0.0, atol=1.0e-12):
            raise ValueError(f"{ledger_path}: group table does not match the first input ledger")

    declared_scales: list[float | None] = []
    for metadata_path, record in zip(metadata_paths, records, strict=True):
        metadata = record[3]
        raw_scale = metadata.get("source_scale_factor")
        if raw_scale is None:
            declared_scales.append(None)
            continue
        try:
            scale = float(raw_scale)
        except (TypeError, ValueError) as error:
            raise ValueError(f"{metadata_path}: source_scale_factor is not numeric") from error
        if not np.isfinite(scale) or scale <= 0.0:
            raise ValueError(f"{metadata_path}: source_scale_factor must be finite and positive")
        if metadata.get("source_scale_factor_uniform") is False and not allow_mixed_epochs:
            raise ValueError(
                f"{metadata_path}: source ledger spans multiple scale factors; "
                "supply --allow-mixed-epochs only for an explicitly non-coeval control"
            )
        declared_scales.append(scale)
    known_scales = np.asarray([scale for scale in declared_scales if scale is not None], dtype=np.float64)
    if (
        known_scales.size > 1
        and np.ptp(known_scales) > 1.0e-12
        and not allow_mixed_epochs
    ):
        raise ValueError(
            "input ledgers declare different source scale factors; "
            "a coeval STAR+AGN merge requires matching epochs; "
            "supply --allow-mixed-epochs only for a labeled integration control"
        )
    mixed_epoch_control = bool(
        allow_mixed_epochs and known_scales.size > 1 and np.ptp(known_scales) > 1.0e-12
    )

    ledgers = [record[0] for record in records]
    source_ids = np.concatenate(
        [ledger.source_id.astype(np.int64) + np.int64(offset) for ledger, offset in zip(ledgers, source_id_offsets, strict=True)]
    )
    if np.unique(source_ids).size != source_ids.size:
        raise ValueError("source_id values collide across input ledgers; remap them explicitly before merging")
    positions = np.concatenate([ledger.position_code for ledger in ledgers], axis=0)
    kinds = np.concatenate([ledger.source_kind for ledger in ledgers])
    luminosities = np.concatenate([ledger.photon_luminosity_s for ledger in ledgers], axis=0)
    combined = PhotonSourceLedger(source_ids, kinds, positions, luminosities)
    totals = combined.photon_luminosity_s.sum(axis=0, dtype=np.float64)
    closure = _combine_closures(
        totals,
        reference_edges,
        [(record[0].photon_luminosity_s.sum(axis=0, dtype=np.float64), record[2]) for record in records],
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source_id",
        "source_kind",
        "x_code",
        "y_code",
        "z_code",
        *[f"q_group_{index}_s" for index in range(reference_edges.size - 1)],
    ]
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for index in range(source_ids.size):
            writer.writerow(
                (
                    int(source_ids[index]),
                    str(kinds[index]),
                    *[f"{value:.17g}" for value in positions[index]],
                    *[f"{value:.17g}" for value in luminosities[index]],
                )
            )

    source_kind = _source_kind(ledgers)
    combined_groups = [
        {
            "index": int(index),
            "energy_interval_ev": [float(reference_edges[index]), float(reference_edges[index + 1])],
            "photon_weighted_mean_energy_ev": float(closure.photon_weighted_energy_ev[index]),
            "total_photon_rate_s": float(totals[index]),
            "closure_status": (
                "combined_spectral_closure" if totals[index] > 0.0 else "empty_combined_group_zero_photons"
            ),
        }
        for index in range(reference_edges.size - 1)
    ]
    metadata: dict[str, object] = {
        "schema": "snrt_combined_photon_source_ledger_v1",
        "status": "noncoeval_integration_control" if mixed_epoch_control else "complete_combined_photon_ledger",
        "source_kind": source_kind,
        "source_count": int(source_ids.size),
        "group_count": int(reference_edges.size - 1),
        "group_edges_ev": reference_edges.tolist(),
        "groups": combined_groups,
        "group_photon_rate_total_s": totals.tolist(),
        "group_spectral_closure": {
            "method": "photon-number-weighted composition of input group closures; absorber-weighted excess recomputed",
            "species_order": ["hydrogen_i", "helium_i", "helium_ii"],
            "cross_sections_cm2": {
                "hydrogen_i": np.asarray(closure.cross_sections.hydrogen_i).tolist(),
                "helium_i": np.asarray(closure.cross_sections.helium_i).tolist(),
                "helium_ii": np.asarray(closure.cross_sections.helium_ii).tolist(),
            },
            "photoelectron_excess_energy_ev": {
                "hydrogen_i": np.asarray(closure.photoelectron_excess_energy_ev[0]).tolist(),
                "helium_i": np.asarray(closure.photoelectron_excess_energy_ev[1]).tolist(),
                "helium_ii": np.asarray(closure.photoelectron_excess_energy_ev[2]).tolist(),
            },
            "group_status": [group["closure_status"] for group in combined_groups],
        },
        "input_ledgers": [
            {
                "ledger": str(ledger_path.resolve()),
                "ledger_sha256": _sha256(ledger_path),
                "metadata": str(metadata_path.resolve()),
                "metadata_sha256": _sha256(metadata_path),
                "source_kind": metadata.get("source_kind", "unspecified"),
                "source_count": metadata.get("source_count", "unspecified"),
                "source_id_offset": int(offset),
            }
            for ledger_path, metadata_path, offset, (_, _, _, metadata) in zip(
                ledger_paths, metadata_paths, source_id_offsets, records, strict=True
            )
        ],
        "source_id_policy": {
            "input_offsets": [int(offset) for offset in source_id_offsets],
            "interpretation": "source IDs are input IDs plus the recorded per-ledger offset",
        },
        "normalization": {
            "operation": "source-ledger concatenation only",
            "source_side_attenuation": "preserved from each input ledger; no new factor applied",
            "dust_attenuation": "not applied; delegated to the RT grid",
        },
        "epoch_policy": {
            "allow_mixed_epochs": allow_mixed_epochs,
            "declared_source_scale_factors": declared_scales,
            "interpretation": (
                "non-coeval integration control; not a science snapshot"
                if mixed_epoch_control
                else "coeval when all input ledgers declare the same scale factor"
            ),
        },
        "limits": [
            "All input ledgers must use the same group boundaries and unique integer source IDs.",
            "The aggregate closure is composed from each input's stored group moments; raw spectra are not reconstructed.",
            "No dust, scattering, IR re-emission, or hydrodynamic feedback is applied by this merger.",
        ],
        "output_csv_sha256": _sha256(output_path),
    }
    metadata_output_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_output_path.open("x", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", nargs="+", type=Path, required=True)
    parser.add_argument("--metadata", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    parser.add_argument(
        "--allow-mixed-epochs",
        action="store_true",
        help="allow a labeled non-coeval integration control; never use it as a science snapshot",
    )
    parser.add_argument(
        "--source-id-offset",
        nargs="+",
        type=int,
        help="one explicit integer ID offset per input ledger; defaults to zero for each",
    )
    args = parser.parse_args()
    try:
        metadata = merge_ledgers(
            args.ledger,
            args.metadata,
            args.output,
            args.metadata_output,
            allow_mixed_epochs=args.allow_mixed_epochs,
            source_id_offsets=args.source_id_offset,
        )
    except (FileExistsError, FileNotFoundError, OSError, ValueError) as error:
        parser.error(str(error))
    print(
        "PHOTON_SOURCE_LEDGER_MERGE_OK "
        f"sources={metadata['source_count']} groups={metadata['group_count']} source_kind={metadata['source_kind']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
