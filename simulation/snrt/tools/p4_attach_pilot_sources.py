"""Attach an audited photon ledger to an existing P4 static gas input."""

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

from snrt_core.snapshot import StaticRTInput, read_static_rt_input, write_static_rt_input
from snrt_core.primordial import group_spectral_closure_from_metadata
from snrt_core.source_ledger import read_photon_source_ledger_csv


DEFAULT_GROUP_EDGES = PROJECT_ROOT / "config" / "p0_photon_group_edges_ev.txt"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _single_source_scale_factor(path: Path) -> float:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    values = {float(row["aexp"]) for row in rows}
    if len(values) != 1:
        raise ValueError("pilot photon ledger must have exactly one source scale factor")
    return values.pop()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gas-input", required=True)
    parser.add_argument("--zoom-manifest", required=True)
    parser.add_argument("--photon-ledger", required=True)
    parser.add_argument("--photon-metadata", required=True)
    parser.add_argument("--gas-metadata")
    parser.add_argument("--group-edges", default=str(DEFAULT_GROUP_EDGES))
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata-output", required=True)
    parser.add_argument("--gas-scale-factor", type=float)
    args = parser.parse_args()

    gas = read_static_rt_input(args.gas_input)
    manifest = json.loads(Path(args.zoom_manifest).read_text())
    final = manifest["final"]
    ledger = read_photon_source_ledger_csv(args.photon_ledger)
    photon_metadata_path = Path(args.photon_metadata)
    photon_metadata = json.loads(photon_metadata_path.read_text(encoding="utf-8"))
    spectral_closure = group_spectral_closure_from_metadata(
        photon_metadata, require_code_manifest=True
    )
    configured_edges = np.loadtxt(args.group_edges, comments="#", dtype=np.float64)
    metadata_edges = np.asarray(photon_metadata["group_edges_ev"], dtype=np.float64)
    if not np.array_equal(metadata_edges, configured_edges):
        raise ValueError(
            "photon-ledger group edges do not exactly match the configured P0 edges"
        )
    if ledger.photon_luminosity_s.shape[1] != len(
        spectral_closure.photon_weighted_energy_ev
    ):
        raise ValueError("photon CSV and spectral metadata group counts differ")
    declared_totals = np.asarray(
        photon_metadata["group_photon_rate_total_s"], dtype=np.float64
    )
    actual_totals = ledger.photon_luminosity_s.sum(axis=0, dtype=np.float64)
    if not np.allclose(actual_totals, declared_totals, rtol=2.0e-15, atol=0.0):
        raise ValueError("photon CSV totals disagree with the spectral metadata")
    sources = ledger.source_catalog_in_cube(final["left_edge_code"], float(final["width_code"]), gas.shape)
    if sources is None:
        raise ValueError("photon ledger has no sources in the static gas cube")
    output = StaticRTInput(
        grid=gas.grid,
        hydrogen_number_density_cm3=gas.hydrogen_number_density_cm3,
        helium_number_density_cm3=gas.helium_number_density_cm3,
        temperature_k=gas.temperature_k,
        dust_relative_abundance=gas.dust_relative_abundance,
        x_hii=gas.x_hii,
        x_heii=gas.x_heii,
        x_heiii=gas.x_heiii,
        sources=sources,
        velocity_cm_s=gas.velocity_cm_s,
        metallicity_solar=gas.metallicity_solar,
        dust_to_metal=gas.dust_to_metal,
        dust_relative_abundance_origin=gas.dust_relative_abundance_origin,
        x_h2=gas.x_h2,
        cell_level=gas.cell_level,
    )
    write_static_rt_input(args.output, output)

    gas_aexp = (
        float(manifest["thermal_initialization"]["scale_factor"])
        if args.gas_scale_factor is None
        else args.gas_scale_factor
    )
    source_aexp = _single_source_scale_factor(Path(args.photon_ledger))
    coeval = abs(source_aexp - gas_aexp) <= 1.0e-12
    gas_input_path = Path(args.gas_input).resolve()
    zoom_manifest_path = Path(args.zoom_manifest).resolve()
    gas_metadata_path = (
        None if args.gas_metadata is None else Path(args.gas_metadata).resolve()
    )
    photon_ledger_path = Path(args.photon_ledger).resolve()
    group_edges_path = Path(args.group_edges).resolve()
    metadata = {
        "schema": "snrt_static_source_rebind_v2",
        "purpose": (
            "coeval static RT source rebind; science status remains inherited from the gas input"
            if coeval
            else "non-coeval transport and chemistry pilot"
        ),
        "gas_input": str(gas_input_path),
        "gas_metadata": None if gas_metadata_path is None else str(gas_metadata_path),
        "zoom_manifest": str(zoom_manifest_path),
        "gas_scale_factor": gas_aexp,
        "photon_ledger": str(photon_ledger_path),
        "photon_metadata": str(photon_metadata_path.resolve()),
        "source_scale_factor": source_aexp,
        "delta_a_source_minus_gas": source_aexp - gas_aexp,
        "coeval_within_1e-12": coeval,
        "source_count": int(len(sources.cell_index)),
        "group_count": int(sources.photon_luminosity_s.shape[1]),
        "group_edges_ev": configured_edges.tolist(),
        "group_interval_convention": photon_metadata["group_interval_convention"],
        "provenance": {
            "source_rebind_tool_sha256": _sha256(Path(__file__).resolve()),
            "gas_input_sha256": _sha256(gas_input_path),
            "gas_metadata_sha256": (
                None
                if gas_metadata_path is None
                else _sha256(gas_metadata_path)
            ),
            "zoom_manifest_sha256": _sha256(zoom_manifest_path),
            "photon_ledger_sha256": _sha256(photon_ledger_path),
            "photon_metadata_sha256": _sha256(photon_metadata_path.resolve()),
            "group_edges_sha256": _sha256(group_edges_path),
        },
        "limits": [
            "Only the source catalogue is replaced; all gas, ionization, grid, and composition fields are copied unchanged.",
            "Scientific eligibility remains limited by the provenance and unresolved fields of the supplied gas input.",
        ],
    }
    metadata_path = Path(args.metadata_output)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(
        "P4_PILOT_ATTACH_OK "
        f"shape={gas.shape} sources={len(sources.cell_index)} groups={sources.photon_luminosity_s.shape[1]}"
    )


if __name__ == "__main__":
    main()
