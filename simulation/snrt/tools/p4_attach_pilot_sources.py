"""Attach an audited photon ledger to an existing P4 static gas input."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.snapshot import StaticRTInput, read_static_rt_input, write_static_rt_input
from snrt_core.source_ledger import read_photon_source_ledger_csv


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
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata-output", required=True)
    parser.add_argument("--gas-scale-factor", type=float)
    args = parser.parse_args()

    gas = read_static_rt_input(args.gas_input)
    manifest = json.loads(Path(args.zoom_manifest).read_text())
    final = manifest["final"]
    ledger = read_photon_source_ledger_csv(args.photon_ledger)
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
    metadata = {
        "purpose": "transport and chemistry pilot only; not a coeval science input",
        "gas_input": str(Path(args.gas_input).resolve()),
        "gas_scale_factor": gas_aexp,
        "photon_ledger": str(Path(args.photon_ledger).resolve()),
        "source_scale_factor": source_aexp,
        "delta_a_source_minus_gas": source_aexp - gas_aexp,
        "source_count": int(len(sources.cell_index)),
        "group_count": int(sources.photon_luminosity_s.shape[1]),
        "limits": [
            "The output_00017 HDF5 snapshot is complete but has adaptive leaves across levels 10-15.",
            "A coeval input requires an adaptive-leaf resampler; this pilot intentionally retains output_00016 gas.",
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
