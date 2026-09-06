"""Stage a density-selected RAMSES RT target without reading particle records.

This P4 driver selects the most massive interior sink only as a host-location
proxy, then selects the maximum-mean-density subcube from gas in a local probe.
It writes no AGN luminosities because the available sink CSV does not document
its final two columns.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.selection import select_high_density_region
from snrt_core.sink_catalog import read_sink_info
from snrt_core.snapshot import RamsesFieldMap, stage_ramses_hydro_only
from snrt_core.source_ledger import read_photon_source_ledger_csv
from snrt_core.thermal_atlas import read_thermal_atlas


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--info", required=True)
    parser.add_argument("--sink-info", required=True)
    parser.add_argument("--thermal-atlas", required=True)
    parser.add_argument("--scale-factor", required=True, type=float)
    parser.add_argument("--output", required=True)
    parser.add_argument("--scratch", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--sink-roster", required=True)
    parser.add_argument("--source-ledger")
    parser.add_argument("--level", type=int, default=5)
    parser.add_argument("--probe-dims", type=int, default=64)
    parser.add_argument("--final-dims", type=int, default=32)
    parser.add_argument("--metallicity-solar", type=float, default=1.0e-6)
    args = parser.parse_args()
    if args.probe_dims < args.final_dims or args.probe_dims % args.final_dims:
        raise ValueError("probe-dims must be an integer multiple of final-dims")
    if args.metallicity_solar <= 0.0:
        raise ValueError("metallicity-solar must be positive")

    root_cells = 1024
    probe_width_code = args.probe_dims / (root_cells * 2**args.level)
    catalog = read_sink_info(args.sink_info)
    host_index = catalog.most_massive_interior(probe_width_code / 2.0)
    probe_left_edge = catalog.position_code[host_index] - probe_width_code / 2.0
    thermal_atlas = read_thermal_atlas(args.thermal_atlas)
    fields = RamsesFieldMap(
        density=("gas", "density"),
        equilibrium_temperature=True,
        velocity=(("gas", "velocity_x"), ("gas", "velocity_y"), ("gas", "velocity_z")),
    )
    probe_path = Path(args.output).with_name("p4_density_probe.h5")
    probe = stage_ramses_hydro_only(
        args.info,
        probe_path,
        level=args.level,
        dimensions=(args.probe_dims,) * 3,
        fields=fields,
        scratch_directory=args.scratch,
        left_edge_code=tuple(float(value) for value in probe_left_edge),
        right_edge_code=tuple(float(value) for value in probe_left_edge + probe_width_code),
        thermal_atlas=thermal_atlas,
        thermal_scale_factor=args.scale_factor,
        grackle_metallicity_solar=args.metallicity_solar,
    )
    selection = select_high_density_region(
        probe.hydrogen_number_density_cm3,
        (args.final_dims,) * 3,
    )
    final_left_edge = probe_left_edge + np.asarray(selection.start_index) * probe_width_code / args.probe_dims
    final_width_code = args.final_dims * probe_width_code / args.probe_dims
    sources = None
    if args.source_ledger is not None:
        ledger = read_photon_source_ledger_csv(args.source_ledger)
        sources = ledger.source_catalog_in_cube(
            final_left_edge,
            final_width_code,
            (args.final_dims,) * 3,
        )
    final = stage_ramses_hydro_only(
        args.info,
        args.output,
        level=args.level,
        dimensions=(args.final_dims,) * 3,
        fields=fields,
        scratch_directory=args.scratch,
        left_edge_code=tuple(float(value) for value in final_left_edge),
        right_edge_code=tuple(float(value) for value in final_left_edge + final_width_code),
        thermal_atlas=thermal_atlas,
        thermal_scale_factor=args.scale_factor,
        grackle_metallicity_solar=args.metallicity_solar,
        sources=sources,
    )
    roster_indices = catalog.indices_in_box(final_left_edge, final_width_code)
    roster_path = Path(args.sink_roster)
    roster_path.parent.mkdir(parents=True, exist_ok=True)
    with roster_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("sink_id", "mass_msun", "x_code", "y_code", "z_code"))
        for index in roster_indices:
            writer.writerow((catalog.sink_id[index], catalog.mass_msun[index], *catalog.position_code[index]))
    manifest = {
        "snapshot_info": str(Path(args.info).resolve()),
        "selection": "maximum mean gas density in a fixed 32^3 subcube of a 64^3 host probe",
        "host_proxy": {
            "sink_id": int(catalog.sink_id[host_index]),
            "mass_msun": float(catalog.mass_msun[host_index]),
            "position_code": catalog.position_code[host_index].tolist(),
        },
        "final": {
            "left_edge_code": final_left_edge.tolist(),
            "width_code": float(final_width_code),
            "shape": final.shape,
            "mean_n_h_cm3": float(final.hydrogen_number_density_cm3.mean()),
            "max_n_h_cm3": float(final.hydrogen_number_density_cm3.max()),
            "sink_count": int(len(roster_indices)),
        },
        "thermal_initialization": {
            "model": "offline Grackle thermal-atlas net-rate-zero fallback",
            "atlas": str(Path(args.thermal_atlas).resolve()),
            "scale_factor": args.scale_factor,
            "metallicity_solar": args.metallicity_solar,
            "initial_rt_ionization": "neutral primordial; Grackle table does not export species fractions",
            "hydro_pressure": "unusable (zero throughout the selected AMR sample)",
        },
        "source_catalog": {
            "ledger": None if args.source_ledger is None else str(Path(args.source_ledger).resolve()),
            "photon_source_count": 0 if sources is None else int(len(sources.cell_index)),
            "status": (
                "not activated: this output has no audited stellar/AGN photon ledger"
                if args.source_ledger is None
                else "external photon ledger has no sources in the final RT cube"
                if sources is None
                else "external photon ledger deposited into the final RT cube"
            ),
        },
    }
    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        "P4_STAGE_OK "
        f"host_sink={manifest['host_proxy']['sink_id']} final_shape={final.shape} sinks={len(roster_indices)}"
    )


if __name__ == "__main__":
    main()
