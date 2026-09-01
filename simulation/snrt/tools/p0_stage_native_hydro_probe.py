#!/usr/bin/env python3
"""Stage a small, explicitly mapped canonical grid from a native RAMSES output.

This is a reproducible interface preflight for the stopped Phase 0 comparison
checkpoint.  It intentionally reads hydro/AMR records only: native particle
headers, sink catalogues, source luminosities, dust, and non-equilibrium
chemistry are not silently inferred.  The resulting HDF5 file is diagnostic
and is never marked production-ready.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.snapshot import RamsesFieldMap, read_static_rt_input, stage_ramses_hydro_only


NATIVE_HYDRO_FIELDS = (
    "Density",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "Pressure",
    "Metallicity",
    "scalar_01",
    "scalar_02",
    "scalar_03",
    "scalar_04",
    "scalar_05",
    "scalar_06",
    "scalar_07",
    "scalar_08",
    "scalar_09",
    "scalar_10",
    "scalar_11",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--info", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--level", type=int, default=9)
    parser.add_argument("--dims", type=int, default=8)
    parser.add_argument("--left-edge-code", type=float, nargs=3, default=(0.0, 0.0, 0.0))
    parser.add_argument("--width-code", type=float, default=0.015625)
    parser.add_argument("--mean-molecular-weight", type=float, default=0.6)
    args = parser.parse_args()

    if args.level < 0:
        raise ValueError("level must be non-negative")
    if args.dims <= 0:
        raise ValueError("dims must be positive")
    if args.width_code <= 0.0:
        raise ValueError("width-code must be positive")
    if args.mean_molecular_weight <= 0.0:
        raise ValueError("mean-molecular-weight must be positive")
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {args.output}")
    if args.manifest.exists():
        raise FileExistsError(f"refusing to overwrite existing manifest: {args.manifest}")

    left_edge = np.asarray(args.left_edge_code, dtype=np.float64)
    right_edge = left_edge + args.width_code
    if np.any(left_edge < 0.0) or np.any(right_edge > 1.0):
        raise ValueError("probe cube must lie within the unit RAMSES domain")

    fields = RamsesFieldMap(
        density=("gas", "density"),
        thermal_pressure=("gas", "pressure"),
        mean_molecular_weight=args.mean_molecular_weight,
        velocity=(("gas", "velocity_x"), ("gas", "velocity_y"), ("gas", "velocity_z")),
    )
    staged = stage_ramses_hydro_only(
        args.info,
        args.output,
        level=args.level,
        dimensions=(args.dims,) * 3,
        fields=fields,
        scratch_directory=args.scratch,
        hydro_fields_in_file=NATIVE_HYDRO_FIELDS,
        left_edge_code=tuple(left_edge),
        right_edge_code=tuple(right_edge),
    )
    reloaded = read_static_rt_input(args.output)
    try:
        reloaded.validate_production_contract(require_sources=False)
    except ValueError as exc:
        contract_error = str(exc)
    else:
        raise AssertionError("native hydro probe unexpectedly passed the production contract")

    manifest = stage_manifest(args, reloaded, args.output, contract_error)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "P0_NATIVE_HYDRO_PROBE_OK "
        f"shape={reloaded.shape} temperature_K={reloaded.temperature_k.min():.6g}..{reloaded.temperature_k.max():.6g}"
    )


def stage_manifest(
    args: argparse.Namespace,
    staged,
    output_path: Path,
    production_contract_error: str,
) -> dict[str, object]:
    density = np.asarray(staged.hydrogen_number_density_cm3)
    temperature = np.asarray(staged.temperature_k)
    velocity = np.asarray(staged.velocity_cm_s)
    return {
        "record_type": "snrt_native_ramses_hydro_probe",
        "schema_version": 1,
        "status": "canonical_hydro_interface_preflight",
        "production_ready": False,
        "production_contract": {
            "passed": False,
            "error": production_contract_error,
        },
        "source_info": str(args.info.resolve()),
        "output": str(output_path.resolve()),
        "output_sha256": _sha256(output_path),
        "geometry": {
            "level": args.level,
            "shape": list(staged.shape),
            "left_edge_code": [float(value) for value in args.left_edge_code],
            "width_code": float(args.width_code),
            "right_edge_code": [
                float(value + args.width_code) for value in args.left_edge_code
            ],
            "cell_width_cm": float(staged.grid.cell_width_cm),
        },
        "native_hydro_field_order": list(NATIVE_HYDRO_FIELDS),
        "thermal_initialization": {
            "source": "native thermal_pressure",
            "mean_molecular_weight": float(args.mean_molecular_weight),
            "hydrogen_mass_fraction": 0.76,
            "thermal_atlas": None,
            "scale_factor": None,
        },
        "field_ranges": {
            "hydrogen_number_density_cm3": [float(density.min()), float(density.max())],
            "temperature_k": [float(temperature.min()), float(temperature.max())],
            "velocity_abs_cm_s": float(np.abs(velocity).max()),
        },
        "not_available_in_probe": [
            "stellar/AGN source catalogue and photon luminosity ledger",
            "dust-to-metal field and dust relative abundance calibration",
            "certified solar-normalized metallicity field",
            "non-equilibrium H/He/H2 chemistry fields",
            "particle and sink records (hydro-only view)",
        ],
        "notes": [
            "The source snapshot is not modified.",
            "This artifact validates only native AMR/hydro reading, pressure mapping, cgs conversion, and canonical HDF5 serialization.",
            "The native descriptor is old/unversioned; the explicit ordered field list prevents yt from applying its NVAR>11 MHD fallback.",
        ],
    }


if __name__ == "__main__":
    main()
