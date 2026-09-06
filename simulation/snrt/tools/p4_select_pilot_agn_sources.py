"""Select audited AGN RT-source candidates inside the P4 zoom-in cube."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def _read_rate_ledger(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("AGN rate ledger requires a header")
        required = {
            "sink_id",
            "aexp",
            "x_code",
            "y_code",
            "z_code",
            "mass_msun",
            "inflow_mdot_msun_per_year",
            "radiative_efficiency",
            "bolometric_luminosity_erg_s",
            "accretion_rate_convention",
        }
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"AGN rate ledger missing columns: {sorted(missing)}")
        return reader.fieldnames, list(reader)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--agn-rate-ledger", required=True)
    parser.add_argument("--zoom-manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata-output", required=True)
    args = parser.parse_args()

    rate_path = Path(args.agn_rate_ledger)
    manifest_path = Path(args.zoom_manifest)
    output_path = Path(args.output)
    metadata_path = Path(args.metadata_output)

    manifest = json.loads(manifest_path.read_text())
    final = manifest["final"]
    left = np.asarray(final["left_edge_code"], dtype=np.float64)
    width = float(final["width_code"])
    if left.shape != (3,) or width <= 0.0 or np.any(left < 0.0) or np.any(left + width > 1.0):
        raise ValueError("P4 zoom cube must be a non-wrapping code-coordinate cube")

    fieldnames, rows = _read_rate_ledger(rate_path)
    selected = []
    for row in rows:
        position = np.asarray([float(row[axis]) for axis in ("x_code", "y_code", "z_code")])
        if np.all(position >= left) and np.all(position < left + width):
            selected.append(row)
    selected.sort(key=lambda row: int(row["sink_id"]))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_fields = ["source_id", "source_kind", *fieldnames]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        for row in selected:
            writer.writerow({"source_id": row["sink_id"], "source_kind": "agn", **row})

    metadata = {
        "purpose": "AGN candidate selection before photon-group SED assignment",
        "agn_rate_ledger": str(rate_path.resolve()),
        "zoom_manifest": str(manifest_path.resolve()),
        "selection": {
            "left_edge_code": left.tolist(),
            "width_code": width,
            "shape": final["shape"],
            "non_wrapping": True,
        },
        "source_count": len(selected),
        "photon_group_luminosity": "not assigned; requires an explicit AGN SED and escape prescription",
        "accretion_rate_convention": (
            selected[0]["accretion_rate_convention"] if selected else "no sources in cube"
        ),
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"P4_PILOT_AGN_SELECTION_OK sources={len(selected)} output={output_path}")


if __name__ == "__main__":
    main()
