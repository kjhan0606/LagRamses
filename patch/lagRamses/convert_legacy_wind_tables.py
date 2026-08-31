#!/usr/bin/env python3
"""Extract legacy AGB and massive-wind tables without inventing time data.

The legacy headers contain mass, metallicity, net elemental yields, and total
wind mass.  They do not contain a release-time axis or actual per-element
ejecta.  This tool therefore writes an intermediate CSV that preserves the
source values.  A later release-history step must add the age axis and convert
net yields to actual ejecta before the Phase 0 canonical table is generated.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ELEMENTS = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
LEGACY_ELEMENTS = ELEMENTS + ("Z_TOTAL",)
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
ENTRY = re.compile(
    rf"\{{\s*({NUMBER})\s*,\s*({NUMBER})\s*,\s*\{{([^{{}}]*)\}}\s*,"
    rf"\s*({NUMBER})\s*,\s*({NUMBER})\s*\}}",
    re.MULTILINE,
)


def parse_header(path: Path):
    text = path.read_text()
    rows = []
    for match in ENTRY.finditer(text):
        mass = float(match.group(1))
        metallicity = float(match.group(2))
        yields = [float(value) for value in re.findall(NUMBER, match.group(3))]
        if len(yields) != len(LEGACY_ELEMENTS):
            raise ValueError(
                f"{path}: entry at character {match.start()} has "
                f"{len(yields)} yields, expected {len(LEGACY_ELEMENTS)}"
            )
        rows.append(
            {
                "initial_mass": mass,
                "birth_metallicity": metallicity,
                "net_yields": yields,
                "wind_mass": float(match.group(4)),
                "wind_fraction": float(match.group(5)),
            }
        )
    if not rows:
        raise ValueError(f"{path}: no YieldTableEntry rows found")
    return rows


def write_intermediate(output: Path, sources):
    fieldnames = [
        "source",
        "channel",
        "initial_mass",
        "birth_metallicity",
        "wind_mass",
        "wind_fraction",
    ]
    fieldnames.extend(f"net_{element}" for element in LEGACY_ELEMENTS)
    with output.open("w", newline="") as stream:
        stream.write(
            "# Intermediate legacy extraction. No age axis is present.\n"
            "# Do not use directly as a time-dependent canonical table.\n"
        )
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for source_name, channel, rows in sources:
            for row in rows:
                record = {
                    "source": source_name,
                    "channel": channel,
                    "initial_mass": f"{row['initial_mass']:.17e}",
                    "birth_metallicity": f"{row['birth_metallicity']:.17e}",
                    "wind_mass": f"{row['wind_mass']:.17e}",
                    "wind_fraction": f"{row['wind_fraction']:.17e}",
                }
                record.update(
                    {
                        f"net_{element}": f"{value:.17e}"
                        for element, value in zip(LEGACY_ELEMENTS, row["net_yields"])
                    }
                )
                writer.writerow(record)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--agb", required=True, type=Path)
    parser.add_argument("--massive", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    sources = [
        (str(args.agb), 2, parse_header(args.agb)),
        (str(args.massive), 1, parse_header(args.massive)),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_intermediate(args.output, sources)


if __name__ == "__main__":
    main()
