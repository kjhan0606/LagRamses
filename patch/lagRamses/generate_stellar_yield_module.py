#!/usr/bin/env python3
"""Generate an embedded Fortran yield-table module from the Phase 0 format.

The input is the common ASCII format documented in stellar_yield_tables.f90.
The generated module contains read-only data and a loader that populates the
runtime stellar_yield_table_t object.  The raw input remains the provenance
source and is never modified.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


N_ELEMENTS = 11
N_CHANNELS = 5
N_FIXED_COLUMNS = 10


def read_rows(path: Path):
    rows = []
    for line_number, raw_line in enumerate(path.read_text().splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        expected = N_FIXED_COLUMNS + 2 * N_ELEMENTS
        if len(fields) != expected:
            raise ValueError(
                f"{path}:{line_number}: expected {expected} columns, "
                f"found {len(fields)}"
            )
        try:
            channel = int(fields[0])
            values = [float(value) for value in fields[1:]]
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid numeric field") from exc
        if not 1 <= channel <= N_CHANNELS:
            raise ValueError(f"{path}:{line_number}: invalid channel {channel}")
        if values[0] <= 0.0 or values[1] < 0.0 or values[2] < 0.0:
            raise ValueError(f"{path}:{line_number}: invalid mass, metallicity, or age")
        if values[3] < 0.0 or values[4] < 0.0:
            raise ValueError(f"{path}:{line_number}: invalid returned or remnant mass")
        rows.append((channel, values))
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return rows


def fortran_real(value: float) -> str:
    return f"{value:.17e}".replace("e", "d")


def emit_vector(name: str, values, type_name: str = "real(stellar_dp)") -> str:
    lines = [f"  {type_name}, parameter :: {name}({len(values)}) = (/ &"]
    chunks = [values[index : index + 4] for index in range(0, len(values), 4)]
    for index, chunk in enumerate(chunks):
        suffix = ", &" if index + 1 < len(chunks) else " /)"
        lines.append("       " + ", ".join(chunk) + suffix)
    return "\n".join(lines)


def emit_matrix(name: str, rows, n_columns: int) -> str:
    """Emit a matrix in Fortran column-major order as a reshape expression."""
    values = []
    for column in range(n_columns):
        values.extend(rows[row][column] for row in range(len(rows)))
    lines = [
        f"  real(stellar_dp), parameter :: {name}({n_columns}, {len(rows)}) = &",
        f"       reshape((/ &",
    ]
    chunks = [values[index : index + 4] for index in range(0, len(values), 4)]
    for index, chunk in enumerate(chunks):
        suffix = ", &" if index + 1 < len(chunks) else " /), (/ "
        lines.append("       " + ", ".join(chunk) + suffix)
    lines[-1] += f"{n_columns}, {len(rows)} /))"
    return "\n".join(lines)


def generate(input_path: Path, output_path: Path) -> None:
    rows = read_rows(input_path)
    digest = hashlib.sha256(input_path.read_bytes()).hexdigest()

    channels = [str(row[0]) for row in rows]
    scalar_rows = [row[1] for row in rows]
    masses = [fortran_real(row[0]) for row in scalar_rows]
    metallicities = [fortran_real(row[1]) for row in scalar_rows]
    # The canonical file stores age_yr.  Embedded data follows the same
    # table-boundary conversion as the runtime ASCII loader.
    ages = [fortran_real(row[2] * 1.0e-9) for row in scalar_rows]
    returned = [fortran_real(row[3]) for row in scalar_rows]
    remnants = [fortran_real(row[4]) for row in scalar_rows]
    energies = [fortran_real(row[5]) for row in scalar_rows]
    momenta = [
        [fortran_real(row[6]), fortran_real(row[7]), fortran_real(row[8])]
        for row in scalar_rows
    ]
    ejected = [
        [fortran_real(value) for value in row[9 : 9 + N_ELEMENTS]]
        for row in scalar_rows
    ]
    net_yields = [
        [fortran_real(value) for value in row[9 + N_ELEMENTS : 9 + 2 * N_ELEMENTS]]
        for row in scalar_rows
    ]

    text = [
        "! Generated file. Do not edit by hand.",
        f"! Source: {input_path}",
        f"! SHA256: {digest}",
        "",
        "module stellar_yield_embedded_data",
        "  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements",
        "  use stellar_yield_tables, only: stellar_yield_table_t, clear_yield_table",
        "  implicit none",
        "",
        "  integer, parameter :: embedded_table_ok = 0",
        "  integer, parameter :: embedded_table_err_alloc = 1",
        f"  integer, parameter :: embedded_n_rows = {len(rows)}",
        f"  character(len=*), parameter :: embedded_source = '{input_path}'",
        f"  character(len=*), parameter :: embedded_sha256 = '{digest}'",
        "",
        emit_vector("embedded_channel", channels, "integer"),
        emit_vector("embedded_initial_mass", masses),
        emit_vector("embedded_birth_metallicity", metallicities),
        emit_vector("embedded_age_gyr", ages),
        emit_vector("embedded_returned_mass", returned),
        emit_vector("embedded_remnant_mass", remnants),
        emit_vector("embedded_energy", energies),
        emit_matrix("embedded_momentum", momenta, 3),
        emit_matrix("embedded_ejected_mass", ejected, N_ELEMENTS),
        emit_matrix("embedded_net_yield", net_yields, N_ELEMENTS),
        "",
        "contains",
        "",
        "  subroutine load_embedded_yield_table(table, ierr)",
        "    type(stellar_yield_table_t), intent(inout) :: table",
        "    integer, intent(out) :: ierr",
        "    integer :: stat",
        "",
        "    call clear_yield_table(table)",
        "    allocate(table%channel(embedded_n_rows), &",
        "         table%initial_mass(embedded_n_rows), &",
        "         table%birth_metallicity(embedded_n_rows), &",
        "         table%age_gyr(embedded_n_rows), &",
        "         table%returned_mass(embedded_n_rows), &",
        "         table%remnant_mass(embedded_n_rows), &",
        "         table%energy(embedded_n_rows), &",
        "         table%momentum(embedded_n_rows,3), &",
        "         table%ejected_mass(embedded_n_rows,n_stellar_elements), &",
        "         table%net_yield(embedded_n_rows,n_stellar_elements), stat=stat)",
        "    if (stat /= 0) then",
        "       ierr = embedded_table_err_alloc",
        "       return",
        "    end if",
        "",
        "    table%channel = embedded_channel",
        "    table%initial_mass = embedded_initial_mass",
        "    table%birth_metallicity = embedded_birth_metallicity",
        "    table%age_gyr = embedded_age_gyr",
        "    table%returned_mass = embedded_returned_mass",
        "    table%remnant_mass = embedded_remnant_mass",
        "    table%energy = embedded_energy",
        "    table%momentum = transpose(embedded_momentum)",
        "    table%ejected_mass = transpose(embedded_ejected_mass)",
        "    table%net_yield = transpose(embedded_net_yield)",
        "    table%n_rows = embedded_n_rows",
        "    table%loaded = .true.",
        "    ierr = embedded_table_ok",
        "  end subroutine load_embedded_yield_table",
        "",
        "end module stellar_yield_embedded_data",
        "",
    ]
    output_path.write_text("\n".join(text))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    generate(args.input, args.output)


if __name__ == "__main__":
    main()
