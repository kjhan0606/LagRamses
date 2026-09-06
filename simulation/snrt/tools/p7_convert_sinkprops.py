#!/usr/bin/env python3
"""Convert RAMSES native sinkprops records into SNRT AGN source JSONL.

The sinkprops record is written at the beginning of AGN_feedback, before the
feedback routine clears the coarse-step accretion accumulators.  Its Bondi and
Eddington rates are therefore the required instantaneous RT source rates.
"""

from __future__ import annotations

import argparse
import json
import math
import struct
from pathlib import Path
from typing import BinaryIO


C_LIGHT_CGS = 2.99792458e10
MSUN_CGS = 1.98847e33
YEAR_S = 365.25 * 24.0 * 3600.0


class FortranSequentialReader:
    """Reader for one compiler-portable RAMSES sequential-unformatted file."""

    def __init__(self, stream: BinaryIO) -> None:
        start = stream.read(16)
        if len(start) < 8:
            raise ValueError("sinkprops file is too short")
        for endian in ("<", ">"):
            for marker_bytes, code in ((4, "i"), (8, "q")):
                size = struct.calcsize(endian + code)
                marker = struct.unpack(endian + code, start[:size])[0]
                if marker not in (4, 8):
                    continue
                payload = start[size : size + marker]
                tail = start[size + marker : size * 2 + marker]
                if len(payload) != marker or len(tail) != size:
                    continue
                if struct.unpack(endian + code, tail)[0] != marker:
                    continue
                nsink = struct.unpack(endian + ("i" if marker == 4 else "q"), payload)[0]
                if 0 <= nsink <= 100000000:
                    self.endian = endian
                    self.marker_bytes = marker_bytes
                    self.marker_code = code
                    stream.seek(0)
                    return
        raise ValueError("cannot determine Fortran record marker size or endianness")

    def read_record(self, stream: BinaryIO) -> bytes:
        marker_raw = stream.read(self.marker_bytes)
        if len(marker_raw) != self.marker_bytes:
            raise ValueError("unexpected end of sinkprops file")
        marker = struct.unpack(self.endian + self.marker_code, marker_raw)[0]
        if marker < 0:
            raise ValueError("negative Fortran record marker is unsupported")
        payload = stream.read(marker)
        tail_raw = stream.read(self.marker_bytes)
        if len(payload) != marker or len(tail_raw) != self.marker_bytes:
            raise ValueError("truncated Fortran record")
        if struct.unpack(self.endian + self.marker_code, tail_raw)[0] != marker:
            raise ValueError("mismatched Fortran record markers")
        return payload

    def integer(self, stream: BinaryIO) -> int:
        payload = self.read_record(stream)
        if len(payload) == 4:
            return struct.unpack(self.endian + "i", payload)[0]
        if len(payload) == 8:
            return struct.unpack(self.endian + "q", payload)[0]
        raise ValueError(f"unexpected integer record size: {len(payload)}")

    def scalar(self, stream: BinaryIO) -> float:
        payload = self.read_record(stream)
        if len(payload) == 8:
            return struct.unpack(self.endian + "d", payload)[0]
        if len(payload) == 4:
            return float(struct.unpack(self.endian + "f", payload)[0])
        raise ValueError(f"unexpected scalar record size: {len(payload)}")

    def float_array(self, stream: BinaryIO, count: int) -> list[float]:
        payload = self.read_record(stream)
        if len(payload) != 8 * count:
            raise ValueError(
                f"expected {count} double precision values, received {len(payload)} bytes"
            )
        return list(struct.unpack(self.endian + f"{count}d", payload))

    def int_array(self, stream: BinaryIO, count: int) -> list[int]:
        payload = self.read_record(stream)
        if len(payload) == 4 * count:
            return list(struct.unpack(self.endian + f"{count}i", payload))
        if len(payload) == 8 * count:
            return list(struct.unpack(self.endian + f"{count}q", payload))
        raise ValueError(f"unexpected idsink record size: {len(payload)}")


def read_sinkprops(
    path: Path, *, nstep_coarse: int
) -> tuple[dict[str, object], list[dict[str, object]]]:
    if nstep_coarse < 0:
        raise ValueError("nstep_coarse must be non-negative")
    with path.open("rb") as stream:
        reader = FortranSequentialReader(stream)
        nsink = reader.integer(stream)
        ndim = reader.integer(stream)
        if nsink < 0 or ndim != 3:
            raise ValueError(f"unsupported sinkprops dimensions: nsink={nsink}, ndim={ndim}")

        aexp = reader.scalar(stream)
        scale_l = reader.scalar(stream)
        scale_d = reader.scalar(stream)
        scale_t = reader.scalar(stream)
        ids = reader.int_array(stream, nsink)
        masses = reader.float_array(stream, nsink)
        positions = [reader.float_array(stream, nsink) for _ in range(ndim)]
        velocities = [reader.float_array(stream, nsink) for _ in range(ndim)]
        angular_momenta = [reader.float_array(stream, nsink) for _ in range(ndim)]
        bondi = reader.float_array(stream, nsink)
        eddington = reader.float_array(stream, nsink)
        accreted = reader.float_array(stream, nsink)
        gas_density = reader.float_array(stream, nsink)
        sound_speed = reader.float_array(stream, nsink)
        relative_speed = reader.float_array(stream, nsink)
        saved_energy = reader.float_array(stream, nsink)
        bh_spin = [reader.float_array(stream, nsink) for _ in range(ndim)]
        spin_magnitude = reader.float_array(stream, nsink)
        radiative_efficiency = reader.float_array(stream, nsink)
        time_code = reader.scalar(stream)

    rate_unit_cgs = scale_d * scale_l**3 / scale_t
    units = {
        "length_cgs": scale_l,
        "density_cgs": scale_d,
        "time_s": scale_t,
        "mass_cgs": scale_d * scale_l**3,
        "mass_rate_cgs": rate_unit_cgs,
    }
    records: list[dict[str, object]] = []
    for index, sink_id in enumerate(ids):
        inflow = min(bondi[index], eddington[index])
        efficiency = radiative_efficiency[index]
        required_fields = {
            "sink_mass_code": masses[index],
            "bondi_rate_code": bondi[index],
            "eddington_rate_code": eddington[index],
            "inflow_rate_code": inflow,
            "effective_radiative_efficiency": efficiency,
            **{f"position_{axis + 1}": positions[axis][index] for axis in range(ndim)},
            **{f"velocity_{axis + 1}": velocities[axis][index] for axis in range(ndim)},
        }
        non_finite = [name for name, value in required_fields.items() if not math.isfinite(value)]
        if non_finite:
            raise ValueError(f"non-finite sink state for id={sink_id}: {','.join(non_finite)}")
        if inflow < 0.0 or not 0.0 <= efficiency < 1.0:
            raise ValueError(f"invalid inflow or radiative efficiency for id={sink_id}")
        inflow_cgs = inflow * rate_unit_cgs
        mass_msun = masses[index] * units["mass_cgs"] / MSUN_CGS
        records.append(
            {
                "record_type": "agn_coarse_state",
                "source_format": "ramses_sinkprops_v1",
                "schema_version": 1,
                "ledger_phase": "pre_feedback_pre_reset",
                "source_interval_kind": "instantaneous_pre_reset_state",
                "julian_year_days": 365.25,
                "nstep_coarse": nstep_coarse,
                "aexp": aexp,
                "t_code": time_code,
                "sink_id": sink_id,
                "mass_code": masses[index],
                "mass_msun": mass_msun,
                "sink_mass_code": masses[index],
                "position_code": [positions[axis][index] for axis in range(ndim)],
                "velocity_code": [velocities[axis][index] for axis in range(ndim)],
                "gas_angular_momentum_code": [angular_momenta[axis][index] for axis in range(ndim)],
                "bondi_rate_code": bondi[index],
                "eddington_rate_code": eddington[index],
                "inflow_rate_code": inflow,
                "inflow_rate_msun_per_yr": inflow_cgs / MSUN_CGS * YEAR_S,
                # sinkprops exposes only eps_sink.  Preserve it in both
                # fields, explicitly marking that no mode-resolved effective
                # efficiency was available in this review-only conversion.
                "radiative_efficiency": efficiency,
                "effective_radiative_efficiency": efficiency,
                "effective_efficiency_status": "sinkprops_raw_equals_effective_mode_not_encoded",
                "bolometric_luminosity_erg_s": efficiency * inflow_cgs * C_LIGHT_CGS**2,
                "accreted_mass_since_feedback_code": accreted[index],
                "ambient_density_code": gas_density[index],
                "ambient_sound_speed_code": sound_speed[index],
                "relative_speed_code": relative_speed[index],
                "saved_feedback_energy_code": saved_energy[index],
                "bh_spin_code": [bh_spin[axis][index] for axis in range(ndim)],
                "bh_spin_magnitude": spin_magnitude[index],
                "units": units,
                "unit_mass_cgs": units["mass_cgs"],
                "unit_time_cgs": units["time_s"],
                "accretion_rate_convention": "instantaneous_inflow_rate_min_bondi_eddington_effective_agn_efficiency",
                "feedback_timing": "written_before_AGN_feedback_accumulator_reset",
            }
        )
    return units, records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="RAMSES sink_XXXXX.dat file")
    parser.add_argument("--output", type=Path, required=True, help="AGN source JSONL output")
    parser.add_argument(
        "--nstep-coarse",
        type=int,
        required=True,
        help="coarse-step identity from the producing RAMSES run (never inferred from the filename)",
    )
    args = parser.parse_args()

    units, records = read_sinkprops(args.input, nstep_coarse=args.nstep_coarse)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    header = {
        "record_type": "agn_coarse_state_header",
        "source_format": "ramses_sinkprops_v1",
        "input": str(args.input),
        "sink_count": len(records),
        "units": units,
        "luminosity_convention": "L_bol=epsilon*min(Mdot_Bondi,Mdot_Edd)*c^2",
        "ledger_phase": "pre_feedback_pre_reset",
        "source_interval_kind": "instantaneous_pre_reset_state",
        "julian_year_days": 365.25,
        "feedback_timing": "written_before_AGN_feedback_accumulator_reset",
    }
    with args.output.open("w", encoding="ascii") as stream:
        stream.write(json.dumps(header, sort_keys=True, allow_nan=False) + "\n")
        for record in sorted(records, key=lambda item: int(item["sink_id"])):
            stream.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
    print(f"SINKPROPS_LEDGER_OK sinks={len(records)} aexp={records[0]['aexp'] if records else 'none'} output={args.output}")


if __name__ == "__main__":
    main()
