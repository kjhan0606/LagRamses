"""Export audited instantaneous AGN rates from a RAMSES sink diagnostic."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.sink_diagnostic import (
    read_agn_coarse_records,
    read_agn_coarse_state,
    read_sink_diagnostic,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--diagnostic")
    source.add_argument("--agn-coarse-json")
    parser.add_argument("--aexp", type=float)
    parser.add_argument("--aexp-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.agn_coarse_json:
        if args.aexp is None:
            parser.error("--aexp is required with --agn-coarse-json")
        coarse_records = read_agn_coarse_records(args.agn_coarse_json)
        selected_records = [
            record
            for record in coarse_records
            if abs(float(record["aexp"]) - args.aexp) <= args.aexp_tolerance
        ]
        non_promotable = [
            int(record["sink_id"])
            for record in selected_records
            if not bool(record["efficiency_contract_ok"])
        ]
        if non_promotable:
            raise ValueError(
                "AGN coarse-state contains non-promotable efficiency contract "
                f"for sink IDs {non_promotable}; refusing P4 artifact"
            )
        diagnostic = read_agn_coarse_state(
            args.agn_coarse_json,
            expansion_factor=args.aexp,
            expansion_factor_tolerance=args.aexp_tolerance,
        )
        rate_convention = "instantaneous_inflow_rate_min_bondi_eddington_effective_agn_efficiency"
        mass_msun = diagnostic.mass_msun
        raw_efficiency = diagnostic.raw_radiative_efficiency
        effective_efficiency = diagnostic.effective_radiative_efficiency
        efficiency_status = diagnostic.efficiency_status
        efficiency_contract_ok = diagnostic.efficiency_contract_ok
        non_promotable = [
            int(sink_id)
            for sink_id, contract_ok in zip(
                diagnostic.sink_id, efficiency_contract_ok, strict=True
            )
            if not bool(contract_ok)
        ]
        if non_promotable:
            raise ValueError(
                "AGN coarse-state contains non-promotable efficiency contract "
                f"for sink IDs {non_promotable}; refusing P4 artifact"
            )
        efficiency_contract_source = "snrt_agn_efficiency_helper"
    else:
        diagnostic = read_sink_diagnostic(args.diagnostic)
        rate_convention = "instantaneous_inflow_rate_min_bondi_eddington"
        mass_msun = diagnostic.mass_code * diagnostic.mass_scale_msun
        raw_efficiency = diagnostic.radiative_efficiency
        # Native sinkprops has no mode-resolved efficiency field, so the
        # review-only diagnostic explicitly records raw == effective rather
        # than making the distinction disappear.
        effective_efficiency = diagnostic.radiative_efficiency
        # The legacy sinkprops record has no mode-resolved status/contract.
        # Keep the CSV schema explicit without claiming helper parity.
        efficiency_status = None
        efficiency_contract_ok = None
        efficiency_contract_source = "legacy_sinkprops_mode_unresolved"
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            (
                "sink_id",
                "aexp",
                "time_code",
                "x_code",
                "y_code",
                "z_code",
                "mass_msun",
                "bondi_mdot_msun_per_year",
                "eddington_mdot_msun_per_year",
                "inflow_mdot_msun_per_year",
                "raw_radiative_efficiency",
                "effective_radiative_efficiency",
                "efficiency_status",
                "efficiency_contract_ok",
                "efficiency_contract_source",
                "bolometric_luminosity_erg_s",
                "accretion_rate_convention",
            )
        )
        for index, sink_id in enumerate(diagnostic.sink_id):
            writer.writerow(
                (
                    int(sink_id),
                    diagnostic.expansion_factor,
                    diagnostic.time_code,
                    *diagnostic.position_code[index],
                    mass_msun[index],
                    diagnostic.bondi_rate_msun_per_year[index],
                    diagnostic.eddington_rate_msun_per_year[index],
                    diagnostic.inflow_rate_msun_per_year[index],
                    raw_efficiency[index],
                    effective_efficiency[index],
                    "" if efficiency_status is None else int(efficiency_status[index]),
                    "" if efficiency_contract_ok is None else str(bool(efficiency_contract_ok[index])).lower(),
                    efficiency_contract_source,
                    diagnostic.bolometric_luminosity_erg_s[index],
                    rate_convention,
                )
            )
    print(f"P4_AGN_RATE_LEDGER_OK sinks={len(diagnostic.sink_id)} output={output}")


if __name__ == "__main__":
    main()
