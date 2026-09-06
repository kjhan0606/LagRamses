"""Synthetic contract test for the native RAMSES output audit."""

from __future__ import annotations

import sys
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.audit_native_ramses_output import audit_output


def main() -> None:
    with TemporaryDirectory() as temporary_directory:
        output_dir = Path(temporary_directory) / "output_00001"
        output_dir.mkdir()
        (output_dir / "COMPLETE").touch()
        (output_dir / "info_00001.txt").write_text(
            "ncpu        =          2\n"
            "ndim        =          3\n"
            "levelmax    =          14\n"
            "nstep_coarse=         10\n"
            "aexp        =  0.2\n",
            encoding="ascii",
        )
        (output_dir / "header_00001.txt").write_text("header\n", encoding="ascii")
        (output_dir / "namelist.txt").write_text("&RUN_PARAMS\n/\n", encoding="ascii")
        (output_dir / "compilation.txt").write_text(
            "compile date = test\nlast commit = abc-dirty\n", encoding="ascii"
        )
        (output_dir / "makefile.txt").write_text(
            "NVAR = 17\nSOLVER = hydro\nPHASE0_STELLAR_ENRICHMENT = 1\nF90 = mpiifx\n",
            encoding="ascii",
        )
        (output_dir / "hydro_file_descriptor.txt").write_text(
            "nvar = 17\n" + "".join(f"variable #{index}: field_{index}\n" for index in range(1, 18)),
            encoding="ascii",
        )
        (output_dir / "resolved_physics_inventory_00001.txt").write_text(
            "sink_info_file = none\nforce_source_ledger_status = unavailable\n",
            encoding="ascii",
        )
        for prefix in ("amr", "hydro", "grav", "part"):
            for rank in (1, 2):
                (output_dir / f"{prefix}_00001.out{rank:05d}").touch()
        run_log = Path(temporary_directory) / "run.log"
        run_log.write_text(
            "Phase 0 stellar enrichment enabled\n"
            "  table rows = 9\n"
            "  total-metal field = 6\n"
            "  first element field = 7\n"
            "  snrt_advance : 0.000 s\n"
            "  snrt_diagnose: 0.000 s\n",
            encoding="ascii",
        )

        record = audit_output(output_dir, run_log)

    assert record["status"] == "complete_native_metadata_audited"
    assert record["expected_ranks"] == 2
    assert record["components"]["hydro"]["rank_count_matches_ncpu"] is True
    assert record["compilation"]["last_commit"] == "abc-dirty"
    assert record["makefile_flags"]["PHASE0_STELLAR_ENRICHMENT"] == "1"
    assert record["source_ledger_available"] is False
    assert record["runtime_log_evidence"]["phase0_enabled_marker"] is True
    assert record["runtime_log_evidence"]["phase0_table_rows"] == 9
    assert record["runtime_log_evidence"]["first_element_field"] == 7
    assert record["runtime_log_evidence"]["snrt_advance_seconds"] == 0.0
    assert record["scientific_readiness"]["direct_snrt_input"] is False
    print("P0_NATIVE_OUTPUT_AUDIT_OK ranks=2 components=4 source_ledger=unavailable")


if __name__ == "__main__":
    main()
