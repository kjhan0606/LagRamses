#!/usr/bin/env python3
"""Regression check for the sink-capture transaction ordering."""

from pathlib import Path
import re
import unittest


SOURCE = Path(__file__).parents[1] / "sink_particle.kjhan.f90"
OUTPUT_SOURCE = Path(__file__).parents[1] / "output_amr.kjhan.f90"
INIT_PART_SOURCE = Path(__file__).parents[1] / "init_part.f90"


def _fortran_call_arguments(source: str, name: str) -> list[str]:
    """Return balanced-parenthesis arguments for calls to one Fortran routine."""

    arguments = []
    pattern = re.compile(rf"\b{re.escape(name)}\s*\(", re.IGNORECASE)
    for match in pattern.finditer(source):
        start = match.end()
        depth = 1
        for index in range(start, len(source)):
            if source[index] == "(":
                depth += 1
            elif source[index] == ")":
                depth -= 1
                if depth == 0:
                    arguments.append(" ".join(source[start:index].split()))
                    break
        else:
            raise AssertionError(f"unterminated Fortran call to {name}")
    return arguments


class SinkMergeLedgerOrderTest(unittest.TestCase):
    def test_attempt_and_checkpoint_hooks_are_at_durable_boundaries(self) -> None:
        init_part = INIT_PART_SOURCE.read_text(encoding="utf-8")
        self.assertLess(
            init_part.index("call init_sink"),
            init_part.index("call initialize_smbh_capture_lineage"),
        )

        output = OUTPUT_SOURCE.read_text(encoding="utf-8")
        lineage = output.index("call write_smbh_capture_checkpoint_lineage")
        complete = output.index("filename='output_'//TRIM(nchar)//'/COMPLETE'", lineage)
        self.assertLess(lineage, complete)
        stale = output.index("'/SMBH_CAPTURE_LINEAGE'")
        self.assertLess(stale, lineage)

        sink_source = SOURCE.read_text(encoding="utf-8")
        initializer = sink_source.split(
            "subroutine initialize_smbh_capture_lineage", 1
        )[1].split("end subroutine initialize_smbh_capture_lineage", 1)[0]
        self.assertNotIn("allocate(commit_seen", initializer)
        self.assertIn('"prepared":true,"committed":false}', initializer)
        self.assertIn('"committed":true}', initializer)
        self.assertIn("parsed_step<=parent_step", initializer)

    def test_ledger_follows_validation_and_precedes_compaction(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        merge = source.split("subroutine merge_sink(ilevel)", 1)[1].split(
            "end subroutine merge_sink", 1
        )[0]
        executable = "\n".join(
            line.split("!", 1)[0] for line in merge.splitlines()
        )

        ledger_calls = list(
            re.finditer(r"(?mi)^\s*call\s+write_smbh_capture_ledger\b", executable)
        )
        self.assertEqual(len(ledger_calls), 1)
        ledger = ledger_calls[0].start()

        validation = executable.index("if(pending_error/=0)then")
        validation_stop = executable.index("call clean_stop", validation)
        pending_commit = executable.index(
            "agn_pending_erg(1:new_sink)=pending_new", ledger
        )
        sink_count_commit = executable.index("nsink=new_sink", pending_commit)

        barriers = list(re.finditer(r"(?mi)^\s*call\s+MPI_BARRIER\b", executable))
        self.assertEqual(len(barriers), 3)
        self.assertLess(ledger, barriers[0].start())
        self.assertLess(barriers[0].start(), pending_commit)
        physical_compaction = executable.rindex("call kjhan_kill_sink")
        batch_commit = executable.index("call commit_smbh_capture_ledger")
        self.assertLess(physical_compaction, barriers[1].start())
        self.assertLess(barriers[1].start(), batch_commit)
        self.assertLess(batch_commit, barriers[2].start())

        original_state = (
            "agn_pending_erg",
            "msink",
            "dMsmbh",
            "Esave",
            "idsink",
            "tsink",
            "xsink",
            "vsink",
            "bhspin",
            "spinmag",
        )
        for name in original_state:
            assignments = list(
                re.finditer(rf"(?mi)^\s*{name}\s*\([^\n=]*\)\s*=", executable)
            )
            self.assertTrue(assignments, f"no concrete assignment found for {name}")
            self.assertTrue(
                all(match.start() > ledger for match in assignments),
                f"{name} is assigned before the ledger transaction",
            )

        self.assertLess(validation_stop, ledger)
        self.assertLess(ledger, pending_commit)
        self.assertLess(pending_commit, sink_count_commit)

    def test_all_group_failures_are_preflighted_before_open(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        match = re.search(
            r"(?mis)^subroutine write_smbh_capture_ledger\([^\n]*&\n.*?"
            r"^end subroutine write_smbh_capture_ledger\s*$",
            source,
        )
        self.assertIsNotNone(match)
        ledger = match.group(0)
        executable = "\n".join(
            line.split("!", 1)[0] for line in ledger.splitlines()
        )
        opened = executable.index("open(newunit=ledger_unit")
        fatal_pattern = re.compile(
            r"call\s+ledger_group_fatal\([^)]*?'([^']+)'\)", re.I | re.S
        )
        before_open = set(fatal_pattern.findall(executable[:opened]))
        after_open = set(fatal_pattern.findall(executable[opened:]))

        self.assertTrue(after_open)
        self.assertTrue(
            after_open.issubset(before_open),
            f"group failures lack a pre-open preflight: {after_open - before_open}",
        )

    def test_fortran_writer_preserves_all_members_and_pairs_for_multiple(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        match = re.search(
            r"(?mis)^subroutine write_smbh_capture_ledger\([^\n]*&\n.*?"
            r"^end subroutine write_smbh_capture_ledger\s*$",
            source,
        )
        self.assertIsNotNone(match)
        writer = match.group(0)
        executable = "\n".join(
            line.split("!", 1)[0] for line in writer.splitlines()
        )

        self.assertIn("count(gsink(1:nsink) == igrp) >= 2", executable)
        self.assertIn("expected_pairs=nmember*(nmember-1)/2", executable)
        self.assertRegex(
            executable,
            r"(?is)if\s*\(nmember\s*==\s*2\)\s*then\s*"
            r"classification\s*=\s*'BINARY'\s*else\s*"
            r"classification\s*=\s*'MULTIPLE'\s*endif",
        )

        member_loop = re.search(
            r"(?is)member_index=0\s*do isink=1,nsink(.*?)\benddo",
            executable,
        )
        self.assertIsNotNone(member_loop)
        self.assertIn("if(gsink(isink) /= igrp) cycle", member_loop.group(0))
        self.assertIn("member_index=member_index+1", member_loop.group(0))
        self.assertIn("json_int(idsink(isink))", member_loop.group(0))

        event_begin_record = executable.split(
            '"record_type":"event_begin"', 1
        )[1].split('"record_type":"member"', 1)[0]
        self.assertIn("trim(json_int(expected_pairs))", event_begin_record)

        pair_records = executable.split("pair_index=0", 1)[1].split(
            '"record_type":"event_end"', 1
        )[0]
        self.assertRegex(pair_records, r"(?is)do isink=1,nsink-1")
        self.assertRegex(pair_records, r"(?is)do jsink_member=isink\+1,nsink")
        self.assertIn("if(gsink(isink) /= igrp) cycle", pair_records)
        self.assertIn("if(gsink(jsink_member) /= igrp) cycle", pair_records)
        self.assertIn("pair_index=pair_index+1", pair_records)
        self.assertIn("json_int(idsink(isink))", pair_records)
        self.assertIn("json_int(idsink(jsink_member))", pair_records)
        self.assertIn("trim(json_int(pair_index))", executable.split(
            '"record_type":"event_end"', 1
        )[1])

    def test_every_serialized_numeric_family_is_preflighted_before_open(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        match = re.search(
            r"(?mis)^subroutine write_smbh_capture_ledger\([^\n]*&\n.*?"
            r"^end subroutine write_smbh_capture_ledger\s*$",
            source,
        )
        self.assertIsNotNone(match)
        ledger = match.group(0)
        executable = "\n".join(
            line.split("!", 1)[0] for line in ledger.splitlines()
        )
        opened = executable.index("open(newunit=ledger_unit")
        preflight = executable[:opened]

        required_checks = (
            "non-finite output units or gravity normalization",
            "non-finite global capture diagnostics",
            "non-finite capture redshift",
            "non-finite or non-positive periodic box extent",
            "non-finite or non-positive member mass",
            "non-finite member kinematics",
            "non-finite member sink state",
            "non-finite last Bondi context",
            "non-finite or non-positive total mass",
            "non-finite group center-of-mass diagnostics",
            "non-finite pair diagnostics",
            "non-finite pair binding diagnostics",
        )
        for diagnostic in required_checks:
            self.assertIn(diagnostic, preflight, f"missing pre-open check: {diagnostic}")

        for expression in (
            "ieee_is_finite(scale_l)",
            "ieee_is_finite(scale_t)",
            "ieee_is_finite(scale_d)",
            "ieee_is_finite(scale_v)",
            "ieee_is_finite(scale_m)",
            "ieee_is_finite(factG)",
            "ieee_is_finite(dx_min)",
            "ieee_is_finite(scale)",
            "ieee_is_finite(rmerge)",
            "ieee_is_finite(rmerge*dx_min)",
            "ieee_is_finite((rmerge*dx_min)**2)",
            "ieee_is_finite(boxlen)",
            "ieee_is_finite(aexp)",
            "ieee_is_finite(t)",
            "ieee_is_finite(texp)",
            "ieee_is_finite(omega_m)",
            "ieee_is_finite(h0)",
            "ieee_is_finite(redshift)",
            "ieee_is_finite(periodic_box_size(idim))",
            "ieee_is_finite(msink(isink))",
            "ieee_is_finite(xsink(isink,1:3))",
            "ieee_is_finite(vsink(isink,1:3))",
            "ieee_is_finite(tsink(isink))",
            "ieee_is_finite(dMsmbh(isink))",
            "ieee_is_finite(dMBH_coarse(isink))",
            "ieee_is_finite(dMEd_coarse(isink))",
            "ieee_is_finite(Esave(isink))",
            "ieee_is_finite(spinmag(isink))",
            "ieee_is_finite(bhspin(isink,1:3))",
            "ieee_is_finite(jsink(isink,1:3))",
            "ieee_is_finite(d_avgptr(isink))",
            "ieee_is_finite(c_avgptr(isink))",
            "ieee_is_finite(v_avgptr(isink))",
            "ieee_is_finite(total_mass)",
            "ieee_is_finite(com_pos(1:3))",
            "ieee_is_finite(com_vel(1:3))",
            "ieee_is_finite(separation2)",
            "ieee_is_finite(separation)",
            "ieee_is_finite(vrel2)",
            "ieee_is_finite(vrel)",
            "ieee_is_finite(reduced_mass)",
            "ieee_is_finite(relative_kinetic)",
            "ieee_is_finite(specific_h)",
            "ieee_is_finite(relative_L)",
            "ieee_is_finite(potential_1overr)",
            "ieee_is_finite(specific_energy)",
            "ieee_is_finite(legacy_binding_proxy)",
        ):
            self.assertIn(expression, preflight, f"not preflighted before append: {expression}")

        writer = executable[opened:].split("\ncontains", 1)[0]
        serialized_reals = set(_fortran_call_arguments(writer, "json_real"))
        serialized_optional_reals = {
            arguments.split(",", 1)[0].strip()
            for arguments in _fortran_call_arguments(writer, "json_optional_real")
        }
        expected_reals = {
            "aexp", "redshift", "t", "texp", "dx_min", "rmerge",
            "rmerge*dx_min", "boxlen", "omega_m", "h0", "factG",
            "scale_l", "scale_t", "scale_d", "scale_v", "scale_m",
            "total_mass", "max_separation", "msink(isink)",
            "tsink(isink)", "dMsmbh(isink)", "dMBH_coarse(isink)",
            "dMEd_coarse(isink)", "Esave(isink)", "spinmag(isink)",
            "separation", "vrel", "reduced_mass", "relative_kinetic",
        }
        for axis in range(1, 4):
            expected_reals.update({
                f"periodic_box_size({axis})",
                f"com_pos({axis})", f"com_vel({axis})",
                f"xsink(isink,{axis})", f"vsink(isink,{axis})",
                f"bhspin(isink,{axis})", f"jsink(isink,{axis})",
                f"delta_pos({axis})", f"delta_vel({axis})",
                f"specific_h({axis})", f"relative_L({axis})",
            })
        self.assertEqual(serialized_reals, expected_reals)
        self.assertEqual(serialized_optional_reals, {
            "last_density", "last_sound_speed", "last_gas_vrel",
            "potential_1overr", "specific_energy", "legacy_binding_proxy",
        })

        finite_arguments = set(_fortran_call_arguments(preflight, "ieee_is_finite"))
        component_checks = {}
        for axis in range(1, 4):
            for name in ("periodic_box_size",):
                component_checks[f"{name}({axis})"] = f"{name}(idim)"
            for name in ("com_pos", "com_vel"):
                component_checks[f"{name}({axis})"] = f"{name}(1:3)"
            for name in ("xsink", "vsink", "bhspin", "jsink"):
                component_checks[f"{name}(isink,{axis})"] = f"{name}(isink,1:3)"
            for name in ("specific_h", "relative_L"):
                component_checks[f"{name}({axis})"] = name
            component_checks[f"delta_pos({axis})"] = "separation2"
            component_checks[f"delta_vel({axis})"] = "vrel2"
        component_checks.update({
            "last_density": "d_avgptr(isink)",
            "last_sound_speed": "c_avgptr(isink)",
            "last_gas_vrel": "v_avgptr(isink)",
            "max_separation": "separation",
        })
        for expression in serialized_reals | serialized_optional_reals:
            check = component_checks.get(expression, expression)
            self.assertIn(
                check, finite_arguments,
                f"serialized numeric expression lacks preflight evidence: {expression} -> {check}",
            )


if __name__ == "__main__":
    unittest.main()
