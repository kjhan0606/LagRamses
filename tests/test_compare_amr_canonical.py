#!/usr/bin/env python3
"""Synthetic regression tests for compare_amr_canonical.py."""

from __future__ import annotations

import pathlib
import struct
import subprocess
import sys
import tempfile
import unittest


HERE = pathlib.Path(__file__).resolve().parent
COMPARATOR = HERE / "compare_amr_canonical.py"


class LegacyWriter:
    def __init__(self, path: pathlib.Path):
        self.stream = path.open("wb")

    def close(self) -> None:
        self.stream.close()

    def record(self, payload: bytes) -> None:
        marker = struct.pack("<i", len(payload))
        self.stream.write(marker)
        self.stream.write(payload)
        self.stream.write(marker)

    def ints(self, *values: int) -> None:
        self.record(struct.pack("<" + "i" * len(values), *values))

    def doubles(self, *values: float) -> None:
        self.record(struct.pack("<" + "d" * len(values), *values))

    def int64s(self, *values: int) -> None:
        self.record(struct.pack("<" + "q" * len(values), *values))


def write_info(path: pathlib.Path, capacity: int) -> None:
    path.write_text(
        "\n".join(
            (
                "ncpu        =          1",
                "ndim        =          1",
                "levelmin    =          1",
                "levelmax    =          2",
                f"ngridmax    = {capacity:10d}",
                "nstep_coarse=          4",
                "",
                "boxlen      =  0.100000000000000E+01",
                "time        =  0.250000000000000E+00",
                "aexp        =  0.500000000000000E+00",
                "H0          =  0.700000000000000E+02",
                "omega_m     =  0.300000000000000E+00",
                "omega_l     =  0.700000000000000E+00",
                "omega_k     =  0.000000000000000E+00",
                "omega_b     =  0.500000000000000E-01",
                "unit_l      =  0.100000000000000E+01",
                "unit_d      =  0.200000000000000E+01",
                "unit_t      =  0.300000000000000E+01",
                "",
                "ordering type=hilbert",
                "   DOMAIN   ind_min                 ind_max",
                "       1   0.0   1.0",
                "",
            )
        ),
        encoding="ascii",
    )


def write_amr(
    path: pathlib.Path,
    *,
    capacity: int,
    level1_id: int,
    level2_id: int,
    flag_delta: int = 0,
    numbtot_rank_noise: int = 0,
    level_total_delta: int = 0,
    nbor_flip: bool = False,
) -> None:
    writer = LegacyWriter(path)
    try:
        writer.ints(1)  # ncpu
        writer.ints(1)  # ndim
        writer.ints(1, 1, 1)
        writer.ints(2)  # nlevelmax
        writer.ints(capacity)
        writer.ints(0)  # nboundary
        writer.ints(2)  # ngrid_current (excluded)
        writer.doubles(1.0)

        writer.ints(1, 1, 2)
        writer.doubles(0.0)
        writer.doubles(0.5)
        writer.doubles(0.25)
        writer.doubles(0.1, 0.05)
        writer.doubles(0.09, 0.04)
        writer.ints(4, 4)
        writer.doubles(1.0, 2.0, 3.0)
        writer.doubles(0.3, 0.7, 0.0, 0.05, 70.0, 0.01, 1.0)
        writer.doubles(0.5, 1.0, 0.49, -2.0, -2.1)
        writer.doubles(4.0)

        writer.ints(level1_id, level2_id)  # headl
        writer.ints(level1_id, level2_id)  # taill
        writer.ints(1, 1)  # numbl
        numbtot = [0] * 20
        numbtot[0] = 1 + level_total_delta
        numbtot[1] = numbtot_rank_noise
        numbtot[2] = numbtot_rank_noise
        numbtot[3] = numbtot_rank_noise
        numbtot[10] = 1
        writer.int64s(*numbtot)  # first index is the global total per level
        writer.ints(0, capacity, capacity - 2, 2, 2)  # free state
        writer.record(b"hilbert".ljust(128, b" "))
        writer.doubles(0.0, 1.0)  # ignored bound_key

        writer.ints(level1_id)  # coarse son
        writer.ints(1)  # coarse flag1
        writer.ints(1)  # coarse cpu_map

        # Level 1 owner block.
        writer.ints(level1_id)
        writer.ints(0)
        writer.ints(0)
        writer.doubles(0.5)
        writer.ints(1)  # father: coarse cell 1
        writer.ints(1)  # left coarse neighbour
        writer.ints(1)  # right coarse neighbour (periodic)
        writer.ints(level2_id)  # lower child is refined
        writer.ints(0)
        writer.ints(1)
        writer.ints(1)
        writer.ints(1 + flag_delta)
        writer.ints(0)

        # Level 2 owner block.  Father is child 0 of the level-1 grid.
        lower_cell = 1 + level1_id
        upper_cell = 1 + capacity + level1_id
        writer.ints(level2_id)
        writer.ints(0)
        writer.ints(0)
        writer.doubles(0.25)
        writer.ints(lower_cell)
        writer.ints(lower_cell if nbor_flip else upper_cell)
        writer.ints(lower_cell)
        writer.ints(0)
        writer.ints(0)
        writer.ints(1)
        writer.ints(1)
        writer.ints(0)
        writer.ints(0)
    finally:
        writer.close()


def write_output(
    directory: pathlib.Path,
    *,
    capacity: int,
    level1_id: int,
    level2_id: int,
    flag_delta: int = 0,
    numbtot_rank_noise: int = 0,
    level_total_delta: int = 0,
    nbor_flip: bool = False,
) -> None:
    directory.mkdir()
    write_amr(
        directory / "amr_00001.out00001",
        capacity=capacity,
        level1_id=level1_id,
        level2_id=level2_id,
        flag_delta=flag_delta,
        numbtot_rank_noise=numbtot_rank_noise,
        level_total_delta=level_total_delta,
        nbor_flip=nbor_flip,
    )
    write_info(directory / "info_00001.txt", capacity)


class ComparatorTest(unittest.TestCase):
    def run_comparator(
        self, left: pathlib.Path, right: pathlib.Path
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            (sys.executable, str(COMPARATOR), str(left), str(right)),
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def test_capacity_and_local_ids_are_canonicalized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            left = root / "left"
            right = root / "right"
            write_output(left, capacity=8, level1_id=2, level2_id=3)
            write_output(right, capacity=16, level1_id=5, level2_id=7)
            result = self.run_comparator(left, right)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("AMR_CANONICAL PASS", result.stdout)

    def test_semantic_flag_corruption_is_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            left = root / "left"
            right = root / "right"
            write_output(left, capacity=8, level1_id=2, level2_id=3)
            write_output(
                right,
                capacity=16,
                level1_id=5,
                level2_id=7,
                flag_delta=1,
            )
            result = self.run_comparator(left, right)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn("fields=['flag1']", result.stderr)

    def test_rank_dependent_numbtot_columns_are_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            left = root / "left"
            right = root / "right"
            write_output(left, capacity=8, level1_id=2, level2_id=3)
            write_output(
                right,
                capacity=16,
                level1_id=5,
                level2_id=7,
                numbtot_rank_noise=99,
            )
            result = self.run_comparator(left, right)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_global_numbtot_must_match_owner_count(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            left = root / "left"
            right = root / "right"
            write_output(
                left,
                capacity=8,
                level1_id=2,
                level2_id=3,
                level_total_delta=1,
            )
            write_output(right, capacity=16, level1_id=5, level2_id=7)
            result = self.run_comparator(left, right)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn("disagrees with numbtot global total", result.stderr)

    def test_valid_neighbour_change_is_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            left = root / "left"
            right = root / "right"
            write_output(left, capacity=8, level1_id=2, level2_id=3)
            write_output(
                right,
                capacity=16,
                level1_id=5,
                level2_id=7,
                nbor_flip=True,
            )
            result = self.run_comparator(left, right)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn("fields=['neighbours']", result.stderr)

    def test_info_capacity_must_match_its_amr_header(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            left = root / "left"
            right = root / "right"
            write_output(left, capacity=8, level1_id=2, level2_id=3)
            write_output(right, capacity=16, level1_id=5, level2_id=7)
            info = right / "info_00001.txt"
            info.write_text(
                info.read_text(encoding="ascii").replace(
                    "ngridmax    =         16", "ngridmax    =         15"
                ),
                encoding="ascii",
            )
            result = self.run_comparator(left, right)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn("info admin fields", result.stderr)

    def test_mixed_output_stem_is_malformed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            left = root / "left"
            right = root / "right"
            write_output(left, capacity=8, level1_id=2, level2_id=3)
            write_output(right, capacity=16, level1_id=5, level2_id=7)
            (right / "amr_00001.out00001").rename(
                right / "amr_00002.out00001"
            )
            result = self.run_comparator(left, right)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn("does not match info_00001.txt", result.stderr)

    def test_truncated_record_is_malformed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            left = root / "left"
            right = root / "right"
            write_output(left, capacity=8, level1_id=2, level2_id=3)
            write_output(right, capacity=16, level1_id=5, level2_id=7)
            amr = right / "amr_00001.out00001"
            payload = amr.read_bytes()
            amr.write_bytes(payload[:-1])
            result = self.run_comparator(left, right)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn("AMR_CANONICAL ERROR", result.stderr)


if __name__ == "__main__":
    unittest.main()
