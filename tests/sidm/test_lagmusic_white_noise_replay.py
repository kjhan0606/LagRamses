#!/usr/bin/env python3
"""Unit checks for the external-white-noise manual replay driver."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from lagmusic_white_noise_replay import CASES, compact_noise_metadata, make_config, run


class WhiteNoiseReplayTests(unittest.TestCase):
    def test_compact_double_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wnoise_0004.bin"
            path.write_bytes((2).to_bytes(4, "little") * 3 + b"\0" * (8 * 8))
            data = compact_noise_metadata(path)
        self.assertEqual(data["dimensions"], [2, 2, 2])
        self.assertEqual(data["precision"], "float64")

    def test_external_anchor_and_explicit_rank_path(self):
        template = "[setup]\nlevelmin = 4\n[random]\ndisk_cached = yes\nseed[4] = 7\n[output]\nfilename = old\n"
        config = make_config(template, Path("/tmp/replay"), Path("/tmp/wnoise_0004.bin"))
        self.assertIn("seed[4] = /tmp/wnoise_0004.bin", config)
        self.assertIn("filename = /tmp/replay/ic", config)
        self.assertIn("slab_rng_refinement = no", config)
        self.assertIn("rng_boost_threads = no", config)

    def test_run_pins_openmp_thread_count_from_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "manifest.json").write_text(
                '{"execution": {"omp_num_threads": 1, "omp_dynamic": "FALSE"}}')
            for label in CASES:
                case = root / label
                case.mkdir()
                (case / "music.conf").write_text("[setup]\n")
            with patch("lagmusic_white_noise_replay.subprocess.run") as mocked:
                run(root, Path("/test/MUSIC"), 30)
            self.assertEqual(mocked.call_count, 2)
            for call in mocked.call_args_list:
                self.assertEqual(call.kwargs["env"]["OMP_NUM_THREADS"], "1")
                self.assertEqual(call.kwargs["env"]["OMP_DYNAMIC"], "FALSE")


if __name__ == "__main__":
    unittest.main()
