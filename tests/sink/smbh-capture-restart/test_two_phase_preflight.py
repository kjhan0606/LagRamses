#!/usr/bin/env python3
"""Focused tests for the sealed prepare/launch boundary."""

from __future__ import annotations

import json
import hashlib
import importlib.util
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock
from pathlib import Path


HERE = Path(__file__).resolve().parent
TOOL = HERE / "prepare_and_validate.py"
SPEC = importlib.util.spec_from_file_location("capture_prepare", TOOL)
assert SPEC is not None and SPEC.loader is not None
PREPARE_MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE_MODULE)


class TwoPhasePreflightTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        base = Path(self.temp.name)
        self.source = base / "source"
        self.source.mkdir()
        subprocess.run(["git", "init", "-q", str(self.source)], check=True)
        subprocess.run(["git", "-C", str(self.source), "config", "user.name", "Test"], check=True)
        subprocess.run(["git", "-C", str(self.source), "config", "user.email", "test@example.invalid"], check=True)
        (self.source / "bin").mkdir()
        (self.source / ".gitignore").write_text("/bin/ramses3d\n", encoding="utf-8")
        (self.source / "tracked.txt").write_text("clean\n", encoding="utf-8")
        self.binary = self.source / "bin" / "ramses3d"
        shutil.copy2("/bin/true", self.binary)
        subprocess.run(["git", "-C", str(self.source), "add", ".gitignore", "tracked.txt"], check=True)
        subprocess.run(["git", "-C", str(self.source), "commit", "-qm", "fixture"], check=True)
        self.commit = subprocess.check_output(
            ["git", "-C", str(self.source), "rev-parse", "HEAD"], text=True
        ).strip()
        self.run_dir = base / "prepared"

    def command(self, *extra: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["python3", str(TOOL), *extra], text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )

    def manifest_hash(self) -> str:
        return hashlib.sha256((self.run_dir / "preflight_manifest.json").read_bytes()).hexdigest()

    def prepare(self) -> subprocess.CompletedProcess[str]:
        return self.command(
            "prepare", str(self.run_dir), "--binary", str(self.binary),
            "--source-tree", str(self.source),
            "--expected-source-commit", self.commit,
            "--bytes-per-output", "1", "--reserve-bytes", "0",
        )

    def test_prepare_then_launch_check(self) -> None:
        result = self.prepare()
        self.assertEqual(result.returncode, 0, result.stdout)
        manifest = json.loads((self.run_dir / "preflight_manifest.json").read_text())
        self.assertEqual(manifest["expected_outputs"], ["output_00001", "output_00002"])
        result = self.command("launch-check", str(self.run_dir), "--binary", str(self.binary),
                              "--manifest-sha256", self.manifest_hash())
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("launch_revalidation=PASS", result.stdout)

    def test_refuses_overwrite_and_detects_tamper(self) -> None:
        self.assertEqual(self.prepare().returncode, 0)
        result = self.prepare()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("refusing existing run directory", result.stdout)
        with (self.run_dir / "stage1.nml").open("a", encoding="utf-8") as stream:
            stream.write("! tampered\n")
        result = self.command("launch-check", str(self.run_dir), "--binary", str(self.binary),
                              "--manifest-sha256", self.manifest_hash())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("hash mismatch", result.stdout)

    def test_malformed_manifest_fails_cleanly(self) -> None:
        self.assertEqual(self.prepare().returncode, 0)
        (self.run_dir / "preflight_manifest.json").write_text("[]\n")
        result = self.command("launch-check", str(self.run_dir), "--binary", str(self.binary),
                              "--manifest-sha256", self.manifest_hash())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("manifest root must be an object", result.stdout)
        self.assertNotIn("Traceback", result.stdout)

    def test_rejects_unapproved_manifest_hash(self) -> None:
        self.assertEqual(self.prepare().returncode, 0)
        result = self.command("launch-check", str(self.run_dir), "--binary", str(self.binary),
                              "--manifest-sha256", "0" * 64)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("operator-approved identity", result.stdout)

    def test_rejects_binary_swap(self) -> None:
        self.assertEqual(self.prepare().returncode, 0)
        shutil.copy2("/bin/false", self.binary)
        result = self.command("launch-check", str(self.run_dir), "--binary", str(self.binary),
                              "--manifest-sha256", self.manifest_hash())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("binary SHA-256 differs", result.stdout)

    def test_rejects_dirty_source(self) -> None:
        self.assertEqual(self.prepare().returncode, 0)
        (self.source / "tracked.txt").write_text("dirty\n", encoding="utf-8")
        result = self.command("launch-check", str(self.run_dir), "--binary", str(self.binary),
                              "--manifest-sha256", self.manifest_hash())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("source worktree is not clean", result.stdout)

    def test_rejects_nonpristine_output_state(self) -> None:
        self.assertEqual(self.prepare().returncode, 0)
        (self.run_dir / "output_00001").mkdir()
        result = self.command("launch-check", str(self.run_dir), "--binary", str(self.binary),
                              "--manifest-sha256", self.manifest_hash())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no longer pristine", result.stdout)

    def test_quota_evidence_uses_stricter_nonzero_ceiling(self) -> None:
        output = """Disk quotas for usr test (uid 1):
Filesystem kbytes quota limit grace files quota limit grace
/scratch 100 500 800 - 1 0 0 -
"""
        with mock.patch.dict(PREPARE_MODULE.os.environ, {"USER": "test"}), \
             mock.patch.object(PREPARE_MODULE.subprocess, "check_output", return_value=output):
            evidence = PREPARE_MODULE.lustre_quota_evidence(Path("/scratch/test/run"))
        self.assertIsNotNone(evidence)
        self.assertEqual(evidence["remaining_bytes_under_stricter_limit"], 400 * 1024)

    def test_zero_quota_ceiling_is_recorded_not_invented(self) -> None:
        output = """Disk quotas for usr test (uid 1):
Filesystem kbytes quota limit grace files quota limit grace
/scratch 100 0 0 - 1 0 0 -
"""
        with mock.patch.dict(PREPARE_MODULE.os.environ, {"USER": "test"}), \
             mock.patch.object(PREPARE_MODULE.subprocess, "check_output", return_value=output):
            evidence = PREPARE_MODULE.lustre_quota_evidence(Path("/scratch/test/run"))
        self.assertIsNotNone(evidence)
        self.assertIsNone(evidence["remaining_bytes_under_stricter_limit"])


if __name__ == "__main__":
    unittest.main()
