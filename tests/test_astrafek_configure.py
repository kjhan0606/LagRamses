"""The additive configure entry point preserves all existing solver sources."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

SOURCE = Path(__file__).resolve().parents[1]
CONFIGURE = SOURCE / "configure"


def invoke(folder, *extra):
    return subprocess.run([sys.executable, str(CONFIGURE), "--build-dir", str(folder), *extra], capture_output=True, text=True)


def test_generates_review_makefile_without_modifying_upstream(tmp_path):
    upstream = SOURCE / "bin/Makefile"
    before = hashlib.sha256(upstream.read_bytes()).hexdigest()
    result = invoke(tmp_path / "build", "--compiler", "intel", "--no-probe")
    assert result.returncode == 0, result.stderr
    generated = tmp_path / "build/bin/Makefile"
    text = generated.read_text()
    assert "override F90 := mpiifx -qopenmp" in text
    assert "override EXEC := ramses_configured" in text
    assert ".NOTPARALLEL:" in text
    assert (tmp_path / "build/amr").resolve() == SOURCE / "amr"
    assert hashlib.sha256(upstream.read_bytes()).hexdigest() == before
    report = json.loads((tmp_path / "build/lagramses-configure.json").read_text())
    assert report["mpi_openmp_probe"] == "skipped"
    assert report["full_solver_build_validated"] is False
    assert report["upstream_makefile_sha256"] == before
    assert subprocess.run(["make", "-C", str(generated.parent), "-n"], capture_output=True).returncode == 0


def test_mpi_openmp_probe_links_without_launching_the_executable(tmp_path):
    compiler = tmp_path / "virtual-mpiifx"
    compiler.write_text("#!" + sys.executable + "\nimport sys\nfrom pathlib import Path\nif '--version' in sys.argv: print('Intel virtual compiler')\nelse:\n p=Path(sys.argv[sys.argv.index('-o')+1]);p.write_text('not an executable; probe must not be run')\n")
    compiler.chmod(0o700)
    result = invoke(tmp_path / "build", "--compiler", "intel", "--fc", str(compiler))
    assert result.returncode == 0, result.stderr
    report = json.loads((tmp_path / "build/lagramses-configure.json").read_text())
    assert report["mpi_openmp_probe"] == "passed"
    assert "virtual compiler" in report["compiler_version"]


def test_unrelated_build_directory_is_preserved(tmp_path):
    (tmp_path / "keep").write_text("user data")
    result = invoke(tmp_path, "--no-probe")
    assert result.returncode != 0 and "Nonempty" in result.stderr
    assert (tmp_path / "keep").read_text() == "user data"
    assert not (tmp_path / "bin").exists()


def test_conflicting_build_link_is_refused(tmp_path):
    root = tmp_path / "build"
    assert invoke(root, "--no-probe").returncode == 0
    (root / "amr").unlink()
    (root / "amr").mkdir()
    (root / "amr/keep").write_text("user data")
    result = invoke(root, "--no-probe")
    assert result.returncode != 0 and "conflicting" in result.stderr
    assert (root / "amr/keep").read_text() == "user data"


def test_experimental_gnu_profile_changes_only_generated_flags(tmp_path):
    before = (SOURCE / "bin/Makefile").read_bytes()
    result = invoke(tmp_path / "gnu", "--compiler", "gnu", "--no-probe")
    assert result.returncode == 0 and "experimental" in result.stdout
    text = (tmp_path / "gnu/bin/Makefile").read_text()
    assert "-fopenmp" in text and "-ffree-line-length-none" in text
    assert (SOURCE / "bin/Makefile").read_bytes() == before
