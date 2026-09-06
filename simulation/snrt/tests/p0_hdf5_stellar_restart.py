"""P0.3 HDF5 stellar restart-state contract test.

This is a bounded schema/bitwise test.  By default it verifies the serialized
contract and legacy fail-closed behavior using a synthetic particle group and
does not mutate durable evidence.  With ``--fortran-runtime`` it additionally
launches the linked HDF5 production binary in temporary fixture directories.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
from tempfile import TemporaryDirectory

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
BACKUP = ROOT / "patch" / "lagRamses" / "backup_hdf5.f90"
RESTORE = ROOT / "patch" / "lagRamses" / "restore_hdf5.f90"
EVIDENCE = ROOT / "simulation" / "snrt" / "data" / "p0_hdf5_stellar_restart_contract.json"

SCHEMA_VERSION = 1
RELEASE_FIELDS = ("tpp", "mp0", "indtab")


def _synthetic_arrays() -> dict[str, np.ndarray]:
    return {
        "birth_epoch": np.array([0.0, -0.0, 0.25, 3.5], dtype=np.float64),
        "mass": np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64),
        "tpp": np.array([0.0, 1.0e-15, 0.125, 2.75], dtype=np.float64),
        "mp0": np.array([8.0, 16.0, 32.0, 64.0], dtype=np.float64),
        "indtab": np.array([0.0, 0.125, 1.0, 2.5], dtype=np.float64),
    }


def _write_synthetic_checkpoint(path: Path, include_release_state: bool) -> None:
    arrays = _synthetic_arrays()
    with h5py.File(path, "w") as handle:
        particles = handle.create_group("particles")
        particles.attrs["stellar_state_schema_version"] = SCHEMA_VERSION
        for name in ("birth_epoch", "mass"):
            particles.create_dataset(name, data=arrays[name])
        if include_release_state:
            for name in RELEASE_FIELDS:
                particles.create_dataset(name, data=arrays[name])


def _read_and_validate(path: Path) -> dict[str, bytes]:
    with h5py.File(path, "r") as handle:
        particles = handle["particles"]
        if "stellar_state_schema_version" not in particles.attrs:
            raise RuntimeError("missing stellar_state_schema_version")
        version = int(particles.attrs["stellar_state_schema_version"])
        if version != SCHEMA_VERSION:
            raise RuntimeError(f"unsupported stellar state schema {version}")
        missing = [name for name in RELEASE_FIELDS if name not in particles]
        if missing:
            raise RuntimeError("missing stellar release fields: " + ",".join(missing))
        return {name: np.asarray(particles[name], dtype=np.float64).tobytes() for name in RELEASE_FIELDS}


def _write_legacy_compatibility_table(path: Path) -> None:
    """Write the tiny legacy table required only by init_part startup."""
    elements = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
    rows = []
    for age in (0.0, 1.0e9):
        rows.append([age, 0.0, 0.0, 0.0, 0.0] + [0.0] * len(elements))
    lines = [
        f"{'nmetal =':<8}{1:10d}",
        f"{'nsteps =':<8}{len(rows):10d}",
        f"{'nelt':<11}{len(elements):10d}",
        f"{'elements':<15}" + "".join(f"{element:>3}" for element in elements),
        " " + f"{0.0:16.5E}",
    ]
    lines.extend("".join(f"{value:19.11E}" for value in row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_ascii_particle_ic(path: Path) -> None:
    """Write a tiny non-empty particle IC for the linked writer fixture."""
    rows = (
        "0.0 0.0 0.0 0.0 0.0 0.0 1.0",
        "1.0 0.0 0.0 0.0 0.0 0.0 2.0",
        "0.0 1.0 0.0 0.0 0.0 0.0 3.0",
    )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _read_fortran_records(path: Path) -> list[bytes]:
    """Read the small native sequential-unformatted particle file."""
    records = []
    with path.open("rb") as handle:
        while True:
            marker = handle.read(4)
            if not marker:
                break
            if len(marker) != 4:
                raise RuntimeError(f"truncated Fortran record marker in {path}")
            size = struct.unpack("<i", marker)[0]
            if size < 0:
                raise RuntimeError(f"unexpected chained Fortran record in {path}")
            payload = handle.read(size)
            if len(payload) != size:
                raise RuntimeError(f"truncated Fortran record payload in {path}")
            closing = handle.read(4)
            if len(closing) != 4 or struct.unpack("<i", closing)[0] != size:
                raise RuntimeError(f"invalid Fortran record trailer in {path}")
            records.append(payload)
    return records


def _run_binary(
    binary: Path,
    namelist: Path,
    run_directory: Path,
    nproc: int = 1,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment.update(
        {
            "I_MPI_FABRICS": "shm",
            "I_MPI_PIN": "0",
            "OMP_NUM_THREADS": "1",
            "OMP_PROC_BIND": "false",
            "OMP_WAIT_POLICY": "passive",
        }
    )
    environment.pop("PHASE0_YIELD_TABLE", None)
    return subprocess.run(
        ["mpirun", "-np", str(nproc), str(binary), str(namelist)],
        cwd=run_directory,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
        timeout=180,
    )


def _read_part_descriptor(path: Path) -> list[tuple[int, str, str]]:
    """Read the native particle descriptor emitted beside a binary output."""
    entries = []
    pattern = re.compile(r"^\s*(\d+)\s*,\s*([^,]+?)\s*,\s*([?bhiqfd])\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            entries.append((int(match.group(1)), match.group(2).strip(), match.group(3)))
    if not entries:
        raise AssertionError(f"particle descriptor has no variable records: {path}")
    return entries


def _binary_release_values(
    binary_output: Path,
    descriptor: Path,
    records: list[bytes],
) -> tuple[dict[str, list[float]], list[str], dict[str, int]]:
    """Map release fields through observed descriptor indices, not constants."""
    entries = _read_part_descriptor(descriptor)
    by_name = {name: (ivar, kind) for ivar, name, kind in entries}
    missing = [name for name in RELEASE_FIELDS if name not in by_name]
    if missing:
        raise AssertionError(f"binary descriptor lacks release fields {missing}: {descriptor}")
    values = {}
    record_indices = {}
    for name in RELEASE_FIELDS:
        ivar, kind = by_name[name]
        if kind != "d":
            raise AssertionError(f"binary release field {name} has unexpected type {kind}")
        # Eight header records precede descriptor ivar=1 in backup_part.
        record_index = 8 + ivar - 1
        if record_index >= len(records):
            raise AssertionError(
                f"descriptor ivar {ivar} for {name} exceeds {len(records)} records in {binary_output}"
            )
        values[name] = np.frombuffer(records[record_index], dtype="<f8").copy().tolist()
        record_indices[name] = record_index
    ordered = [name for _, name, _ in sorted(entries) if name in RELEASE_FIELDS]
    if ordered != list(RELEASE_FIELDS):
        raise AssertionError(f"binary release descriptor order={ordered}, expected={list(RELEASE_FIELDS)}")
    return values, ordered, record_indices


def _binary_provenance(binary: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    with binary.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    stat = binary.stat()
    return {
        "path": str(binary),
        "sha256": digest.hexdigest(),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "source_inputs": _source_provenance(),
    }


def _source_provenance() -> dict[str, object]:
    """Bind the linked binary evidence to the restart source files and git state."""
    source_paths = (
        ROOT / "patch" / "cuRamses" / "ramses_hdf5_io.f90",
        ROOT / "patch" / "lagRamses" / "restore_hdf5.f90",
        ROOT / "patch" / "lagRamses" / "backup_hdf5.f90",
        ROOT / "patch" / "cuRamses" / "output_part.f90",
        ROOT / "bin" / "Makefile",
    )
    source_hashes = {}
    for path in source_paths:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        source_hashes[str(path.relative_to(ROOT))] = digest.hexdigest()
    head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    status = subprocess.check_output(
        ["git", "status", "--short"], cwd=ROOT, text=True
    ).splitlines()
    return {
        "git_head": head,
        "worktree_dirty": bool(status),
        "git_status_short": status,
        "sha256": source_hashes,
    }


def _rewrite_checkpoint_for_same_ncpu(path: Path, ncpu: int) -> None:
    """Make a bounded metadata-only ncpu=2 view of a one-rank fixture."""
    with h5py.File(path, "r+") as handle:
        header = handle["header"]
        header.attrs["ncpu"] = ncpu
        particles = handle["particles"]
        old_counts = np.asarray(particles["npart_per_cpu"][...], dtype=np.int32)
        if old_counts.size != 1:
            raise AssertionError(f"fixture expected one original particle-count entry, got {old_counts}")
        del particles["npart_per_cpu"]
        particles.create_dataset(
            "npart_per_cpu", data=np.array([2, 1], dtype=np.int32)
        )
        grid_datasets = []
        def collect(name: str, item: object) -> None:
            if isinstance(item, h5py.Dataset) and name.rsplit("/", 1)[-1] == "ngrid_per_cpu":
                grid_datasets.append(name)
        handle.visititems(collect)
        for name in grid_datasets:
            values = np.asarray(handle[name][...], dtype=np.int32)
            if values.size != 1:
                raise AssertionError(f"fixture expected one original grid-count entry for {name}")
            del handle[name]
            handle.create_dataset(
                name, data=np.array([int(values[0]), 0], dtype=np.int32)
            )


def _materialize_sink_payload(path: Path, levelmin: int, nlevelmax: int) -> None:
    """Install one self-consistent sink row for the linked restart fixture."""
    with h5py.File(path, "r+") as handle:
        sinks = handle.require_group("sinks")
        sinks.attrs["nsink"] = 1
        sinks.attrs["nindsink"] = 1
        sinks.attrs["levelmin"] = levelmin
        sinks.attrs["nlevelmax"] = nlevelmax
        integer_fields = {"idsink": np.array([1], dtype=np.int32)}
        double_fields = {
            "msink": 1.0,
            "tsink": 0.0,
            "dMsmbh": 0.0,
            "dMBH_coarse": 0.0,
            "dMEd_coarse": 0.0,
            "Esave": 0.0,
            "spinmag": 0.0,
            "eps_sink": 0.1,
        }
        for dimension in range(1, 4):
            double_fields[f"xsink_{dimension}"] = 0.5
            double_fields[f"vsink_{dimension}"] = 0.0
            double_fields[f"jsink_{dimension}"] = 0.0
            double_fields[f"bhspin_{dimension}"] = 0.0
        for name, values in integer_fields.items():
            if name in sinks:
                del sinks[name]
            sinks.create_dataset(name, data=values)
        for name, value in double_fields.items():
            if name in sinks:
                del sinks[name]
            sinks.create_dataset(name, data=np.array([value], dtype=np.float64))
        for statistic in range(1, 2 * 3 + 2):
            for level in range(levelmin, nlevelmax + 1):
                name = f"sink_stat_{statistic}_{level}"
                if name in sinks:
                    del sinks[name]
                sinks.create_dataset(name, data=np.zeros(1, dtype=np.float64))


def _first_dataset_named(path: Path, basename: str) -> str:
    matches: list[str] = []
    with h5py.File(path, "r") as handle:
        def collect(name: str, item: object) -> None:
            if isinstance(item, h5py.Dataset) and name.rsplit("/", 1)[-1] == basename:
                matches.append(name)
        handle.visititems(collect)
    if not matches:
        raise AssertionError(f"checkpoint contains no {basename} dataset: {path}")
    return sorted(matches)[0]


def _run_fortran_hdf5_roundtrip(root: Path) -> dict[str, object]:
    """Exercise the linked HDF5 writer and reader on a tiny sink fixture."""
    binary = ROOT / "bin" / "ramses_final3d"
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise RuntimeError(f"HDF5-enabled production binary is unavailable: {binary}")
    if shutil.which("mpirun") is None:
        raise RuntimeError("mpirun is unavailable for the linked HDF5 runtime check")

    source_nml = ROOT / "tests" / "sink" / "stellar-HII" / "stellar-HII.nml"
    source_ic_sink = ROOT / "tests" / "sink" / "stellar-HII" / "ic_sink"
    with TemporaryDirectory(prefix="p0h5_", dir=ROOT / "build") as temporary_directory:
        run_root = Path(temporary_directory)
        continuous = run_root / "c"
        restart = run_root / "r"
        continuous.mkdir()
        restart.mkdir()
        shutil.copy2(source_ic_sink, continuous / "ic_sink")
        _write_ascii_particle_ic(continuous / "ic_part")

        table_path = continuous / "legacy_compatibility_yields.dat"
        _write_legacy_compatibility_table(table_path)
        namelist_text = source_nml.read_text(encoding="utf-8")
        # ``stellar`` is not a RUN_PARAMS variable in the active lagRamses
        # namelist; sink=True is sufficient to allocate the release arrays
        # for this serialization fixture.  Keep RT disabled so this remains
        # a small HDF5 particle test.
        namelist_text = namelist_text.split("&RT_PARAMS", 1)[0]
        namelist_text = namelist_text.replace(
            "&INIT_PARAMS\n", "&INIT_PARAMS\nfiletype='ascii'\ninitfile(1)='.'\n", 1
        )
        namelist_text = namelist_text.replace("stellar=.true.\n", "", 1)
        namelist_text = namelist_text.replace("rt=.true.\n", "rt=.false.\n", 1)
        namelist_text = namelist_text.replace("nremap=10", "nremap=0", 1)
        namelist_text = namelist_text.replace("levelmin=6", "levelmin=4", 1)
        namelist_text = namelist_text.replace("levelmax=8", "levelmax=4", 1)
        namelist_text = namelist_text.replace("ngridmax=350000", "ngridmax=1000", 1)
        namelist_text = namelist_text.replace("npartmax=50000", "npartmax=100", 1)
        namelist_text = namelist_text.replace("nlevelmax_sink=7", "nlevelmax_sink=4", 1)
        namelist_text = namelist_text.replace(
            "noutput=1\ntout=0.18",
            "noutput=1\ntout=0.18\ninformat='hdf5'\noutformat='hdf5'\nfoutput=1000000\nfbackup=1000000",
            1,
        )
        namelist_text = namelist_text.replace(
            "&COOLING_PARAMS",
            f"&PHYSICS_PARAMS\nyieldtablefilename='{table_path}'\n/\n\n&COOLING_PARAMS",
            1,
        )
        continuous_nml = continuous / "run.nml"
        continuous_nml.write_text(namelist_text, encoding="utf-8")
        continuous_result = _run_binary(binary, continuous_nml, continuous)
        (continuous / "run.log").write_text(continuous_result.stdout, encoding="utf-8")
        if continuous_result.returncode != 0:
            raise RuntimeError(
                f"linked HDF5 writer run failed (rc={continuous_result.returncode}):\n"
                + continuous_result.stdout[-12000:]
            )

        outputs = sorted(continuous.glob("output_*/data_*.h5"))
        if len(outputs) != 1:
            raise AssertionError(
                f"expected one bounded HDF5 output, found {outputs}; "
                f"files={[str(path.relative_to(run_root)) for path in run_root.rglob('*') if path.is_file()]}\n"
                + continuous_result.stdout[-16000:]
            )
        output_path = outputs[0]
        with h5py.File(output_path, "r") as handle:
            if "particles" not in handle:
                raise RuntimeError(
                    "linked HDF5 output lacks /particles; groups="
                    + ",".join(sorted(handle.keys()))
                )
            particles = handle["particles"]
            version = int(np.asarray(particles.attrs["stellar_state_schema_version"]).reshape(-1)[0])
            missing = [name for name in RELEASE_FIELDS if name not in particles]
            release_lengths = {name: int(particles[name].shape[0]) for name in RELEASE_FIELDS}
            writer_release_values = {name: particles[name][...].tolist() for name in RELEASE_FIELDS}
            particle_types = np.asarray(particles["ptypep"][...], dtype=np.int32)
            release_nonempty = all(release_lengths[name] > 0 for name in RELEASE_FIELDS)
            if not release_nonempty:
                raise AssertionError("linked HDF5 writer produced empty release datasets")
            writer_release_zero = all(
                np.all(np.asarray(particles[name][...], dtype=np.float64) == 0.0)
                for name in RELEASE_FIELDS
            )
            writer_zero_expected = {
                name: np.zeros(release_lengths[name], dtype=np.float64)
                for name in RELEASE_FIELDS
            }
            writer_zero_match = all(
                np.array_equal(
                    np.asarray(particles[name][...], dtype=np.float64),
                    writer_zero_expected[name],
                )
                for name in RELEASE_FIELDS
            )
            writer_release_values = {
                name: np.asarray(particles[name][...], dtype=np.float64).copy().tolist()
                for name in RELEASE_FIELDS
            }
        # Exercise the same-ncpu ksection checkpoint-tree branch with a real
        # two-rank writer output.  This is separate from the Hilbert fixture
        # below, whose ncpu metadata is intentionally rewritten only for the
        # cross-rank reader contract.
        ksection_initial = run_root / "ks"
        ksection_initial.mkdir()
        shutil.copy2(source_ic_sink, ksection_initial / "ic_sink")
        _write_ascii_particle_ic(ksection_initial / "ic_part")
        ksection_table = ksection_initial / "legacy_compatibility_yields.dat"
        _write_legacy_compatibility_table(ksection_table)
        ksection_text = namelist_text.replace(
            f"yieldtablefilename='{table_path}'",
            f"yieldtablefilename='{ksection_table}'",
        ).replace("pic=.true.\n", "pic=.true.\nordering='ksection'\n", 1)
        ksection_nml = ksection_initial / "run.nml"
        ksection_nml.write_text(ksection_text, encoding="utf-8")
        ksection_initial_result = _run_binary(
            binary, ksection_nml, ksection_initial, nproc=2
        )
        (ksection_initial / "run.log").write_text(
            ksection_initial_result.stdout, encoding="utf-8"
        )
        if ksection_initial_result.returncode != 0:
            raise RuntimeError(
                f"linked ksection writer run failed "
                f"(rc={ksection_initial_result.returncode}):\n"
                + ksection_initial_result.stdout[-12000:]
            )
        ksection_outputs = sorted(ksection_initial.glob("output_*/data_*.h5"))
        if len(ksection_outputs) != 1:
            raise AssertionError(
                f"expected one ksection HDF5 output, found {ksection_outputs}"
            )
        ksection_output = ksection_outputs[0]
        with h5py.File(ksection_output, "r") as handle:
            domain = handle["domain"]
            if "ksec_kmax" not in domain.attrs or "ksec_nbinodes" not in domain.attrs:
                raise AssertionError("ksection writer omitted tree dimension attributes")
            ksection_dimensions = (
                int(np.asarray(domain.attrs["ksec_kmax"]).reshape(-1)[0]),
                int(np.asarray(domain.attrs["ksec_nbinodes"]).reshape(-1)[0]),
            )
            if ksection_dimensions[0] < 2 or ksection_dimensions[1] < 1:
                raise AssertionError(f"invalid ksection writer dimensions {ksection_dimensions}")

        ksection_restart = run_root / "kr"
        ksection_restart.mkdir()
        shutil.copytree(ksection_output.parent, ksection_restart / ksection_output.parent.name)
        ksection_restart_text = ksection_text.replace("nrestart=0", "nrestart=1", 1).replace(
            "tout=0.18", "tout=0.25", 1
        )
        ksection_restart_nml = ksection_restart / "run.nml"
        ksection_restart_nml.write_text(ksection_restart_text, encoding="utf-8")
        ksection_restart_result = _run_binary(
            binary, ksection_restart_nml, ksection_restart, nproc=2
        )
        (ksection_restart / "run.log").write_text(
            ksection_restart_result.stdout, encoding="utf-8"
        )
        if ksection_restart_result.returncode != 0:
            raise RuntimeError(
                f"linked ksection restart failed "
                f"(rc={ksection_restart_result.returncode}):\n"
                + ksection_restart_result.stdout[-12000:]
            )
        if "HDF5: ksection tree restored" not in ksection_restart_result.stdout:
            raise AssertionError("ksection restart did not use the checkpoint tree branch")

        ksection_negative = run_root / "kn"
        ksection_negative.mkdir()
        shutil.copytree(ksection_output.parent, ksection_negative / ksection_output.parent.name)
        ksection_negative_file = ksection_negative / ksection_output.parent.name / ksection_output.name
        with h5py.File(ksection_negative_file, "r+") as handle:
            handle["domain"].attrs["ksec_nbinodes"] = 1000000
        ksection_negative_nml = ksection_negative / "run.nml"
        ksection_negative_nml.write_text(ksection_restart_text, encoding="utf-8")
        ksection_negative_result = _run_binary(
            binary, ksection_negative_nml, ksection_negative, nproc=2
        )
        (ksection_negative / "run.log").write_text(
            ksection_negative_result.stdout, encoding="utf-8"
        )
        if ksection_negative_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted oversized HDF5 ksection dimensions; output:\n"
                + ksection_negative_result.stdout[-12000:]
            )
        if "invalid HDF5 ksection dimensions" not in ksection_negative_result.stdout:
            raise AssertionError(
                "oversized-ksection negative run did not expose dimension diagnostic"
            )
        if version != SCHEMA_VERSION or missing:
            raise AssertionError(
                f"invalid linked HDF5 stellar payload: version={version} missing={missing} "
                f"release_lengths={release_lengths} writer_values={writer_release_values}"
            )

        # Materialize a nonzero PTYPE_STAR state only in this temporary
        # checkpoint.  The linked restart reader restores it into the real
        # particle arrays, and the subsequent linked HDF5 writer must emit the
        # same payload without Python rewriting the writer output.
        restart_seed_types = np.ones(release_lengths["tpp"], dtype=np.int32)
        restart_seed_release = {
            "tpp": np.array([0.125, 0.25, 0.5], dtype=np.float64),
            "mp0": np.array([8.0, 16.0, 32.0], dtype=np.float64),
            "indtab": np.array([1.0, 2.0, 3.0], dtype=np.float64),
        }
        if release_lengths["tpp"] != len(restart_seed_types):
            raise AssertionError("restart seed assumes three active fixture particles")
        with h5py.File(output_path, "r+") as handle:
            particles = handle["particles"]
            particles["ptypep"][...] = restart_seed_types
            for name, values in restart_seed_release.items():
                particles[name][...] = values
        with h5py.File(output_path, "r") as handle:
            restart_seed_types = np.asarray(handle["particles"]["ptypep"][...], dtype=np.int32)
            expected_release = {
                name: np.asarray(handle["particles"][name][...], dtype=np.float64).copy()
                for name in RELEASE_FIELDS
            }
        if not np.all(restart_seed_types == 1):
            raise AssertionError("temporary restart seed did not materialize PTYPE_STAR")
        expected_star_sample = tuple(
            float(expected_release[name][0]) for name in RELEASE_FIELDS
        )
        # The original tiny stellar-HII evolution can finish with no active
        # sinks in its first output.  Materialize one complete, physically
        # bounded sink row in the temporary checkpoint so the linked restart
        # exercises every sink reader and the nindsink consistency guard.
        _materialize_sink_payload(output_path, levelmin=4, nlevelmax=4)
        checkpoint = output_path.parent
        shutil.copytree(checkpoint, restart / checkpoint.name)
        restart_text = namelist_text.replace("nrestart=0", "nrestart=1", 1).replace(
            "tout=0.18", "tout=0.25", 1
        )
        restart_nml = restart / "run.nml"
        restart_nml.write_text(restart_text, encoding="utf-8")
        restart_result = _run_binary(binary, restart_nml, restart)
        (restart / "run.log").write_text(restart_result.stdout, encoding="utf-8")
        if restart_result.returncode != 0:
            raise RuntimeError(
                f"linked HDF5 reader run failed (rc={restart_result.returncode}):\n"
                + restart_result.stdout[-12000:]
            )
        if "HDF5 particle restore done." not in restart_result.stdout:
            raise AssertionError("linked HDF5 reader did not report particle restoration")
        sample_match = re.search(
            r"HDF5 stellar release state sample:\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)",
            restart_result.stdout,
        )
        if sample_match is None:
            raise AssertionError("linked HDF5 reader did not expose a release-state sample")
        restored_sample = tuple(float(value) for value in sample_match.groups())
        expected_sample = expected_star_sample
        if restored_sample != expected_sample:
            raise AssertionError(
                f"linked HDF5 reader restored wrong release-state sample: "
                f"actual={restored_sample} expected={expected_sample}"
            )

        restart_outputs = sorted(
            path
            for path in restart.glob("output_*/data_*.h5")
            if path.parent.name != checkpoint.name
        )
        if len(restart_outputs) != 1:
            raise AssertionError(
                f"expected one post-restart HDF5 writer output, found {restart_outputs};\n"
                + restart_result.stdout[-16000:]
            )
        restart_writer_output = restart_outputs[0]
        with h5py.File(restart_writer_output, "r") as handle:
            restart_writer_particles = handle["particles"]
            restart_writer_types = np.asarray(
                restart_writer_particles["ptypep"][...], dtype=np.int32
            )
            restart_writer_release_values = {
                name: np.asarray(restart_writer_particles[name][...], dtype=np.float64).copy()
                for name in RELEASE_FIELDS
            }
            if "sinks" not in handle:
                raise AssertionError("post-restart HDF5 writer omitted /sinks")
            restart_writer_sinks = handle["sinks"]
            restart_writer_nsink = int(
                np.asarray(restart_writer_sinks.attrs["nsink"]).reshape(-1)[0]
            )
            sink_writer_required = [
                "idsink",
                "msink",
                "xsink_1",
                "xsink_2",
                "xsink_3",
                "vsink_1",
                "vsink_2",
                "vsink_3",
                "tsink",
                "dMsmbh",
                "dMBH_coarse",
                "dMEd_coarse",
                "Esave",
                "jsink_1",
                "jsink_2",
                "jsink_3",
                "bhspin_1",
                "bhspin_2",
                "bhspin_3",
                "spinmag",
                "eps_sink",
            ]
            missing_sink_writer = [
                name for name in sink_writer_required if name not in restart_writer_sinks
            ]
            sink_writer_lengths = {
                name: int(restart_writer_sinks[name].shape[0])
                for name in sink_writer_required
                if name in restart_writer_sinks
            }
            sink_writer_roundtrip = (
                restart_writer_nsink > 0
                and not missing_sink_writer
                and all(
                    length == restart_writer_nsink
                    for length in sink_writer_lengths.values()
                )
            )
            if not sink_writer_roundtrip:
                raise AssertionError(
                    "post-restart HDF5 writer sink payload is incomplete: "
                    f"nsink={restart_writer_nsink} missing={missing_sink_writer} "
                    f"lengths={sink_writer_lengths}"
                )
        restart_writer_nonzero = all(
            np.any(np.abs(values) > 0.0)
            for values in restart_writer_release_values.values()
        )
        restart_writer_payload_match = (
            np.all(restart_writer_types == 1)
            and all(
                np.array_equal(restart_writer_release_values[name], expected_release[name])
                for name in RELEASE_FIELDS
            )
        )
        if not restart_writer_nonzero or not restart_writer_payload_match:
            raise AssertionError(
                "post-restart linked HDF5 writer did not preserve the nonzero "
                f"PTYPE_STAR payload: types={restart_writer_types.tolist()} "
                f"values={ {name: values.tolist() for name, values in restart_writer_release_values.items()} }"
            )

        # The same restored, nonzero PTYPE_STAR checkpoint is also written in
        # native binary format.  Compare that observed stream with the HDF5
        # output using the descriptor emitted by backup_part, so the test does
        # not carry a second implicit record-order contract.
        binary_restart = run_root / "b"
        binary_restart.mkdir()
        shutil.copytree(checkpoint, binary_restart / checkpoint.name)
        binary_restart_nml = binary_restart / "run.nml"
        binary_restart_text = restart_text.replace("outformat='hdf5'", "outformat='binary'", 1)
        binary_restart_nml.write_text(binary_restart_text, encoding="utf-8")
        binary_restart_result = _run_binary(binary, binary_restart_nml, binary_restart)
        (binary_restart / "run.log").write_text(
            binary_restart_result.stdout, encoding="utf-8"
        )
        if binary_restart_result.returncode != 0:
            raise RuntimeError(
                f"linked binary output from HDF5 restart failed "
                f"(rc={binary_restart_result.returncode}):\n"
                + binary_restart_result.stdout[-12000:]
            )
        binary_outputs = sorted(binary_restart.glob("output_*/part_*.out*"))
        if len(binary_outputs) != 1:
            raise AssertionError(
                f"expected one bounded binary particle output, found {binary_outputs}"
            )
        binary_output = binary_outputs[0]
        binary_descriptor = binary_output.parent / "part_file_descriptor.txt"
        if not binary_descriptor.is_file():
            raise AssertionError(f"linked binary writer did not emit {binary_descriptor}")
        binary_records = _read_fortran_records(binary_output)
        binary_release_values, binary_record_order, binary_record_indices = _binary_release_values(
            binary_output, binary_descriptor, binary_records
        )
        binary_value_match = all(
            np.asarray(binary_release_values[name], dtype=np.float64).tobytes()
            == np.asarray(restart_writer_release_values[name], dtype=np.float64).tobytes()
            for name in RELEASE_FIELDS
        )
        binary_nonzero = all(
            np.any(np.abs(np.asarray(binary_release_values[name], dtype=np.float64)) > 0.0)
            for name in RELEASE_FIELDS
        )
        if not binary_nonzero or not binary_value_match:
            raise AssertionError(
                f"binary/HDF5 nonzero release-state mismatch: "
                f"binary={binary_release_values} hdf5={ {name: values.tolist() for name, values in restart_writer_release_values.items()} }"
            )

        cross_ncpu = run_root / "v"
        cross_ncpu.mkdir()
        shutil.copytree(checkpoint, cross_ncpu / checkpoint.name)
        cross_ncpu_nml = cross_ncpu / "run.nml"
        cross_ncpu_nml.write_text(restart_text, encoding="utf-8")
        cross_ncpu_result = _run_binary(
            binary, cross_ncpu_nml, cross_ncpu, nproc=4
        )
        (cross_ncpu / "run.log").write_text(
            cross_ncpu_result.stdout, encoding="utf-8"
        )
        if cross_ncpu_result.returncode != 0:
            raise RuntimeError(
                f"linked HDF5 ncpu-file-to-ncpu=4 reader failed "
                f"(rc={cross_ncpu_result.returncode}):\n"
                + cross_ncpu_result.stdout[-12000:]
            )
        if "HDF5 particle restore done." not in cross_ncpu_result.stdout:
            raise AssertionError("ncpu-file-to-ncpu=4 reader did not restore particles")
        cross_local_count_markers = cross_ncpu_result.stdout.count(
            "HDF5 stellar release local count rank="
        )
        if cross_local_count_markers < 4:
            raise AssertionError(
                "ncpu-file-to-ncpu=4 reader did not report all rank-local counts"
            )
        cross_state_lines = re.findall(
            r"HDF5 stellar release local state rank=\s*(\d+)\s*:\s*(\d+)\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)",
            cross_ncpu_result.stdout,
        )
        cross_empty_ranks = re.findall(
            r"HDF5 stellar release local state rank=\s*(\d+)\s*:\s*0\s+EMPTY",
            cross_ncpu_result.stdout,
        )
        expected_rank_samples = {
            (1, 1): expected_star_sample,
            (2, 1): tuple(float(expected_release[name][1]) for name in RELEASE_FIELDS),
            (3, 1): tuple(float(expected_release[name][2]) for name in RELEASE_FIELDS),
        }
        cross_rank_value_match = True
        for rank_text, count_text, *sample_text in cross_state_lines:
            key = (int(rank_text), int(count_text))
            if key not in expected_rank_samples:
                cross_rank_value_match = False
                break
            cross_rank_value_match = cross_rank_value_match and tuple(
                float(value) for value in sample_text
            ) == expected_rank_samples[key]
        cross_rank_value_match = cross_rank_value_match and len(cross_state_lines) == 3
        cross_zero_rank_match = len(cross_empty_ranks) == 1
        if not cross_rank_value_match or not cross_zero_rank_match:
            raise AssertionError(
                f"cross-ncpu per-rank release state mismatch: states={cross_state_lines} "
                f"empty={cross_empty_ranks}"
            )

        same_ncpu = run_root / "s"
        same_ncpu.mkdir()
        shutil.copytree(checkpoint, same_ncpu / checkpoint.name)
        same_ncpu_file = same_ncpu / checkpoint.name / output_path.name
        _rewrite_checkpoint_for_same_ncpu(same_ncpu_file, ncpu=2)
        same_ncpu_nml = same_ncpu / "run.nml"
        same_ncpu_nml.write_text(restart_text, encoding="utf-8")
        same_ncpu_result = _run_binary(
            binary, same_ncpu_nml, same_ncpu, nproc=2
        )
        (same_ncpu / "run.log").write_text(
            same_ncpu_result.stdout, encoding="utf-8"
        )
        if same_ncpu_result.returncode != 0:
            raise RuntimeError(
                f"linked same-ncpu=2 HDF5 reader failed "
                f"(rc={same_ncpu_result.returncode}):\n"
                + same_ncpu_result.stdout[-12000:]
            )
        same_ncpu_state_lines = re.findall(
            r"HDF5 stellar release local state rank=\s*(\d+)\s*:\s*(\d+)\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)\s+"
            r"([+-]?[0-9.]+[Ee][+-][0-9]+)",
            same_ncpu_result.stdout,
        )
        same_ncpu_expected = {
            (1, 2): tuple(float(expected_release[name][0]) for name in RELEASE_FIELDS),
            (2, 1): tuple(float(expected_release[name][2]) for name in RELEASE_FIELDS),
        }
        same_ncpu_value_match = all(
            (int(rank), int(count)) in same_ncpu_expected
            and tuple(float(value) for value in sample)
            == same_ncpu_expected[(int(rank), int(count))]
            for rank, count, *sample in same_ncpu_state_lines
        ) and len(same_ncpu_state_lines) == 2
        if (
            "HDF5 particle restore done." not in same_ncpu_result.stdout
            or not same_ncpu_value_match
        ):
            raise AssertionError(
                f"same-ncpu=2 restore mismatch: states={same_ncpu_state_lines}\n"
                + same_ncpu_result.stdout[-12000:]
            )

        ngrid_dataset_name = _first_dataset_named(output_path, "ngrid_per_cpu")
        negative_ngrid = run_root / "gn"
        negative_ngrid.mkdir()
        shutil.copytree(checkpoint, negative_ngrid / checkpoint.name)
        negative_ngrid_file = negative_ngrid / checkpoint.name / output_path.name
        with h5py.File(negative_ngrid_file, "r+") as handle:
            del handle[ngrid_dataset_name]
            handle.create_dataset(ngrid_dataset_name, data=np.array([-1], dtype=np.int32))
        negative_ngrid_nml = negative_ngrid / "run.nml"
        negative_ngrid_nml.write_text(restart_text, encoding="utf-8")
        negative_ngrid_result = _run_binary(
            binary, negative_ngrid_nml, negative_ngrid, nproc=2
        )
        (negative_ngrid / "run.log").write_text(
            negative_ngrid_result.stdout, encoding="utf-8"
        )
        if negative_ngrid_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted negative ngrid_per_cpu entry; output:\n"
                + negative_ngrid_result.stdout[-12000:]
            )
        if "invalid ngrid_per_cpu: negative entry" not in negative_ngrid_result.stdout:
            raise AssertionError(
                "negative-ngrid_per_cpu run did not expose its entry diagnostic"
            )

        wrong_length_ngrid = run_root / "gl"
        wrong_length_ngrid.mkdir()
        shutil.copytree(checkpoint, wrong_length_ngrid / checkpoint.name)
        wrong_length_ngrid_file = wrong_length_ngrid / checkpoint.name / output_path.name
        with h5py.File(wrong_length_ngrid_file, "r+") as handle:
            values = np.asarray(handle[ngrid_dataset_name][...], dtype=np.int32)
            del handle[ngrid_dataset_name]
            handle.create_dataset(
                ngrid_dataset_name,
                data=np.concatenate([values, np.array([0], dtype=np.int32)]),
            )
        wrong_length_ngrid_nml = wrong_length_ngrid / "run.nml"
        wrong_length_ngrid_nml.write_text(restart_text, encoding="utf-8")
        wrong_length_ngrid_result = _run_binary(
            binary, wrong_length_ngrid_nml, wrong_length_ngrid, nproc=2
        )
        (wrong_length_ngrid / "run.log").write_text(
            wrong_length_ngrid_result.stdout, encoding="utf-8"
        )
        if wrong_length_ngrid_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted wrong-length ngrid_per_cpu dataset; output:\n"
                + wrong_length_ngrid_result.stdout[-12000:]
            )
        if "failed checked HDF5 read for ngrid_per_cpu" not in wrong_length_ngrid_result.stdout:
            raise AssertionError(
                "wrong-length ngrid_per_cpu run did not expose its extent diagnostic"
            )

        son_flag_dataset_name = _first_dataset_named(output_path, "son_flag")
        negative_amr_payload = run_root / "ap"
        negative_amr_payload.mkdir()
        shutil.copytree(checkpoint, negative_amr_payload / checkpoint.name)
        negative_amr_payload_file = (
            negative_amr_payload / checkpoint.name / output_path.name
        )
        with h5py.File(negative_amr_payload_file, "r+") as handle:
            values = np.asarray(handle[son_flag_dataset_name][...], dtype=np.int32)
            if values.size < 2:
                raise AssertionError(
                    f"fixture {son_flag_dataset_name} is too small for an extent negative"
                )
            del handle[son_flag_dataset_name]
            handle.create_dataset(son_flag_dataset_name, data=values[:-1])
        negative_amr_payload_nml = negative_amr_payload / "run.nml"
        negative_amr_payload_nml.write_text(restart_text, encoding="utf-8")
        negative_amr_payload_result = _run_binary(
            binary, negative_amr_payload_nml, negative_amr_payload, nproc=4
        )
        (negative_amr_payload / "run.log").write_text(
            negative_amr_payload_result.stdout, encoding="utf-8"
        )
        if negative_amr_payload_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted truncated AMR son_flag payload; output:\n"
                + negative_amr_payload_result.stdout[-12000:]
            )
        if "checked AMR son_flag read failed" not in negative_amr_payload_result.stdout:
            raise AssertionError(
                "truncated-AMR son_flag negative run did not expose its extent diagnostic"
            )

        negative_nsink = run_root / "sk"
        negative_nsink.mkdir()
        shutil.copytree(checkpoint, negative_nsink / checkpoint.name)
        negative_nsink_file = negative_nsink / checkpoint.name / output_path.name
        with h5py.File(negative_nsink_file, "r+") as handle:
            if "sinks" not in handle:
                raise AssertionError("fixture lacks /sinks group for nsink guard")
            handle["sinks"].attrs["nsink"] = 1000000
        negative_nsink_nml = negative_nsink / "run.nml"
        negative_nsink_nml.write_text(restart_text, encoding="utf-8")
        negative_nsink_result = _run_binary(
            binary, negative_nsink_nml, negative_nsink, nproc=2
        )
        (negative_nsink / "run.log").write_text(
            negative_nsink_result.stdout, encoding="utf-8"
        )
        if negative_nsink_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted nsink beyond nsinkmax; output:\n"
                + negative_nsink_result.stdout[-12000:]
            )
        if "invalid HDF5 nsink=" not in negative_nsink_result.stdout:
            raise AssertionError(
                "oversized-nsink negative run did not expose its bound diagnostic"
            )

        negative_nindsink = run_root / "si"
        negative_nindsink.mkdir()
        shutil.copytree(checkpoint, negative_nindsink / checkpoint.name)
        negative_nindsink_file = negative_nindsink / checkpoint.name / output_path.name
        with h5py.File(negative_nindsink_file, "r+") as handle:
            sinks = handle["sinks"]
            if int(np.asarray(sinks.attrs["nsink"]).reshape(-1)[0]) > 0:
                sinks.attrs["nindsink"] = 0
            else:
                raise AssertionError(
                    "fixture has no active sinks; nindsink consistency negative is not executable"
                )
        negative_nindsink_nml = negative_nindsink / "run.nml"
        negative_nindsink_nml.write_text(restart_text, encoding="utf-8")
        negative_nindsink_result = _run_binary(
            binary, negative_nindsink_nml, negative_nindsink, nproc=2
        )
        (negative_nindsink / "run.log").write_text(
            negative_nindsink_result.stdout, encoding="utf-8"
        )
        if negative_nindsink_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted nindsink below the idsink maximum; output:\n"
                + negative_nindsink_result.stdout[-12000:]
            )
        if "nindsink is below idsink maximum" not in negative_nindsink_result.stdout:
            raise AssertionError(
                "inconsistent-nindsink negative run did not expose its consistency diagnostic"
            )

        negative_stellar_ledger = run_root / "sl"
        negative_stellar_ledger.mkdir()
        shutil.copytree(checkpoint, negative_stellar_ledger / checkpoint.name)
        negative_stellar_ledger_file = (
            negative_stellar_ledger / checkpoint.name / output_path.name
        )
        with h5py.File(negative_stellar_ledger_file, "r+") as handle:
            del handle["particles"].attrs["nstar_tot"]
        negative_stellar_ledger_nml = negative_stellar_ledger / "run.nml"
        negative_stellar_ledger_nml.write_text(restart_text, encoding="utf-8")
        negative_stellar_ledger_result = _run_binary(
            binary, negative_stellar_ledger_nml, negative_stellar_ledger, nproc=2
        )
        (negative_stellar_ledger / "run.log").write_text(
            negative_stellar_ledger_result.stdout, encoding="utf-8"
        )
        if negative_stellar_ledger_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted HDF5 restart without nstar_tot; output:\n"
                + negative_stellar_ledger_result.stdout[-12000:]
            )
        if "failed checked HDF5 read for particles/nstar_tot" not in negative_stellar_ledger_result.stdout:
            raise AssertionError(
                "missing-nstar_tot negative run did not expose its ledger diagnostic"
            )

        negative_header = run_root / "ha"
        negative_header.mkdir()
        shutil.copytree(checkpoint, negative_header / checkpoint.name)
        negative_header_file = negative_header / checkpoint.name / output_path.name
        with h5py.File(negative_header_file, "r+") as handle:
            del handle["header"].attrs["aexp"]
        negative_header_nml = negative_header / "run.nml"
        negative_header_nml.write_text(restart_text, encoding="utf-8")
        negative_header_result = _run_binary(
            binary, negative_header_nml, negative_header, nproc=2
        )
        (negative_header / "run.log").write_text(
            negative_header_result.stdout, encoding="utf-8"
        )
        if negative_header_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted HDF5 restart without header/aexp; output:\n"
                + negative_header_result.stdout[-12000:]
            )
        if "failed checked HDF5 read for header/aexp" not in negative_header_result.stdout:
            raise AssertionError(
                "missing-header/aexp negative run did not expose its clock diagnostic"
            )

        negative_npart_negative = run_root / "nneg"
        negative_npart_negative.mkdir()
        shutil.copytree(checkpoint, negative_npart_negative / checkpoint.name)
        negative_npart_negative_file = (
            negative_npart_negative / checkpoint.name / output_path.name
        )
        with h5py.File(negative_npart_negative_file, "r+") as handle:
            particles = handle["particles"]
            del particles["npart_per_cpu"]
            particles.create_dataset("npart_per_cpu", data=np.array([-1], dtype=np.int32))
        negative_npart_negative_nml = negative_npart_negative / "run.nml"
        negative_npart_negative_nml.write_text(restart_text, encoding="utf-8")
        negative_npart_negative_result = _run_binary(
            binary, negative_npart_negative_nml, negative_npart_negative, nproc=2
        )
        (negative_npart_negative / "run.log").write_text(
            negative_npart_negative_result.stdout, encoding="utf-8"
        )
        if negative_npart_negative_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted negative npart_per_cpu entry; output:\n"
                + negative_npart_negative_result.stdout[-12000:]
            )
        if "invalid npart_per_cpu: negative entry" not in negative_npart_negative_result.stdout:
            raise AssertionError(
                "negative-npart_per_cpu run did not expose its entry diagnostic"
            )

        negative_npart_length = run_root / "nlen"
        negative_npart_length.mkdir()
        shutil.copytree(checkpoint, negative_npart_length / checkpoint.name)
        negative_npart_length_file = negative_npart_length / checkpoint.name / output_path.name
        with h5py.File(negative_npart_length_file, "r+") as handle:
            particles = handle["particles"]
            del particles["npart_per_cpu"]
            particles.create_dataset(
                "npart_per_cpu", data=np.array([3, 0], dtype=np.int32)
            )
        negative_npart_length_nml = negative_npart_length / "run.nml"
        negative_npart_length_nml.write_text(restart_text, encoding="utf-8")
        negative_npart_length_result = _run_binary(
            binary, negative_npart_length_nml, negative_npart_length, nproc=2
        )
        (negative_npart_length / "run.log").write_text(
            negative_npart_length_result.stdout, encoding="utf-8"
        )
        if negative_npart_length_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted wrong-length npart_per_cpu dataset; output:\n"
                + negative_npart_length_result.stdout[-12000:]
            )
        if "failed checked HDF5 read for npart_per_cpu" not in negative_npart_length_result.stdout:
            raise AssertionError(
                "wrong-length npart_per_cpu run did not expose its extent diagnostic"
            )

        negative_npart = run_root / "nn"
        negative_npart.mkdir()
        shutil.copytree(checkpoint, negative_npart / checkpoint.name)
        negative_npart_file = negative_npart / checkpoint.name / output_path.name
        with h5py.File(negative_npart_file, "r+") as handle:
            del handle["particles"]["npart_per_cpu"]
        negative_npart_nml = negative_npart / "run.nml"
        negative_npart_nml.write_text(restart_text, encoding="utf-8")
        negative_npart_result = _run_binary(
            binary, negative_npart_nml, negative_npart, nproc=2
        )
        (negative_npart / "run.log").write_text(
            negative_npart_result.stdout, encoding="utf-8"
        )
        if negative_npart_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted HDF5 restart without npart_per_cpu; output:\n"
                + negative_npart_result.stdout[-12000:]
            )
        if "failed checked HDF5 read for npart_per_cpu" not in negative_npart_result.stdout:
            raise AssertionError(
                "missing-npart_per_cpu negative run did not expose its specific diagnostic"
            )

        negative_count = run_root / "nc"
        negative_count.mkdir()
        shutil.copytree(checkpoint, negative_count / checkpoint.name)
        negative_count_file = negative_count / checkpoint.name / output_path.name
        with h5py.File(negative_count_file, "r+") as handle:
            particles = handle["particles"]
            del particles["npart_per_cpu"]
            particles.create_dataset("npart_per_cpu", data=np.array([2], dtype=np.int32))
        negative_count_nml = negative_count / "run.nml"
        negative_count_nml.write_text(restart_text, encoding="utf-8")
        negative_count_result = _run_binary(
            binary, negative_count_nml, negative_count, nproc=2
        )
        (negative_count / "run.log").write_text(
            negative_count_result.stdout, encoding="utf-8"
        )
        if negative_count_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted HDF5 restart with inconsistent npart_per_cpu sum; output:\n"
                + negative_count_result.stdout[-12000:]
            )
        if "invalid npart_per_cpu: sum=" not in negative_count_result.stdout:
            raise AssertionError(
                "inconsistent-npart_per_cpu negative run did not expose its sum diagnostic"
            )

        negative_schema = run_root / "ns"
        negative_schema.mkdir()
        shutil.copytree(checkpoint, negative_schema / checkpoint.name)
        negative_schema_file = negative_schema / checkpoint.name / output_path.name
        with h5py.File(negative_schema_file, "r+") as handle:
            del handle["particles"].attrs["stellar_state_schema_version"]
        negative_schema_nml = negative_schema / "run.nml"
        negative_schema_nml.write_text(restart_text, encoding="utf-8")
        negative_schema_result = _run_binary(
            binary, negative_schema_nml, negative_schema
        )
        (negative_schema / "run.log").write_text(
            negative_schema_result.stdout, encoding="utf-8"
        )
        if negative_schema_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted HDF5 restart without schema marker; output:\n"
                + negative_schema_result.stdout[-12000:]
            )
        if "stellar_state_schema_version" not in negative_schema_result.stdout:
            raise AssertionError(
                "missing-schema negative run did not expose the schema-specific diagnostic"
            )

        negative_schema_extent = run_root / "nse"
        negative_schema_extent.mkdir()
        shutil.copytree(checkpoint, negative_schema_extent / checkpoint.name)
        negative_schema_extent_file = (
            negative_schema_extent / checkpoint.name / output_path.name
        )
        with h5py.File(negative_schema_extent_file, "r+") as handle:
            handle["particles"].attrs["stellar_state_schema_version"] = np.array(
                [SCHEMA_VERSION, SCHEMA_VERSION], dtype=np.int32
            )
        negative_schema_extent_nml = negative_schema_extent / "run.nml"
        negative_schema_extent_nml.write_text(restart_text, encoding="utf-8")
        negative_schema_extent_result = _run_binary(
            binary, negative_schema_extent_nml, negative_schema_extent, nproc=2
        )
        (negative_schema_extent / "run.log").write_text(
            negative_schema_extent_result.stdout, encoding="utf-8"
        )
        if negative_schema_extent_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted non-scalar stellar schema attribute; output:\n"
                + negative_schema_extent_result.stdout[-12000:]
            )
        if "failed checked HDF5 read for particles/stellar_state_schema_version" not in negative_schema_extent_result.stdout:
            raise AssertionError(
                "non-scalar schema negative run did not expose its extent diagnostic"
            )

        negative_extent = run_root / "ne"
        negative_extent.mkdir()
        shutil.copytree(checkpoint, negative_extent / checkpoint.name)
        negative_extent_file = negative_extent / checkpoint.name / output_path.name
        with h5py.File(negative_extent_file, "r+") as handle:
            particles = handle["particles"]
            del particles["indtab"]
            particles.create_dataset("indtab", data=np.array([0.0], dtype=np.float64))
        negative_extent_nml = negative_extent / "run.nml"
        negative_extent_nml.write_text(restart_text, encoding="utf-8")
        negative_extent_result = _run_binary(
            binary, negative_extent_nml, negative_extent
        )
        (negative_extent / "run.log").write_text(
            negative_extent_result.stdout, encoding="utf-8"
        )
        if negative_extent_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted HDF5 release dataset with bad extent; output:\n"
                + negative_extent_result.stdout[-12000:]
            )
        if "failed checked HDF5 read for indtab" not in negative_extent_result.stdout:
            raise AssertionError(
                "bad-extent negative run did not expose the indtab-specific diagnostic"
            )

        negative_ptype = run_root / "np"
        negative_ptype.mkdir()
        shutil.copytree(checkpoint, negative_ptype / checkpoint.name)
        negative_ptype_file = negative_ptype / checkpoint.name / output_path.name
        with h5py.File(negative_ptype_file, "r+") as handle:
            particles = handle["particles"]
            del particles["ptypep"]
            particles.create_dataset("ptypep", data=np.array([1], dtype=np.int32))
        negative_ptype_nml = negative_ptype / "run.nml"
        negative_ptype_nml.write_text(restart_text, encoding="utf-8")
        negative_ptype_result = _run_binary(binary, negative_ptype_nml, negative_ptype)
        (negative_ptype / "run.log").write_text(
            negative_ptype_result.stdout, encoding="utf-8"
        )
        if negative_ptype_result.returncode == 0:
            raise AssertionError(
                "Fortran accepted HDF5 restart with truncated ptypep; output:\n"
                + negative_ptype_result.stdout[-12000:]
            )
        if "failed checked HDF5 read for ptypep" not in negative_ptype_result.stdout:
            raise AssertionError(
                "truncated-ptypep negative run did not expose the ptypep-specific diagnostic"
            )

    return {
        "writer_output": "linked HDF5 production binary",
        "reader_output": "linked HDF5 production binary",
        "schema_version": version,
        "active_particle_count": release_lengths["tpp"],
        "release_dataset_lengths": release_lengths,
        "writer_release_values": writer_release_values,
        "writer_particle_type_counts": {
            "PTYPE_DM": int(np.count_nonzero(particle_types == 0)),
            "PTYPE_STAR": int(np.count_nonzero(particle_types == 1)),
            "PTYPE_SINK": int(np.count_nonzero(particle_types == 2)),
        },
        "writer_release_state_zero": writer_release_zero,
        "writer_release_zero_match": writer_zero_match,
        "restart_seed_particle_type": "PTYPE_STAR",
        "restart_seed_release_values": {
            name: values.tolist() for name, values in restart_seed_release.items()
        },
        "restart_writer_nonzero_payload": {
            name: values.tolist() for name, values in restart_writer_release_values.items()
        },
        "restart_writer_nonzero": restart_writer_nonzero,
        "restart_writer_nonzero_payload_match": restart_writer_payload_match,
        "binary_release_values": binary_release_values,
        "binary_release_record_indices": binary_record_indices,
        "binary_nonzero_payload": binary_nonzero,
        "binary_hdf5_writer_bitwise_match": binary_value_match,
        "release_state_value_match": restart_writer_payload_match,
        "binary_release_record_order": binary_record_order,
        "restored_release_state_sample": list(restored_sample),
        "restored_release_state_value_match": restored_sample == expected_sample,
        "restart_particle_restore_marker": "HDF5 particle restore done." in restart_result.stdout,
        "cross_ncpu_reader_nproc": 4,
        "cross_ncpu_particle_restore_marker": "HDF5 particle restore done." in cross_ncpu_result.stdout,
        "cross_ncpu_local_count_markers": cross_local_count_markers,
        "cross_ncpu_rank_value_match": cross_rank_value_match,
        "cross_ncpu_zero_particle_rank_covered": cross_zero_rank_match,
        "same_ncpu_reader_nproc": 2,
        "same_ncpu_rank_value_match": same_ncpu_value_match,
        "linked_sink_payload_exercised": (
            "HDF5 particle restore done." in restart_result.stdout
            and sink_writer_roundtrip
        ),
        "linked_sink_writer_nsink": restart_writer_nsink,
        "linked_sink_writer_roundtrip": sink_writer_roundtrip,
        "fortran_missing_schema_rejected": (
            negative_schema_result.returncode != 0
            and "stellar_state_schema_version" in negative_schema_result.stdout
        ),
        "fortran_non_scalar_schema_attribute_rejected": (
            negative_schema_extent_result.returncode != 0
            and "failed checked HDF5 read for particles/stellar_state_schema_version"
            in negative_schema_extent_result.stdout
        ),
        "fortran_bad_extent_rejected": (
            negative_extent_result.returncode != 0
            and "failed checked HDF5 read for indtab" in negative_extent_result.stdout
        ),
        "fortran_truncated_ptypep_rejected": (
            negative_ptype_result.returncode != 0
            and "failed checked HDF5 read for ptypep" in negative_ptype_result.stdout
        ),
        "fortran_missing_npart_per_cpu_rejected": (
            negative_npart_result.returncode != 0
            and "failed checked HDF5 read for npart_per_cpu" in negative_npart_result.stdout
        ),
        "fortran_inconsistent_npart_per_cpu_rejected": (
            negative_count_result.returncode != 0
            and "invalid npart_per_cpu: sum=" in negative_count_result.stdout
        ),
        "fortran_negative_npart_per_cpu_rejected": (
            negative_npart_negative_result.returncode != 0
            and "invalid npart_per_cpu: negative entry" in negative_npart_negative_result.stdout
        ),
        "fortran_wrong_length_npart_per_cpu_rejected": (
            negative_npart_length_result.returncode != 0
            and "failed checked HDF5 read for npart_per_cpu" in negative_npart_length_result.stdout
        ),
        "fortran_negative_ngrid_per_cpu_rejected": (
            negative_ngrid_result.returncode != 0
            and "invalid ngrid_per_cpu: negative entry" in negative_ngrid_result.stdout
        ),
        "fortran_wrong_length_ngrid_per_cpu_rejected": (
            wrong_length_ngrid_result.returncode != 0
            and "failed checked HDF5 read for ngrid_per_cpu" in wrong_length_ngrid_result.stdout
        ),
        "fortran_truncated_amr_son_flag_rejected": (
            negative_amr_payload_result.returncode != 0
            and "checked AMR son_flag read failed" in negative_amr_payload_result.stdout
        ),
        "fortran_oversized_nsink_rejected": (
            negative_nsink_result.returncode != 0
            and "invalid HDF5 nsink=" in negative_nsink_result.stdout
        ),
        "fortran_inconsistent_nindsink_rejected": (
            negative_nindsink_result.returncode != 0
            and "nindsink is below idsink maximum" in negative_nindsink_result.stdout
        ),
        "fortran_missing_stellar_ledger_attribute_rejected": (
            negative_stellar_ledger_result.returncode != 0
            and "failed checked HDF5 read for particles/nstar_tot"
            in negative_stellar_ledger_result.stdout
        ),
        "fortran_missing_header_aexp_rejected": (
            negative_header_result.returncode != 0
            and "failed checked HDF5 read for header/aexp" in negative_header_result.stdout
        ),
        "ksection_runtime_writer_dimensions": list(ksection_dimensions),
        "ksection_runtime_restart_marker": (
            "HDF5: ksection tree restored" in ksection_restart_result.stdout
        ),
        "ksection_oversized_dimensions_rejected": (
            ksection_negative_result.returncode != 0
            and "invalid HDF5 ksection dimensions" in ksection_negative_result.stdout
        ),
        "same_ncpu_metadata_injected": True,
        "injected_restart_payload": True,
        "true_ptype_star_continuation_equivalence": False,
        "binary_build": _binary_provenance(binary),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fortran-runtime",
        action="store_true",
        help="also run the bounded one-rank linked HDF5 writer/reader fixture",
    )
    args = parser.parse_args()

    backup_text = BACKUP.read_text(encoding="utf-8")
    restore_text = RESTORE.read_text(encoding="utf-8")
    source_checks = {
        "backup_schema_marker": "stellar_state_schema_version" in backup_text,
        "backup_tpp": "'tpp'" in backup_text,
        "backup_mp0": "'mp0'" in backup_text,
        "backup_indtab": "'indtab'" in backup_text,
        "restore_schema_guard": "h5aexists_f" in restore_text
        and "unsupported HDF5 stellar state schema" in restore_text,
        "restore_checked_schema_attribute": "hdf5_read_attr_int_checked" in restore_text
        and "stellar_state_schema_version" in restore_text
        and "failed checked HDF5 read for particles/stellar_state_schema_version" in restore_text,
        "restore_tpp": "hdf5_read_dataset_1d_dp_checked(grp_id, 'tpp'" in restore_text,
        "restore_checked_release_reads": "hdf5_read_dataset_1d_dp_checked" in restore_text
        and all(f"'{name}'" in restore_text for name in RELEASE_FIELDS),
        "restore_checked_particle_reads": all(
            name in restore_text
            for name in (
                "hdf5_read_dataset_1d_dp_checked(grp_id, 'x_",
                "hdf5_read_dataset_1d_dp_checked(grp_id, 'v_",
                "hdf5_read_dataset_1d_dp_checked(grp_id, 'mass'",
                "hdf5_read_dataset_1d_int8_checked(grp_id, 'identity'",
                "hdf5_read_dataset_1d_int_checked(grp_id, 'levelp'",
                "hdf5_read_dataset_1d_int_checked(grp_id, 'ptypep'",
            )
        ),
        "restore_checked_npart_per_cpu": "hdf5_read_dataset_all_int_checked" in restore_text
        and "failed checked HDF5 read for npart_per_cpu" in restore_text
        and "invalid npart_per_cpu: sum=" in restore_text
        and "invalid npart_per_cpu: negative entry" in restore_text,
        "restore_checked_ngrid_per_cpu": "failed checked HDF5 read for ngrid_per_cpu" in restore_text
        and "invalid ngrid_per_cpu: negative entry" in restore_text,
        "restore_checked_sink_count": "invalid HDF5 nsink=" in restore_text
        and "nindsink is below idsink maximum" in restore_text
        and "nsink > nsinkmax" in restore_text,
        "restore_checked_sink_payload": "hdf5_restore_sink_dp_checked" in restore_text
        and "hdf5_restore_sink_int_checked" in restore_text
        and "failed checked HDF5 read for sink dataset" in restore_text,
        "restore_checked_adm_payload": "dark_energy_int" in restore_text
        and "dark_h2_frac" in restore_text
        and "hdf5_read_dataset_1d_dp_checked" in restore_text,
        "restore_checked_stellar_ledger": "hdf5_read_attr_int8_checked" in restore_text
        and "hdf5_read_attr_dp_checked" in restore_text
        and "failed checked HDF5 read for particles/nstar_tot" in restore_text
        and "failed checked HDF5 read for particles/mstar_tot" in restore_text
        and "failed checked HDF5 read for particles/mstar_lost" in restore_text,
        "restore_checked_header_clock": "hdf5_restore_header_dp_checked" in restore_text
        and "hdf5_restore_header_int_checked" in restore_text
        and "hdf5_restore_header_dp_array_checked" in restore_text
        and "hdf5_restore_header_string_checked" in restore_text
        and "hdf5_restore_header_dp_checked(grp_id, 'aexp', aexp)" in restore_text,
        "restore_checked_amr_payload": "hdf5_read_dataset_collective_dp_checked" in restore_text
        and "hdf5_read_dataset_collective_int_checked" in restore_text
        and "checked AMR son_flag read failed" in restore_text
        and "hdf5_read_dataset_all_int_checked" in restore_text,
        "restore_checked_ksection": "hdf5_read_dataset_all_dp_checked" in restore_text
        and "invalid HDF5 ksection dimensions" in restore_text
        and "failed checked HDF5 read for ksec_wall" in restore_text,
        "restore_restart_abort": "hdf5_restart_abort" in restore_text,
        "restore_local_count_diagnostic": "HDF5 stellar release local count rank=" in restore_text,
    }
    if not all(source_checks.values()):
        raise AssertionError(f"incomplete production HDF5 source contract: {source_checks}")

    build_root = ROOT / "build"
    build_root.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="p0_hdf5_stellar_restart_", dir=build_root) as temporary_directory:
        root = Path(temporary_directory)
        valid = root / "valid.h5"
        legacy = root / "legacy_missing_release_state.h5"
        _write_synthetic_checkpoint(valid, include_release_state=True)
        _write_synthetic_checkpoint(legacy, include_release_state=False)
        valid_state = _read_and_validate(valid)
        expected_state = {
            name: _synthetic_arrays()[name].tobytes() for name in RELEASE_FIELDS
        }
        if valid_state != expected_state:
            raise AssertionError("synthetic HDF5 state was not bitwise preserved")
        try:
            _read_and_validate(legacy)
        except RuntimeError as error:
            legacy_rejected = "missing stellar release fields" in str(error)
        else:
            legacy_rejected = False
        if not legacy_rejected:
            raise AssertionError("legacy/incomplete HDF5 state was not rejected")

    evidence = {
        "schema": "lagramses-hdf5-stellar-restart-state",
        "schema_version": SCHEMA_VERSION,
        "release_fields": list(RELEASE_FIELDS),
        "bitwise_synthetic_roundtrip": valid_state == expected_state,
        "legacy_missing_state_rejected": legacy_rejected,
        "production_fortran_source_guard": all(source_checks.values()),
        "production_fortran_runtime_roundtrip": False,
        "scope": "bounded synthetic schema contract; no RAMSES output opened",
    }
    if args.fortran_runtime:
        evidence["production_fortran_runtime_roundtrip"] = _run_fortran_hdf5_roundtrip(ROOT)
        evidence["scope"] = (
            "bounded one-rank HDF5/binary writer plus one-rank and four-rank "
            "HDF5 reader fixtures, two-rank same-ncpu and negative fixtures; "
            "post-restart writer uses a temporary "
            "PTYPE_STAR seed; no production data used; fixture runtime outputs "
            "were opened"
        )
    if not evidence["bitwise_synthetic_roundtrip"]:
        raise AssertionError("unexpected synthetic payload size")
    # The default synthetic check is intentionally non-mutating.  Only the
    # explicit linked-runtime run may update durable provenance, preventing a
    # quick smoke test from overwriting a passing runtime artifact.
    if args.fortran_runtime:
        EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    runtime_status = bool(args.fortran_runtime)
    print(
        "P0_HDF5_STELLAR_RESTART_CONTRACT_OK "
        f"bitwise_synthetic=True legacy_fail_closed=True "
        f"fortran_runtime_roundtrip={str(runtime_status).lower()}"
    )


if __name__ == "__main__":
    main()
