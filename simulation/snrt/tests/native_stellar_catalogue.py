"""Synthetic native cuRAMSES stellar-catalogue reader test."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import struct
import sys
from tempfile import TemporaryDirectory

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
READER_PATH = PROJECT_ROOT / "tools" / "read_native_stellar_catalogue.py"


def _load_reader():
    spec = importlib.util.spec_from_file_location("read_native_stellar_catalogue", READER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("unable to load native stellar reader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _record(payload: bytes) -> bytes:
    marker = struct.pack("<i", len(payload))
    return marker + payload + marker


def _part_file(path: Path, rank: int, ptype: np.ndarray, source_ids: np.ndarray, current_proper: float) -> None:
    npart = len(ptype)
    star_mask = ptype == 1
    positions = np.arange(npart, dtype=np.float64)[:, None] * 0.1 + rank * 0.01
    position_x = positions[:, 0]
    position_y = position_x + 0.1
    position_z = position_x + 0.2
    velocity_x = np.full(npart, 0.01)
    velocity_y = np.full(npart, -0.02)
    velocity_z = np.full(npart, 0.03)
    mass = np.where(star_mask, 2.0e-4, 1.0e-3)
    level = np.full(npart, 9, dtype=np.int32)
    birth_epoch = np.where(star_mask, -2.0 - np.arange(npart), 0.0)
    metallicity = np.where(star_mask, 0.01 + np.arange(npart) * 0.001, 0.0)
    birth_proper = np.where(star_mask, current_proper - 0.05 - np.arange(npart) * 0.01, 0.0)
    initial_mass = np.where(star_mask, 2.5e-4, 0.0)
    yield_index = np.where(star_mask, 0.02, 0.0)
    fields = (
        position_x,
        position_y,
        position_z,
        velocity_x,
        velocity_y,
        velocity_z,
        mass,
        source_ids.astype(np.int64),
        level,
        ptype.astype(np.int8),
        np.zeros(npart, dtype=np.float64),
        birth_epoch,
        metallicity,
        birth_proper,
        initial_mass,
        yield_index,
    )
    dtypes = ("<f8", "<f8", "<f8", "<f8", "<f8", "<f8", "<f8", "<i8", "<i4", "<i1", "<f8", "<f8", "<f8", "<f8", "<f8", "<f8")
    headers = (
        _record(struct.pack("<i", 2)),
        _record(struct.pack("<i", 3)),
        _record(struct.pack("<i", npart)),
        _record(bytes(16)),
        _record(struct.pack("<q", 3)),
        _record(struct.pack("<d", 0.0)),
        _record(struct.pack("<d", 0.0)),
        _record(struct.pack("<i", 0)),
    )
    with path.open("wb") as handle:
        for payload in headers:
            handle.write(payload)
        for values, dtype in zip(fields, dtypes, strict=True):
            handle.write(_record(np.asarray(values, dtype=dtype).tobytes()))


def main() -> None:
    reader = _load_reader()
    with TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        output_dir = root / "output_00001"
        output_dir.mkdir()
        info = output_dir / "info_00001.txt"
        info.write_text(
            "aexp = 0.5\n"
            "time = -1.0\n"
            "H0 = 68.0\n"
            "omega_m = 0.3\n"
            "omega_l = 0.7\n"
            "omega_k = 0.0\n"
            "unit_l = 1.0e24\n"
            "unit_d = 1.0e-27\n"
            "unit_t = 1.0e16\n"
        )
        info_values = reader._read_info(info)
        current_proper = reader._proper_time_code_from_aexp(info_values, -1.0, 0.0)
        _part_file(
            output_dir / "part_00001.out00001",
            1,
            np.array([0, 1, 1], dtype=np.int8),
            np.array([11, 12, 13], dtype=np.int64),
            current_proper,
        )
        _part_file(
            output_dir / "part_00001.out00002",
            2,
            np.array([1, 0], dtype=np.int8),
            np.array([14, 15], dtype=np.int64),
            current_proper,
        )
        rows, manifest = reader.extract_catalogue(output_dir, info, include_velocity=True)
        assert rows["source_id"].tolist() == [12, 13, 14]
        assert np.all(rows["age_myr"] > 0.0)
        assert manifest["star_count"] == 3
        assert manifest["header"]["nstar_tot_record_bytes"] == 8
        assert manifest["source_catalogue"]["photon_luminosity_assigned"] is False
        output = root / "stars.csv"
        metadata = root / "stars.json"
        reader.write_catalogue(output, metadata, rows, manifest, include_velocity=True)
        payload = json.loads(metadata.read_text())
        assert payload["catalogue_csv_sha256"] == __import__("hashlib").sha256(output.read_bytes()).hexdigest()
        assert len(output.read_text().splitlines()) == 4
    print("NATIVE_STELLAR_CATALOGUE_TEST_OK stars=3 velocity=1")


if __name__ == "__main__":
    main()
