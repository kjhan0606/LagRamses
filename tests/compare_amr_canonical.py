#!/usr/bin/env python3
"""Capacity- and decomposition-independent comparator for legacy RAMSES AMR.

Legacy AMR files encode cell references using ``ngridmax`` as a plane stride.
Consequently, byte comparison fails after a harmless capacity change.  This
tool decodes each file with its own stride and re-keys grids/cells by exact
dyadic geometry before comparing topology.

The implementation uses only the Python standard library.  SQLite holds the
global owner-grid index so memory use is bounded when an output contains
millions of grids.
"""

from __future__ import annotations

import argparse
import dataclasses
import decimal
import hashlib
import json
import math
import os
import pathlib
import re
import sqlite3
import struct
import sys
import tempfile
from collections.abc import Iterable, Sequence


SCHEMA_VERSION = "lagRamses-amr-canonical-v1"
EXCLUDED_FIELDS = (
    "ngridmax, ngrid_current, ncpu across compared outputs",
    "headl/taill/numbl, headb/tailb/numbb, free-list and used-memory state",
    "numbtot rank min/max/average/reserved columns (only global totals retained)",
    "ordering/domain decomposition records and all cpu_map values",
    "raw grid ids and next/prev linked-list ids",
    "info ncpu/ngridmax/ordering/DOMAIN rows",
    "physical boundary owner topology (boundary maps are reference-only)",
)


class AmrFormatError(RuntimeError):
    """Malformed or unsupported AMR input (CLI exit status 2)."""


class CanonicalMismatch(RuntimeError):
    """Well-formed inputs have different semantics (CLI exit status 1)."""


GridKey = tuple[int, int, int, int]
CellKey = tuple[int, int, int, int] | None


def _float_bits(value: float) -> str:
    if not math.isfinite(value):
        raise AmrFormatError(f"non-finite real(dp) value {value!r}")
    return struct.pack(">d", value).hex()


def _float_tuple(values: Sequence[float]) -> tuple[str, ...]:
    return tuple(_float_bits(value) for value in values)


@dataclasses.dataclass(frozen=True)
class SemanticHeader:
    ndim: int
    coarse_shape: tuple[int, int, int]
    nlevelmax: int
    nboundary: int
    boxlen: str
    output_control: tuple[int, int, int]
    tout: tuple[str, ...]
    aout: tuple[str, ...]
    time: str
    dtold: tuple[str, ...]
    dtnew: tuple[str, ...]
    steps: tuple[int, int]
    mass_state: tuple[str, ...]
    cosmology: tuple[str, ...]
    expansion_energy: tuple[str, ...]
    mass_sph: str
    level_grid_totals: tuple[int, ...]

    def digest(self) -> str:
        return hashlib.sha256(_json_bytes(dataclasses.asdict(self))).hexdigest()


@dataclasses.dataclass
class FileLayout:
    ncpu: int
    ndim: int
    coarse_shape: tuple[int, int, int]
    nlevelmax: int
    ngridmax: int
    nboundary: int
    numbl: tuple[int, ...]
    numbb: tuple[int, ...]
    coarse_sons: tuple[int, ...]
    coarse_flags: tuple[int, ...]
    header: SemanticHeader

    @property
    def ncoarse(self) -> int:
        return math.prod(self.coarse_shape)


@dataclasses.dataclass(frozen=True)
class FieldDigests:
    father: bytes
    neighbours: bytes
    sons: bytes
    flags: bytes

    @property
    def combined(self) -> bytes:
        digest = hashlib.sha256()
        digest.update(b"grid-fields-v1\0")
        digest.update(self.father)
        digest.update(self.neighbours)
        digest.update(self.sons)
        digest.update(self.flags)
        return digest.digest()


@dataclasses.dataclass
class OutputIndex:
    output_dir: pathlib.Path
    db_path: pathlib.Path
    header: SemanticHeader
    info: tuple[tuple[str, str], ...]
    ncpu: int
    level_counts: dict[int, int]
    coarse_digest: str
    topology_digest: str
    local_layout_digests: dict[int, str]
    input_files: list[dict[str, object]]

    def report(self) -> dict[str, object]:
        return {
            "path": str(self.output_dir.resolve()),
            "ncpu": self.ncpu,
            "header_sha256": self.header.digest(),
            "info": dict(self.info),
            "level_counts": self.level_counts,
            "coarse_sha256": self.coarse_digest,
            "topology_sha256": self.topology_digest,
            "local_layout_sha256": self.local_layout_digests,
            "inputs": self.input_files,
        }


def _json_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


class FortranSequentialReader:
    """Fail-closed reader for 4-byte legacy Fortran sequential records."""

    def __init__(self, path: pathlib.Path, *, calculate_sha256: bool = False):
        self.path = path
        self._stream = path.open("rb")
        self.endian = self._detect_endian()
        self._marker = struct.Struct(self.endian + "i")
        self._hasher = hashlib.sha256() if calculate_sha256 else None
        self.record_number = 0

    def _detect_endian(self) -> str:
        prefix = self._stream.read(12)
        self._stream.seek(0)
        for endian in ("<", ">"):
            if len(prefix) < 12:
                continue
            if struct.unpack(endian + "i", prefix[:4])[0] != 4:
                continue
            if struct.unpack(endian + "i", prefix[8:12])[0] == 4:
                return endian
        raise AmrFormatError(
            f"{self.path}: expected 4-byte Fortran markers and first payload length 4"
        )

    def __enter__(self) -> "FortranSequentialReader":
        return self

    def __exit__(self, *_: object) -> None:
        self._stream.close()

    @property
    def sha256(self) -> str:
        if self._hasher is None:
            raise RuntimeError("SHA256 was not requested")
        return self._hasher.hexdigest()

    def record(self) -> bytes:
        leading = self._stream.read(4)
        if not leading:
            raise AmrFormatError(
                f"{self.path}: unexpected EOF before record {self.record_number}"
            )
        if len(leading) != 4:
            raise AmrFormatError(f"{self.path}: truncated leading marker")
        length = self._marker.unpack(leading)[0]
        if length < 0 or length > 2**31 - 1:
            raise AmrFormatError(
                f"{self.path}: invalid length {length} at record {self.record_number}"
            )
        payload = self._stream.read(length)
        trailing = self._stream.read(4)
        if len(payload) != length or len(trailing) != 4:
            raise AmrFormatError(
                f"{self.path}: truncated record {self.record_number}"
            )
        if self._marker.unpack(trailing)[0] != length:
            raise AmrFormatError(
                f"{self.path}: marker mismatch at record {self.record_number}"
            )
        if self._hasher is not None:
            self._hasher.update(leading)
            self._hasher.update(payload)
            self._hasher.update(trailing)
        self.record_number += 1
        return payload

    def ints(self, count: int | None = None) -> tuple[int, ...]:
        payload = self.record()
        if len(payload) % 4:
            raise AmrFormatError(
                f"{self.path}: integer record length {len(payload)} is not divisible by 4"
            )
        values = tuple(v[0] for v in struct.iter_unpack(self.endian + "i", payload))
        if count is not None and len(values) != count:
            raise AmrFormatError(
                f"{self.path}: expected {count} integers, got {len(values)}"
            )
        return values

    def doubles(self, count: int | None = None) -> tuple[float, ...]:
        payload = self.record()
        if len(payload) % 8:
            raise AmrFormatError(
                f"{self.path}: real(dp) record length {len(payload)} is not divisible by 8"
            )
        values = tuple(v[0] for v in struct.iter_unpack(self.endian + "d", payload))
        if count is not None and len(values) != count:
            raise AmrFormatError(
                f"{self.path}: expected {count} reals, got {len(values)}"
            )
        if not all(math.isfinite(value) for value in values):
            raise AmrFormatError(f"{self.path}: non-finite real(dp) payload")
        return values

    def int64s(self, count: int | None = None) -> tuple[int, ...]:
        payload = self.record()
        if len(payload) % 8:
            raise AmrFormatError(
                f"{self.path}: integer(i8b) record length {len(payload)} is not "
                "divisible by 8"
            )
        values = tuple(v[0] for v in struct.iter_unpack(self.endian + "q", payload))
        if count is not None and len(values) != count:
            raise AmrFormatError(
                f"{self.path}: expected {count} integer(i8b), got {len(values)}"
            )
        return values

    def expect_eof(self) -> None:
        if self._stream.read(1):
            raise AmrFormatError(
                f"{self.path}: trailing data after record {self.record_number}"
            )


def _one(values: Sequence[int] | Sequence[float], label: str) -> int | float:
    if len(values) != 1:
        raise AmrFormatError(f"{label}: expected one value, got {len(values)}")
    return values[0]


def _read_layout(reader: FortranSequentialReader) -> FileLayout:
    ncpu = int(_one(reader.ints(), "ncpu"))
    ndim = int(_one(reader.ints(), "ndim"))
    coarse_shape = tuple(reader.ints(3))
    nlevelmax = int(_one(reader.ints(), "nlevelmax"))
    ngridmax = int(_one(reader.ints(), "ngridmax"))
    nboundary = int(_one(reader.ints(), "nboundary"))
    _ngrid_current = int(_one(reader.ints(), "ngrid_current"))
    boxlen_values = reader.doubles(1)

    if ncpu <= 0 or ndim not in (1, 2, 3) or nlevelmax <= 0:
        raise AmrFormatError(
            f"{reader.path}: invalid ncpu/ndim/nlevelmax: {ncpu}/{ndim}/{nlevelmax}"
        )
    if ngridmax <= 0 or nboundary < 0 or any(value <= 0 for value in coarse_shape):
        raise AmrFormatError(f"{reader.path}: invalid capacity/coarse domain")

    output_control = tuple(reader.ints(3))
    noutput = output_control[0]
    if noutput < 0:
        raise AmrFormatError(f"{reader.path}: negative noutput")
    tout_values = reader.doubles(noutput)
    aout_values = reader.doubles(noutput)
    time_values = reader.doubles(1)
    dtold_values = reader.doubles(nlevelmax)
    dtnew_values = reader.doubles(nlevelmax)
    steps = tuple(reader.ints(2))
    mass_state_values = reader.doubles(3)
    cosmology_values = reader.doubles(7)
    expansion_values = reader.doubles(5)
    mass_sph_values = reader.doubles(1)

    matrix_size = ncpu * nlevelmax
    reader.ints(matrix_size)  # linked-list heads
    reader.ints(matrix_size)  # linked-list tails
    numbl = reader.ints(matrix_size)
    numbtot = tuple(reader.int64s(10 * nlevelmax))
    level_grid_totals = tuple(numbtot[10 * level] for level in range(nlevelmax))

    numbb: tuple[int, ...] = ()
    if nboundary:
        boundary_size = nboundary * nlevelmax
        reader.ints(boundary_size)
        reader.ints(boundary_size)
        numbb = reader.ints(boundary_size)
    reader.ints(5)  # headf, tailf, numbf, used_mem, used_mem_tot

    ordering = reader.record().decode("ascii", errors="strict").strip()
    if ordering == "bisection":
        ordering_records = 5
    elif ordering == "ksection":
        ordering_records = 10
    else:
        ordering_records = 1
    for _ in range(ordering_records):
        reader.record()

    ncoarse = math.prod(coarse_shape)
    coarse_sons = tuple(reader.ints(ncoarse))
    coarse_flags = tuple(reader.ints(ncoarse))
    reader.ints(ncoarse)  # coarse cpu_map

    if any(count < 0 for count in (*numbl, *numbb)):
        raise AmrFormatError(f"{reader.path}: negative per-block grid count")

    header = SemanticHeader(
        ndim=ndim,
        coarse_shape=coarse_shape,
        nlevelmax=nlevelmax,
        nboundary=nboundary,
        boxlen=_float_bits(boxlen_values[0]),
        output_control=output_control,
        tout=_float_tuple(tout_values),
        aout=_float_tuple(aout_values),
        time=_float_bits(time_values[0]),
        dtold=_float_tuple(dtold_values),
        dtnew=_float_tuple(dtnew_values),
        steps=steps,
        mass_state=_float_tuple(mass_state_values),
        cosmology=_float_tuple(cosmology_values),
        expansion_energy=_float_tuple(expansion_values),
        mass_sph=_float_bits(mass_sph_values[0]),
        level_grid_totals=level_grid_totals,
    )
    return FileLayout(
        ncpu=ncpu,
        ndim=ndim,
        coarse_shape=coarse_shape,
        nlevelmax=nlevelmax,
        ngridmax=ngridmax,
        nboundary=nboundary,
        numbl=numbl,
        numbb=numbb,
        coarse_sons=coarse_sons,
        coarse_flags=coarse_flags,
        header=header,
    )


def _block_count(layout: FileLayout, level: int, owner: int) -> int:
    if owner <= layout.ncpu:
        return layout.numbl[(level - 1) * layout.ncpu + owner - 1]
    boundary = owner - layout.ncpu
    return layout.numbb[(level - 1) * layout.nboundary + boundary - 1]


def _canonical_grid_key(
    level: int,
    coordinates: Sequence[float],
    layout: FileLayout,
    path: pathlib.Path,
) -> GridKey:
    canonical = [level, 0, 0, 0]
    for axis in range(layout.ndim):
        numerator, denominator = coordinates[axis].as_integer_ratio()
        numerator <<= level
        quotient, remainder = divmod(numerator, denominator)
        if remainder:
            raise AmrFormatError(
                f"{path}: level {level} xg={coordinates[axis]:.17g} is not "
                "an exact dyadic grid centre"
            )
        if quotient % 2 != 1:
            raise AmrFormatError(
                f"{path}: level {level} grid coordinate {quotient} is not odd"
            )
        period = layout.coarse_shape[axis] * (1 << level)
        canonical[axis + 1] = quotient % period
    return tuple(canonical)  # type: ignore[return-value]


def _skip_topology_records(reader: FortranSequentialReader, layout: FileLayout) -> None:
    count = 1 + 2 * layout.ndim + 3 * (2**layout.ndim)
    for _ in range(count):
        reader.record()


def _build_local_map(
    path: pathlib.Path,
) -> tuple[FileLayout, dict[int, GridKey], str, str]:
    local_map: dict[int, GridKey] = {}
    layout_digest = hashlib.sha256(b"local-block-order-v1\0")
    with FortranSequentialReader(path, calculate_sha256=True) as reader:
        layout = _read_layout(reader)
        layout_digest.update(
            _pack_i64((layout.ncpu, layout.nboundary, layout.nlevelmax))
        )
        for level in range(1, layout.nlevelmax + 1):
            for owner in range(1, layout.ncpu + layout.nboundary + 1):
                count = _block_count(layout, level, owner)
                layout_digest.update(_pack_i64((level, owner, count)))
                if not count:
                    continue
                indices = reader.ints(count)
                reader.ints(count)
                reader.ints(count)
                coordinates = [reader.doubles(count) for _ in range(layout.ndim)]
                for offset, local_id in enumerate(indices):
                    if local_id <= 0 or local_id > layout.ngridmax:
                        raise AmrFormatError(
                            f"{path}: local grid {local_id} outside 1..{layout.ngridmax}"
                        )
                    key = _canonical_grid_key(
                        level,
                        tuple(axis[offset] for axis in coordinates),
                        layout,
                        path,
                    )
                    previous = local_map.setdefault(local_id, key)
                    if previous != key:
                        raise AmrFormatError(
                            f"{path}: local grid {local_id} has conflicting geometry"
                        )
                    # Gravity payloads use this exact level/owner/linked-list order.
                    # Canonical keys make the proof independent of local grid IDs.
                    layout_digest.update(_serialize_grid(key))
                _skip_topology_records(reader, layout)
        reader.expect_eof()
        file_sha256 = reader.sha256
    return layout, local_map, file_sha256, layout_digest.hexdigest()


def _child_bits(child0: int, ndim: int) -> tuple[int, int, int]:
    return tuple((child0 >> axis) & 1 if axis < ndim else 0 for axis in range(3))


def _coarse_cell(index: int, shape: tuple[int, int, int]) -> CellKey:
    if index <= 0 or index > math.prod(shape):
        raise AmrFormatError(f"coarse cell index {index} is out of range")
    linear = index - 1
    x = linear % shape[0]
    linear //= shape[0]
    y = linear % shape[1]
    z = linear // shape[1]
    return 0, x, y, z


def _decode_cell(
    value: int,
    layout: FileLayout,
    local_map: dict[int, GridKey],
    path: pathlib.Path,
) -> CellKey:
    if value == 0:
        return None
    if value < 0:
        raise AmrFormatError(f"{path}: negative cell reference {value}")
    if value <= layout.ncoarse:
        return _coarse_cell(value, layout.coarse_shape)
    quotient = value - layout.ncoarse - 1
    child0, grid0 = divmod(quotient, layout.ngridmax)
    children = 2**layout.ndim
    if child0 < 0 or child0 >= children:
        raise AmrFormatError(
            f"{path}: legacy cell {value} has child plane {child0} outside 0..{children-1}"
        )
    grid_id = grid0 + 1
    parent = local_map.get(grid_id)
    if parent is None:
        raise AmrFormatError(
            f"{path}: legacy cell {value} references absent grid {grid_id}"
        )
    bits = _child_bits(child0, layout.ndim)
    result = [parent[0], 0, 0, 0]
    for axis in range(layout.ndim):
        period = layout.coarse_shape[axis] * (1 << parent[0])
        result[axis + 1] = (parent[axis + 1] + bits[axis] - 1) % period
    return tuple(result)  # type: ignore[return-value]


def _decode_grid(
    value: int, local_map: dict[int, GridKey], path: pathlib.Path
) -> GridKey | None:
    if value == 0:
        return None
    if value < 0:
        raise AmrFormatError(f"{path}: negative grid reference {value}")
    result = local_map.get(value)
    if result is None:
        raise AmrFormatError(f"{path}: reference names absent local grid {value}")
    return result


def _expected_father(grid: GridKey, layout: FileLayout) -> CellKey:
    if grid[0] < 1:
        raise AmrFormatError(f"invalid grid level in key {grid}")
    return (grid[0] - 1, grid[1] // 2, grid[2] // 2, grid[3] // 2)


def _expected_son(
    grid: GridKey, child0: int, layout: FileLayout
) -> GridKey:
    bits = _child_bits(child0, layout.ndim)
    result = [grid[0] + 1, 0, 0, 0]
    for axis in range(layout.ndim):
        cell_period = layout.coarse_shape[axis] * (1 << grid[0])
        cell_index = (grid[axis + 1] + bits[axis] - 1) % cell_period
        result[axis + 1] = 2 * cell_index + 1
    return tuple(result)  # type: ignore[return-value]


def _pack_i64(values: Sequence[int]) -> bytes:
    try:
        return struct.pack(">" + "q" * len(values), *values)
    except struct.error as error:
        raise AmrFormatError(f"canonical integer exceeds signed int64: {values}") from error


def _serialize_cell(cell: CellKey) -> bytes:
    return _pack_i64((-1, 0, 0, 0) if cell is None else cell)


def _serialize_grid(grid: GridKey | None) -> bytes:
    return _pack_i64((-1, 0, 0, 0) if grid is None else grid)


def _field_digests(
    key: GridKey,
    father_value: int,
    neighbour_values: Sequence[int],
    son_values: Sequence[int],
    flag_values: Sequence[int],
    layout: FileLayout,
    local_map: dict[int, GridKey],
    path: pathlib.Path,
) -> FieldDigests:
    father = _decode_cell(father_value, layout, local_map, path)
    expected_father = _expected_father(key, layout)
    if father != expected_father:
        raise AmrFormatError(
            f"{path}: grid {key} father {father} disagrees with geometry "
            f"{expected_father}"
        )
    neighbours = tuple(
        _decode_cell(value, layout, local_map, path) for value in neighbour_values
    )
    sons = tuple(_decode_grid(value, local_map, path) for value in son_values)
    for child0, son in enumerate(sons):
        if son is not None:
            expected = _expected_son(key, child0, layout)
            if son != expected:
                raise AmrFormatError(
                    f"{path}: grid {key} child {child0} son {son} disagrees "
                    f"with geometry {expected}"
                )

    father_digest = hashlib.sha256(b"father\0" + _serialize_cell(father)).digest()
    neighbour_digest = hashlib.sha256(
        b"neighbours\0" + b"".join(_serialize_cell(cell) for cell in neighbours)
    ).digest()
    son_digest = hashlib.sha256(
        b"sons\0" + b"".join(_serialize_grid(son) for son in sons)
    ).digest()
    flag_digest = hashlib.sha256(
        b"flags\0" + _pack_i64(tuple(flag_values))
    ).digest()
    return FieldDigests(father_digest, neighbour_digest, son_digest, flag_digest)


def _canonical_coarse(
    layout: FileLayout, local_map: dict[int, GridKey], path: pathlib.Path
) -> tuple[tuple[GridKey, GridKey | None, int], ...]:
    rows: list[tuple[GridKey, GridKey | None, int]] = []
    for index, (son_value, flag) in enumerate(
        zip(layout.coarse_sons, layout.coarse_flags, strict=True), start=1
    ):
        cell = _coarse_cell(index, layout.coarse_shape)
        assert cell is not None
        son = _decode_grid(son_value, local_map, path)
        if son is not None:
            expected_values = [1, 0, 0, 0]
            for axis in range(layout.ndim):
                expected_values[axis + 1] = 2 * cell[axis + 1] + 1
            expected = tuple(expected_values)
            if son != expected:
                raise AmrFormatError(
                    f"{path}: coarse cell {cell} son {son} disagrees with {expected}"
                )
        rows.append((cell, son, flag))
    return tuple(rows)


def _coarse_digest(rows: Sequence[tuple[GridKey, GridKey | None, int]]) -> str:
    digest = hashlib.sha256(b"coarse-topology-v1\0")
    for cell, son, flag in rows:
        digest.update(_serialize_grid(cell))
        digest.update(_serialize_grid(son))
        digest.update(_pack_i64((flag,)))
    return digest.hexdigest()


def _create_database(path: pathlib.Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=FILE;
        CREATE TABLE grids (
          level INTEGER NOT NULL,
          gx INTEGER NOT NULL,
          gy INTEGER NOT NULL,
          gz INTEGER NOT NULL,
          digest BLOB NOT NULL,
          father BLOB NOT NULL,
          neighbours BLOB NOT NULL,
          sons BLOB NOT NULL,
          flags BLOB NOT NULL,
          src_cpu INTEGER NOT NULL,
          src_grid INTEGER NOT NULL,
          PRIMARY KEY(level,gx,gy,gz)
        ) WITHOUT ROWID;
        """
    )
    return connection


def _emit_owner_rows(
    path: pathlib.Path,
    source_cpu: int,
    expected_layout: FileLayout,
    local_map: dict[int, GridKey],
    connection: sqlite3.Connection,
) -> tuple[int, dict[int, int], str]:
    inserted = 0
    inserted_by_level: dict[int, int] = {}
    with FortranSequentialReader(path, calculate_sha256=True) as reader:
        layout = _read_layout(reader)
        if layout != expected_layout:
            raise AmrFormatError(f"{path}: header/layout changed between parser passes")
        with connection:
            for level in range(1, layout.nlevelmax + 1):
                for owner in range(1, layout.ncpu + layout.nboundary + 1):
                    count = _block_count(layout, level, owner)
                    if not count:
                        continue
                    indices = reader.ints(count)
                    reader.ints(count)
                    reader.ints(count)
                    coordinates = [reader.doubles(count) for _ in range(layout.ndim)]
                    fathers = reader.ints(count)
                    neighbours = [reader.ints(count) for _ in range(2 * layout.ndim)]
                    sons = [reader.ints(count) for _ in range(2**layout.ndim)]
                    for _ in range(2**layout.ndim):
                        reader.ints(count)  # cpu_map
                    flags = [reader.ints(count) for _ in range(2**layout.ndim)]
                    if owner != source_cpu:
                        continue
                    rows = []
                    for offset, local_id in enumerate(indices):
                        key = local_map.get(local_id)
                        if key is None:
                            raise AmrFormatError(
                                f"{path}: owner grid {local_id} is absent from pass-1 map"
                            )
                        repeated_key = _canonical_grid_key(
                            level,
                            tuple(axis[offset] for axis in coordinates),
                            layout,
                            path,
                        )
                        if repeated_key != key:
                            raise AmrFormatError(
                                f"{path}: grid {local_id} changed between parser passes"
                            )
                        fields = _field_digests(
                            key,
                            fathers[offset],
                            tuple(axis[offset] for axis in neighbours),
                            tuple(axis[offset] for axis in sons),
                            tuple(axis[offset] for axis in flags),
                            layout,
                            local_map,
                            path,
                        )
                        rows.append(
                            (*key, fields.combined, fields.father, fields.neighbours,
                             fields.sons, fields.flags, source_cpu, local_id)
                        )
                    try:
                        connection.executemany(
                            "INSERT INTO grids VALUES (?,?,?,?,?,?,?,?,?,?,?)", rows
                        )
                    except sqlite3.IntegrityError as error:
                        raise AmrFormatError(
                            f"{path}: duplicate canonical owner-grid key"
                        ) from error
                    inserted += len(rows)
                    inserted_by_level[level] = inserted_by_level.get(level, 0) + len(rows)
        reader.expect_eof()
        second_pass_sha256 = reader.sha256
    return inserted, inserted_by_level, second_pass_sha256


_AMR_NAME = re.compile(r"^amr_([0-9]{5})\.out([0-9]{5})$")
_INFO_ASSIGNMENT = re.compile(r"^\s*([A-Za-z0-9_]+)\s*=\s*(.*?)\s*$")
_INFO_FIELDS = (
    "ndim", "levelmin", "levelmax", "nstep_coarse", "boxlen", "time",
    "aexp", "h0", "omega_m", "omega_l", "omega_k", "omega_b", "unit_l",
    "unit_d", "unit_t",
)


def _source_cpu(path: pathlib.Path) -> int:
    match = _AMR_NAME.match(path.name)
    if not match:
        raise AmrFormatError(f"{path}: cannot identify source CPU suffix")
    return int(match.group(2))


def _output_stem(path: pathlib.Path) -> str:
    match = _AMR_NAME.match(path.name)
    if not match:
        raise AmrFormatError(f"{path}: invalid AMR output filename")
    return match.group(1)


def _amr_files(output_dir: pathlib.Path) -> list[pathlib.Path]:
    files = sorted(output_dir.glob("amr_*.out[0-9][0-9][0-9][0-9][0-9]"))
    if not files:
        raise AmrFormatError(f"{output_dir}: no amr_*.outNNNNN files")
    return files


def _info_file(output_dir: pathlib.Path, output_stem: str) -> pathlib.Path:
    files = sorted(output_dir.glob("info_*.txt"))
    if len(files) != 1:
        raise AmrFormatError(
            f"{output_dir}: expected one info_*.txt, found {len(files)}"
        )
    expected_name = f"info_{output_stem}.txt"
    if files[0].name != expected_name:
        raise AmrFormatError(
            f"{output_dir}: AMR stem {output_stem} does not match {files[0].name}"
        )
    return files[0]


def _parse_info(
    path: pathlib.Path,
) -> tuple[tuple[tuple[str, str], ...], dict[str, int], str]:
    found: dict[str, decimal.Decimal] = {}
    payload = path.read_bytes()
    for line in payload.decode("ascii").splitlines():
        if line.strip().startswith("ordering type"):
            break
        match = _INFO_ASSIGNMENT.match(line)
        if not match:
            continue
        name, value = match.groups()
        name = name.lower()
        if name not in (*_INFO_FIELDS, "ncpu", "ngridmax"):
            continue
        if name in found:
            raise AmrFormatError(f"{path}: duplicate info assignment for {name}")
        try:
            found[name] = decimal.Decimal(value)
        except decimal.InvalidOperation as error:
            raise AmrFormatError(f"{path}: invalid Decimal for {name}: {value}") from error
        if not found[name].is_finite():
            raise AmrFormatError(f"{path}: non-finite Decimal for {name}: {value}")
    missing = set((*_INFO_FIELDS, "ncpu", "ngridmax")) - set(found)
    if missing:
        raise AmrFormatError(f"{path}: missing info fields {sorted(missing)}")
    admin: dict[str, int] = {}
    for name in ("ncpu", "ngridmax"):
        integral = found[name].to_integral_value()
        if found[name] != integral:
            raise AmrFormatError(f"{path}: non-integral {name}={found[name]}")
        admin[name] = int(integral)
    semantic = tuple((name, str(found[name].normalize())) for name in _INFO_FIELDS)
    return semantic, admin, hashlib.sha256(payload).hexdigest()


def _check_info_header(
    info: tuple[tuple[str, str], ...],
    admin: dict[str, int],
    header: SemanticHeader,
    ncpu: int,
    ngridmax: int,
    path: pathlib.Path,
) -> None:
    values = dict(info)
    integer_pairs = {
        "ndim": header.ndim,
        "levelmax": header.nlevelmax,
        "nstep_coarse": header.steps[1],
    }
    for name, expected in integer_pairs.items():
        if decimal.Decimal(values[name]) != expected:
            raise AmrFormatError(
                f"{path}: info {name}={values[name]} disagrees with AMR {expected}"
            )
    if admin != {"ncpu": ncpu, "ngridmax": ngridmax}:
        raise AmrFormatError(
            f"{path}: info admin fields {admin} disagree with AMR "
            f"ncpu={ncpu}, ngridmax={ngridmax}"
        )
    binary_pairs = {
        "time": header.time,
        "aexp": header.expansion_energy[0],
        "h0": header.cosmology[4],
        "omega_m": header.cosmology[0],
        "omega_l": header.cosmology[1],
        "omega_k": header.cosmology[2],
        "omega_b": header.cosmology[3],
    }
    for name, bits in binary_pairs.items():
        expected = struct.unpack(">d", bytes.fromhex(bits))[0]
        actual = float(decimal.Decimal(values[name]))
        if not math.isclose(actual, expected, rel_tol=5.0e-15, abs_tol=5.0e-300):
            raise AmrFormatError(
                f"{path}: info {name}={values[name]} disagrees with AMR {expected:.17g}"
            )


def _topology_digest(connection: sqlite3.Connection) -> str:
    digest = hashlib.sha256(b"owner-topology-v1\0")
    rows = connection.execute(
        "SELECT level,gx,gy,gz,digest FROM grids ORDER BY level,gx,gy,gz"
    )
    for level, gx, gy, gz, payload_digest in rows:
        digest.update(_pack_i64((level, gx, gy, gz)))
        digest.update(payload_digest)
    return digest.hexdigest()


def _level_counts(connection: sqlite3.Connection) -> dict[int, int]:
    return {
        int(level): int(count)
        for level, count in connection.execute(
            "SELECT level,count(*) FROM grids GROUP BY level ORDER BY level"
        )
    }


def build_output_index(
    output_dir: pathlib.Path, db_path: pathlib.Path
) -> OutputIndex:
    output_dir = output_dir.resolve()
    files = _amr_files(output_dir)
    source_cpus = [_source_cpu(path) for path in files]
    output_stems = {_output_stem(path) for path in files}
    if len(output_stems) != 1:
        raise AmrFormatError(
            f"{output_dir}: mixed AMR output stems {sorted(output_stems)}"
        )
    output_stem = next(iter(output_stems))
    if len(set(source_cpus)) != len(source_cpus):
        raise AmrFormatError(f"{output_dir}: duplicate AMR CPU suffix")

    connection = _create_database(db_path)
    first_header: SemanticHeader | None = None
    first_ncpu: int | None = None
    first_ngridmax: int | None = None
    first_coarse: tuple[tuple[GridKey, GridKey | None, int], ...] | None = None
    input_files: list[dict[str, object]] = []
    local_layout_digests: dict[int, str] = {}
    try:
        for path, source_cpu in zip(files, source_cpus, strict=True):
            layout, local_map, file_sha256, local_layout_digest = _build_local_map(path)
            if source_cpu < 1 or source_cpu > layout.ncpu:
                raise AmrFormatError(
                    f"{path}: source CPU {source_cpu} outside 1..{layout.ncpu}"
                )
            if first_ncpu is None:
                first_ncpu = layout.ncpu
                first_ngridmax = layout.ngridmax
                expected = set(range(1, first_ncpu + 1))
                if set(source_cpus) != expected:
                    raise AmrFormatError(
                        f"{output_dir}: AMR file CPU set is incomplete; expected "
                        f"1..{first_ncpu}"
                    )
            elif layout.ncpu != first_ncpu:
                raise AmrFormatError(f"{path}: inconsistent ncpu header")
            elif layout.ngridmax != first_ngridmax:
                raise AmrFormatError(f"{path}: inconsistent ngridmax header")
            if first_header is None:
                first_header = layout.header
            elif layout.header != first_header:
                raise AmrFormatError(f"{path}: semantic header differs within output")

            coarse = _canonical_coarse(layout, local_map, path)
            if first_coarse is None:
                first_coarse = coarse
            elif coarse != first_coarse:
                raise AmrFormatError(f"{path}: coarse topology differs within output")

            owner_count, owner_level_counts, second_pass_sha256 = _emit_owner_rows(
                path, source_cpu, layout, local_map, connection
            )
            if second_pass_sha256 != file_sha256:
                raise AmrFormatError(
                    f"{path}: file changed between parser passes "
                    f"({file_sha256} != {second_pass_sha256})"
                )
            input_files.append(
                {
                    "path": str(path.resolve()),
                    "sha256": file_sha256,
                    "source_cpu": source_cpu,
                    "local_map_grids": len(local_map),
                    "owner_grids": owner_count,
                    "owner_level_counts": owner_level_counts,
                }
            )
            local_layout_digests[source_cpu] = local_layout_digest
        assert (
            first_header is not None
            and first_ncpu is not None
            and first_ngridmax is not None
            and first_coarse is not None
        )
        info_path = _info_file(output_dir, output_stem)
        info, info_admin, info_sha256 = _parse_info(info_path)
        _check_info_header(
            info, info_admin, first_header, first_ncpu, first_ngridmax, info_path
        )
        input_files.append(
            {
                "path": str(info_path.resolve()),
                "sha256": info_sha256,
                "kind": "info",
            }
        )
        level_counts = _level_counts(connection)
        for level, expected in enumerate(first_header.level_grid_totals, start=1):
            actual = level_counts.get(level, 0)
            if actual != expected:
                raise AmrFormatError(
                    f"{output_dir}: level {level} canonical owner count {actual} "
                    f"disagrees with numbtot global total {expected}"
                )
        return OutputIndex(
            output_dir=output_dir,
            db_path=db_path,
            header=first_header,
            info=info,
            ncpu=first_ncpu,
            level_counts=level_counts,
            coarse_digest=_coarse_digest(first_coarse),
            topology_digest=_topology_digest(connection),
            local_layout_digests=local_layout_digests,
            input_files=input_files,
        )
    finally:
        connection.close()


def _compare_databases(left: OutputIndex, right: OutputIndex) -> None:
    connection = sqlite3.connect(left.db_path)
    try:
        connection.execute("ATTACH DATABASE ? AS rhs", (str(right.db_path),))
        missing = connection.execute(
            """
            SELECT l.level,l.gx,l.gy,l.gz FROM main.grids l
            LEFT JOIN rhs.grids r USING(level,gx,gy,gz)
            WHERE r.level IS NULL LIMIT 1
            """
        ).fetchone()
        if missing:
            raise CanonicalMismatch(f"owner grid exists only on left: {missing}")
        extra = connection.execute(
            """
            SELECT r.level,r.gx,r.gy,r.gz FROM rhs.grids r
            LEFT JOIN main.grids l USING(level,gx,gy,gz)
            WHERE l.level IS NULL LIMIT 1
            """
        ).fetchone()
        if extra:
            raise CanonicalMismatch(f"owner grid exists only on right: {extra}")
        mismatch = connection.execute(
            """
            SELECT l.level,l.gx,l.gy,l.gz,
                   l.father=r.father,l.neighbours=r.neighbours,
                   l.sons=r.sons,l.flags=r.flags,
                   l.src_cpu,l.src_grid,r.src_cpu,r.src_grid
            FROM main.grids l JOIN rhs.grids r USING(level,gx,gy,gz)
            WHERE l.digest != r.digest LIMIT 1
            """
        ).fetchone()
        if mismatch:
            key = tuple(mismatch[:4])
            labels = ("father", "neighbours", "sons", "flag1")
            changed = [label for label, equal in zip(labels, mismatch[4:8]) if not equal]
            raise CanonicalMismatch(
                f"topology differs at grid {key}, fields={changed}, "
                f"left_source=cpu{mismatch[8]}/grid{mismatch[9]}, "
                f"right_source=cpu{mismatch[10]}/grid{mismatch[11]}"
            )
    finally:
        connection.close()


def compare_indices(left: OutputIndex, right: OutputIndex) -> None:
    if left.header != right.header:
        raise CanonicalMismatch(
            f"semantic header differs: left_sha={left.header.digest()} "
            f"right_sha={right.header.digest()}"
        )
    if left.info != right.info:
        raise CanonicalMismatch(
            f"semantic info fields differ: left={dict(left.info)} right={dict(right.info)}"
        )
    if left.coarse_digest != right.coarse_digest:
        raise CanonicalMismatch(
            f"coarse topology differs: left={left.coarse_digest} "
            f"right={right.coarse_digest}"
        )
    if left.level_counts != right.level_counts:
        raise CanonicalMismatch(
            f"level owner-grid counts differ: left={left.level_counts} "
            f"right={right.level_counts}"
        )
    _compare_databases(left, right)


def compare_topology_indices(left: OutputIndex, right: OutputIndex) -> None:
    """Compare topology while allowing physical time/state headers to differ."""
    left_geometry = (
        left.header.ndim,
        left.header.coarse_shape,
        left.header.nlevelmax,
        left.header.nboundary,
        left.header.boxlen,
    )
    right_geometry = (
        right.header.ndim,
        right.header.coarse_shape,
        right.header.nlevelmax,
        right.header.nboundary,
        right.header.boxlen,
    )
    if left_geometry != right_geometry:
        raise CanonicalMismatch(
            f"AMR geometry differs: left={left_geometry} right={right_geometry}"
        )
    if left.coarse_digest != right.coarse_digest:
        raise CanonicalMismatch(
            f"coarse topology differs: left={left.coarse_digest} "
            f"right={right.coarse_digest}"
        )
    if left.level_counts != right.level_counts:
        raise CanonicalMismatch(
            f"level owner-grid counts differ: left={left.level_counts} "
            f"right={right.level_counts}"
        )
    _compare_databases(left, right)


def compare_local_layout(left: OutputIndex, right: OutputIndex) -> None:
    """Require identical canonical grid order in every rank/owner block.

    Legacy grav files contain no grid IDs and follow the AMR linked-list order.
    This optional condition proves that a record-by-record grav comparison is
    topology keyed rather than an accidental comparison of different cells.
    """
    if left.ncpu != right.ncpu:
        raise CanonicalMismatch(
            f"local layout requires equal ncpu: left={left.ncpu} right={right.ncpu}"
        )
    if left.local_layout_digests != right.local_layout_digests:
        ranks = sorted(
            rank
            for rank in set(left.local_layout_digests) | set(right.local_layout_digests)
            if left.local_layout_digests.get(rank)
            != right.local_layout_digests.get(rank)
        )
        raise CanonicalMismatch(
            f"canonical per-rank/owner grid order differs for ranks {ranks[:8]}"
        )


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", type=pathlib.Path, help="left output_NNNNN directory")
    parser.add_argument("right", type=pathlib.Path, help="right output_NNNNN directory")
    parser.add_argument(
        "--scratch", type=pathlib.Path, help="parent for temporary SQLite databases"
    )
    parser.add_argument("--json", type=pathlib.Path, help="write a JSON report")
    parser.add_argument(
        "--topology-only",
        action="store_true",
        help="compare geometry/topology but not time, dt, cosmology, or info state",
    )
    parser.add_argument(
        "--require-same-local-layout",
        action="store_true",
        help="also require identical canonical grid order for legacy grav payloads",
    )
    return parser


def _run(args: argparse.Namespace) -> dict[str, object]:
    scratch_parent: str | None = None
    if args.scratch is not None:
        args.scratch.mkdir(parents=True, exist_ok=True)
        scratch_parent = str(args.scratch.resolve())
    with tempfile.TemporaryDirectory(
        prefix="amr-canonical-", dir=scratch_parent
    ) as temporary:
        temporary_path = pathlib.Path(temporary)
        left = build_output_index(args.left, temporary_path / "left.sqlite3")
        right = build_output_index(args.right, temporary_path / "right.sqlite3")
        if args.topology_only:
            compare_topology_indices(left, right)
        else:
            compare_indices(left, right)
        if args.require_same_local_layout:
            compare_local_layout(left, right)
        return {
            "schema": SCHEMA_VERSION,
            "mode": "topology-only" if args.topology_only else "semantic-exact",
            "same_local_layout": bool(args.require_same_local_layout),
            "status": "PASS",
            "left": left.report(),
            "right": right.report(),
            "excluded": EXCLUDED_FIELDS,
        }


def main(argv: Iterable[str] | None = None) -> int:
    args = _argument_parser().parse_args(argv)
    try:
        report = _run(args)
    except CanonicalMismatch as error:
        print(f"AMR_CANONICAL FAIL: {error}", file=sys.stderr)
        return 1
    except (AmrFormatError, OSError, sqlite3.Error, UnicodeError) as error:
        print(f"AMR_CANONICAL ERROR: {error}", file=sys.stderr)
        return 2

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        temporary_report = args.json.with_name(args.json.name + f".tmp.{os.getpid()}")
        temporary_report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary_report.replace(args.json)
    left_report = report["left"]
    assert isinstance(left_report, dict)
    print(
        "AMR_CANONICAL PASS "
        f"levels={left_report['level_counts']} "
        f"topology_sha256={left_report['topology_sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
