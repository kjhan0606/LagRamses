"""Audited external photon-source ledgers for post-processed RT."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from snrt_core.snapshot import SourceCatalog


_GROUP_COLUMN = re.compile(r"q_group_(\d+)_s")


@dataclass(frozen=True)
class PhotonSourceLedger:
    """Photon-number source luminosities with positions in RAMSES code units.

    The ledger is intentionally downstream of stellar-population and AGN SED
    modelling.  It never infers luminosity from a sink mass or checkpoint
    accumulator.
    """

    source_id: np.ndarray
    source_kind: np.ndarray
    position_code: np.ndarray
    photon_luminosity_s: np.ndarray

    def __post_init__(self) -> None:
        source_id = np.asarray(self.source_id, dtype=np.int64)
        source_kind = np.asarray(self.source_kind, dtype=str)
        position = np.asarray(self.position_code, dtype=np.float64)
        luminosity = np.asarray(self.photon_luminosity_s, dtype=np.float64)
        if source_id.ndim != 1 or source_kind.shape != source_id.shape:
            raise ValueError("source IDs and kinds must be matching one-dimensional arrays")
        if position.shape != (len(source_id), 3):
            raise ValueError("source positions must have shape (n_source, 3)")
        if luminosity.ndim != 2 or luminosity.shape[0] != len(source_id) or luminosity.shape[1] == 0:
            raise ValueError("photon luminosity must have shape (n_source, n_group)")
        if not np.isfinite(position).all() or np.any(position < 0.0) or np.any(position > 1.0):
            raise ValueError("source code positions must be finite and lie in [0, 1]")
        if not np.isfinite(luminosity).all() or np.any(luminosity < 0.0):
            raise ValueError("photon luminosities must be finite and non-negative")
        object.__setattr__(self, "source_id", source_id)
        object.__setattr__(self, "source_kind", source_kind)
        object.__setattr__(self, "position_code", position)
        object.__setattr__(self, "photon_luminosity_s", luminosity)

    def source_catalog_in_cube(
        self,
        left_edge_code: np.ndarray,
        width_code: float,
        dimensions: tuple[int, int, int],
    ) -> SourceCatalog | None:
        """Return ledger sources in a non-wrapping cube as cell-centered RT sources."""

        left = np.asarray(left_edge_code, dtype=np.float64)
        shape = np.asarray(dimensions, dtype=np.int64)
        if left.shape != (3,) or np.any(left < 0.0) or width_code <= 0.0 or np.any(left + width_code > 1.0):
            raise ValueError("source deposition requires a non-wrapping code-coordinate cube")
        if shape.shape != (3,) or np.any(shape <= 0):
            raise ValueError("dimensions must contain three positive integers")
        inside = np.all((self.position_code >= left) & (self.position_code < left + width_code), axis=1)
        if not np.any(inside):
            return None
        local_coordinate = (self.position_code[inside] - left) / width_code
        cell_index = np.floor(local_coordinate * shape).astype(np.int64)
        cell_index = np.minimum(cell_index, shape - 1)
        return SourceCatalog(cell_index=cell_index, photon_luminosity_s=self.photon_luminosity_s[inside])


def read_photon_source_ledger_csv(path: str | Path) -> PhotonSourceLedger:
    """Read a pre-audited photon ledger without applying any SED conversion.

    Required columns are ``source_id``, ``source_kind``, ``x_code``,
    ``y_code``, ``z_code``, and contiguous ``q_group_N_s`` columns beginning
    at zero.  The source-model provenance belongs in the accompanying ledger
    metadata described in ``P4_SOURCE_LEDGER.md``.
    """

    ledger_path = Path(path)
    with ledger_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("source ledger CSV requires a header")
        group_indices = sorted(
            int(match.group(1))
            for name in reader.fieldnames
            if (match := _GROUP_COLUMN.fullmatch(name)) is not None
        )
        if group_indices != list(range(len(group_indices))):
            raise ValueError("photon group columns must be contiguous from q_group_0_s")
        required = {"source_id", "source_kind", "x_code", "y_code", "z_code"}
        missing = required - set(reader.fieldnames)
        if missing or not group_indices:
            raise ValueError(f"source ledger missing required columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("source ledger contains no sources")
    return PhotonSourceLedger(
        source_id=np.asarray([int(row["source_id"]) for row in rows]),
        source_kind=np.asarray([row["source_kind"] for row in rows]),
        position_code=np.asarray(
            [[float(row[axis]) for axis in ("x_code", "y_code", "z_code")] for row in rows], dtype=np.float64
        ),
        photon_luminosity_s=np.asarray(
            [[float(row[f"q_group_{group}_s"]) for group in group_indices] for row in rows], dtype=np.float64
        ),
    )
