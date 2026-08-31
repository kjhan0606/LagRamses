"""Audited reader for the documented columns of RAMSES sink_XXXXX.info."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SinkMassCatalog:
    sink_id: np.ndarray
    mass_msun: np.ndarray
    position_code: np.ndarray
    velocity_code: np.ndarray

    def __post_init__(self) -> None:
        sink_id = np.asarray(self.sink_id, dtype=np.int64)
        mass = np.asarray(self.mass_msun, dtype=np.float64)
        position = np.asarray(self.position_code, dtype=np.float64)
        velocity = np.asarray(self.velocity_code, dtype=np.float64)
        if sink_id.ndim != 1 or mass.shape != sink_id.shape:
            raise ValueError("sink IDs and masses must be matching one-dimensional arrays")
        if position.shape != (len(sink_id), 3) or velocity.shape != position.shape:
            raise ValueError("sink positions and velocities must have shape (n_sink, 3)")
        if np.any(mass <= 0.0) or not np.isfinite(mass).all() or not np.isfinite(position).all():
            raise ValueError("sink masses and positions must be finite, with positive mass")
        object.__setattr__(self, "sink_id", sink_id)
        object.__setattr__(self, "mass_msun", mass)
        object.__setattr__(self, "position_code", position)
        object.__setattr__(self, "velocity_code", velocity)

    def most_massive_interior(self, half_width_code: float) -> int:
        """Return a sink whose centered cube does not cross a periodic boundary."""

        if not 0.0 < half_width_code < 0.5:
            raise ValueError("half_width_code must be between zero and one half")
        interior = np.all(
            (self.position_code >= half_width_code) & (self.position_code <= 1.0 - half_width_code), axis=1
        )
        if not np.any(interior):
            raise ValueError("no sink can host a non-wrapping target cube")
        indices = np.flatnonzero(interior)
        return int(indices[np.argmax(self.mass_msun[indices])])

    def indices_in_box(self, left_edge_code: np.ndarray, width_code: float) -> np.ndarray:
        """Return sinks inside a non-wrapping code-coordinate cube."""

        left_edge = np.asarray(left_edge_code, dtype=np.float64)
        if left_edge.shape != (3,) or width_code <= 0.0 or np.any(left_edge < 0.0) or np.any(left_edge + width_code > 1.0):
            raise ValueError("only non-wrapping sink boxes are supported")
        return np.flatnonzero(np.all((self.position_code >= left_edge) & (self.position_code < left_edge + width_code), axis=1))


def read_sink_info(path: str | Path) -> SinkMassCatalog:
    """Read only the eight columns documented by sink_XXXXX.info.

    The companion CSV has additional undocumented columns and is deliberately
    not used as an accretion-rate source.
    """

    rows = []
    with Path(path).open() as handle:
        for line in handle:
            columns = line.split()
            if len(columns) < 8:
                continue
            try:
                rows.append([float(value) for value in columns[:8]])
            except ValueError:
                continue
    if not rows:
        raise ValueError("sink info contains no documented numeric rows")
    data = np.asarray(rows, dtype=np.float64)
    return SinkMassCatalog(
        sink_id=data[:, 0].astype(np.int64),
        mass_msun=data[:, 1],
        position_code=data[:, 2:5],
        velocity_code=data[:, 5:8],
    )
