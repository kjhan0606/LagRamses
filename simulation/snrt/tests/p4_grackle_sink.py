"""P4 Grackle binary and documented sink-info parsing checks."""

from __future__ import annotations

from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.grackle import read_grackle_equilibrium_table
from snrt_core.sink_catalog import read_sink_info


def main() -> None:
    with TemporaryDirectory() as temporary_directory:
        directory = Path(temporary_directory)
        table_path = directory / "grackle.bin"
        with table_path.open("wb") as handle:
            np.asarray([2, 2, 2], dtype=np.int32).tofile(handle)
            np.asarray([3.8, -6.0, 2.0, -6.0, 0.0, 1.0, 9.0], dtype=np.float64).tofile(handle)
            np.asarray([1.0, 9.0], dtype=np.float64).tofile(handle)
            np.zeros(8, dtype=np.float64).tofile(handle)
            np.full(8, 1.1, dtype=np.float64).tofile(handle)
        table = read_grackle_equilibrium_table(table_path)
        assert np.isclose(table.mean_mu(1.0e5, 1.0e-2, 1.0e-3), 1.1)

        sink_path = directory / "sink.info"
        sink_path.write_text(
            " Number of sink =       2\n===\n Id Mass x y z vx vy vz\n===\n"
            " 1 1.0e6 0.2 0.3 0.4 0.0 0.0 0.0\n 2 2.0e6 0.8 0.7 0.6 0.0 0.0 0.0\n"
        )
        sinks = read_sink_info(sink_path)
        assert sinks.sink_id.tolist() == [1, 2]
        assert sinks.most_massive_interior(0.1) == 1
    print("P4_GRACKLE_SINK_OK redshift=3.8 sinks=2")


if __name__ == "__main__":
    main()
