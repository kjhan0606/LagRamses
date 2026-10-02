"""Regression gates for rank-local HDF5 particle-order equivalence."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from verify_multiple_hdf5 import compare_particle_records


def _write_particles(path, order, *, changed_x=False):
    values = np.asarray([10.0, 20.0, 30.0, 40.0])
    with h5py.File(path, "w") as stream:
        group = stream.create_group("particles")
        group.create_dataset("npart_per_cpu", data=[2, 2])
        group.create_dataset("indtab", data=np.zeros(4))
        group.create_dataset("identity", data=np.asarray([1, 2, 3, 4])[order])
        group.create_dataset("levelp", data=np.asarray([4, 4, 5, 5])[order])
        group.create_dataset("ptypep", data=np.zeros(4, dtype=int))
        group.create_dataset("mass", data=values[order])
        for axis in range(1, 4):
            position = values[order].copy() + axis
            if changed_x and axis == 1:
                position[0] += 1.0
            group.create_dataset(f"x_{axis}", data=position)
            group.create_dataset(f"v_{axis}", data=values[order] - axis)


def test_rank_local_particle_permutation_is_equivalent(tmp_path):
    reference = tmp_path / "reference.h5"
    reordered = tmp_path / "reordered.h5"
    _write_particles(reference, [0, 1, 2, 3])
    _write_particles(reordered, [1, 0, 3, 2])
    with h5py.File(reference) as left, h5py.File(reordered) as right:
        fields = compare_particle_records(left["particles"], right["particles"])
    assert {"identity", "mass", "x_1", "v_3"}.issubset(fields)


@pytest.mark.parametrize(
    ("order", "changed_x"),
    [([1, 0, 3, 2], True), ([2, 1, 0, 3], False)],
)
def test_changed_particle_or_cross_rank_move_is_rejected(tmp_path, order, changed_x):
    reference = tmp_path / "reference.h5"
    changed = tmp_path / "changed.h5"
    _write_particles(reference, [0, 1, 2, 3])
    _write_particles(changed, order, changed_x=changed_x)
    with h5py.File(reference) as left, h5py.File(changed) as right:
        with pytest.raises(ValueError, match="particle records differ on rank"):
            compare_particle_records(left["particles"], right["particles"])
