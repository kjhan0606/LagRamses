#!/usr/bin/env python3
"""Small in-memory regression checks for particle-order-independent replay."""

import unittest

import h5py
import numpy as np

from compare_evolved_hdf5 import particle_rows_equal


class ParticleRowsTest(unittest.TestCase):
    def setUp(self):
        self.reference = h5py.File("reference", "w", driver="core", backing_store=False)
        self.replay = h5py.File("replay", "w", driver="core", backing_store=False)
        self.reference_group = self.reference.create_group("particles")
        self.replay_group = self.replay.create_group("particles")
        self.rows = {
            "identity": np.array([-1, -1, -1], dtype=np.int64),
            "mass": np.array([2.0, 2.0, 2.0]),
            "x_1": np.array([0.1, 0.2, 0.3]),
            "v_1": np.array([1.0, 2.0, 3.0]),
        }
        self.fill(self.reference_group, self.rows, np.array([0, 1, 2]))
        self.fill(self.replay_group, self.rows, np.array([2, 0, 1]))

    def tearDown(self):
        self.reference.close()
        self.replay.close()

    @staticmethod
    def fill(group, rows, order):
        for name, values in rows.items():
            group.create_dataset(name, data=values[order])
        group.create_dataset("npart_per_cpu", data=np.array([len(order)], dtype=np.int32))

    def test_output_permutation_preserves_particle_records(self):
        self.assertTrue(particle_rows_equal(self.reference, self.replay))

    def test_changed_velocity_is_rejected_even_with_duplicate_identity(self):
        self.replay_group["v_1"][0] = 4.0
        self.assertFalse(particle_rows_equal(self.reference, self.replay))

    def test_missing_particle_field_is_rejected(self):
        del self.replay_group["mass"]
        self.assertFalse(particle_rows_equal(self.reference, self.replay))


if __name__ == "__main__":
    unittest.main()
