import unittest

import numpy as np
from lagmusic_cause_probe import metrics, periodic_last_write


class CauseProbeTests(unittest.TestCase):
    def test_relative_l2(self):
        self.assertEqual(metrics([1, 2], [1, 2])["relative_l2"], 0)
        self.assertAlmostEqual(metrics([1, 2], [2, 4])["relative_l2"], 1)

    def test_projection_matches_serial_overwrite(self):
        for shape, period in (((5, 5, 5), 4), ((10, 10, 10), 4), ((3, 3, 3), 4)):
            src = np.arange(np.prod(shape)).reshape(shape)
            target = np.zeros((period,)*3, dtype=src.dtype)
            for i in range(shape[0]):
                for j in range(shape[1]):
                    for k in range(shape[2]):
                        target[i % period, j % period, k % period] = src[i, j, k]
            expected = target[np.ix_(*(np.arange(n) % period for n in shape))]
            np.testing.assert_array_equal(periodic_last_write(src, period), expected)


if __name__ == "__main__":
    unittest.main()
