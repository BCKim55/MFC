import tempfile
import unittest
from pathlib import Path

import numpy as np

from mfc import remap


class TestRestartRemap(unittest.TestCase):
    def test_remap_fields_linear_1d(self):
        old_grid = (np.array([0.0, 0.25, 0.5, 0.75, 1.0]), np.array([-0.5, 0.5]), np.array([-0.5, 0.5]))
        new_grid = (np.array([0.125, 0.375, 0.625, 0.875]), np.array([-0.5, 0.5]), np.array([-0.5, 0.5]))
        old_x = remap._cell_centers(old_grid[0])

        fields = np.zeros((2, old_x.size, 1, 1))
        fields[0, :, 0, 0] = 2.0 * old_x + 1.0
        fields[1, :, 0, 0] = -old_x

        out = remap.remap_fields(fields, old_grid, new_grid)
        new_x = remap._cell_centers(new_grid[0])

        np.testing.assert_allclose(out[0, :, 0, 0], 2.0 * new_x + 1.0)
        np.testing.assert_allclose(out[1, :, 0, 0], -new_x)

    def test_remap_allows_refined_same_domain_edges(self):
        old_grid = (np.array([0.0, 0.5, 1.0]), np.array([-0.5, 0.5]), np.array([-0.5, 0.5]))
        new_grid = (np.array([0.0, 0.25, 0.5, 0.75, 1.0]), np.array([-0.5, 0.5]), np.array([-0.5, 0.5]))
        fields = np.array([[[[1.0]], [[3.0]]]])

        out = remap.remap_fields(fields, old_grid, new_grid)

        np.testing.assert_allclose(out[0, :, 0, 0], [1.0, 1.5, 2.5, 3.0])

    def test_restart_roundtrip_fortran_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "7.dat"
            fields = np.arange(2 * 3 * 2, dtype=np.float64).reshape((2, 3, 2, 1))

            remap._write_restart(path, fields, np.float64)
            loaded = remap._read_restart(path, (3, 2, 1), np.float64)

            np.testing.assert_allclose(loaded, fields)


if __name__ == "__main__":
    unittest.main()
