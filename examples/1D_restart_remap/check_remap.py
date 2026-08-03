#!/usr/bin/env python3
from pathlib import Path

import numpy as np

here = Path(__file__).resolve().parent
source_dir = here / "source_restart_data"
output_dir = here / "remapped_restart_data"

old_x_cb = np.fromfile(source_dir / "x_cb.dat", dtype=np.float64)
new_x_cb = np.fromfile(output_dir / "x_cb.dat", dtype=np.float64)
old_x_cc = 0.5 * (old_x_cb[:-1] + old_x_cb[1:])
new_x_cc = 0.5 * (new_x_cb[:-1] + new_x_cb[1:])

source_fields = [
    1.0 + 2.0 * old_x_cc,
    2.0 - old_x_cc,
    1.0 + old_x_cc**2,
    np.full_like(old_x_cc, 0.25),
]

expected = np.vstack([np.interp(new_x_cc, old_x_cc, field) for field in source_fields])
actual = np.fromfile(output_dir / "0.dat", dtype=np.float64).reshape((4, 8), order="C")

np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1.0e-14)

print("Remap example passed.")
print(f"max abs error: {np.max(np.abs(actual - expected)):.3e}")
