#!/usr/bin/env python3
from pathlib import Path

import numpy as np

here = Path(__file__).resolve().parent
restart_dir = here / "source_restart_data"
restart_dir.mkdir(exist_ok=True)

x_cb = np.linspace(0.0, 1.0, 5)
x_cc = 0.5 * (x_cb[:-1] + x_cb[1:])

fields = [
    1.0 + 2.0 * x_cc,
    2.0 - x_cc,
    1.0 + x_cc**2,
    np.full_like(x_cc, 0.25),
]

x_cb.astype(np.float64).tofile(restart_dir / "x_cb.dat")

with open(restart_dir / "100.dat", "wb") as restart_file:
    for field in fields:
        field.reshape((4, 1, 1), order="F").astype(np.float64).ravel(order="F").tofile(restart_file)

print(f"Wrote {restart_dir / '100.dat'}")
