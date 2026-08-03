# 1D Restart Remap Example

This example validates `./mfc.sh remap` without running the Fortran solver.
It creates a synthetic single-file restart on a coarse 1D grid, remaps it to a
finer 1D grid, and checks the remapped conservative-variable fields against the
known analytic values.

Run from the MFC repository root:

```bash
python3 examples/1D_restart_remap/make_restart.py
./mfc.sh remap \
  examples/1D_restart_remap/source_case.py \
  examples/1D_restart_remap/target_case.py \
  --step 100 \
  --source-restart-dir examples/1D_restart_remap/source_restart_data \
  --output-dir examples/1D_restart_remap/remapped_restart_data \
  --force
python3 examples/1D_restart_remap/check_remap.py
```

The synthetic restart contains four conservative fields:

```text
q1 = 1 + 2x
q2 = 2 - x
q3 = 1 + x^2
q4 = 0.25
```

The current remap implementation performs linear cell-center interpolation. For
the quadratic field, the check compares against the linear interpolation of the
coarse cell-center samples, not the exact quadratic function.
