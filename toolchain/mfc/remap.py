"""
Restart-field remapping for rectilinear MFC grids.

This command intentionally handles the narrow restart format used by
parallel_io=T, file_per_process=F, down_sample=F: one raw Lustre file with
conservative-variable blocks written in Fortran order.
"""

import os
from pathlib import Path

import numpy as np

from .common import MFCException, create_directory, delete_directory
from .printer import cons
from .run import input as run_input
from .state import ARG


def _case_params(case_file):
    if not os.path.isfile(case_file):
        raise MFCException(f"Case file not found: {case_file}")
    return run_input.load(case_file, do_print=False).params


def _active_shape(params):
    m = int(params["m"])
    n = int(params.get("n", 0))
    p = int(params.get("p", 0))
    return (m + 1, max(n, 0) + 1, max(p, 0) + 1)


def _uniform_boundaries(params, axis, count):
    beg_key = f"{axis}_domain%beg"
    end_key = f"{axis}_domain%end"
    stretch_key = f"stretch_{axis}"

    if params.get(stretch_key, "F") == "T":
        raise MFCException(f"Remap currently requires unstretched target grids; found {stretch_key}=T.")
    if beg_key not in params or end_key not in params:
        raise MFCException(f"Target case must define {beg_key} and {end_key}.")

    beg = float(params[beg_key])
    end = float(params[end_key])
    return np.linspace(beg, end, count + 1, dtype=np.float64)


def _target_grid(params):
    nx, ny, nz = _active_shape(params)
    x_cb = _uniform_boundaries(params, "x", nx)
    y_cb = _uniform_boundaries(params, "y", ny) if ny > 1 else np.array([-0.5, 0.5], dtype=np.float64)
    z_cb = _uniform_boundaries(params, "z", nz) if nz > 1 else np.array([-0.5, 0.5], dtype=np.float64)
    return x_cb, y_cb, z_cb


def _cell_centers(boundaries):
    return 0.5 * (boundaries[:-1] + boundaries[1:])


def _read_grid(restart_dir, shape):
    axes = ["x", "y", "z"]
    grid = []
    for axis, count in zip(axes, shape):
        path = restart_dir / f"{axis}_cb.dat"
        if count == 1 and axis != "x" and not path.exists():
            grid.append(np.array([-0.5, 0.5], dtype=np.float64))
            continue
        if not path.exists():
            raise MFCException(f"Missing restart grid file: {path}")
        values = np.fromfile(path, dtype=np.float64)
        expected = count + 1
        if values.size != expected:
            raise MFCException(f"{path} has {values.size} values; expected {expected}.")
        grid.append(values)
    return tuple(grid)


def _dtype_from_precision(params):
    precision = int(params.get("precision", 2))
    if precision == 1:
        return np.float32
    if precision == 2:
        return np.float64
    raise MFCException(f"Unsupported precision={precision}; expected 1 (single) or 2 (double).")


def _read_restart(path, shape, dtype):
    if not path.exists():
        raise MFCException(f"Missing restart file: {path}")

    cells = int(np.prod(shape))
    values = np.fromfile(path, dtype=dtype)
    if values.size % cells != 0:
        raise MFCException(f"{path} has {values.size} values, not divisible by {cells} cells.")

    nvars = values.size // cells
    fields = np.empty((nvars, *shape), dtype=np.float64)
    for var_id in range(nvars):
        beg = var_id * cells
        end = beg + cells
        fields[var_id] = values[beg:end].reshape(shape, order="F")
    return fields


def _write_restart(path, fields, dtype):
    create_directory(str(path.parent))
    with open(path, "wb") as restart_file:
        for field in fields:
            np.asarray(field, dtype=dtype).ravel(order="F").tofile(restart_file)


def _write_grid(restart_dir, grid):
    create_directory(str(restart_dir))
    for axis, values in zip(["x", "y", "z"], grid):
        if axis != "x" and values.size == 2:
            continue
        np.asarray(values, dtype=np.float64).tofile(restart_dir / f"{axis}_cb.dat")


def _interp_axis(values, old_coords, new_coords, axis):
    moved = np.moveaxis(values, axis, 0)
    out = np.empty((new_coords.size, *moved.shape[1:]), dtype=np.float64)
    flat_in = moved.reshape((old_coords.size, -1))
    flat_out = out.reshape((new_coords.size, -1))

    for col in range(flat_in.shape[1]):
        flat_out[:, col] = np.interp(new_coords, old_coords, flat_in[:, col])

    return np.moveaxis(out.reshape(out.shape), 0, axis)


def remap_fields(fields, old_grid, new_grid):
    old_centers = [_cell_centers(axis) for axis in old_grid]
    new_centers = [_cell_centers(axis) for axis in new_grid]

    out = fields
    for axis, (old_axis, new_axis, old_bounds) in enumerate(zip(old_centers, new_centers, old_grid), start=1):
        if old_axis.size == 1 and new_axis.size == 1:
            continue
        if new_axis[0] < old_bounds[0] or new_axis[-1] > old_bounds[-1]:
            raise MFCException("Target grid cell centers must lie inside the source grid extent.")
        out = _interp_axis(out, old_axis, new_axis, axis)
    return out


def remap():
    source_case = ARG("source_case")
    target_case = ARG("target_case")
    step = int(ARG("step"))

    source_params = _case_params(source_case)
    target_params = _case_params(target_case)
    dtype = _dtype_from_precision(source_params)

    source_dir = Path(ARG("source_restart_dir") or Path(source_case).parent / "restart_data")
    output_dir = Path(ARG("output_dir") or Path(target_case).parent / "restart_data")
    output_step = int(ARG("output_step"))

    source_shape = _active_shape(source_params)
    source_grid = _read_grid(source_dir, source_shape)
    target_grid = _target_grid(target_params)
    source_file = source_dir / f"{step}.dat"

    fields = _read_restart(source_file, source_shape, dtype)
    remapped = remap_fields(fields, source_grid, target_grid)

    if output_dir.exists() and not ARG("force"):
        raise MFCException(f"Output restart directory already exists: {output_dir}. Use --force to replace it.")
    if output_dir.exists():
        delete_directory(str(output_dir))

    _write_grid(output_dir, target_grid)
    _write_restart(output_dir / f"{output_step}.dat", remapped, dtype)

    cons.print(f"Remapped [bold magenta]{fields.shape[0]}[/bold magenta] conservative fields")
    cons.print(f"  source: {source_file}")
    cons.print(f"  output: {output_dir / f'{output_step}.dat'}")
