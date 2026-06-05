#!/usr/bin/env python3
"""Compare dumped DAzi TI depth-kernel arrays before G assembly."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


ARRAYS = {
    "pv_ti": ("pv_ti.bin", "float64", ("nxy", "kmax")),
    "lsen_gsc_ti": ("lsen_gsc_ti.bin", "float32", ("nxy", "kmax", "nz_minus_1")),
}


def read_array(path: Path, dtype: str, shape: tuple[int, ...]) -> np.ndarray:
    data = np.fromfile(path, dtype=np.dtype(dtype))
    expected = int(np.prod(shape))
    if data.size != expected:
        raise ValueError(f"{path} has {data.size} values, expected {expected}")
    return data.reshape(shape, order="F")


def metrics(name: str, base: np.ndarray, test: np.ndarray) -> dict[str, float | int | str]:
    a = base.astype(np.float64, copy=False)
    b = test.astype(np.float64, copy=False)
    diff = b - a
    abs_diff = np.abs(diff)
    norm_a = float(np.linalg.norm(a.ravel()))
    row: dict[str, float | int | str] = {
        "array": name,
        "shape": "x".join(map(str, a.shape)),
        "n": int(a.size),
        "base_min": float(np.min(a)),
        "base_max": float(np.max(a)),
        "test_min": float(np.min(b)),
        "test_max": float(np.max(b)),
        "max_abs": float(np.max(abs_diff)),
        "mean_abs": float(np.mean(abs_diff)),
        "rms_abs": float(np.sqrt(np.mean(diff * diff))),
        "rel_l2": float(np.linalg.norm(diff.ravel()) / max(norm_a, 1.0e-30)),
    }
    for level in (1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3):
        row[f"base_abs_lt_{level:g}"] = int(np.count_nonzero(np.abs(a) < level))
        row[f"cross_abs_{level:g}"] = int(np.count_nonzero((np.abs(a) <= level) != (np.abs(b) <= level)))
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--test", required=True, type=Path)
    parser.add_argument("--nx", required=True, type=int)
    parser.add_argument("--ny", required=True, type=int)
    parser.add_argument("--nz", required=True, type=int)
    parser.add_argument("--kmax", required=True, type=int)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    nxy = args.nx * args.ny
    shapes = {
        ("nxy", "kmax"): (nxy, args.kmax),
        ("nxy", "kmax", "nz_minus_1"): (nxy, args.kmax, args.nz - 1),
    }

    rows = []
    for name, (filename, dtype, shape_key) in ARRAYS.items():
        shape = shapes[shape_key]
        base = read_array(args.base / filename, dtype, shape)
        test = read_array(args.test / filename, dtype, shape)
        rows.append(metrics(name, base, test))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['array']}: max_abs={row['max_abs']:.6g} "
            f"rms_abs={row['rms_abs']:.6g} rel_l2={row['rel_l2']:.6g} "
            f"cross_abs_1e-4={row.get('cross_abs_0.0001')}"
        )


if __name__ == "__main__":
    main()
