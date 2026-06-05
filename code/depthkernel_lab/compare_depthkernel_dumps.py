#!/usr/bin/env python3
"""Compare dumped DAzi depth-kernel arrays before G assembly."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


ARRAYS = {
    "pv": ("pv_depthkernel.bin", ("nxy", "kmax")),
    "sen_vs": ("sen_vs_depthkernel.bin", ("nxy", "kmax", "nz")),
    "sen_vp": ("sen_vp_depthkernel.bin", ("nxy", "kmax", "nz")),
    "sen_rho": ("sen_rho_depthkernel.bin", ("nxy", "kmax", "nz")),
}


def read_f64(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    data = np.fromfile(path, dtype=np.float64)
    expected = int(np.prod(shape))
    if data.size != expected:
        raise ValueError(f"{path} has {data.size} values, expected {expected}")
    return data.reshape(shape, order="F")


def empirical_coefficients(vs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coe_a = 2.0947 - 2.0 * 0.8206 * vs + 3.0 * 0.2683 * vs**2 - 4.0 * 0.0251 * vs**3
    vp = 0.9409 + 2.0947 * vs - 0.8206 * vs**2 + 0.2683 * vs**3 - 0.0251 * vs**4
    coe_rho = coe_a * (1.6612 - 2.0 * 0.4721 * vp + 3.0 * 0.0671 * vp**2 - 4.0 * 0.0043 * vp**3 + 5.0 * 0.000106 * vp**4)
    return coe_a, coe_rho


def metrics(name: str, a: np.ndarray, b: np.ndarray) -> dict[str, float | int | str]:
    diff = b - a
    abs_diff = np.abs(diff)
    norm_a = float(np.linalg.norm(a.ravel()))
    near_levels = [1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3]
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
    for level in near_levels:
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
    parser.add_argument("--vel-bin", type=Path, help="Optional float32 vel_senK.bin for effective dC/dVs comparison.")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    nxy = args.nx * args.ny
    shapes = {
        ("nxy", "kmax"): (nxy, args.kmax),
        ("nxy", "kmax", "nz"): (nxy, args.kmax, args.nz),
    }

    base_arrays: dict[str, np.ndarray] = {}
    test_arrays: dict[str, np.ndarray] = {}
    rows = []
    for name, (filename, shape_key) in ARRAYS.items():
        shape = shapes[shape_key]
        base_arrays[name] = read_f64(args.base / filename, shape)
        test_arrays[name] = read_f64(args.test / filename, shape)
        rows.append(metrics(name, base_arrays[name], test_arrays[name]))

    if args.vel_bin:
        vel = np.fromfile(args.vel_bin, dtype=np.float32).astype(np.float64).reshape((args.nx, args.ny, args.nz), order="F")
        vel_nxy = vel.reshape((nxy, args.nz), order="F")[:, None, :]
        coe_a, coe_rho = empirical_coefficients(vel_nxy)
        base_eff = base_arrays["sen_vs"] + base_arrays["sen_vp"] * coe_a + base_arrays["sen_rho"] * coe_rho
        test_eff = test_arrays["sen_vs"] + test_arrays["sen_vp"] * coe_a + test_arrays["sen_rho"] * coe_rho
        rows.append(metrics("effective_dcdvs", base_eff, test_eff))

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
