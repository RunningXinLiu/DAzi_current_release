#!/usr/bin/env python3
"""Compare dumped sparse G triplets keyed by (data row, model column)."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def load_triplets(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rw = np.fromfile(path / "rw_triplets.bin", dtype=np.float32)
    iw = np.fromfile(path / "iw_triplets.bin", dtype=np.int32)
    col = np.fromfile(path / "col_triplets.bin", dtype=np.int32)
    if not (rw.size == iw.size == col.size):
        raise ValueError(f"Triplet size mismatch under {path}: rw={rw.size} iw={iw.size} col={col.size}")
    return iw, col, rw


def structured_keys(iw: np.ndarray, col: np.ndarray) -> np.ndarray:
    keys = np.empty(iw.size, dtype=[("iw", np.int32), ("col", np.int32)])
    keys["iw"] = iw
    keys["col"] = col
    return keys


def value_map(iw: np.ndarray, col: np.ndarray, rw: np.ndarray) -> dict[tuple[int, int], float]:
    return {(int(i), int(c)): float(v) for i, c, v in zip(iw, col, rw)}


def abs_stats(values: np.ndarray) -> dict[str, float | int]:
    if values.size == 0:
        return {"n": 0, "min_abs": 0.0, "max_abs": 0.0, "mean_abs": 0.0, "rms": 0.0}
    av = np.abs(values.astype(np.float64))
    return {
        "n": int(values.size),
        "min_abs": float(np.min(av)),
        "max_abs": float(np.max(av)),
        "mean_abs": float(np.mean(av)),
        "rms": float(np.sqrt(np.mean(values.astype(np.float64) ** 2))),
    }


def column_family(col: np.ndarray, nparpi: int) -> dict[str, int]:
    return {
        "Vs": int(np.count_nonzero(col <= nparpi)),
        "Gc": int(np.count_nonzero((col > nparpi) & (col <= 2 * nparpi))),
        "Gs": int(np.count_nonzero(col > 2 * nparpi)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--test", required=True, type=Path)
    parser.add_argument("--nparpi", required=True, type=int)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    biw, bcol, brw = load_triplets(args.base)
    tiw, tcol, trw = load_triplets(args.test)
    bkeys = structured_keys(biw, bcol)
    tkeys = structured_keys(tiw, tcol)

    common_keys = np.intersect1d(bkeys, tkeys)
    base_only = np.setdiff1d(bkeys, tkeys)
    test_only = np.setdiff1d(tkeys, bkeys)

    bmap = value_map(biw, bcol, brw)
    tmap = value_map(tiw, tcol, trw)
    common_diff = np.array([tmap[(int(k["iw"]), int(k["col"]))] - bmap[(int(k["iw"]), int(k["col"]))] for k in common_keys], dtype=np.float64)
    common_base = np.array([bmap[(int(k["iw"]), int(k["col"]))] for k in common_keys], dtype=np.float64)
    base_only_values = np.array([bmap[(int(k["iw"]), int(k["col"]))] for k in base_only], dtype=np.float64)
    test_only_values = np.array([tmap[(int(k["iw"]), int(k["col"]))] for k in test_only], dtype=np.float64)

    rows = []
    rows.append({"set": "base_all", **abs_stats(brw), **column_family(bcol, args.nparpi)})
    rows.append({"set": "test_all", **abs_stats(trw), **column_family(tcol, args.nparpi)})
    rows.append({"set": "common_diff", **abs_stats(common_diff), "rel_l2": float(np.linalg.norm(common_diff) / max(np.linalg.norm(common_base), 1.0e-30))})
    common_cols = common_keys["col"].astype(np.int32)
    for label, mask in [
        ("common_diff_Vs", common_cols <= args.nparpi),
        ("common_diff_Gc", (common_cols > args.nparpi) & (common_cols <= 2 * args.nparpi)),
        ("common_diff_Gs", common_cols > 2 * args.nparpi),
    ]:
        part_diff = common_diff[mask]
        part_base = common_base[mask]
        rows.append({
            "set": label,
            **abs_stats(part_diff),
            "rel_l2": float(np.linalg.norm(part_diff) / max(np.linalg.norm(part_base), 1.0e-30)),
        })
    rows.append({"set": "base_only", **abs_stats(base_only_values), **column_family(base_only["col"], args.nparpi)})
    rows.append({"set": "test_only", **abs_stats(test_only_values), **column_family(test_only["col"], args.nparpi)})

    near_levels = [1.0e-4, 1.1e-4, 1.5e-4, 2.0e-4, 5.0e-4, 1.0e-3]
    for row, vals in [(rows[-2], base_only_values), (rows[-1], test_only_values)]:
        av = np.abs(vals)
        for level in near_levels:
            row[f"abs_le_{level:g}"] = int(np.count_nonzero(av <= level))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"base_n={brw.size} test_n={trw.size} common={common_keys.size} base_only={base_only.size} test_only={test_only.size}")
    for row in rows:
        msg = f"{row['set']}: n={row['n']} min_abs={row['min_abs']:.6g} max_abs={row['max_abs']:.6g} mean_abs={row['mean_abs']:.6g}"
        if "rel_l2" in row:
            msg += f" rel_l2={row['rel_l2']:.6g}"
        if row["set"] in ("base_only", "test_only"):
            msg += f" Vs={row.get('Vs', 0)} Gc={row.get('Gc', 0)} Gs={row.get('Gs', 0)} abs<=1.5e-4={row.get('abs_le_0.00015', 0)}"
        print(msg)


if __name__ == "__main__":
    main()
