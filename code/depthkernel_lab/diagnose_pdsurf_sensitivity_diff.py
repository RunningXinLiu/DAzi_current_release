#!/usr/bin/env python3
"""Diagnose where pDSurf sensitivity differences enter depth kernels and G."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

import numpy as np


DEPTH_ARRAYS = {
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


def read_triplets(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rw = np.fromfile(path / "rw_triplets.bin", dtype=np.float32)
    iw = np.fromfile(path / "iw_triplets.bin", dtype=np.int32)
    col = np.fromfile(path / "col_triplets.bin", dtype=np.int32)
    if not (rw.size == iw.size == col.size):
        raise ValueError(f"Triplet size mismatch under {path}: rw={rw.size} iw={iw.size} col={col.size}")
    return iw, col, rw.astype(np.float64)


def empirical_coefficients(vs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coe_a = 2.0947 - 2.0 * 0.8206 * vs + 3.0 * 0.2683 * vs**2 - 4.0 * 0.0251 * vs**3
    vp = 0.9409 + 2.0947 * vs - 0.8206 * vs**2 + 0.2683 * vs**3 - 0.0251 * vs**4
    coe_rho = coe_a * (
        1.6612
        - 2.0 * 0.4721 * vp
        + 3.0 * 0.0671 * vp**2
        - 4.0 * 0.0043 * vp**3
        + 5.0 * 0.000106 * vp**4
    )
    return coe_a, coe_rho


def abs_metrics(values: np.ndarray, base: np.ndarray | None = None) -> dict[str, float | int]:
    if values.size == 0:
        return {"n": 0, "max_abs": 0.0, "mean_abs": 0.0, "rms_abs": 0.0, "rel_l2": 0.0}
    av = np.abs(values.astype(np.float64))
    out: dict[str, float | int] = {
        "n": int(values.size),
        "max_abs": float(np.max(av)),
        "mean_abs": float(np.mean(av)),
        "rms_abs": float(np.sqrt(np.mean(values.astype(np.float64) ** 2))),
    }
    if base is not None:
        out["rel_l2"] = float(np.linalg.norm(values.ravel()) / max(np.linalg.norm(base.ravel()), 1.0e-30))
    else:
        out["rel_l2"] = 0.0
    return out


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def depth_rows(
    name: str,
    base: np.ndarray,
    test: np.ndarray,
    periods: np.ndarray | None,
    depths: np.ndarray | None,
) -> list[dict[str, object]]:
    diff = test - base
    rows: list[dict[str, object]] = []
    rows.append({"array": name, "bucket": "all", **abs_metrics(diff, base)})
    if diff.ndim >= 2:
        for ip in range(diff.shape[1]):
            period = float(periods[ip]) if periods is not None and ip < periods.size else ip + 1
            rows.append({
                "array": name,
                "bucket": "period",
                "period_index": ip + 1,
                "period": period,
                **abs_metrics(diff[:, ip, ...], base[:, ip, ...]),
            })
    if diff.ndim == 3:
        for iz in range(diff.shape[2]):
            depth = float(depths[iz]) if depths is not None and iz < depths.size else iz + 1
            rows.append({
                "array": name,
                "bucket": "depth",
                "depth_index": iz + 1,
                "depth": depth,
                **abs_metrics(diff[:, :, iz], base[:, :, iz]),
            })
        for ip in range(diff.shape[1]):
            period = float(periods[ip]) if periods is not None and ip < periods.size else ip + 1
            for iz in range(diff.shape[2]):
                depth = float(depths[iz]) if depths is not None and iz < depths.size else iz + 1
                rows.append({
                    "array": name,
                    "bucket": "period_depth",
                    "period_index": ip + 1,
                    "period": period,
                    "depth_index": iz + 1,
                    "depth": depth,
                    **abs_metrics(diff[:, ip, iz], base[:, ip, iz]),
                })
    return rows


def structured_keys(iw: np.ndarray, col: np.ndarray) -> np.ndarray:
    keys = np.empty(iw.size, dtype=[("iw", np.int32), ("col", np.int32)])
    keys["iw"] = iw
    keys["col"] = col
    return keys


def block_depth_codes(col: np.ndarray, nparpi: int, cells_per_depth: int) -> tuple[np.ndarray, np.ndarray]:
    local0 = col.astype(np.int64) - 1
    block = local0 // nparpi
    within = local0 % nparpi
    depth0 = within // cells_per_depth
    return block.astype(np.int16), depth0.astype(np.int16)


def block_name(block: int) -> str:
    return ("Vs", "Gc", "Gs")[block] if 0 <= block <= 2 else f"block_{block}"


def grouped_triplet_rows(
    label: str,
    values: np.ndarray,
    col: np.ndarray,
    nparpi: int,
    cells_per_depth: int,
    depths: np.ndarray | None,
    ftol: float,
    base: np.ndarray | None = None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if values.size == 0:
        rows.append({"set": label, "block": "all", "bucket": "all", **abs_metrics(values, base)})
        return rows

    block, depth0 = block_depth_codes(col, nparpi, cells_per_depth)
    rows.append({
        "set": label,
        "block": "all",
        "bucket": "all",
        **abs_metrics(values, base),
        "near_ftol_count": int(np.count_nonzero(np.abs(values) <= 1.5 * ftol)),
    })
    for blk in sorted(np.unique(block).tolist()):
        mask = block == blk
        bbase = base[mask] if base is not None else None
        vals = values[mask]
        rows.append({
            "set": label,
            "block": block_name(int(blk)),
            "bucket": "block",
            **abs_metrics(vals, bbase),
            "near_ftol_count": int(np.count_nonzero(np.abs(vals) <= 1.5 * ftol)),
        })
        for iz in sorted(np.unique(depth0[mask]).tolist()):
            zmask = mask & (depth0 == iz)
            zbase = base[zmask] if base is not None else None
            zvals = values[zmask]
            depth = float(depths[iz]) if depths is not None and iz < depths.size else int(iz + 1)
            rows.append({
                "set": label,
                "block": block_name(int(blk)),
                "bucket": "block_depth",
                "depth_index": int(iz + 1),
                "depth": depth,
                **abs_metrics(zvals, zbase),
                "near_ftol_count": int(np.count_nonzero(np.abs(zvals) <= 1.5 * ftol)),
            })
    return rows


def top_rows(rows: Iterable[dict[str, object]], metric: str, n: int) -> list[dict[str, object]]:
    return sorted(rows, key=lambda row: float(row.get(metric, 0.0) or 0.0), reverse=True)[:n]


def write_markdown_summary(
    path: Path,
    depth_summary: list[dict[str, object]],
    triplet_summary: list[dict[str, object]],
    top_n: int,
) -> None:
    lines = [
        "# pDSurf Sensitivity Difference Diagnostics",
        "",
        "## Top Depth-Kernel Buckets By rel_l2",
        "",
        "| array | bucket | period | depth | rel_l2 | max_abs | rms_abs | n |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in top_rows((r for r in depth_summary if r.get("bucket") == "period_depth"), "rel_l2", top_n):
        lines.append(
            f"| {row.get('array','')} | {row.get('bucket','')} | {row.get('period','')} | "
            f"{row.get('depth','')} | {float(row.get('rel_l2', 0.0)):.6g} | "
            f"{float(row.get('max_abs', 0.0)):.6g} | {float(row.get('rms_abs', 0.0)):.6g} | {row.get('n', 0)} |"
        )
    lines.extend([
        "",
        "## G Triplet Block Summary",
        "",
        "| set | block | bucket | depth | rel_l2 | max_abs | rms_abs | n | near_ftol_count |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in triplet_summary:
        if row.get("bucket") not in {"all", "block"}:
            continue
        lines.append(
            f"| {row.get('set','')} | {row.get('block','')} | {row.get('bucket','')} | {row.get('depth','')} | "
            f"{float(row.get('rel_l2', 0.0)):.6g} | {float(row.get('max_abs', 0.0)):.6g} | "
            f"{float(row.get('rms_abs', 0.0)):.6g} | {row.get('n', 0)} | {row.get('near_ftol_count', 0)} |"
        )
    lines.extend([
        "",
        "## Top G Common-Diff Block/Depth Buckets By rel_l2",
        "",
        "| set | block | depth | rel_l2 | max_abs | rms_abs | n |",
        "|---|---|---:|---:|---:|---:|---:|",
    ])
    common_depth = [
        row for row in triplet_summary
        if row.get("set") == "common_diff" and row.get("bucket") == "block_depth"
    ]
    for row in top_rows(common_depth, "rel_l2", top_n):
        lines.append(
            f"| {row.get('set','')} | {row.get('block','')} | {row.get('depth','')} | "
            f"{float(row.get('rel_l2', 0.0)):.6g} | {float(row.get('max_abs', 0.0)):.6g} | "
            f"{float(row.get('rms_abs', 0.0)):.6g} | {row.get('n', 0)} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-depth", required=True, type=Path)
    parser.add_argument("--test-depth", required=True, type=Path)
    parser.add_argument("--base-triplets", type=Path)
    parser.add_argument("--test-triplets", type=Path)
    parser.add_argument("--vel-bin", type=Path)
    parser.add_argument("--depz-bin", type=Path)
    parser.add_argument("--periods-bin", type=Path)
    parser.add_argument("--nx", required=True, type=int)
    parser.add_argument("--ny", required=True, type=int)
    parser.add_argument("--nz", required=True, type=int)
    parser.add_argument("--kmax", required=True, type=int)
    parser.add_argument("--ftol", type=float, default=1.0e-4)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    nxy = args.nx * args.ny
    periods = np.fromfile(args.periods_bin, dtype=np.float64) if args.periods_bin else None
    depths = np.fromfile(args.depz_bin, dtype=np.float32).astype(np.float64) if args.depz_bin else None
    shapes = {
        ("nxy", "kmax"): (nxy, args.kmax),
        ("nxy", "kmax", "nz"): (nxy, args.kmax, args.nz),
    }

    base_arrays: dict[str, np.ndarray] = {}
    test_arrays: dict[str, np.ndarray] = {}
    depth_summary: list[dict[str, object]] = []
    for name, (filename, shape_key) in DEPTH_ARRAYS.items():
        shape = shapes[shape_key]
        base_arrays[name] = read_f64(args.base_depth / filename, shape)
        test_arrays[name] = read_f64(args.test_depth / filename, shape)
        depth_summary.extend(depth_rows(name, base_arrays[name], test_arrays[name], periods, depths))

    if args.vel_bin:
        vel = np.fromfile(args.vel_bin, dtype=np.float32).astype(np.float64).reshape((args.nx, args.ny, args.nz), order="F")
        vel_nxy = vel.reshape((nxy, args.nz), order="F")[:, None, :]
        coe_a, coe_rho = empirical_coefficients(vel_nxy)
        base_eff = base_arrays["sen_vs"] + base_arrays["sen_vp"] * coe_a + base_arrays["sen_rho"] * coe_rho
        test_eff = test_arrays["sen_vs"] + test_arrays["sen_vp"] * coe_a + test_arrays["sen_rho"] * coe_rho
        depth_summary.extend(depth_rows("effective_dcdvs", base_eff, test_eff, periods, depths))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "depthkernel_bucket_summary.csv", depth_summary)

    triplet_summary: list[dict[str, object]] = []
    if args.base_triplets and args.test_triplets:
        biw, bcol, brw = read_triplets(args.base_triplets)
        tiw, tcol, trw = read_triplets(args.test_triplets)
        bkeys = structured_keys(biw, bcol)
        tkeys = structured_keys(tiw, tcol)
        common_keys, bidx, tidx = np.intersect1d(bkeys, tkeys, return_indices=True)
        base_only_mask = ~np.isin(bkeys, tkeys)
        test_only_mask = ~np.isin(tkeys, bkeys)
        base_only_idx = np.nonzero(base_only_mask)[0]
        test_only_idx = np.nonzero(test_only_mask)[0]
        base_only_keys = bkeys[base_only_idx]
        test_only_keys = tkeys[test_only_idx]

        cells_per_depth = (args.nx - 2) * (args.ny - 2)
        common_diff = trw[tidx] - brw[bidx]
        common_base = brw[bidx]
        triplet_summary.extend(grouped_triplet_rows(
            "common_diff",
            common_diff,
            common_keys["col"],
            cells_per_depth * (args.nz - 1),
            cells_per_depth,
            depths,
            args.ftol,
            base=common_base,
        ))
        triplet_summary.extend(grouped_triplet_rows(
            "base_only",
            brw[base_only_idx],
            base_only_keys["col"],
            cells_per_depth * (args.nz - 1),
            cells_per_depth,
            depths,
            args.ftol,
        ))
        triplet_summary.extend(grouped_triplet_rows(
            "test_only",
            trw[test_only_idx],
            test_only_keys["col"],
            cells_per_depth * (args.nz - 1),
            cells_per_depth,
            depths,
            args.ftol,
        ))
        write_csv(args.out_dir / "g_triplet_bucket_summary.csv", triplet_summary)

    write_markdown_summary(args.out_dir / "PDSURF_SENSITIVITY_DIAGNOSTICS.md", depth_summary, triplet_summary, args.top_n)
    print(f"wrote {args.out_dir / 'depthkernel_bucket_summary.csv'}")
    if triplet_summary:
        print(f"wrote {args.out_dir / 'g_triplet_bucket_summary.csv'}")
    print(f"wrote {args.out_dir / 'PDSURF_SENSITIVITY_DIAGNOSTICS.md'}")


if __name__ == "__main__":
    main()
