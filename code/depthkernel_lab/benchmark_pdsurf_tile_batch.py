#!/usr/bin/env python3
"""Benchmark memory-bounded pDSurfTomo Numba surf96 tile batches."""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

import numpy as np

from batch_disba_fd import make_jobs
from batch_pdsurf_numba_fd import run_pdsurf_batch
from sweep_real_columns import parse_case_param, read_mod, representative_periods


ROOT = Path(__file__).resolve().parent


def parse_int_list(text: str) -> list[int]:
    vals = []
    for item in text.split(","):
        item = item.strip()
        if item:
            vals.append(int(item))
    if not vals:
        raise ValueError("Expected at least one integer")
    return vals


def estimate_batch_mib(tile_columns: int, nz: int, nk: int, nlayer: int) -> float:
    jobs = tile_columns * (1 + 6 * nz)
    # t_batch, d_batch, a_batch, b_batch, rho_batch, curves. This is the
    # dominant Numba input/output footprint; Python job objects add overhead.
    bytes_total = jobs * (nk + 4 * nlayer + nk) * 8
    return bytes_total / (1024.0 * 1024.0)


def make_tile_jobs(dep, columns, minthk, dln):
    all_jobs = []
    meta = []
    for col_idx, col_vs in columns:
        start = len(all_jobs)
        jobs, vp, rho = make_jobs(dep, col_vs, minthk, dln, dln, dln)
        all_jobs.extend(jobs)
        meta.append((col_idx, start, len(jobs), vp, rho))
    return all_jobs, meta


def run_tile_benchmark(dep, vs_columns, periods, minthk, iwave, igr, tile_size, dln):
    nz = len(dep)
    nk = len(periods)
    total_forward_jobs = 0
    total_batch_seconds = 0.0
    total_wall = 0.0
    checksum = 0.0

    for start in range(0, len(vs_columns), tile_size):
        tile = vs_columns[start:start + tile_size]
        t0 = time.perf_counter()
        jobs, meta = make_tile_jobs(dep, tile, minthk, dln)
        curves, batch_seconds = run_pdsurf_batch(jobs, periods, iwave, igr)

        # Assemble a small checksum so this benchmark exercises the same data
        # movement pattern without retaining the whole 3D kernel in memory.
        for local_idx, (_col_idx, job_start, njobs, vp, rho) in enumerate(meta):
            col_vs = tile[local_idx][1]
            pairs = {
                ("vs", -1): np.zeros((nz, nk)),
                ("vs", 1): np.zeros((nz, nk)),
                ("vp", -1): np.zeros((nz, nk)),
                ("vp", 1): np.zeros((nz, nk)),
                ("rho", -1): np.zeros((nz, nk)),
                ("rho", 1): np.zeros((nz, nk)),
            }
            for job, curve in zip(jobs[job_start:job_start + njobs], curves[job_start:job_start + njobs]):
                param, iz, side = job[:3]
                if param != "pv":
                    pairs[(param, side)][iz, :] = curve
            checksum += float(np.sum((pairs[("vs", 1)] - pairs[("vs", -1)]) / (dln * col_vs[:, None])))
            checksum += float(np.sum((pairs[("vp", 1)] - pairs[("vp", -1)]) / (dln * vp[:, None])))
            checksum += float(np.sum((pairs[("rho", 1)] - pairs[("rho", -1)]) / (dln * rho[:, None])))

        total_forward_jobs += len(jobs)
        total_batch_seconds += batch_seconds
        total_wall += time.perf_counter() - t0

    return {
        "tile_columns": tile_size,
        "columns": len(vs_columns),
        "forward_jobs": total_forward_jobs,
        "batch_seconds": total_batch_seconds,
        "wall_seconds": total_wall,
        "jobs_per_second": total_forward_jobs / total_wall if total_wall > 0 else np.nan,
        "checksum": checksum,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case-dir",
        default="/Users/liuxin/Desktop/DSurf_test/DAzi_large_data_260128/example/anning_real_SRL_ani_cell30_15_nreal050_iaratio0p5",
    )
    parser.add_argument("--param", default=None)
    parser.add_argument(
        "--mod",
        default="/Users/liuxin/Desktop/DSurf_test/DAzi_large_data_260128/example/anning_real_SRL_ani_cell30_15_nreal050_iaratio0p5/MOD_Ref",
    )
    parser.add_argument("--outdir", default="pdsurf_tile_benchmark")
    parser.add_argument("--period-count", type=int, default=5)
    parser.add_argument("--max-columns", type=int, default=128)
    parser.add_argument("--tile-sizes", default="1,2,4,8,16,32")
    parser.add_argument("--dln", type=float, default=0.01)
    args = parser.parse_args()

    case_dir = Path(args.case_dir)
    param = Path(args.param) if args.param else case_dir / "para.in"
    mod = Path(args.mod)
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    nx, ny, nz, minthk, periods, _ = parse_case_param(param)
    dep, vs_grid = read_mod(mod, nx, ny, nz)
    selected_periods = representative_periods(periods, args.period_count)

    flat_columns = []
    for j in range(ny):
        for i in range(nx):
            flat_columns.append((j * nx + i, vs_grid[:, j, i].copy()))
    if args.max_columns > 0:
        flat_columns = flat_columns[:args.max_columns]

    # Warm up numba compilation outside timing.
    warm_jobs, _ = make_tile_jobs(dep, flat_columns[:1], minthk, args.dln)
    run_pdsurf_batch(warm_jobs, selected_periods, iwave=2, igr=0)

    rows = []
    for tile_size in parse_int_list(args.tile_sizes):
        stats = run_tile_benchmark(dep, flat_columns, selected_periods, minthk, iwave=2, igr=0, tile_size=tile_size, dln=args.dln)
        nlayer = len(warm_jobs[0][3])
        stats["estimated_batch_mib"] = estimate_batch_mib(tile_size, nz, len(selected_periods), nlayer)
        stats["nz"] = nz
        stats["period_count"] = len(selected_periods)
        stats["nlayer"] = nlayer
        rows.append(stats)
        print(
            f"tile={tile_size:4d} columns={stats['columns']:4d} jobs={stats['forward_jobs']:6d} "
            f"wall={stats['wall_seconds']:.4f}s jobs/s={stats['jobs_per_second']:.1f} "
            f"est_batch={stats['estimated_batch_mib']:.2f} MiB"
        )

    summary = outdir / "tile_benchmark_summary.csv"
    with summary.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
