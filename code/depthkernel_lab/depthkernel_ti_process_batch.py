#!/usr/bin/env python3
"""Process-level batch driver for DAzi TI depth kernels."""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np


def partition(ncols: int, workers: int, tile_columns: int) -> list[tuple[int, int]]:
    if ncols <= 0:
        raise ValueError("ncols must be positive")
    tile_columns = max(1, tile_columns)
    starts = list(range(1, ncols + 1, tile_columns))
    tasks = [(start, min(tile_columns, ncols - start + 1)) for start in starts]
    if len(tasks) < workers and ncols >= workers:
        chunk = int(math.ceil(ncols / workers))
        starts = list(range(1, ncols + 1, chunk))
        tasks = [(start, min(chunk, ncols - start + 1)) for start in starts]
    return tasks


def run_worker(
    worker: Path,
    vel_bin: Path,
    depz_bin: Path,
    periods_bin: Path,
    nx: int,
    ny: int,
    nz: int,
    kmax: int,
    minthk: float,
    iwave: int,
    igr: int,
    start: int,
    count: int,
    outdir: Path,
) -> tuple[int, int, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(worker),
        str(vel_bin),
        str(depz_bin),
        str(periods_bin),
        str(nx),
        str(ny),
        str(nz),
        str(kmax),
        f"{minthk:.8f}",
        str(iwave),
        str(igr),
        str(start),
        str(count),
        str(outdir),
    ]
    worker_env = os.environ.copy()
    worker_env["OMP_NUM_THREADS"] = "1"
    worker_env["DAZI_OMP_THREADS"] = "1"
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=worker_env)
    (outdir / "worker.log").write_text(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(f"worker failed start={start} count={count} rc={proc.returncode}\n{proc.stdout}")
    return start, count, outdir


def read_part(part_dir: Path, count: int, kmax: int, nz: int) -> tuple[np.ndarray, np.ndarray]:
    pv_path = part_dir / "pv_part.bin"
    lsen_path = part_dir / "lsen_gsc_part.bin"
    if not pv_path.exists():
        raise FileNotFoundError(pv_path)
    if not lsen_path.exists():
        raise FileNotFoundError(lsen_path)
    pv = np.fromfile(pv_path, dtype=np.float64)
    expected_pv = count * kmax
    if pv.size != expected_pv:
        raise ValueError(f"{pv_path} has {pv.size} values, expected {expected_pv}")
    pv = pv.reshape((count, kmax), order="F")
    lsen = np.fromfile(lsen_path, dtype=np.float32)
    expected_lsen = count * kmax * (nz - 1)
    if lsen.size != expected_lsen:
        raise ValueError(f"{lsen_path} has {lsen.size} values, expected {expected_lsen}")
    lsen = lsen.reshape((count, kmax, nz - 1), order="F")
    return pv, lsen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vel-bin", required=True, type=Path)
    parser.add_argument("--depz-bin", required=True, type=Path)
    parser.add_argument("--periods-bin", required=True, type=Path)
    parser.add_argument("--nx", required=True, type=int)
    parser.add_argument("--ny", required=True, type=int)
    parser.add_argument("--nz", required=True, type=int)
    parser.add_argument("--kmax", required=True, type=int)
    parser.add_argument("--minthk", required=True, type=float)
    parser.add_argument("--iwave", required=True, type=int)
    parser.add_argument("--igr", required=True, type=int)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--tile-columns", type=int, default=64)
    parser.add_argument("--worker", type=Path, default=Path(__file__).resolve().parent / "depthkernel_ti_worker")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    t0 = time.perf_counter()
    args.outdir.mkdir(parents=True, exist_ok=True)
    ncols = args.nx * args.ny
    tasks = partition(ncols, max(1, args.workers), args.tile_columns)

    completed: list[tuple[int, int, Path]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = []
        for index, (start, count) in enumerate(tasks):
            part_dir = args.outdir / f"part_{index:05d}_{start:08d}_{count:08d}"
            futures.append(
                pool.submit(
                    run_worker,
                    args.worker,
                    args.vel_bin,
                    args.depz_bin,
                    args.periods_bin,
                    args.nx,
                    args.ny,
                    args.nz,
                    args.kmax,
                    args.minthk,
                    args.iwave,
                    args.igr,
                    start,
                    count,
                    part_dir,
                )
            )
        for future in as_completed(futures):
            completed.append(future.result())

    pv_full = np.zeros((ncols, args.kmax), dtype=np.float64, order="F")
    lsen_full = np.zeros((ncols, args.kmax, args.nz - 1), dtype=np.float32, order="F")
    for start, count, part_dir in sorted(completed):
        pv, lsen = read_part(part_dir, count, args.kmax, args.nz)
        lo = start - 1
        hi = lo + count
        pv_full[lo:hi, :] = pv
        lsen_full[lo:hi, :, :] = lsen

    np.asfortranarray(pv_full).ravel(order="F").tofile(args.outdir / "pv_ti_batch.bin")
    np.asfortranarray(lsen_full).ravel(order="F").tofile(args.outdir / "lsen_gsc_ti_batch.bin")
    wall = time.perf_counter() - t0
    meta = (
        f"columns={ncols} periods={args.kmax} nz={args.nz} "
        f"workers={max(1, args.workers)} tasks={len(tasks)} wall_seconds={wall:.6f}\n"
    )
    (args.outdir / "batch_meta.txt").write_text(meta)
    print(f"wrote {args.outdir}")
    print(meta.strip())


if __name__ == "__main__":
    main()
