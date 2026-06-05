#!/usr/bin/env python3
"""Solve DSurfTomo LSMR systems exported in the native triplet layout.

Fortran stores sparse matrices as two integer headers and one value array:
row = iw(1:nar), col = iw(nar+1:2*nar), data = rw(1:nar).
Indices are one-based.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SciPy LSMR for DSurfTomo triplets")
    parser.add_argument("--row-bin", required=True)
    parser.add_argument("--col-bin", required=True)
    parser.add_argument("--rw-bin", required=True)
    parser.add_argument("--b-bin", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--m", type=int, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--nar", type=int, required=True)
    parser.add_argument("--damp", type=float, required=True)
    parser.add_argument("--atol", type=float, required=True)
    parser.add_argument("--btol", type=float, required=True)
    parser.add_argument("--conlim", type=float, required=True)
    parser.add_argument("--maxiter", type=int, required=True)
    parser.add_argument("--index-base", type=int, default=1, choices=(0, 1))
    parser.add_argument("--float-dtype", default="float32", choices=("float32", "float64"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    row = np.memmap(args.row_bin, dtype=np.int32, mode="r", shape=(args.nar,))
    col = np.memmap(args.col_bin, dtype=np.int32, mode="r", shape=(args.nar,))
    data = np.memmap(args.rw_bin, dtype=np.float32, mode="r", shape=(args.nar,))
    b = np.memmap(args.b_bin, dtype=np.float32, mode="r", shape=(args.m,))

    row_i = np.asarray(row, dtype=np.int64) - args.index_base
    col_i = np.asarray(col, dtype=np.int64) - args.index_base
    if row_i.size and (row_i.min() < 0 or row_i.max() >= args.m):
        raise ValueError("row indices are outside matrix shape")
    if col_i.size and (col_i.min() < 0 or col_i.max() >= args.n):
        raise ValueError("col indices are outside matrix shape")

    dtype = np.float64 if args.float_dtype == "float64" else np.float32
    matrix = coo_matrix(
        (np.asarray(data, dtype=dtype), (row_i, col_i)),
        shape=(args.m, args.n),
    ).tocsr()
    rhs = np.asarray(b, dtype=dtype)
    t_build = time.perf_counter()

    result = lsmr(
        matrix,
        rhs,
        damp=args.damp,
        atol=args.atol,
        btol=args.btol,
        conlim=args.conlim,
        maxiter=args.maxiter,
        show=False,
    )
    t_solve = time.perf_counter()

    x, istop, itn, normr, normar, norma, conda, normx = result[:8]
    np.asarray(x, dtype=np.float32).tofile(outdir / "x_lsmr.bin")
    np.asarray([istop, itn], dtype=np.int32).tofile(outdir / "meta_i.bin")
    np.asarray([norma, conda, normr, normar, normx], dtype=np.float64).tofile(outdir / "meta_r.bin")

    meta = {
        "m": args.m,
        "n": args.n,
        "nar": args.nar,
        "damp": args.damp,
        "atol": args.atol,
        "btol": args.btol,
        "conlim": args.conlim,
        "maxiter": args.maxiter,
        "istop": int(istop),
        "itn": int(itn),
        "anorm": float(norma),
        "acond": float(conda),
        "rnorm": float(normr),
        "arnorm": float(normar),
        "xnorm": float(normx),
        "build_seconds": t_build - t0,
        "solve_seconds": t_solve - t_build,
        "total_seconds": t_solve - t0,
    }
    (outdir / "meta_lsmr.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(
        "SciPy LSMR: "
        f"m={args.m} n={args.n} nnz={args.nar} "
        f"itn={itn} istop={istop} total={meta['total_seconds']:.3f}s"
    )


if __name__ == "__main__":
    main()
