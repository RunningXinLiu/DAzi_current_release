#!/usr/bin/env python3
"""pDSurfTomo-style Numba batch finite-difference depth kernel.

pDSurfTomo calls this backend "disba", but the local implementation is a
Numba-vectorized translation of surf96 in pDSurfTomo-main/bin/surf96_disba.py,
not the pip package named disba.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from batch_disba_fd import make_jobs
from compare_disba import read_fortran_baseline, read_input, write_compare


ROOT = Path(__file__).resolve().parent
PDSURF_BIN = ROOT.parent / "pDSurfTomo-main" / "bin"
sys.path.insert(0, str(PDSURF_BIN))

from surf96_disba import surf96_batch  # noqa: E402


def run_pdsurf_batch(jobs, periods, iwave, igr, workers=None):
    if iwave not in (1, 2):
        raise ValueError(f"Unsupported iwave={iwave}")
    if igr not in (0, 1):
        raise ValueError(f"Unsupported igr={igr}")

    batch_size = len(jobs)
    n_periods = len(periods)
    n_layers = len(jobs[0][3])

    t_batch = np.repeat(np.asarray(periods, dtype=np.float64)[None, :], batch_size, axis=0)
    d_batch = np.empty((batch_size, n_layers), dtype=np.float64)
    a_batch = np.empty_like(d_batch)
    b_batch = np.empty_like(d_batch)
    rho_batch = np.empty_like(d_batch)

    for idx, job in enumerate(jobs):
        _, _, _, thk, rvp, rvs, rrho = job
        d_batch[idx, :] = thk
        a_batch[idx, :] = rvp
        b_batch[idx, :] = rvs
        rho_batch[idx, :] = rrho

    ifunc = 1 if iwave == 1 else 2
    itype = 1 if igr > 0 else 0
    t0 = time.perf_counter()
    curves = surf96_batch(
        t_batch,
        d_batch,
        a_batch,
        b_batch,
        rho_batch,
        iflsph=1,
        mode=0,
        itype=itype,
        ifunc=ifunc,
        dc=0.005,
        dt=0.005,
    )
    elapsed = time.perf_counter() - t0
    if curves.shape != (batch_size, n_periods):
        raise RuntimeError(f"Unexpected curve shape {curves.shape}")
    return curves, elapsed


def batch_node_fd_pdsurf(dep, vs, periods, minthk, iwave, igr, dln_vs=0.01, dln_vp=0.01, dln_rho=0.01):
    jobs, vp, rho = make_jobs(dep, vs, minthk, dln_vs, dln_vp, dln_rho)
    curves, elapsed = run_pdsurf_batch(jobs, periods, iwave, igr)

    nz = len(vs)
    nk = len(periods)
    pairs = {
        ("vs", -1): np.zeros((nz, nk)),
        ("vs", 1): np.zeros((nz, nk)),
        ("vp", -1): np.zeros((nz, nk)),
        ("vp", 1): np.zeros((nz, nk)),
        ("rho", -1): np.zeros((nz, nk)),
        ("rho", 1): np.zeros((nz, nk)),
    }
    pv = None
    for job, curve in zip(jobs, curves):
        param, iz, side = job[:3]
        if param == "pv":
            pv = curve
        else:
            pairs[(param, side)][iz, :] = curve

    if pv is None:
        raise RuntimeError("Baseline dispersion result was not computed.")

    k_vs = (pairs[("vs", 1)] - pairs[("vs", -1)]) / (dln_vs * vs[:, None])
    k_vp = (pairs[("vp", 1)] - pairs[("vp", -1)]) / (dln_vp * vp[:, None])
    k_rho = (pairs[("rho", 1)] - pairs[("rho", -1)]) / (dln_rho * rho[:, None])
    return pv, k_vs, k_vp, k_rho, elapsed, len(jobs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="sample_column.in")
    parser.add_argument("--fortran", default="baseline_kernel.out")
    parser.add_argument("--compare-out", default="pdsurf_numba_compare.out")
    parser.add_argument("--dln-vs", type=float, default=0.01)
    parser.add_argument("--dln-vp", type=float, default=0.01)
    parser.add_argument("--dln-rho", type=float, default=0.01)
    parser.add_argument("--minthk", type=float, default=None)
    args = parser.parse_args()

    dep, vs, periods, minthk, iwave, igr = read_input(Path(args.input))
    if args.minthk is not None:
        minthk = args.minthk
    f_periods, f_pv, f_kernels = read_fortran_baseline(Path(args.fortran))
    if not np.allclose(periods, f_periods):
        raise SystemExit("Input periods and Fortran output periods do not match.")

    d_pv, d_kvs, d_kvp, d_krho, elapsed, njobs = batch_node_fd_pdsurf(
        dep,
        vs,
        periods,
        minthk,
        iwave,
        igr,
        dln_vs=args.dln_vs,
        dln_vp=args.dln_vp,
        dln_rho=args.dln_rho,
    )
    write_compare(Path(args.compare_out), dep, periods, f_pv, d_pv, f_kernels, d_kvs, d_kvp, d_krho)
    print(f"Wrote {args.compare_out}")
    print(f"jobs: {njobs}")
    print(f"pDSurf Numba surf96 batch seconds: {elapsed:.6f}")
    print(f"max |velocity_pdsurf_numba - velocity_fortran| = {float(np.max(np.abs(d_pv - f_pv))):.6e} km/s")


if __name__ == "__main__":
    main()
