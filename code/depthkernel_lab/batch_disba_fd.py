#!/usr/bin/env python3
"""CPU batch finite-difference depth kernel using disba.

This is a local prototype of the pDSurfTomo optimization idea:
build all +/- perturbation models first, run the independent dispersion
forward jobs in parallel, then assemble the finite-difference kernels.
It keeps the same node-perturbation parameterization as DSurfTomo.
"""

from __future__ import annotations

import argparse
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from compare_disba import (
    disba_forward,
    empirical_vp_rho,
    read_fortran_baseline,
    read_input,
    refine_grid_to_layer,
    write_compare,
)


def make_jobs(dep, vs, minthk, dln_vs, dln_vp, dln_rho):
    vp, rho = empirical_vp_rho(vs)
    jobs = []

    thk, rvp, rvs, rrho, _, _ = refine_grid_to_layer(minthk, dep, vp, vs, rho)
    jobs.append(("pv", -1, 0, thk, rvp, rvs, rrho))

    for iz in range(len(vs)):
        for param, side in (("vs", -1), ("vs", 1), ("vp", -1), ("vp", 1), ("rho", -1), ("rho", 1)):
            vsp = vs.copy()
            vpp = vp.copy()
            rhop = rho.copy()
            if param == "vs":
                vsp[iz] *= 1.0 + side * 0.5 * dln_vs
            elif param == "vp":
                vpp[iz] *= 1.0 + side * 0.5 * dln_vp
            else:
                rhop[iz] *= 1.0 + side * 0.5 * dln_rho
            thk, rvp, rvs, rrho, _, _ = refine_grid_to_layer(minthk, dep, vpp, vsp, rhop)
            jobs.append((param, iz, side, thk, rvp, rvs, rrho))
    return jobs, vp, rho


def run_job(job, periods, iwave, igr):
    param, iz, side, thk, rvp, rvs, rrho = job
    curve = disba_forward(periods, thk, rvp, rvs, rrho, iwave, igr)
    return param, iz, side, curve


def batch_node_fd(dep, vs, periods, minthk, iwave, igr, dln_vs=0.01, dln_vp=0.01, dln_rho=0.01, workers=1):
    jobs, vp, rho = make_jobs(dep, vs, minthk, dln_vs, dln_vp, dln_rho)
    t0 = time.perf_counter()
    if workers <= 1:
        results = [run_job(job, periods, iwave, igr) for job in jobs]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(lambda job: run_job(job, periods, iwave, igr), jobs))
    elapsed = time.perf_counter() - t0

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
    for param, iz, side, curve in results:
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
    parser.add_argument("--compare-out", default="batch_disba_compare.out")
    parser.add_argument("--dln-vs", type=float, default=0.01)
    parser.add_argument("--dln-vp", type=float, default=0.01)
    parser.add_argument("--dln-rho", type=float, default=0.01)
    parser.add_argument("--minthk", type=float, default=None)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    dep, vs, periods, minthk, iwave, igr = read_input(Path(args.input))
    if args.minthk is not None:
        minthk = args.minthk
    f_periods, f_pv, f_kernels = read_fortran_baseline(Path(args.fortran))
    if not np.allclose(periods, f_periods):
        raise SystemExit("Input periods and Fortran output periods do not match.")

    d_pv, d_kvs, d_kvp, d_krho, elapsed, njobs = batch_node_fd(
        dep,
        vs,
        periods,
        minthk,
        iwave,
        igr,
        dln_vs=args.dln_vs,
        dln_vp=args.dln_vp,
        dln_rho=args.dln_rho,
        workers=args.workers,
    )
    write_compare(Path(args.compare_out), dep, periods, f_pv, d_pv, f_kernels, d_kvs, d_kvp, d_krho)
    print(f"Wrote {args.compare_out}")
    print(f"jobs: {njobs}")
    print(f"workers: {args.workers}")
    print(f"batch disba seconds: {elapsed:.6f}")
    print(f"max |velocity_disba - velocity_fortran| = {float(np.max(np.abs(d_pv - f_pv))):.6e} km/s")


if __name__ == "__main__":
    main()
