#!/usr/bin/env python3
"""Compare the isolated Fortran depth-kernel baseline with disba.

Two disba products are useful here:

1. A directly comparable finite-difference kernel.  This reproduces the
   DSurfTomo node-perturbation logic, but uses disba for the dispersion
   forward calculation.
2. disba's native layer sensitivity kernels.  These are not the same
   parameterization as DSurfTomo's node kernels, but they are useful
   diagnostics for future kernel refactoring.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from disba import GroupDispersion, GroupSensitivity, PhaseDispersion, PhaseSensitivity


def read_input(path: Path):
    rows: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        rows.append(line)
    nz, kmax = map(int, rows[0].split()[:2])
    dep = np.array([float(x) for x in rows[1].split()[:nz]], dtype=float)
    vs = np.array([float(x) for x in rows[2].split()[:nz]], dtype=float)
    periods = np.array([float(x) for x in rows[3].split()[:kmax]], dtype=float)
    minthk, iwave, igr = rows[4].split()[:3]
    return dep, vs, periods, float(minthk), int(iwave), int(igr)


def empirical_vp_rho(vs: np.ndarray):
    vp = 0.9409 + 2.0947 * vs - 0.8206 * vs**2 + 0.2683 * vs**3 - 0.0251 * vs**4
    rho = (
        1.6612 * vp
        - 0.4721 * vp**2
        + 0.0671 * vp**3
        - 0.0043 * vp**4
        + 0.000106 * vp**5
    )
    return vp, rho


def refine_grid_to_layer(minthk0: float, dep: np.ndarray, vp: np.ndarray, vs: np.ndarray, rho: np.ndarray):
    thickness: list[float] = []
    rdep: list[float] = []
    rvp: list[float] = []
    rvs: list[float] = []
    rrho: list[float] = []
    sublayer_counts: list[int] = []
    initdep = 0.0
    for i in range(len(dep) - 1):
        thk = dep[i + 1] - dep[i]
        minthk = thk / minthk0
        nsublay = int((thk + 1.0e-4) / minthk) + 1
        sublayer_counts.append(nsublay)
        newthk = thk / nsublay
        for j in range(1, nsublay + 1):
            thickness.append(newthk)
            initdep += newthk
            rdep.append(initdep)
            frac = (2 * j - 1) / (2 * nsublay)
            rvp.append(vp[i] + frac * (vp[i + 1] - vp[i]))
            rvs.append(vs[i] + frac * (vs[i + 1] - vs[i]))
            rrho.append(rho[i] + frac * (rho[i + 1] - rho[i]))
    thickness.append(0.0)
    rdep.append(dep[-1])
    rvp.append(vp[-1])
    rvs.append(vs[-1])
    rrho.append(rho[-1])
    return (
        np.array(thickness, dtype=float),
        np.array(rvp, dtype=float),
        np.array(rvs, dtype=float),
        np.array(rrho, dtype=float),
        np.array(rdep, dtype=float),
        sublayer_counts,
    )


def dispersion_class(igr: int):
    return GroupDispersion if igr > 0 else PhaseDispersion


def sensitivity_class(igr: int):
    return GroupSensitivity if igr > 0 else PhaseSensitivity


def wave_name(iwave: int) -> str:
    if iwave == 1:
        return "love"
    if iwave == 2:
        return "rayleigh"
    raise ValueError(f"Unsupported iwave={iwave}")


def disba_forward(periods: np.ndarray, thk: np.ndarray, vp: np.ndarray, vs: np.ndarray, rho: np.ndarray, iwave: int, igr: int):
    disp = dispersion_class(igr)(thk, vp, vs, rho, algorithm="dunkin")
    curve = disp(periods, mode=0, wave=wave_name(iwave))
    return np.asarray(curve.velocity, dtype=float)


def disba_node_fd(dep, vs, periods, minthk, iwave, igr, dln_vs=0.01, dln_vp=0.01, dln_rho=0.01):
    vp, rho = empirical_vp_rho(vs)
    thk, rvp, rvs, rrho, _, _ = refine_grid_to_layer(minthk, dep, vp, vs, rho)
    pv = disba_forward(periods, thk, rvp, rvs, rrho, iwave, igr)
    k_vs = np.zeros((len(vs), len(periods)))
    k_vp = np.zeros_like(k_vs)
    k_rho = np.zeros_like(k_vs)

    for i in range(len(vs)):
        vsm = vs.copy()
        vpm = vp.copy()
        rhom = rho.copy()

        vsm[i] = vs[i] - 0.5 * dln_vs * vs[i]
        thk1, rvp1, rvs1, rrho1, _, _ = refine_grid_to_layer(minthk, dep, vpm, vsm, rhom)
        c1 = disba_forward(periods, thk1, rvp1, rvs1, rrho1, iwave, igr)
        vsm[i] = vs[i] + 0.5 * dln_vs * vs[i]
        thk2, rvp2, rvs2, rrho2, _, _ = refine_grid_to_layer(minthk, dep, vpm, vsm, rhom)
        c2 = disba_forward(periods, thk2, rvp2, rvs2, rrho2, iwave, igr)
        k_vs[i, :] = (c2 - c1) / (dln_vs * vs[i])

        vsm = vs.copy()
        vpm = vp.copy()
        rhom = rho.copy()
        vpm[i] = vp[i] - 0.5 * dln_vp * vp[i]
        thk1, rvp1, rvs1, rrho1, _, _ = refine_grid_to_layer(minthk, dep, vpm, vsm, rhom)
        c1 = disba_forward(periods, thk1, rvp1, rvs1, rrho1, iwave, igr)
        vpm[i] = vp[i] + 0.5 * dln_vp * vp[i]
        thk2, rvp2, rvs2, rrho2, _, _ = refine_grid_to_layer(minthk, dep, vpm, vsm, rhom)
        c2 = disba_forward(periods, thk2, rvp2, rvs2, rrho2, iwave, igr)
        k_vp[i, :] = (c2 - c1) / (dln_vp * vp[i])

        vsm = vs.copy()
        vpm = vp.copy()
        rhom = rho.copy()
        rhom[i] = rho[i] - 0.5 * dln_rho * rho[i]
        thk1, rvp1, rvs1, rrho1, _, _ = refine_grid_to_layer(minthk, dep, vpm, vsm, rhom)
        c1 = disba_forward(periods, thk1, rvp1, rvs1, rrho1, iwave, igr)
        rhom[i] = rho[i] + 0.5 * dln_rho * rho[i]
        thk2, rvp2, rvs2, rrho2, _, _ = refine_grid_to_layer(minthk, dep, vpm, vsm, rhom)
        c2 = disba_forward(periods, thk2, rvp2, rvs2, rrho2, iwave, igr)
        k_rho[i, :] = (c2 - c1) / (dln_rho * rho[i])

    return pv, k_vs, k_vp, k_rho


def read_fortran_baseline(path: Path):
    periods: list[float] = []
    velocities: list[float] = []
    kernels: dict[tuple[int, float], tuple[float, float, float]] = {}
    section = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("# periods_and"):
            section = "pv"
            continue
        if line.startswith("# kernels"):
            section = "kernels"
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if section == "pv" and len(parts) == 2:
            periods.append(float(parts[0]))
            velocities.append(float(parts[1]))
        elif section == "kernels" and len(parts) == 6:
            iz = int(parts[0])
            per = float(parts[2])
            kernels[(iz, per)] = (float(parts[3]), float(parts[4]), float(parts[5]))
    return np.array(periods), np.array(velocities), kernels


def write_compare(path: Path, dep, periods, f_pv, d_pv, f_kernels, d_kvs, d_kvp, d_krho):
    with path.open("w") as f:
        f.write("# Fortran surfdisp96 finite-difference baseline vs disba finite-difference baseline\n")
        f.write("# velocity table: period fortran_velocity disba_velocity abs_diff rel_diff\n")
        for p, fc, dc in zip(periods, f_pv, d_pv):
            diff = dc - fc
            rel = diff / fc if fc else math.nan
            f.write(f"VEL {p:12.6f} {fc:16.8f} {dc:16.8f} {diff:16.8e} {rel:16.8e}\n")
        f.write("# kernel table: iz depth period param fortran disba_fd abs_diff rel_diff\n")
        for iz, depth in enumerate(dep, start=1):
            for ip, p in enumerate(periods):
                fk = f_kernels[(iz, float(p))]
                vals = [("Vs", fk[0], d_kvs[iz - 1, ip]), ("Vp", fk[1], d_kvp[iz - 1, ip]), ("rho", fk[2], d_krho[iz - 1, ip])]
                for name, fv, dv in vals:
                    diff = dv - fv
                    rel = diff / fv if fv else math.nan
                    f.write(f"KER {iz:4d} {depth:12.6f} {p:12.6f} {name:>4s} {fv:16.8e} {dv:16.8e} {diff:16.8e} {rel:16.8e}\n")


def write_disba_layer_kernels(path: Path, periods, thk, vp, vs, rho, iwave, igr):
    sens = sensitivity_class(igr)(thk, vp, vs, rho, algorithm="dunkin")
    params = ["velocity_s", "velocity_p", "density", "thickness"]
    with path.open("w") as f:
        f.write("# disba native layer sensitivity kernels\n")
        f.write("# parameter period velocity layer_index depth_top_or_layer_depth kernel\n")
        for p in periods:
            for param in params:
                ker = sens(float(p), mode=0, wave=wave_name(iwave), parameter=param)
                for i, (depth, kval) in enumerate(zip(ker.depth, ker.kernel), start=1):
                    f.write(f"{param:>12s} {p:12.6f} {ker.velocity:16.8f} {i:5d} {depth:16.8f} {kval:16.8e}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="sample_column.in")
    parser.add_argument("--fortran", default="baseline_kernel.out")
    parser.add_argument("--compare-out", default="disba_compare.out")
    parser.add_argument("--layer-out", default="disba_layer_kernels.out")
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

    d_pv, d_kvs, d_kvp, d_krho = disba_node_fd(
        dep, vs, periods, minthk, iwave, igr,
        dln_vs=args.dln_vs, dln_vp=args.dln_vp, dln_rho=args.dln_rho,
    )
    write_compare(Path(args.compare_out), dep, periods, f_pv, d_pv, f_kernels, d_kvs, d_kvp, d_krho)

    vp, rho = empirical_vp_rho(vs)
    thk, rvp, rvs, rrho, _, _ = refine_grid_to_layer(minthk, dep, vp, vs, rho)
    write_disba_layer_kernels(Path(args.layer_out), periods, thk, rvp, rvs, rrho, iwave, igr)

    max_vdiff = float(np.max(np.abs(d_pv - f_pv)))
    print(f"Wrote {args.compare_out}")
    print(f"Wrote {args.layer_out}")
    print(f"max |velocity_disba - velocity_fortran| = {max_vdiff:.6e} km/s")


if __name__ == "__main__":
    main()
