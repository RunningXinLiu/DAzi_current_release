#!/usr/bin/env python3
"""Memory-bounded pDSurfTomo Numba SenK solver with Fortran binary I/O."""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PDSURF_BIN_CANDIDATES = [
    ROOT.parent / "pDSurfTomo-main" / "bin",
    ROOT.parent.parent / "pDSurfTomo-main" / "bin",
]
PDSURF_BIN = next((path for path in PDSURF_BIN_CANDIDATES if path.exists()), None)
if PDSURF_BIN is None:
    raise RuntimeError("could not locate pDSurfTomo-main/bin")
sys.path.insert(0, str(PDSURF_BIN))

from surf96_disba import surf96_batch  # noqa: E402

ENGINE_TRANSFORMS = {
    "default": (),
    "strict_fastmath": (('"fastmath": True', '"fastmath": False'),),
    "tight_root": (("1.0e-6 * c1", "1.0e-8 * c1"),),
    "strict_fastmath_tight_root": (
        ('"fastmath": True', '"fastmath": False'),
        ("1.0e-6 * c1", "1.0e-8 * c1"),
    ),
}


def load_transformed_surf96(out_dir: Path, engine: str):
    """Load a diagnostic pDSurf module copy with source-level transforms."""
    if engine == "default":
        return surf96_batch
    if engine not in ENGINE_TRANSFORMS:
        raise ValueError(f"unsupported pDSurf engine={engine}")
    if PDSURF_BIN is None:
        raise RuntimeError("could not locate pDSurfTomo-main/bin")

    src = PDSURF_BIN / "surf96_disba.py"
    text = src.read_text(encoding="utf-8")
    for old, new in ENGINE_TRANSFORMS[engine]:
        if old not in text:
            raise RuntimeError(f"did not find {old} in {src}")
        text = text.replace(old, new)

    module_dir = out_dir / f"_pdsurf_{engine}"
    module_dir.mkdir(parents=True, exist_ok=True)
    module_path = module_dir / f"surf96_disba_{engine}.py"
    module_path.write_text(text, encoding="utf-8")

    module_name = f"surf96_disba_{engine}"
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load module spec for {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module.surf96_batch


def empirical_vp_rho(vs: np.ndarray):
    vp = 0.9409 + vs * (2.0947 + vs * (-0.8206 + vs * (0.2683 - 0.0251 * vs)))
    rho = vp * (1.6612 + vp * (-0.4721 + vp * (0.0671 + vp * (-0.0043 + 0.000106 * vp))))
    return vp, rho


def refine_grid_to_layer(minthk0: float, dep: np.ndarray, vp: np.ndarray, vs: np.ndarray, rho: np.ndarray):
    thickness = []
    rvp = []
    rvs = []
    rrho = []
    initdep = 0.0
    for i in range(len(dep) - 1):
        thk = dep[i + 1] - dep[i]
        minthk = thk / minthk0
        nsublay = int((thk + 1.0e-4) / minthk) + 1
        newthk = thk / nsublay
        for j in range(1, nsublay + 1):
            thickness.append(newthk)
            initdep += newthk
            frac = (2 * j - 1) / (2 * nsublay)
            rvp.append(vp[i] + frac * (vp[i + 1] - vp[i]))
            rvs.append(vs[i] + frac * (vs[i + 1] - vs[i]))
            rrho.append(rho[i] + frac * (rho[i + 1] - rho[i]))
    thickness.append(0.0)
    rvp.append(vp[-1])
    rvs.append(vs[-1])
    rrho.append(rho[-1])
    return (
        np.asarray(thickness, dtype=np.float64),
        np.asarray(rvp, dtype=np.float64),
        np.asarray(rvs, dtype=np.float64),
        np.asarray(rrho, dtype=np.float64),
    )


def make_jobs(dep, vs, minthk, dln_vs, dln_vp, dln_rho):
    vp, rho = empirical_vp_rho(vs)
    jobs = []
    thk, rvp, rvs, rrho = refine_grid_to_layer(minthk, dep, vp, vs, rho)
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
            thk, rvp, rvs, rrho = refine_grid_to_layer(minthk, dep, vpp, vsp, rhop)
            jobs.append((param, iz, side, thk, rvp, rvs, rrho))
    return jobs, vp, rho


def run_pdsurf_batch(
    jobs,
    periods,
    iwave,
    igr,
    dc,
    dt,
    model_precision="float64",
    output_precision="float64",
    iflsph=1,
    surf96_batch_func=surf96_batch,
):
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
    if model_precision == "float32":
        d_batch = d_batch.astype(np.float32).astype(np.float64)
        a_batch = a_batch.astype(np.float32).astype(np.float64)
        b_batch = b_batch.astype(np.float32).astype(np.float64)
        rho_batch = rho_batch.astype(np.float32).astype(np.float64)
    elif model_precision != "float64":
        raise ValueError(f"unsupported model_precision={model_precision}")
    if surf96_batch_func is None:
        surf96_batch_func = surf96_batch
    ifunc = 1 if iwave == 1 else 2
    itype = 1 if igr > 0 else 0
    t0 = time.perf_counter()
    curves = surf96_batch_func(
        t_batch, d_batch, a_batch, b_batch, rho_batch,
        iflsph=iflsph, mode=0, itype=itype, ifunc=ifunc, dc=dc, dt=dt,
    )
    if output_precision == "float32":
        curves = curves.astype(np.float32).astype(np.float64)
    elif output_precision != "float64":
        raise ValueError(f"unsupported output_precision={output_precision}")
    return curves, time.perf_counter() - t0


def load_vector(path: Path, dtype) -> np.ndarray:
    return np.fromfile(path, dtype=dtype)


def save_float64_fortran(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as f:
        f.write(np.asarray(array, dtype=np.float64).tobytes(order="F"))


def process_tile(
    tile_items,
    depz,
    periods,
    minthk,
    iwave,
    igr,
    dln,
    dc,
    dt,
    model_precision="float64",
    output_precision="float64",
    iflsph=1,
    surf96_batch_func=surf96_batch,
):
    nz = len(depz)
    nk = len(periods)
    jobs = []
    meta = []
    for lin_idx, col_vs in tile_items:
        start = len(jobs)
        col_jobs, vp, rho = make_jobs(depz, col_vs, minthk, dln, dln, dln)
        jobs.extend(col_jobs)
        meta.append((lin_idx, start, len(col_jobs), col_vs, vp, rho))

    curves, batch_seconds = run_pdsurf_batch(
        jobs,
        periods,
        iwave,
        igr,
        dc,
        dt,
        model_precision=model_precision,
        output_precision=output_precision,
        iflsph=iflsph,
        surf96_batch_func=surf96_batch_func,
    )
    pv_tile = np.empty((len(tile_items), nk), dtype=np.float64)
    sen_vs_tile = np.empty((len(tile_items), nk, nz), dtype=np.float64)
    sen_vp_tile = np.empty_like(sen_vs_tile)
    sen_rho_tile = np.empty_like(sen_vs_tile)

    for out_idx, (_lin_idx, start, njobs, col_vs, vp, rho) in enumerate(meta):
        pairs = {
            ("vs", -1): np.zeros((nz, nk)),
            ("vs", 1): np.zeros((nz, nk)),
            ("vp", -1): np.zeros((nz, nk)),
            ("vp", 1): np.zeros((nz, nk)),
            ("rho", -1): np.zeros((nz, nk)),
            ("rho", 1): np.zeros((nz, nk)),
        }
        for job, curve in zip(jobs[start:start + njobs], curves[start:start + njobs]):
            param, iz, side = job[:3]
            if param == "pv":
                pv_tile[out_idx, :] = curve
            else:
                pairs[(param, side)][iz, :] = curve
        sen_vs_tile[out_idx, :, :] = ((pairs[("vs", 1)] - pairs[("vs", -1)]) / (dln * col_vs[:, None])).T
        sen_vp_tile[out_idx, :, :] = ((pairs[("vp", 1)] - pairs[("vp", -1)]) / (dln * vp[:, None])).T
        sen_rho_tile[out_idx, :, :] = ((pairs[("rho", 1)] - pairs[("rho", -1)]) / (dln * rho[:, None])).T

    return pv_tile, sen_vs_tile, sen_vp_tile, sen_rho_tile, batch_seconds, len(jobs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vel-bin", required=True)
    parser.add_argument("--depz-bin", required=True)
    parser.add_argument("--periods-bin", required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--minthk", type=float, required=True)
    parser.add_argument("--iwave", type=int, default=2)
    parser.add_argument("--igr", type=int, default=0)
    parser.add_argument("--iflsph", type=int, choices=(0, 1), default=1)
    parser.add_argument("--dln", type=float, default=0.01)
    parser.add_argument("--dc", type=float, default=0.005, help="pDSurf surf96 root-search velocity increment")
    parser.add_argument("--dt", type=float, default=0.005, help="pDSurf surf96 group-velocity frequency increment")
    parser.add_argument("--model-precision", choices=("float64", "float32"), default="float64")
    parser.add_argument("--output-precision", choices=("float64", "float32"), default="float64")
    parser.add_argument("--engine", choices=tuple(ENGINE_TRANSFORMS), default="default")
    parser.add_argument("--tile-columns", type=int, default=64)
    parser.add_argument("--outdir", default="senK_tile")
    args = parser.parse_args()

    nx, ny, nz = args.nx, args.ny, args.nz
    outdir = Path(args.outdir)
    vel = load_vector(Path(args.vel_bin), np.float32).astype(np.float64).reshape((nx, ny, nz), order="F")
    depz = load_vector(Path(args.depz_bin), np.float32).astype(np.float64)
    periods = load_vector(Path(args.periods_bin), np.float64)
    if depz.size != nz:
        raise ValueError(f"depz size {depz.size} != nz {nz}")
    surf96_batch_func = load_transformed_surf96(outdir, args.engine)

    nxy = nx * ny
    nk = len(periods)
    pv = np.empty((nxy, nk), dtype=np.float64)
    sen_vs = np.empty((nxy, nk, nz), dtype=np.float64)
    sen_vp = np.empty_like(sen_vs)
    sen_rho = np.empty_like(sen_vs)

    flat_items = []
    for jj in range(ny):
        for ii in range(nx):
            lin_idx = jj * nx + ii
            flat_items.append((lin_idx, vel[ii, jj, :].copy()))

    total_jobs = 0
    total_batch = 0.0
    t0 = time.perf_counter()
    for start in range(0, nxy, args.tile_columns):
        tile_items = flat_items[start:start + args.tile_columns]
        pv_tile, sv_tile, sp_tile, sr_tile, batch_seconds, njobs = process_tile(
            tile_items,
            depz,
            periods,
            args.minthk,
            args.iwave,
            args.igr,
            args.dln,
            args.dc,
            args.dt,
            model_precision=args.model_precision,
            output_precision=args.output_precision,
            iflsph=args.iflsph,
            surf96_batch_func=surf96_batch_func,
        )
        for local_idx, (lin_idx, _col_vs) in enumerate(tile_items):
            pv[lin_idx, :] = pv_tile[local_idx, :]
            sen_vs[lin_idx, :, :] = sv_tile[local_idx, :, :]
            sen_vp[lin_idx, :, :] = sp_tile[local_idx, :, :]
            sen_rho[lin_idx, :, :] = sr_tile[local_idx, :, :]
        total_jobs += njobs
        total_batch += batch_seconds
    wall = time.perf_counter() - t0

    save_float64_fortran(outdir / "pv_pdsurf_numba.bin", pv)
    save_float64_fortran(outdir / "sen_vs_pdsurf_numba.bin", sen_vs)
    save_float64_fortran(outdir / "sen_vp_pdsurf_numba.bin", sen_vp)
    save_float64_fortran(outdir / "sen_rho_pdsurf_numba.bin", sen_rho)
    print(f"wrote {outdir}")
    print(f"columns={nxy} periods={nk} nz={nz} jobs={total_jobs}")
    print(f"iflsph={args.iflsph}")
    print(f"engine={args.engine}")
    print(f"model_precision={args.model_precision} output_precision={args.output_precision}")
    print(f"batch_seconds={total_batch:.6f} wall_seconds={wall:.6f}")


if __name__ == "__main__":
    main()
