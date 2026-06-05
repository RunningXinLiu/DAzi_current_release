#!/usr/bin/env python3
"""Compare pDSurf surf96_batch forward curves against Fortran surfdisp96."""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from pdsurf_tile_senk_solver import load_transformed_surf96, make_jobs, run_pdsurf_batch


ROOT = Path(__file__).resolve().parent
DEFAULT_FORTRAN_EXE = ROOT / "surfdisp96_forward_batch"
VARIANT_SPECS = {
    "default": ("float64", "float64", "default"),
    "model_f32": ("float32", "float64", "default"),
    "output_f32": ("float64", "float32", "default"),
    "model_output_f32": ("float32", "float32", "default"),
    "strict_fastmath": ("float64", "float64", "strict_fastmath"),
    "tight_root": ("float64", "float64", "tight_root"),
    "strict_fastmath_tight_root": ("float64", "float64", "strict_fastmath_tight_root"),
}


def read_f32(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    data = np.fromfile(path, dtype=np.float32)
    expected = int(np.prod(shape))
    if data.size != expected:
        raise ValueError(f"{path} has {data.size} values, expected {expected}")
    return data.reshape(shape, order="F").astype(np.float64)


def parse_indices(text: str | None, nxy: int) -> list[int]:
    if text is None or not text.strip():
        candidates = [1, (nxy + 1) // 2, nxy]
    else:
        candidates = [int(item.strip()) for item in text.split(",") if item.strip()]
    out: list[int] = []
    seen = set()
    for idx1 in candidates:
        if idx1 < 1 or idx1 > nxy:
            raise ValueError(f"column index {idx1} outside 1..{nxy}")
        idx0 = idx1 - 1
        if idx0 in seen:
            continue
        seen.add(idx0)
        out.append(idx0)
    return out


def write_fortran_input(path: Path, jobs: list[tuple], periods: np.ndarray, iflsph: int, iwave: int, igr: int) -> None:
    nmodels = len(jobs)
    nlay = len(jobs[0][3])
    with path.open("w", encoding="utf-8") as f:
        f.write("# nmodels nlay kmax iflsph iwave igr\n")
        f.write(f"{nmodels} {nlay} {len(periods)} {iflsph} {iwave} {igr}\n")
        f.write("# periods\n")
        f.write(" ".join(f"{x:.16g}" for x in periods) + "\n")
        for i, job in enumerate(jobs, 1):
            _param, _iz, _side, thk, rvp, rvs, rrho = job
            f.write(f"# model {i}\n")
            f.write(" ".join(f"{x:.16g}" for x in thk) + "\n")
            f.write(" ".join(f"{x:.16g}" for x in rvp) + "\n")
            f.write(" ".join(f"{x:.16g}" for x in rvs) + "\n")
            f.write(" ".join(f"{x:.16g}" for x in rrho) + "\n")


def read_fortran_output(path: Path, nmodels: int, kmax: int) -> np.ndarray:
    curves = np.empty((nmodels, kmax), dtype=np.float64)
    seen = np.zeros((nmodels, kmax), dtype=bool)
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        model_idx, period_idx, _period, velocity = line.split()[:4]
        im = int(model_idx) - 1
        ip = int(period_idx) - 1
        curves[im, ip] = float(velocity)
        seen[im, ip] = True
    if not np.all(seen):
        raise ValueError(f"{path} missing {np.size(seen) - np.count_nonzero(seen)} curve values")
    return curves


def curve_metrics(diff: np.ndarray, base: np.ndarray) -> dict[str, float]:
    return {
        "max_abs": float(np.max(np.abs(diff))),
        "mean_abs": float(np.mean(np.abs(diff))),
        "rms_abs": float(np.sqrt(np.mean(diff * diff))),
        "rel_l2": float(np.linalg.norm(diff.ravel()) / max(np.linalg.norm(base.ravel()), 1.0e-30)),
    }


def job_rows(
    metadata: list[dict[str, object]],
    fortran_curves: np.ndarray,
    pdsurf_curves: np.ndarray,
    variant: str = "default",
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for i, meta in enumerate(metadata):
        diff = pdsurf_curves[i] - fortran_curves[i]
        rows.append({
            **meta,
            "variant": variant,
            **curve_metrics(diff, fortran_curves[i]),
            "fortran_min": float(np.min(fortran_curves[i])),
            "fortran_max": float(np.max(fortran_curves[i])),
            "pdsurf_min": float(np.min(pdsurf_curves[i])),
            "pdsurf_max": float(np.max(pdsurf_curves[i])),
        })
    return rows


def grouped_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        variant = row.get("variant", "default")
        groups[(variant, "all", "", "")].append(row)
        groups[(variant, "param", row["param"], "")].append(row)
        groups[(variant, "param_side", row["param"], row["side"])].append(row)

    out: list[dict[str, object]] = []
    for key, members in groups.items():
        variant, bucket, param, side = key
        max_abs = max(float(row["max_abs"]) for row in members)
        rms_abs = float(np.sqrt(np.mean([float(row["rms_abs"]) ** 2 for row in members])))
        mean_abs = float(np.mean([float(row["mean_abs"]) for row in members]))
        rel_l2 = float(np.sqrt(np.mean([float(row["rel_l2"]) ** 2 for row in members])))
        out.append({
            "variant": variant,
            "bucket": bucket,
            "param": param,
            "side": side,
            "n_jobs": len(members),
            "max_abs": max_abs,
            "mean_abs": mean_abs,
            "rms_abs": rms_abs,
            "rel_l2_rms": rel_l2,
        })
    return sorted(out, key=lambda row: (str(row["variant"]), str(row["bucket"]), str(row["param"]), str(row["side"])))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row})
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    path: Path,
    rows: list[dict[str, object]],
    grouped: list[dict[str, object]],
    variant_timings: dict[str, float],
    args: argparse.Namespace,
) -> None:
    top = sorted(rows, key=lambda row: float(row["max_abs"]), reverse=True)[:20]
    variant_all = [row for row in grouped if row["bucket"] == "all"]
    variant_all = sorted(variant_all, key=lambda row: float(row["max_abs"]))
    lines = [
        "# pDSurf Forward Engine Parity",
        "",
        "## Settings",
        "",
        f"- `dc`: `{args.dc}`",
        f"- `dt`: `{args.dt}`",
        f"- `dln`: `{args.dln}`",
        f"- `minthk`: `{args.minthk}`",
        f"- `iflsph`: `{args.iflsph}`",
        "",
        "## Variant Summary",
        "",
        "| variant | engine | model_precision | output_precision | pDSurf seconds | max_abs | rms_abs | mean_abs | rel_l2_rms |",
        "|---|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in variant_all:
        model_precision, output_precision, engine = VARIANT_SPECS[str(row["variant"])]
        lines.append(
            f"| {row['variant']} | {engine} | {model_precision} | {output_precision} | "
            f"{variant_timings.get(str(row['variant']), float('nan')):.6g} | "
            f"{float(row['max_abs']):.6g} | {float(row['rms_abs']):.6g} | "
            f"{float(row['mean_abs']):.6g} | {float(row['rel_l2_rms']):.6g} |"
        )
    lines.extend([
        "",
        "## Grouped Summary",
        "",
        "| variant | bucket | param | side | n_jobs | max_abs | rms_abs | mean_abs | rel_l2_rms |",
        "|---|---|---|---|---:|---:|---:|---:|---:|",
    ])
    for row in grouped:
        lines.append(
            f"| {row['variant']} | {row['bucket']} | {row['param']} | {row['side']} | {row['n_jobs']} | "
            f"{float(row['max_abs']):.6g} | {float(row['rms_abs']):.6g} | "
            f"{float(row['mean_abs']):.6g} | {float(row['rel_l2_rms']):.6g} |"
        )
    lines.extend([
        "",
        "## Top Job Differences",
        "",
        "| variant | column | param | depth_index | side | max_abs | rms_abs | mean_abs |",
        "|---|---:|---|---:|---:|---:|---:|---:|",
    ])
    for row in top:
        lines.append(
            f"| {row['variant']} | {row['column_index']} | {row['param']} | {row['depth_index']} | {row['side']} | "
            f"{float(row['max_abs']):.6g} | {float(row['rms_abs']):.6g} | {float(row['mean_abs']):.6g} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vel-bin", required=True, type=Path)
    parser.add_argument("--depz-bin", required=True, type=Path)
    parser.add_argument("--periods-bin", required=True, type=Path)
    parser.add_argument("--nx", required=True, type=int)
    parser.add_argument("--ny", required=True, type=int)
    parser.add_argument("--nz", required=True, type=int)
    parser.add_argument("--columns", default=None, help="1-based Fortran-order linear column indices, comma-separated")
    parser.add_argument("--minthk", required=True, type=float)
    parser.add_argument("--iflsph", type=int, choices=(0, 1), default=1)
    parser.add_argument("--iwave", type=int, default=2)
    parser.add_argument("--igr", type=int, default=0)
    parser.add_argument("--dln", type=float, default=0.01)
    parser.add_argument("--dc", type=float, default=0.005)
    parser.add_argument("--dt", type=float, default=0.005)
    parser.add_argument(
        "--variants",
        default="default,model_f32,output_f32,model_output_f32",
        help=f"comma-separated pDSurf variants. Available: {','.join(VARIANT_SPECS)}",
    )
    parser.add_argument("--fortran-exe", type=Path, default=DEFAULT_FORTRAN_EXE)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    vel = read_f32(args.vel_bin, (args.nx, args.ny, args.nz))
    depz = np.fromfile(args.depz_bin, dtype=np.float32).astype(np.float64)
    periods = np.fromfile(args.periods_bin, dtype=np.float64)
    if depz.size != args.nz:
        raise ValueError(f"depz size {depz.size} != nz {args.nz}")

    column_indices = parse_indices(args.columns, args.nx * args.ny)
    all_jobs: list[tuple] = []
    metadata: list[dict[str, object]] = []
    for col0 in column_indices:
        jj = col0 // args.nx
        ii = col0 % args.nx
        jobs, _vp, _rho = make_jobs(depz, vel[ii, jj, :].copy(), args.minthk, args.dln, args.dln, args.dln)
        for job in jobs:
            param, iz, side = job[:3]
            all_jobs.append(job)
            metadata.append({
                "column_index": col0 + 1,
                "ii": ii + 1,
                "jj": jj + 1,
                "param": param,
                "depth_index": int(iz + 1) if iz >= 0 else 0,
                "side": int(side),
            })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fortran_input = args.out_dir / "surfdisp96_forward_batch.in"
    fortran_output = args.out_dir / "surfdisp96_forward_batch.out"
    write_fortran_input(fortran_input, all_jobs, periods, args.iflsph, args.iwave, args.igr)

    fortran_exe = args.fortran_exe.expanduser().resolve()
    t0 = time.perf_counter()
    subprocess.run([str(fortran_exe), str(fortran_input), str(fortran_output)], check=True, cwd=ROOT)
    fortran_seconds = time.perf_counter() - t0
    fortran_curves = read_fortran_output(fortran_output, len(all_jobs), len(periods))

    variants = [item.strip() for item in args.variants.split(",") if item.strip()]
    unknown = [variant for variant in variants if variant not in VARIANT_SPECS]
    if unknown:
        raise ValueError(f"unknown variants: {','.join(unknown)}")

    rows: list[dict[str, object]] = []
    variant_timings: dict[str, float] = {}
    default_curves = None
    transformed_funcs = {}
    for variant in variants:
        model_precision, output_precision, engine = VARIANT_SPECS[variant]
        surf96_batch_func = None
        if engine != "default":
            if engine not in transformed_funcs:
                transformed_funcs[engine] = load_transformed_surf96(args.out_dir, engine)
            surf96_batch_func = transformed_funcs[engine]
        pdsurf_curves, pdsurf_seconds = run_pdsurf_batch(
            all_jobs,
            periods,
            args.iwave,
            args.igr,
            args.dc,
            args.dt,
            model_precision=model_precision,
            output_precision=output_precision,
            iflsph=args.iflsph,
            surf96_batch_func=surf96_batch_func,
        )
        variant_timings[variant] = pdsurf_seconds
        rows.extend(job_rows(metadata, fortran_curves, pdsurf_curves, variant))
        np.save(args.out_dir / f"pdsurf_curves_{variant}.npy", pdsurf_curves)
        if variant == "default":
            default_curves = pdsurf_curves

    grouped = grouped_rows(rows)

    write_csv(args.out_dir / "forward_job_diff.csv", rows)
    write_csv(args.out_dir / "forward_group_summary.csv", grouped)
    np.save(args.out_dir / "fortran_curves.npy", fortran_curves)
    if default_curves is not None:
        np.save(args.out_dir / "pdsurf_curves.npy", default_curves)
    write_summary(args.out_dir / "FORWARD_ENGINE_PARITY.md", rows, grouped, variant_timings, args)

    all_group = [row for row in grouped if row["bucket"] == "all"]
    best = min(all_group, key=lambda row: float(row["max_abs"]))
    print(f"columns={','.join(str(i + 1) for i in column_indices)} jobs={len(all_jobs)} periods={len(periods)}")
    print(f"fortran_seconds={fortran_seconds:.6f}")
    for variant in variants:
        variant_max = max(float(row["max_abs"]) for row in rows if row["variant"] == variant)
        print(f"{variant}_pdsurf_seconds={variant_timings[variant]:.6f} max_abs_forward_diff={variant_max:.6g}")
    print(f"best_variant={best['variant']} best_max_abs_forward_diff={float(best['max_abs']):.6g}")
    print(f"summary={args.out_dir / 'FORWARD_ENGINE_PARITY.md'}")


if __name__ == "__main__":
    main()
