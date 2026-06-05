#!/usr/bin/env python3
"""Extract representative MOD columns from a real case and run kernel sweeps."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from pathlib import Path

import numpy as np

from sweep_kernel_params import parse_compare, summarize_kernel


ROOT = Path(__file__).resolve().parent
PYTHON = os.environ.get("PYTHON_DISBA", "/opt/miniconda3/envs/dispa/bin/python")
EXE = ROOT / "depthkernel_baseline"


def clean_data_lines(path: Path):
    lines = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.lower().startswith("c"):
            continue
        if " c:" in line:
            line = line.split(" c:", 1)[0].strip()
        lines.append(line)
    return lines


def parse_case_param(path: Path):
    lines = clean_data_lines(path)
    grid_idx = None
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) >= 3:
            try:
                vals = [int(float(x)) for x in parts[:3]]
            except ValueError:
                continue
            if all(v > 1 for v in vals):
                grid_idx = i
                nx, ny, nz = vals
                break
    if grid_idx is None:
        raise ValueError(f"Could not parse grid size from {path}")

    # DSurfTomo.in has minthk/nsublayer after weight+damp; DAzi para.in has it
    # directly after grid origin and spacing. In both common files here, it is
    # the first scalar integer/float line after the grid spacing line.
    scalar_after_spacing = lines[grid_idx + 3].split()[0]
    minthk = float(scalar_after_spacing)

    kmax_idx = None
    periods = None
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 1:
            try:
                k = int(float(parts[0]))
            except ValueError:
                continue
            if k > 0 and i + 1 < len(lines):
                next_parts = lines[i + 1].split()
                if len(next_parts) >= k:
                    try:
                        cand = [float(x) for x in next_parts[:k]]
                    except ValueError:
                        continue
                    if all(x > 0 for x in cand):
                        kmax_idx = i
                        periods = np.array(cand, dtype=float)
                        break
    if periods is None:
        raise ValueError(f"Could not parse periods from {path}")
    return nx, ny, nz, minthk, periods, kmax_idx


def read_mod(path: Path, nx: int, ny: int, nz: int):
    nums = np.fromstring(path.read_text(), sep=" ")
    expected = nz + nx * ny * nz
    if nums.size < expected:
        raise ValueError(f"{path} has {nums.size} numbers; expected at least {expected}")
    dep = nums[:nz]
    data = nums[nz:nz + nx * ny * nz]
    # File order follows Fortran: k depth, j row, i column.
    vs = data.reshape((nz, ny, nx))
    return dep, vs


def representative_columns(vs: np.ndarray):
    nz, ny, nx = vs.shape
    mean_vs = vs.mean(axis=0)
    shallow_vs = vs[0]
    grad = np.max(np.abs(np.diff(vs, axis=0)), axis=0)
    candidates = {
        "center": (ny // 2, nx // 2),
        "min_mean_vs": tuple(np.unravel_index(np.argmin(mean_vs), mean_vs.shape)),
        "max_mean_vs": tuple(np.unravel_index(np.argmax(mean_vs), mean_vs.shape)),
        "min_shallow_vs": tuple(np.unravel_index(np.argmin(shallow_vs), shallow_vs.shape)),
        "max_vertical_gradient": tuple(np.unravel_index(np.argmax(grad), grad.shape)),
    }
    seen = set()
    unique = []
    for label, ij in candidates.items():
        if ij in seen:
            continue
        seen.add(ij)
        unique.append((label, ij[1] + 1, ij[0] + 1))  # 1-based i,j
    return unique


def representative_periods(periods: np.ndarray, n: int = 5):
    idx = np.linspace(0, len(periods) - 1, min(n, len(periods))).round().astype(int)
    idx = np.unique(idx)
    return periods[idx]


def write_column_input(path: Path, dep, column_vs, periods, minthk, iwave=2, igr=0):
    with path.open("w") as f:
        f.write("# nz kmax\n")
        f.write(f"{len(dep)} {len(periods)}\n")
        f.write("# depz(km)\n")
        f.write(" ".join(f"{x:g}" for x in dep) + "\n")
        f.write("# vs(km/s)\n")
        f.write(" ".join(f"{x:g}" for x in column_vs) + "\n")
        f.write("# periods(s)\n")
        f.write(" ".join(f"{x:g}" for x in periods) + "\n")
        f.write("# minthk iwave igr\n")
        f.write(f"{minthk:g} {iwave} {igr}\n")


def run_one(input_path: Path, outdir: Path, label: str, dln: float, minthk: float):
    stem = f"{label}_dln{dln:g}_minthk{minthk:g}".replace(".", "p")
    baseline = outdir / f"baseline_{stem}.out"
    compare = outdir / f"compare_{stem}.out"
    layer = outdir / f"layer_{stem}.out"

    subprocess.run(
        [str(EXE), str(input_path), str(baseline), str(dln), str(minthk)],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    env = os.environ.copy()
    env.setdefault("NUMBA_CACHE_DIR", "/private/tmp/numba_cache")
    subprocess.run(
        [
            PYTHON,
            str(ROOT / "compare_disba.py"),
            "--input", str(input_path),
            "--fortran", str(baseline),
            "--compare-out", str(compare),
            "--layer-out", str(layer),
            "--dln-vs", str(dln),
            "--dln-vp", str(dln),
            "--dln-rho", str(dln),
            "--minthk", str(minthk),
        ],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    vdiffs, rows = parse_compare(compare)
    all_stats = summarize_kernel(rows)
    vs_stats = summarize_kernel(rows, "Vs")
    vp_stats = summarize_kernel(rows, "Vp")
    rho_stats = summarize_kernel(rows, "rho")
    return {
        "column": label,
        "dln": dln,
        "minthk": minthk,
        "max_abs_velocity_diff": float(np.max(np.abs(vdiffs))),
        "rms_velocity_diff": float(np.sqrt(np.mean(vdiffs**2))),
        "kernel_rel_l2_all": all_stats["rel_l2"],
        "kernel_rel_l2_vs": vs_stats["rel_l2"],
        "kernel_rel_l2_vp": vp_stats["rel_l2"],
        "kernel_rel_l2_rho": rho_stats["rel_l2"],
        "kernel_sign_mismatch_all": all_stats["sign_mismatch"],
        "compare_file": compare.name,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", default="/Users/liuxin/Desktop/DSurf_test/DAzi_large_data_260128/example/anning_real_SRL_ani_cell30_15_nreal050_iaratio0p5")
    parser.add_argument("--param", default=None)
    parser.add_argument("--mod", default=None)
    parser.add_argument("--outdir", default="real_column_sweep")
    parser.add_argument("--period-count", type=int, default=5)
    args = parser.parse_args()

    case_dir = Path(args.case_dir)
    param = Path(args.param) if args.param else case_dir / "para.in"
    mod = Path(args.mod) if args.mod else case_dir / "MOD"
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    nx, ny, nz, minthk, periods, _ = parse_case_param(param)
    dep, vs = read_mod(mod, nx, ny, nz)
    selected_periods = representative_periods(periods, args.period_count)
    columns = representative_columns(vs)

    dln_values = [0.02, 0.01, 0.005, 0.0025]
    minthk_values = sorted(set([float(minthk), 1.0, 2.0, 4.0]))

    rows = []
    manifest = outdir / "columns_manifest.csv"
    with manifest.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["label", "i", "j", "min_vs", "max_vs", "mean_vs", "input_file"])
        writer.writeheader()
        for label, i, j in columns:
            col = vs[:, j - 1, i - 1]
            input_path = outdir / f"column_{label}_i{i}_j{j}.in"
            write_column_input(input_path, dep, col, selected_periods, minthk)
            writer.writerow({
                "label": label,
                "i": i,
                "j": j,
                "min_vs": float(np.min(col)),
                "max_vs": float(np.max(col)),
                "mean_vs": float(np.mean(col)),
                "input_file": input_path.name,
            })
            for mt in minthk_values:
                for dln in dln_values:
                    run_label = f"{label}_i{i}_j{j}"
                    print(f"running {run_label} dln={dln:g} minthk={mt:g}")
                    rows.append(run_one(input_path, outdir, run_label, dln, mt))

    summary = outdir / "real_columns_sweep_summary.csv"
    with summary.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {manifest}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
