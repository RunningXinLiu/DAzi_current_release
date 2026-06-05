#!/usr/bin/env python3
"""Regression-test production 1D depth kernels on representative real columns."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import numpy as np

from sweep_real_columns import (
    parse_case_param,
    read_mod,
    representative_columns,
    representative_periods,
    write_column_input,
)


ROOT = Path(__file__).resolve().parent
BASELINE_EXE = ROOT / "depthkernel_baseline"
PROD_EXE = ROOT / "depthkernel_prod1d_driver"


def parse_numeric_output(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"No numeric rows found in {path}")
    width = max(len(row) for row in rows)
    arr = np.full((len(rows), width), np.nan, dtype=float)
    for i, row in enumerate(rows):
        arr[i, : len(row)] = row
    return arr


def compare_outputs(baseline: Path, prod: Path) -> dict[str, float]:
    a = parse_numeric_output(baseline)
    b = parse_numeric_output(prod)
    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: {baseline} {a.shape}, {prod} {b.shape}")
    diff = np.abs(a - b)
    scale = np.maximum(np.abs(a), 1.0)
    return {
        "max_abs": float(np.nanmax(diff)),
        "rms_abs": float(np.sqrt(np.nanmean(diff**2))),
        "max_rel": float(np.nanmax(diff / scale)),
    }


def run_kernel(exe: Path, input_path: Path, output_path: Path, dln: float, minthk: float) -> None:
    subprocess.run(
        [str(exe), str(input_path), str(output_path), str(dln), str(minthk)],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def parse_float_list(text: str) -> list[float]:
    vals = []
    for item in text.split(","):
        item = item.strip()
        if item:
            vals.append(float(item))
    if not vals:
        raise ValueError("Expected at least one float value")
    return vals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case-dir",
        default="/Users/liuxin/Desktop/DSurf_test/DAzi_large_data_260128/example/anning_real_SRL_ani_cell30_15_nreal050_iaratio0p5",
    )
    parser.add_argument("--param", default=None)
    parser.add_argument("--mod", default=None)
    parser.add_argument("--outdir", default="real_column_prod1d_regression")
    parser.add_argument("--period-count", type=int, default=5)
    parser.add_argument("--dln-values", default="0.01")
    parser.add_argument("--minthk-values", default=None)
    parser.add_argument("--abs-tol", type=float, default=1.0e-10)
    args = parser.parse_args()

    case_dir = Path(args.case_dir)
    param = Path(args.param) if args.param else case_dir / "para.in"
    mod = Path(args.mod) if args.mod else case_dir / "MOD"
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    nx, ny, nz, default_minthk, periods, _ = parse_case_param(param)
    dep, vs = read_mod(mod, nx, ny, nz)
    selected_periods = representative_periods(periods, args.period_count)
    columns = representative_columns(vs)
    dln_values = parse_float_list(args.dln_values)
    minthk_values = parse_float_list(args.minthk_values) if args.minthk_values else [float(default_minthk)]

    rows = []
    manifest_path = outdir / "columns_manifest.csv"
    with manifest_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["label", "i", "j", "min_vs", "max_vs", "mean_vs", "input_file"])
        writer.writeheader()
        for label, i, j in columns:
            col = vs[:, j - 1, i - 1]
            input_path = outdir / f"column_{label}_i{i}_j{j}.in"
            write_column_input(input_path, dep, col, selected_periods, default_minthk)
            writer.writerow({
                "label": label,
                "i": i,
                "j": j,
                "min_vs": float(np.min(col)),
                "max_vs": float(np.max(col)),
                "mean_vs": float(np.mean(col)),
                "input_file": input_path.name,
            })

            for minthk in minthk_values:
                for dln in dln_values:
                    stem = f"{label}_i{i}_j{j}_dln{dln:g}_minthk{minthk:g}".replace(".", "p")
                    baseline_out = outdir / f"baseline_{stem}.out"
                    prod_out = outdir / f"prod1d_{stem}.out"
                    run_kernel(BASELINE_EXE, input_path, baseline_out, dln, minthk)
                    run_kernel(PROD_EXE, input_path, prod_out, dln, minthk)
                    stats = compare_outputs(baseline_out, prod_out)
                    rows.append({
                        "column": label,
                        "i": i,
                        "j": j,
                        "dln": dln,
                        "minthk": minthk,
                        "max_abs": stats["max_abs"],
                        "rms_abs": stats["rms_abs"],
                        "max_rel": stats["max_rel"],
                        "baseline_file": baseline_out.name,
                        "prod_file": prod_out.name,
                    })

    summary_path = outdir / "prod1d_real_columns_summary.csv"
    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    max_abs = max(row["max_abs"] for row in rows)
    print(f"wrote {manifest_path}")
    print(f"wrote {summary_path}")
    print(f"tested_columns={len(columns)} cases={len(rows)} max_abs={max_abs:.6e}")
    if max_abs > args.abs_tol:
        raise SystemExit(f"FAILED: max_abs {max_abs:.6e} > abs_tol {args.abs_tol:.6e}")
    print("Production 1D real-column regression passed.")


if __name__ == "__main__":
    main()
