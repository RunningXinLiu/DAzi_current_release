#!/usr/bin/env python3
"""Sweep finite-difference step and sublayer refinement parameters."""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
PYTHON = os.environ.get("PYTHON_DISBA", "/opt/miniconda3/envs/dispa/bin/python")
INPUT = ROOT / "sample_column.in"
EXE = ROOT / "depthkernel_baseline"


def parse_compare(path: Path):
    velocity_diffs = []
    kernel_rows = []
    for raw in path.read_text().splitlines():
        parts = raw.split()
        if not parts or parts[0].startswith("#"):
            continue
        if parts[0] == "VEL":
            velocity_diffs.append(float(parts[4]))
        elif parts[0] == "KER":
            param = parts[4]
            fval = float(parts[5])
            dval = float(parts[6])
            diff = float(parts[7])
            kernel_rows.append((param, fval, dval, diff))
    return np.array(velocity_diffs, dtype=float), kernel_rows


def summarize_kernel(rows, param: str | None = None):
    if param is not None:
        rows = [r for r in rows if r[0] == param]
    fvals = np.array([r[1] for r in rows], dtype=float)
    dvals = np.array([r[2] for r in rows], dtype=float)
    diffs = np.array([r[3] for r in rows], dtype=float)
    if len(diffs) == 0:
        return dict(max_abs=np.nan, rms=np.nan, rel_l2=np.nan, sign_mismatch=np.nan)
    denom = np.linalg.norm(fvals)
    sign_mask = (np.abs(fvals) > 1.0e-6) & (np.abs(dvals) > 1.0e-6)
    sign_mismatch = np.mean(np.sign(fvals[sign_mask]) != np.sign(dvals[sign_mask])) if np.any(sign_mask) else 0.0
    return {
        "max_abs": float(np.max(np.abs(diffs))),
        "rms": float(np.sqrt(np.mean(diffs**2))),
        "rel_l2": float(np.linalg.norm(diffs) / denom) if denom else np.nan,
        "sign_mismatch": float(sign_mismatch),
    }


def run_case(dln: float, minthk: float):
    stem = f"dln{dln:g}_minthk{minthk:g}".replace(".", "p")
    baseline = ROOT / f"baseline_{stem}.out"
    compare = ROOT / f"compare_{stem}.out"
    layer = ROOT / f"layer_{stem}.out"

    subprocess.run(
        [str(EXE), str(INPUT), str(baseline), str(dln), str(minthk)],
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
            "--input",
            str(INPUT),
            "--fortran",
            str(baseline),
            "--compare-out",
            str(compare),
            "--layer-out",
            str(layer),
            "--dln-vs",
            str(dln),
            "--dln-vp",
            str(dln),
            "--dln-rho",
            str(dln),
            "--minthk",
            str(minthk),
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
        "dln": dln,
        "minthk": minthk,
        "max_abs_velocity_diff": float(np.max(np.abs(vdiffs))),
        "rms_velocity_diff": float(np.sqrt(np.mean(vdiffs**2))),
        "kernel_max_abs_all": all_stats["max_abs"],
        "kernel_rms_all": all_stats["rms"],
        "kernel_rel_l2_all": all_stats["rel_l2"],
        "kernel_sign_mismatch_all": all_stats["sign_mismatch"],
        "kernel_rel_l2_vs": vs_stats["rel_l2"],
        "kernel_rel_l2_vp": vp_stats["rel_l2"],
        "kernel_rel_l2_rho": rho_stats["rel_l2"],
        "kernel_sign_mismatch_vs": vs_stats["sign_mismatch"],
        "kernel_sign_mismatch_vp": vp_stats["sign_mismatch"],
        "kernel_sign_mismatch_rho": rho_stats["sign_mismatch"],
        "baseline_file": baseline.name,
        "compare_file": compare.name,
    }


def main():
    dln_values = [0.02, 0.01, 0.005, 0.0025]
    minthk_values = [1.0, 2.0, 4.0, 8.0]
    rows = []
    for minthk in minthk_values:
        for dln in dln_values:
            print(f"running dln={dln:g} minthk={minthk:g}")
            rows.append(run_case(dln, minthk))

    out = ROOT / "sweep_summary.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
