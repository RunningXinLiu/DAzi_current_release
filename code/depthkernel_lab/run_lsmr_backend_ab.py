#!/usr/bin/env python3
"""Run a small A/B check for the optional SciPy LSMR backend."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path("/Users/liuxin/Desktop/DSurf_test")
DEFAULT_BINARY = ROOT / "DAzi_large_data_260128/bin/DAzimSurfTomo"
DEFAULT_PYTHON = Path("/opt/miniconda3/envs/dispa/bin/python")

OUTPUT_FILES = [
    "DSurfTomo.inv",
    "Gc_Gs_model.inv",
    "period_Azm_tomo.inv",
    "period_phaseVMOD.dat",
    "phaseV_FWD.dat",
    "Traveltime_statis_01th.dat",
    "Traveltime_statis_05th.dat",
    "Traveltime_use_01th.dat",
    "Traveltime_use_05th.dat",
    "cost_terms.dat",
    "costVSiter.dat",
    "misfit_sorted.dat",
    "misfit_voro_raw_sanitized.dat",
    "surf_tomo.log",
    "para.in_inv.log",
    "lsmr.txt",
]

CLEAN_DIRS = [
    ".dazi_lsmr",
    ".dazi_senk_tile",
    ".dazi_ti_batch",
    "depthkernel_dump",
    "g_triplet_dump",
    "figs",
    "figs_srl_main",
]

COMPARE_FILES = [
    "DSurfTomo.inv",
    "Gc_Gs_model.inv",
    "phaseV_FWD.dat",
    "cost_terms.dat",
    "misfit_sorted.dat",
    "misfit_voro_raw_sanitized.dat",
]


def run(
    cmd: list[str],
    cwd: Path,
    env: dict[str, str],
    stdout: Path,
) -> float:
    merged_env = os.environ.copy()
    merged_env.update(env)
    t0 = time.perf_counter()
    with stdout.open("w") as f:
        subprocess.run(cmd, cwd=cwd, env=merged_env, stdout=f, stderr=subprocess.STDOUT, check=True)
    return time.perf_counter() - t0


def prepare_case(src: Path, dst: Path, force: bool) -> None:
    if dst.exists():
        if not force:
            raise FileExistsError(f"{dst} exists; pass --force or use a new --out-root")
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    for path in dst.glob("run_*.log"):
        path.unlink()
    for name in OUTPUT_FILES:
        path = dst / name
        if path.exists():
            path.unlink()
    for name in CLEAN_DIRS:
        path = dst / name
        if path.exists():
            shutil.rmtree(path)


def parse_log(path: Path) -> dict[str, Any]:
    text = path.read_text(errors="replace") if path.exists() else ""
    patterns = {
        "finished": r"Program finishes successfully",
        "measurements": r"Number of all measurements\s+([0-9]+)",
        "nar": r"Number of non-zero element \(nar\) of G is\s+([0-9]+)",
        "solver_s": r"Solver/projection/ensemble =\s+([0-9.Ee+-]+)",
        "iteration_total_s": r"Iteration total =\s+([0-9.Ee+-]+)",
        "all_time_s": r"All time cost=\s+([0-9.Ee+-]+)s",
        "after_abs_mean_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+([0-9.Ee+-]+)",
        "after_rms_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+[0-9.Ee+-]+\s+s\s+[0-9.Ee+-]+\s+s\s+([0-9.Ee+-]+)",
    }
    out: dict[str, Any] = {"log": str(path)}
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if key == "finished":
            out[key] = match is not None
        elif not match:
            out[key] = None
        elif key in {"measurements", "nar"}:
            out[key] = int(match.group(1))
        else:
            out[key] = float(match.group(1))
    out["scipy_backend_calls"] = len(re.findall(r"LSMR backend: scipy", text))
    timer_rows = []
    timer_pattern = re.compile(
        r"Voro timer\s+(\w+)\s+iter=\s*([0-9]+)\s+real=\s*([0-9]+)"
        r"\s+cells=\s*([0-9.Ee+-]+)s\s+map=\s*([0-9.Ee+-]+)s"
        r"\s+gp=\s*([0-9.Ee+-]+)s\s+lsmr=\s*([0-9.Ee+-]+)s"
        r"\s+back=\s*([0-9.Ee+-]+)s"
    )
    for match in timer_pattern.finditer(text):
        timer_rows.append(
            {
                "kind": match.group(1),
                "iter": int(match.group(2)),
                "real": int(match.group(3)),
                "cells_s": float(match.group(4)),
                "map_s": float(match.group(5)),
                "gp_s": float(match.group(6)),
                "lsmr_s": float(match.group(7)),
                "back_s": float(match.group(8)),
            }
        )
    out["voro_timers"] = timer_rows
    return out


def collect_scipy_meta(run_dir: Path) -> dict[str, Any]:
    meta_files = sorted((run_dir / ".dazi_lsmr").glob("*/meta_lsmr.json"))
    metas = []
    for path in meta_files:
        try:
            metas.append(json.loads(path.read_text()))
        except json.JSONDecodeError:
            continue
    return {
        "calls": len(metas),
        "total_seconds": sum(float(item.get("total_seconds", 0.0)) for item in metas),
        "solve_seconds": sum(float(item.get("solve_seconds", 0.0)) for item in metas),
        "build_seconds": sum(float(item.get("build_seconds", 0.0)) for item in metas),
        "itn": [int(item.get("itn", 0)) for item in metas],
        "istop": [int(item.get("istop", 0)) for item in metas],
    }


def load_numeric(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None
    rows: list[list[float]] = []
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split()
        if not parts:
            continue
        try:
            rows.append([float(part) for part in parts])
        except ValueError:
            continue
    if not rows:
        return None
    width = max(len(row) for row in rows)
    if any(len(row) != width for row in rows):
        return np.asarray([value for row in rows for value in row], dtype=np.float64)
    return np.asarray(rows, dtype=np.float64)


def compare_outputs(baseline: Path, scipy: Path) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name in COMPARE_FILES:
        a = load_numeric(baseline / name)
        b = load_numeric(scipy / name)
        if a is None or b is None:
            out[name] = {"status": "missing_or_non_numeric"}
            continue
        if a.shape != b.shape:
            out[name] = {"status": "shape_mismatch", "baseline_shape": a.shape, "scipy_shape": b.shape}
            continue
        diff = b - a
        out[name] = {
            "status": "ok",
            "shape": a.shape,
            "max_abs": float(np.max(np.abs(diff))) if diff.size else 0.0,
            "rms": float(np.sqrt(np.mean(diff * diff))) if diff.size else 0.0,
            "mean_abs": float(np.mean(np.abs(diff))) if diff.size else 0.0,
        }
    return out


def write_summary(out_root: Path, summary: dict[str, Any]) -> None:
    (out_root / "summary_lsmr_backend_ab.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    base = summary["baseline"]
    scipy = summary["scipy"]
    lines = [
        "# LSMR Backend A/B",
        "",
        f"- case: `{summary['case']}`",
        f"- backend scope: `{summary['backend_scope']}`",
        f"- baseline wall: `{base.get('wall_s'):.3f} s`",
        f"- scipy wall: `{scipy.get('wall_s'):.3f} s`",
        f"- baseline solver/projection: `{base.get('solver_s')} s`",
        f"- scipy solver/projection: `{scipy.get('solver_s')} s`",
        f"- SciPy backend calls: `{summary['scipy_meta']['calls']}`",
        f"- SciPy backend total subprocess time: `{summary['scipy_meta']['total_seconds']:.3f} s`",
        "",
        "## Output Differences",
        "",
        "| file | status | max_abs | rms | mean_abs |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, row in summary["diffs"].items():
        lines.append(
            f"| `{name}` | {row.get('status')} | "
            f"{row.get('max_abs', '')} | {row.get('rms', '')} | {row.get('mean_abs', '')} |"
        )
    if base.get("voro_timers") or scipy.get("voro_timers"):
        lines.extend(["", "## Vorotomo Timer Rows", ""])
        lines.append("| run | kind | iter | real | cells | map | gp | lsmr | back |")
        lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for label, run_data in (("baseline", base), ("scipy", scipy)):
            for row in run_data.get("voro_timers", []):
                lines.append(
                    f"| {label} | {row['kind']} | {row['iter']} | {row['real']} | "
                    f"{row['cells_s']:.4f} | {row['map_s']:.4f} | {row['gp_s']:.4f} | "
                    f"{row['lsmr_s']:.4f} | {row['back_s']:.4f} |"
                )
    (out_root / "summary_lsmr_backend_ab.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="A/B test DAzi LSMR backend")
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--case-label", default=None)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--omp-threads", default="10")
    parser.add_argument("--ti-workers", default="4")
    parser.add_argument("--tile-columns", default="64")
    parser.add_argument("--backend-scope", choices=("voro", "main", "all"), default="voro")
    parser.add_argument("--voro-timer", action="store_true", help="Enable DAZI_VORO_TIMER=1 in both runs.")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.out_root.mkdir(parents=True, exist_ok=True)
    baseline_dir = args.out_root / "baseline_fortran"
    scipy_dir = args.out_root / "scipy_lsmr"
    prepare_case(args.case, baseline_dir, args.force)
    prepare_case(args.case, scipy_dir, args.force)

    common_env = {
        "OMP_NUM_THREADS": str(args.omp_threads),
        "DAZI_PYTHON": str(args.python),
        "DAZI_DEPTHKERNEL_BACKEND": "pdsurf_numba_fortran_pv",
        "DAZI_DEPTHKERNEL_TI_BACKEND": "process_batch",
        "DAZI_DEPTHKERNEL_TI_WORKERS": str(args.ti_workers),
        "DAZI_PDSURF_TILE_COLUMNS": str(args.tile_columns),
        "DAZI_DEPTHKERNEL_TI_TILE_COLUMNS": str(args.tile_columns),
    }
    if args.voro_timer:
        common_env["DAZI_VORO_TIMER"] = "1"

    baseline_wall = run(
        [str(args.binary), "para.in"],
        cwd=baseline_dir,
        env={**common_env},
        stdout=baseline_dir / "run_lsmr_baseline.log",
    )

    scipy_env = {
        **common_env,
        "DAZI_LSMR_DIR": str(scipy_dir / ".dazi_lsmr"),
    }
    if args.backend_scope == "voro":
        scipy_env["DAZI_VORO_LSMR_BACKEND"] = "scipy"
    elif args.backend_scope == "main":
        scipy_env["DAZI_MAIN_LSMR_BACKEND"] = "scipy"
    else:
        scipy_env["DAZI_LSMR_BACKEND"] = "scipy"

    scipy_wall = run(
        [str(args.binary), "para.in"],
        cwd=scipy_dir,
        env=scipy_env,
        stdout=scipy_dir / "run_lsmr_scipy.log",
    )

    baseline = parse_log(baseline_dir / "run_lsmr_baseline.log")
    scipy = parse_log(scipy_dir / "run_lsmr_scipy.log")
    baseline["wall_s"] = baseline_wall
    scipy["wall_s"] = scipy_wall
    summary = {
        "case": args.case_label or args.case.name,
        "backend_scope": args.backend_scope,
        "baseline": baseline,
        "scipy": scipy,
        "scipy_meta": collect_scipy_meta(scipy_dir),
        "diffs": compare_outputs(baseline_dir, scipy_dir),
    }
    write_summary(args.out_root, summary)

    print(f"case={summary['case']}")
    print(f"baseline_wall={baseline_wall:.3f}s scipy_wall={scipy_wall:.3f}s")
    print(f"baseline_solver={baseline.get('solver_s')} scipy_solver={scipy.get('solver_s')}")
    print(f"scipy_calls={summary['scipy_meta']['calls']}")
    print(f"summary={args.out_root / 'summary_lsmr_backend_ab.md'}")


if __name__ == "__main__":
    main()
