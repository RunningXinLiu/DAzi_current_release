#!/usr/bin/env python3
"""Run DAzi TI depth-kernel A/B validation."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any


ROOT = Path("/Users/liuxin/Desktop/DSurf_test")
LAB = ROOT / "depthkernel_lab"
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
]


def run(cmd: list[str], cwd: Path | None = None, env: dict[str, str] | None = None, stdout: Path | None = None) -> None:
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    if stdout:
        with stdout.open("w") as f:
            subprocess.run(cmd, cwd=cwd, env=merged_env, stdout=f, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(cmd, cwd=cwd, env=merged_env, check=True)


def prepare_case(src: Path, dst: Path, force: bool) -> None:
    if dst.exists():
        if not force:
            raise FileExistsError(f"{dst} exists; pass --force or choose a new --out-root")
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    for name in OUTPUT_FILES:
        path = dst / name
        if path.exists():
            path.unlink()
    for name in [".dazi_senk_tile", ".dazi_ti_batch", "ti_dump", "g_triplet_dump"]:
        path = dst / name
        if path.exists():
            shutil.rmtree(path)


def parse_log(path: Path) -> dict[str, Any]:
    text = path.read_text(errors="replace")
    out: dict[str, Any] = {
        "finished": "Program finishes successfully" in text,
        "text_path": str(path),
    }
    patterns = {
        "nx": r"model dimension:nx,ny,nz\s+([0-9]+)\s+[0-9]+\s+[0-9]+",
        "ny": r"model dimension:nx,ny,nz\s+[0-9]+\s+([0-9]+)\s+[0-9]+",
        "nz": r"model dimension:nx,ny,nz\s+[0-9]+\s+[0-9]+\s+([0-9]+)",
        "kmax": r"Depth kernel TI grid/periods: nx=[0-9]+ ny=[0-9]+ nz=[0-9]+ kmax=([0-9]+)",
        "measurements": r"Number of all measurements\s+([0-9]+)",
        "nar": r"Number of non-zero element \(nar\) of G is\s+([0-9]+)",
        "depthkernel_ti_s": r"depthkernelTI anisotropic kernel =\s+([0-9.Ee+-]+)",
        "depthkernel_backend_s": r"depthkernel backend =\s+([0-9.Ee+-]+)",
        "g_build_s": r"G build incl\. kernels/FMM/ray/G rows =\s+([0-9.Ee+-]+)",
        "iteration_total_s": r"Iteration total =\s+([0-9.Ee+-]+)",
        "all_time_s": r"All time cost=\s+([0-9.Ee+-]+)s",
        "ti_batch_wall_s": r"workers=[0-9]+ tasks=[0-9]+ wall_seconds=([0-9.Ee+-]+)",
        "after_abs_mean_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+([0-9.Ee+-]+)",
        "after_rms_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+[0-9.Ee+-]+\s+s\s+[0-9.Ee+-]+\s+s\s+([0-9.Ee+-]+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if not match:
            out[key] = None
            continue
        value = match.group(1)
        out[key] = int(value) if key in {"nx", "ny", "nz", "kmax", "measurements", "nar"} else float(value)
    return out


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def as_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def evaluate(compare_dir: Path, thresholds: dict[str, float]) -> tuple[bool, list[str]]:
    failures: list[str] = []
    ti_rows = {row["array"]: row for row in read_rows(compare_dir / "ti_depthkernel_diff_summary.csv")}
    for name in ["pv_ti", "lsen_gsc_ti"]:
        rel = as_float(ti_rows[name], "rel_l2")
        max_abs = as_float(ti_rows[name], "max_abs")
        if rel > thresholds[f"{name}_rel_l2"]:
            failures.append(f"{name} rel_l2 {rel:.6g} > {thresholds[f'{name}_rel_l2']:.6g}")
        if max_abs > thresholds[f"{name}_max_abs"]:
            failures.append(f"{name} max_abs {max_abs:.6g} > {thresholds[f'{name}_max_abs']:.6g}")

    if (compare_dir / "g_triplet_diff_summary.csv").exists():
        triplet_rows = {row["set"]: row for row in read_rows(compare_dir / "g_triplet_diff_summary.csv")}
        for set_name in ["common_diff_Gc", "common_diff_Gs"]:
            if as_float(triplet_rows[set_name], "rel_l2") > thresholds["g_common_aniso_rel_l2"]:
                failures.append(f"{set_name} rel_l2 exceeds threshold")
        if as_float(triplet_rows["common_diff"], "max_abs") > thresholds["g_common_max_abs"]:
            failures.append("G common_diff max_abs exceeds threshold")
    return not failures, failures


def write_summary(out_root: Path, label: str, baseline: dict[str, Any], batch: dict[str, Any], thresholds: dict[str, float], failures: list[str]) -> None:
    compare_dir = out_root / "compare"
    speedup_ti = None
    if baseline.get("depthkernel_ti_s") and batch.get("depthkernel_ti_s"):
        speedup_ti = baseline["depthkernel_ti_s"] / batch["depthkernel_ti_s"]
    payload = {
        "case": label,
        "passed": not failures,
        "failures": failures,
        "thresholds": thresholds,
        "speedup_ti": speedup_ti,
        "baseline": baseline,
        "process_batch": batch,
    }
    (compare_dir / "ti_ab_gate_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    lines = [
        f"# TI depth-kernel A/B: {label}",
        "",
        f"Status: {'PASS' if not failures else 'FAIL'}",
        "",
        "| metric | baseline | process_batch | ratio |",
        "|---|---:|---:|---:|",
    ]
    for key in ["depthkernel_ti_s", "depthkernel_backend_s", "g_build_s", "iteration_total_s", "all_time_s", "nar", "after_abs_mean_s", "after_rms_s"]:
        b = baseline.get(key)
        p = batch.get(key)
        ratio = ""
        if isinstance(b, (int, float)) and isinstance(p, (int, float)) and p:
            ratio = f"{b / p:.3f}"
        lines.append(f"| {key} | {b} | {p} | {ratio} |")
    lines.extend(["", "Failures:"])
    if failures:
        lines.extend(f"- {failure}" for failure in failures)
    else:
        lines.append("- none")
    lines.extend([
        "",
        "Artifacts:",
        f"- `{compare_dir / 'ti_depthkernel_diff_summary.csv'}`",
        f"- `{compare_dir / 'g_triplet_diff_summary.csv'}`",
        f"- `{compare_dir / 'file_diff_summary.csv'}`",
        f"- `{compare_dir / 'log_summary.csv'}`",
    ])
    (compare_dir / "TI_AB_SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--case-label", default=None)
    parser.add_argument("--binary", default=DEFAULT_BINARY, type=Path)
    parser.add_argument("--python", default=DEFAULT_PYTHON, type=Path)
    parser.add_argument("--omp-threads", default="10")
    parser.add_argument("--iso-tile-columns", default="64")
    parser.add_argument("--ti-workers", default="4")
    parser.add_argument("--ti-tile-columns", default="64")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    args = parser.parse_args()

    label = args.case_label or args.case.name
    baseline_dir = args.out_root / "baseline"
    batch_dir = args.out_root / "process_batch"
    compare_dir = args.out_root / "compare"
    compare_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_run:
        prepare_case(args.case, baseline_dir, args.force)
        prepare_case(args.case, batch_dir, args.force)
        common_env = {
            "OMP_NUM_THREADS": str(args.omp_threads),
            "DAZI_DEPTHKERNEL_BACKEND": "pdsurf_numba_fortran_pv",
            "DAZI_PDSURF_TILE_COLUMNS": str(args.iso_tile_columns),
            "DAZI_PYTHON": str(args.python),
        }
        run(
            [str(args.binary), "para.in"],
            cwd=baseline_dir,
            env={
                **common_env,
                "DAZI_DEPTHKERNEL_TI_DUMP_DIR": str(baseline_dir / "ti_dump"),
                "DAZI_G_TRIPLET_DUMP_DIR": str(baseline_dir / "g_triplet_dump"),
            },
            stdout=baseline_dir / "run_baseline.log",
        )
        run(
            [str(args.binary), "para.in"],
            cwd=batch_dir,
            env={
                **common_env,
                "DAZI_DEPTHKERNEL_TI_BACKEND": "process_batch",
                "DAZI_DEPTHKERNEL_TI_DUMP_DIR": str(batch_dir / "ti_dump"),
                "DAZI_DEPTHKERNEL_TI_BATCH_DIR": str(batch_dir / ".dazi_ti_batch"),
                "DAZI_DEPTHKERNEL_TI_WORKERS": str(args.ti_workers),
                "DAZI_DEPTHKERNEL_TI_TILE_COLUMNS": str(args.ti_tile_columns),
                "DAZI_G_TRIPLET_DUMP_DIR": str(batch_dir / "g_triplet_dump"),
            },
            stdout=batch_dir / "run_process_batch.log",
        )

    baseline_log = parse_log(baseline_dir / "run_baseline.log")
    batch_log = parse_log(batch_dir / "run_process_batch.log")
    if not baseline_log["finished"]:
        raise RuntimeError(f"baseline run did not finish: {baseline_dir / 'run_baseline.log'}")
    if not batch_log["finished"]:
        raise RuntimeError(f"process batch run did not finish: {batch_dir / 'run_process_batch.log'}")

    nx, ny, nz = int(baseline_log["nx"]), int(baseline_log["ny"]), int(baseline_log["nz"])
    kmax = int(baseline_log["kmax"])
    nparpi = (nx - 2) * (ny - 2) * (nz - 1)

    run([
        str(args.python),
        str(LAB / "compare_depthkernel_ti_dumps.py"),
        "--base",
        str(baseline_dir / "ti_dump"),
        "--test",
        str(batch_dir / "ti_dump"),
        "--nx",
        str(nx),
        "--ny",
        str(ny),
        "--nz",
        str(nz),
        "--kmax",
        str(kmax),
        "--out",
        str(compare_dir / "ti_depthkernel_diff_summary.csv"),
    ])
    run([
        str(args.python),
        str(LAB / "compare_g_triplet_dumps.py"),
        "--base",
        str(baseline_dir / "g_triplet_dump"),
        "--test",
        str(batch_dir / "g_triplet_dump"),
        "--nparpi",
        str(nparpi),
        "--out",
        str(compare_dir / "g_triplet_diff_summary.csv"),
    ])
    run([
        str(args.python),
        str(LAB / "compare_dazi_backend_outputs.py"),
        "--base",
        str(baseline_dir),
        "--test",
        str(batch_dir),
        "--base-log",
        "run_baseline.log",
        "--test-log",
        "run_process_batch.log",
        "--base-label",
        "default_fortran_tregn96",
        "--test-label",
        "process_batch",
        "--outdir",
        str(compare_dir),
    ])

    thresholds = {
        "pv_ti_rel_l2": 1.0e-6,
        "pv_ti_max_abs": 1.0e-5,
        "lsen_gsc_ti_rel_l2": 1.0e-5,
        "lsen_gsc_ti_max_abs": 1.0e-4,
        "g_common_aniso_rel_l2": 1.0e-4,
        "g_common_max_abs": 1.0e-4,
    }
    passed, failures = evaluate(compare_dir, thresholds)
    write_summary(args.out_root, label, baseline_log, batch_log, thresholds, failures)
    print(f"case={label} status={'PASS' if passed else 'FAIL'}")
    print(f"summary={compare_dir / 'TI_AB_SUMMARY.md'}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)


if __name__ == "__main__":
    start = time.time()
    main()
    print(f"elapsed_s={time.time() - start:.3f}")
