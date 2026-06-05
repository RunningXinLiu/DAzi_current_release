#!/usr/bin/env python3
"""Run a complete DAzi iso depth-kernel A/B validation.

The runner executes the default Fortran finite-difference backend and the
hybrid pDSurf-Numba + Fortran phase-velocity backend from identical case
templates, dumps iso depth-kernel arrays and G triplets, compares them, and
emits a pass/fail summary.
"""

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
    for name in [".dazi_senk_tile", "depthkernel_dump", "g_triplet_dump"]:
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
        "kmax": r"Depth kernel(?: TI)? grid/periods: nx=[0-9]+ ny=[0-9]+ nz=[0-9]+ kmax=([0-9]+)",
        "measurements": r"Number of all measurements\s+([0-9]+)",
        "nar": r"Number of non-zero element \(nar\) of G is\s+([0-9]+)",
        "depthkernel_ti_s": r"depthkernelTI anisotropic kernel =\s+([0-9.Ee+-]+)",
        "depthkernel_backend_s": r"depthkernel backend =\s+([0-9.Ee+-]+)",
        "trace_g_loop_s": r"trace/G loop total =\s+([0-9.Ee+-]+)",
        "fmm_s": r"travel/FMM only =\s+([0-9.Ee+-]+)",
        "g_build_s": r"G build incl\. kernels/FMM/ray/G rows =\s+([0-9.Ee+-]+)",
        "solver_s": r"Solver/projection/ensemble =\s+([0-9.Ee+-]+)",
        "iteration_total_s": r"Iteration total =\s+([0-9.Ee+-]+)",
        "all_time_s": r"All time cost=\s+([0-9.Ee+-]+)s",
        "pdsurf_batch_s": r"batch_seconds=([0-9.Ee+-]+)",
        "pdsurf_wall_s": r"wall_seconds=([0-9.Ee+-]+)",
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


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def csv_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    val = row.get(key, "")
    if val == "":
        return default
    return float(val)


def evaluate(compare_dir: Path, thresholds: dict[str, float]) -> tuple[bool, list[str]]:
    failures: list[str] = []

    depth_rows = {row["array"]: row for row in read_csv_rows(compare_dir / "depthkernel_dump_diff_summary.csv")}
    for name in ["pv", "sen_vs", "sen_vp", "sen_rho", "effective_dcdvs"]:
        rel = csv_float(depth_rows[name], "rel_l2")
        max_abs = csv_float(depth_rows[name], "max_abs")
        if rel > thresholds[f"depth_{name}_rel_l2"]:
            failures.append(f"depth {name} rel_l2 {rel:.6g} > {thresholds[f'depth_{name}_rel_l2']:.6g}")
        if max_abs > thresholds[f"depth_{name}_max_abs"]:
            failures.append(f"depth {name} max_abs {max_abs:.6g} > {thresholds[f'depth_{name}_max_abs']:.6g}")

    triplet_rows = {row["set"]: row for row in read_csv_rows(compare_dir / "g_triplet_diff_summary.csv")}
    common_vs = triplet_rows["common_diff_Vs"]
    common_all = triplet_rows["common_diff"]
    if csv_float(common_vs, "rel_l2") > thresholds["g_common_vs_rel_l2"]:
        failures.append("G common_diff_Vs rel_l2 exceeds threshold")
    if csv_float(common_all, "max_abs") > thresholds["g_common_max_abs"]:
        failures.append("G common_diff max_abs exceeds threshold")

    base_only_n = int(float(triplet_rows["base_only"]["n"]))
    test_only_n = int(float(triplet_rows["test_only"]["n"]))
    base_n = int(float(triplet_rows["base_all"]["n"]))
    test_n = int(float(triplet_rows["test_all"]["n"]))
    unmatched_frac = (base_only_n + test_only_n) / max(base_n + test_n, 1)
    if unmatched_frac > thresholds["g_unmatched_frac"]:
        failures.append(f"G unmatched fraction {unmatched_frac:.6g} > {thresholds['g_unmatched_frac']:.6g}")

    return not failures, failures


def write_case_summary(out_root: Path, summary: dict[str, Any], thresholds: dict[str, float], failures: list[str]) -> None:
    compare_dir = out_root / "compare"
    default_log = summary["default"]
    hybrid_log = summary["hybrid"]
    speedup_backend = None
    if default_log.get("depthkernel_backend_s") and hybrid_log.get("depthkernel_backend_s"):
        speedup_backend = default_log["depthkernel_backend_s"] / hybrid_log["depthkernel_backend_s"]
    speedup_iter = None
    if default_log.get("iteration_total_s") and hybrid_log.get("iteration_total_s"):
        speedup_iter = default_log["iteration_total_s"] / hybrid_log["iteration_total_s"]

    payload = {
        "case": summary["case"],
        "passed": not failures,
        "failures": failures,
        "thresholds": thresholds,
        "speedup_backend": speedup_backend,
        "speedup_iteration": speedup_iter,
        "default": default_log,
        "hybrid": hybrid_log,
    }
    (compare_dir / "ab_gate_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    lines = [
        f"# ISO depth-kernel A/B: {summary['case']}",
        "",
        f"Status: {'PASS' if not failures else 'FAIL'}",
        "",
        "| metric | default | hybrid | ratio |",
        "|---|---:|---:|---:|",
    ]
    for key in ["depthkernel_backend_s", "g_build_s", "iteration_total_s", "all_time_s", "nar", "after_abs_mean_s", "after_rms_s"]:
        d = default_log.get(key)
        h = hybrid_log.get(key)
        ratio = ""
        if isinstance(d, (int, float)) and isinstance(h, (int, float)) and h:
            ratio = f"{d / h:.3f}"
        lines.append(f"| {key} | {d} | {h} | {ratio} |")
    lines.extend(["", "Failures:"])
    if failures:
        lines.extend(f"- {failure}" for failure in failures)
    else:
        lines.append("- none")
    lines.extend([
        "",
        "Artifacts:",
        f"- `{compare_dir / 'depthkernel_dump_diff_summary.csv'}`",
        f"- `{compare_dir / 'g_triplet_diff_summary.csv'}`",
        f"- `{compare_dir / 'file_diff_summary.csv'}`",
        f"- `{compare_dir / 'log_summary.csv'}`",
    ])
    (compare_dir / "AB_SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, type=Path, help="Template case directory containing para.in and data files.")
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--case-label", default=None)
    parser.add_argument("--binary", default=DEFAULT_BINARY, type=Path)
    parser.add_argument("--python", default=DEFAULT_PYTHON, type=Path)
    parser.add_argument("--omp-threads", default="10")
    parser.add_argument("--tile-columns", default="64")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-run", action="store_true", help="Only compare existing default/hybrid directories.")
    args = parser.parse_args()

    case_label = args.case_label or args.case.name
    default_dir = args.out_root / "default"
    hybrid_dir = args.out_root / "hybrid"
    compare_dir = args.out_root / "compare"
    compare_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_run:
        prepare_case(args.case, default_dir, args.force)
        prepare_case(args.case, hybrid_dir, args.force)

        common_env = {"OMP_NUM_THREADS": str(args.omp_threads)}
        run(
            [str(args.binary), "para.in"],
            cwd=default_dir,
            env={
                **common_env,
                "DAZI_DEPTHKERNEL_DUMP_DIR": str(default_dir / "depthkernel_dump"),
                "DAZI_G_TRIPLET_DUMP_DIR": str(default_dir / "g_triplet_dump"),
            },
            stdout=default_dir / "run_default.log",
        )
        run(
            [str(args.binary), "para.in"],
            cwd=hybrid_dir,
            env={
                **common_env,
                "DAZI_DEPTHKERNEL_BACKEND": "pdsurf_numba_fortran_pv",
                "DAZI_PDSURF_TILE_COLUMNS": str(args.tile_columns),
                "DAZI_PYTHON": str(args.python),
                "DAZI_DEPTHKERNEL_DUMP_DIR": str(hybrid_dir / "depthkernel_dump"),
                "DAZI_G_TRIPLET_DUMP_DIR": str(hybrid_dir / "g_triplet_dump"),
            },
            stdout=hybrid_dir / "run_hybrid.log",
        )

    default_log = parse_log(default_dir / "run_default.log")
    hybrid_log = parse_log(hybrid_dir / "run_hybrid.log")
    if not default_log["finished"]:
        raise RuntimeError(f"default run did not finish: {default_dir / 'run_default.log'}")
    if not hybrid_log["finished"]:
        raise RuntimeError(f"hybrid run did not finish: {hybrid_dir / 'run_hybrid.log'}")

    nx, ny, nz = int(default_log["nx"]), int(default_log["ny"]), int(default_log["nz"])
    kmax = int(default_log["kmax"])
    nparpi = (nx - 2) * (ny - 2) * (nz - 1)

    run(
        [
            str(args.python),
            str(LAB / "compare_depthkernel_dumps.py"),
            "--base",
            str(default_dir / "depthkernel_dump"),
            "--test",
            str(hybrid_dir / "depthkernel_dump"),
            "--nx",
            str(nx),
            "--ny",
            str(ny),
            "--nz",
            str(nz),
            "--kmax",
            str(kmax),
            "--vel-bin",
            str(hybrid_dir / ".dazi_senk_tile/vel_senK.bin"),
            "--out",
            str(compare_dir / "depthkernel_dump_diff_summary.csv"),
        ]
    )
    run(
        [
            str(args.python),
            str(LAB / "compare_g_triplet_dumps.py"),
            "--base",
            str(default_dir / "g_triplet_dump"),
            "--test",
            str(hybrid_dir / "g_triplet_dump"),
            "--nparpi",
            str(nparpi),
            "--out",
            str(compare_dir / "g_triplet_diff_summary.csv"),
        ]
    )
    run(
        [
            str(args.python),
            str(LAB / "compare_dazi_backend_outputs.py"),
            "--base",
            str(default_dir),
            "--test",
            str(hybrid_dir),
            "--base-log",
            "run_default.log",
            "--test-log",
            "run_hybrid.log",
            "--base-label",
            "default_fortran_fd",
            "--test-label",
            "pdsurf_numba_fortran_pv",
            "--outdir",
            str(compare_dir),
        ]
    )

    thresholds = {
        "depth_pv_rel_l2": 1.0e-5,
        "depth_pv_max_abs": 1.0e-4,
        "depth_sen_vs_rel_l2": 1.0e-3,
        "depth_sen_vs_max_abs": 2.0e-4,
        "depth_sen_vp_rel_l2": 2.0e-3,
        "depth_sen_vp_max_abs": 2.0e-4,
        "depth_sen_rho_rel_l2": 2.0e-3,
        "depth_sen_rho_max_abs": 2.0e-4,
        "depth_effective_dcdvs_rel_l2": 1.0e-3,
        "depth_effective_dcdvs_max_abs": 3.0e-4,
        "g_common_vs_rel_l2": 2.0e-3,
        "g_common_max_abs": 2.0e-4,
        "g_unmatched_frac": 2.0e-3,
    }
    passed, failures = evaluate(compare_dir, thresholds)
    summary = {"case": case_label, "default": default_log, "hybrid": hybrid_log}
    write_case_summary(args.out_root, summary, thresholds, failures)

    print(f"case={case_label} status={'PASS' if passed else 'FAIL'}")
    print(f"summary={compare_dir / 'AB_SUMMARY.md'}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)


if __name__ == "__main__":
    start = time.time()
    main()
    print(f"elapsed_s={time.time() - start:.3f}")
