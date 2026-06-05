#!/usr/bin/env python3
"""Run iso/ani x voro/novoro acceptance jobs with unified DAzi settings."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from compare_anisotropy_model import compare_anisotropy_files
from run_dazi_configured import DEFAULT_BINARY, DEFAULT_PYTHON, parse_config


COMBOS = ["iso_novoro", "ani_novoro", "iso_voro", "ani_voro"]
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


def parse_log(path: Path) -> dict[str, Any]:
    text = path.read_text(errors="replace") if path.exists() else ""
    patterns = {
        "finished": r"Program finishes successfully",
        "measurements": r"Number of all measurements\s+([0-9]+)",
        "nar": r"Number of non-zero element \(nar\) of G is\s+([0-9]+)",
        "solver_s": r"Solver/projection/ensemble =\s+([0-9.Ee+-]+)",
        "g_matrix_s": r"G build incl\. kernels/FMM/ray/G rows =\s+([0-9.Ee+-]+)",
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
        elif match is None:
            out[key] = None
        elif key in {"measurements", "nar"}:
            out[key] = int(match.group(1))
        else:
            out[key] = float(match.group(1))
    out["depthkernel_backend"] = first_match(text, r"Depth kernel backend:\s+(\S+)")
    out["ti_backend"] = first_match(text, r"TI depth kernel backend:\s+(\S+)")
    out["voro_gp_backend"] = first_match(text, r"Voro GP backend:\s+(\S+)")
    return out


def first_match(text: str, pattern: str) -> str | None:
    match = re.search(pattern, text)
    return match.group(1) if match else None


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


def load_numeric(path: Path) -> list[list[float]] | None:
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
    return rows or None


def compare_numeric_files(base: Path, test: Path) -> dict[str, Any]:
    a = load_numeric(base)
    b = load_numeric(test)
    if a is None or b is None:
        return {"status": "missing_or_non_numeric"}
    if len(a) != len(b) or [len(row) for row in a] != [len(row) for row in b]:
        return {
            "status": "shape_mismatch",
            "baseline_rows": len(a),
            "test_rows": len(b),
        }
    count = 0
    sum_sq = 0.0
    sum_abs = 0.0
    max_abs = 0.0
    for row_a, row_b in zip(a, b):
        for value_a, value_b in zip(row_a, row_b):
            diff = value_b - value_a
            abs_diff = abs(diff)
            max_abs = max(max_abs, abs_diff)
            sum_abs += abs_diff
            sum_sq += diff * diff
            count += 1
    return {
        "status": "ok",
        "rows": len(a),
        "values": count,
        "max_abs": max_abs,
        "rms": math.sqrt(sum_sq / count) if count else 0.0,
        "mean_abs": sum_abs / count if count else 0.0,
    }


def compare_outputs(base: Path, test: Path) -> dict[str, Any]:
    return {name: compare_numeric_files(base / name, test / name) for name in COMPARE_FILES}


def compare_anisotropy_outputs(base: Path, test: Path) -> dict[str, Any]:
    return compare_anisotropy_files(base / "Gc_Gs_model.inv", test / "Gc_Gs_model.inv")


def fmt_summary_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def run_configured(case_dir: Path, label: str, args: argparse.Namespace, optimized: bool) -> dict[str, Any]:
    wrapper = Path(__file__).with_name("run_dazi_configured.py")
    log = case_dir / f"run_{label}.log"
    meta = case_dir / f"run_{label}_metadata.json"
    cmd = [
        sys.executable,
        str(wrapper),
        str(case_dir),
        "--binary",
        str(args.binary),
        "--python",
        str(args.python),
        "--infile",
        args.infile,
        "--threads",
        str(args.threads),
        "--omp-threads",
        str(args.omp_threads),
        "--ti-workers",
        str(args.ti_workers),
        "--pdsurf-tile-columns",
        str(args.tile_columns),
        "--ti-tile-columns",
        str(args.ti_tile_columns),
        "--log",
        str(log),
        "--metadata",
        str(meta),
    ]
    if optimized:
        cmd.extend(
            [
                "--depthkernel-backend",
                "pdsurf_numba_fortran_pv",
                "--depthkernel-ti-backend",
                "process_batch",
                "--voro-gp-backend",
                "direct",
            ]
        )
    else:
        cmd.extend(
            [
                "--depthkernel-backend",
                "default",
                "--depthkernel-ti-backend",
                "default",
                "--voro-gp-backend",
                "legacy",
            ]
        )
    if args.voro_timer:
        cmd.append("--voro-timer")

    completed = subprocess.run(cmd, check=False)
    metadata = json.loads(meta.read_text(encoding="utf-8")) if meta.exists() else {}
    return {
        "returncode": completed.returncode,
        "metadata": metadata,
        "log_summary": parse_log(log),
    }


def resolve_cases(args: argparse.Namespace) -> dict[str, Path]:
    config = parse_config(args.matrix_config) if args.matrix_config is not None else {}
    raw_cases = {
        "iso_novoro": args.iso_novoro or config.get("iso_novoro"),
        "ani_novoro": args.ani_novoro or config.get("ani_novoro"),
        "iso_voro": args.iso_voro or config.get("iso_voro"),
        "ani_voro": args.ani_voro or config.get("ani_voro"),
    }
    missing = [name for name, value in raw_cases.items() if value is None]
    if missing:
        raise ValueError(f"missing case paths: {', '.join(missing)}")
    return {name: Path(value).expanduser().resolve() for name, value in raw_cases.items() if value is not None}


def apply_matrix_defaults(args: argparse.Namespace) -> None:
    if args.matrix_config is None:
        return
    config = parse_config(args.matrix_config)
    keys = {
        "threads": ["threads"],
        "omp_threads": ["omp_threads"],
        "ti_workers": ["ti_workers"],
        "tile_columns": ["tile_columns", "pdsurf_tile_columns"],
        "ti_tile_columns": ["ti_tile_columns"],
        "infile": ["infile"],
    }
    for attr, names in keys.items():
        if getattr(args, attr) is not None:
            continue
        for name in names:
            if name in config:
                setattr(args, attr, config[name])
                break


def write_summary(out_root: Path, summary: dict[str, Any]) -> None:
    (out_root / "summary_four_combo_acceptance.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# Four-Combo DAzi Acceptance Matrix",
        "",
        f"- started: `{summary['started_at']}`",
        f"- optimized only: `{not summary['with_baseline']}`",
        f"- max RSS gate: `{summary.get('max_rss_gb')}` GB",
        "",
        "| combo | mode | pass | return | finished | wall_s | max_rss_gb | nar | G_s | solver_s |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for combo, combo_data in summary["combos"].items():
        for mode in ("baseline", "optimized"):
            run_data = combo_data.get(mode)
            if not run_data:
                continue
            meta = run_data.get("metadata", {})
            log = run_data.get("log_summary", {})
            lines.append(
                f"| `{combo}` | `{mode}` | {run_data.get('pass')} | "
                f"{run_data.get('returncode')} | {log.get('finished')} | "
                f"{float(meta.get('wall_s', 0.0)):.3f} | {float(meta.get('max_rss_gb', 0.0)):.3f} | "
                f"{log.get('nar')} | {log.get('g_matrix_s')} | {log.get('solver_s')} |"
            )
    if summary["with_baseline"]:
        lines.extend(["", "## Baseline vs Optimized Differences", ""])
        lines.append("| combo | file | status | max_abs | rms | mean_abs |")
        lines.append("|---|---|---:|---:|---:|---:|")
        for combo, combo_data in summary["combos"].items():
            for name, diff in combo_data.get("diffs", {}).items():
                lines.append(
                    f"| `{combo}` | `{name}` | {diff.get('status')} | "
                    f"{diff.get('max_abs', '')} | {diff.get('rms', '')} | {diff.get('mean_abs', '')} |"
                )
        lines.extend(["", "## Anisotropy-Aware Differences", ""])
        lines.append(
            "| combo | status | vector_rel_l2 | angle file max/p95 | "
            "angle GcGs max/p95 | Gc max/rms | Gs max/rms | amp pct max/rms |"
        )
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
        for combo, combo_data in summary["combos"].items():
            report = combo_data.get("anisotropy_diffs", {})
            if report.get("status") != "ok":
                lines.append(f"| `{combo}` | {report.get('status')} |  |  |  |  |  |  |")
                continue
            metrics = report.get("metrics", {})
            angle_file = metrics.get("azim_angle_file_deg", {})
            angle_gcgs = metrics.get("fast_axis_angle_from_gcgs_deg", {})
            gc = metrics.get("gc_pct", {})
            gs = metrics.get("gs_pct", {})
            amp = metrics.get("amp_from_gcgs_pct", {})
            lines.append(
                f"| `{combo}` | ok | {fmt_summary_value(report.get('vector_rel_l2'))} | "
                f"{fmt_summary_value(angle_file.get('max_abs'))}/{fmt_summary_value(angle_file.get('p95_abs'))} | "
                f"{fmt_summary_value(angle_gcgs.get('max_abs'))}/{fmt_summary_value(angle_gcgs.get('p95_abs'))} | "
                f"{fmt_summary_value(gc.get('max_abs'))}/{fmt_summary_value(gc.get('rms'))} | "
                f"{fmt_summary_value(gs.get('max_abs'))}/{fmt_summary_value(gs.get('rms'))} | "
                f"{fmt_summary_value(amp.get('max_abs'))}/{fmt_summary_value(amp.get('rms'))} |"
            )
    (out_root / "summary_four_combo_acceptance.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def pass_gate(run_data: dict[str, Any], max_rss_gb: float | None) -> bool:
    meta = run_data.get("metadata", {})
    log = run_data.get("log_summary", {})
    if run_data.get("returncode") != 0:
        return False
    if not log.get("finished"):
        return False
    if max_rss_gb is not None and float(meta.get("max_rss_gb", 0.0)) > max_rss_gb:
        return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Run iso/ani x voro/novoro DAzi acceptance matrix")
    parser.add_argument("--matrix-config", type=Path, default=None)
    parser.add_argument("--iso-novoro", default=None)
    parser.add_argument("--ani-novoro", default=None)
    parser.add_argument("--iso-voro", default=None)
    parser.add_argument("--ani-voro", default=None)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--infile", default=None)
    parser.add_argument("--threads", default=None)
    parser.add_argument("--omp-threads", default=None)
    parser.add_argument("--ti-workers", default=None)
    parser.add_argument("--tile-columns", default=None)
    parser.add_argument("--ti-tile-columns", default=None)
    parser.add_argument("--max-rss-gb", type=float, default=None)
    parser.add_argument("--with-baseline", action="store_true")
    parser.add_argument("--voro-timer", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    apply_matrix_defaults(args)
    args.infile = args.infile or "para.in"
    args.threads = args.threads or "4"
    args.omp_threads = args.omp_threads or args.threads
    args.ti_workers = args.ti_workers or args.threads
    args.tile_columns = args.tile_columns or "64"
    args.ti_tile_columns = args.ti_tile_columns or args.tile_columns

    cases = resolve_cases(args)
    out_root = args.out_root.expanduser().resolve()
    if args.dry_run:
        print(
            json.dumps(
                {
                    "cases": {name: str(path) for name, path in cases.items()},
                    "out_root": str(out_root),
                    "with_baseline": args.with_baseline,
                    "settings": {
                        "threads": args.threads,
                        "omp_threads": args.omp_threads,
                        "ti_workers": args.ti_workers,
                        "tile_columns": args.tile_columns,
                        "ti_tile_columns": args.ti_tile_columns,
                        "infile": args.infile,
                    },
                },
                indent=2,
                ensure_ascii=False,
            )
        )
        return
    out_root.mkdir(parents=True, exist_ok=True)

    summary: dict[str, Any] = {
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "with_baseline": args.with_baseline,
        "max_rss_gb": args.max_rss_gb,
        "settings": {
            "threads": args.threads,
            "omp_threads": args.omp_threads,
            "ti_workers": args.ti_workers,
            "tile_columns": args.tile_columns,
            "ti_tile_columns": args.ti_tile_columns,
        },
        "combos": {},
    }

    for combo in COMBOS:
        src = cases[combo]
        combo_dir = out_root / combo
        combo_summary: dict[str, Any] = {"source": str(src)}
        if args.with_baseline:
            baseline_dir = combo_dir / "baseline"
            prepare_case(src, baseline_dir, args.force)
            baseline = run_configured(baseline_dir, "baseline", args, optimized=False)
            baseline["pass"] = pass_gate(baseline, args.max_rss_gb)
            combo_summary["baseline"] = baseline
        optimized_dir = combo_dir / "optimized"
        prepare_case(src, optimized_dir, args.force)
        optimized = run_configured(optimized_dir, "optimized", args, optimized=True)
        optimized["pass"] = pass_gate(optimized, args.max_rss_gb)
        combo_summary["optimized"] = optimized
        if args.with_baseline:
            combo_summary["diffs"] = compare_outputs(combo_dir / "baseline", combo_dir / "optimized")
            combo_summary["anisotropy_diffs"] = compare_anisotropy_outputs(
                combo_dir / "baseline",
                combo_dir / "optimized",
            )
        summary["combos"][combo] = combo_summary
        write_summary(out_root, summary)

    write_summary(out_root, summary)
    print(f"summary={out_root / 'summary_four_combo_acceptance.md'}")


if __name__ == "__main__":
    main()
