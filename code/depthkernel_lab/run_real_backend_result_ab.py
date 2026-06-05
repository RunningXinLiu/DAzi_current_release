#!/usr/bin/env python3
"""Run real-data strict-vs-fast backend A/B and compare inversion products."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from compare_anisotropy_model import compare_anisotropy_files, json_safe


REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_CONFIGURED = Path(__file__).resolve().parent / "run_dazi_configured.py"
DEFAULT_BINARY = REPO_ROOT / "bin" / "DAzimSurfTomo"
DEFAULT_PYTHON = Path("/opt/miniconda3/envs/dispa/bin/python")

CLEAN_FILE_PATTERNS = [
    "DSurfTomo.inv",
    "Gc_Gs_model.inv",
    "period_Azm_tomo.inv",
    "period_phaseVMOD.dat",
    "phaseV_FWD.dat",
    "Traveltime_statis_*th.dat",
    "Traveltime_use_*th.dat",
    "cost_terms.dat",
    "costVSiter.dat",
    "misfit_sorted.dat",
    "misfit_voro_raw_sanitized.dat",
    "para.in_inv.log",
    "surf_tomo.log",
    "lsmr.txt",
    "IterVel.out",
    "col.txt",
    "pvpoints_acell.dat",
    "vorosites_latlon*.txt",
    "run_*.log",
    "runner_*.log",
    "metadata_*.json",
]

CLEAN_DIR_NAMES = [
    ".dazi_lsmr",
    ".dazi_senk_tile",
    ".dazi_ti_batch",
    "depthkernel_dump",
    "g_triplet_dump",
    "figs",
    "figs_srl_main",
]

GENERIC_COMPARE_FILES = [
    "DSurfTomo.inv",
    "Gc_Gs_model.inv",
    "period_Azm_tomo.inv",
    "period_phaseVMOD.dat",
    "phaseV_FWD.dat",
    "Traveltime_statis_01th.dat",
    "Traveltime_use_01th.dat",
    "Traveltime_statis_05th.dat",
    "Traveltime_use_05th.dat",
    "cost_terms.dat",
    "costVSiter.dat",
    "misfit_sorted.dat",
    "misfit_voro_raw_sanitized.dat",
]


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.6g}"
    return str(value)


def replace_numeric_prefix(line: str, value: int | float) -> str:
    newline = "\n" if line.endswith("\n") else ""
    body = line[:-1] if newline else line
    stripped = body.lstrip()
    indent = body[: len(body) - len(stripped)]
    parts = stripped.split(None, 1)
    suffix = "" if len(parts) == 1 else " " + parts[1]
    if isinstance(value, int):
        text = str(value)
    else:
        text = f"{value:g}"
    return f"{indent}{text}{suffix}{newline}"


def parse_log(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
    patterns = {
        "finished": r"Program finishes successfully",
        "measurements": r"Number of all measurements\s+([0-9]+)",
        "nar": r"Number of non-zero element \(nar\) of G is\s+([0-9]+)",
        "g_build_s": r"G build incl\. kernels/FMM/ray/G rows =\s+([0-9.Ee+-]+)",
        "solver_s": r"Solver/projection/ensemble =\s+([0-9.Ee+-]+)",
        "iteration_total_s": r"Iteration total =\s+([0-9.Ee+-]+)",
        "all_time_s": r"All time cost=\s+([0-9.Ee+-]+)s",
        "after_abs_mean_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+([0-9.Ee+-]+)",
        "after_rms_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+[0-9.Ee+-]+\s+s\s+[0-9.Ee+-]+\s+s\s+([0-9.Ee+-]+)",
        "pdsurf_wall_s": r"batch_seconds=[0-9.Ee+-]+\s+wall_seconds=([0-9.Ee+-]+)",
        "ti_batch_wall_s": r"workers=[0-9]+\s+tasks=[0-9]+\s+wall_seconds=([0-9.Ee+-]+)",
    }
    out: dict[str, Any] = {"log": str(path)}
    for key, pattern in patterns.items():
        matches = list(re.finditer(pattern, text))
        match = matches[-1] if matches else None
        if key == "finished":
            out[key] = match is not None
        elif not match:
            out[key] = None
        elif key in {"measurements", "nar"}:
            out[key] = int(match.group(1))
        else:
            out[key] = float(match.group(1))
    return out


def rewrite_para_in(
    para_path: Path,
    *,
    maxiter: int | None = None,
    nrealizations: int | None = None,
    sparse_fraction: float | None = None,
    ncell: int | None = None,
    acell: int | None = None,
    iaratio: float | None = None,
) -> dict[str, int | float]:
    lines = para_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    changed: dict[str, int | float] = {}
    for idx, line in enumerate(lines):
        lowered = line.lower()
        if maxiter is not None and "maximum of interation" in lowered:
            lines[idx] = replace_numeric_prefix(line, maxiter)
            changed["maxiter"] = maxiter
        elif sparse_fraction is not None and "sparsity fraction" in lowered:
            lines[idx] = replace_numeric_prefix(line, sparse_fraction)
            changed["sparse_fraction"] = sparse_fraction
        elif ncell is not None and "normal voro cell number" in lowered:
            lines[idx] = replace_numeric_prefix(line, ncell)
            changed["ncell"] = ncell
        elif acell is not None and "adaptive voro cell number" in lowered:
            lines[idx] = replace_numeric_prefix(line, acell)
            changed["acell"] = acell
        elif nrealizations is not None and "number of realizations" in lowered:
            lines[idx] = replace_numeric_prefix(line, nrealizations)
            changed["nrealizations"] = nrealizations
        elif iaratio is not None and "iaratio" in lowered:
            lines[idx] = replace_numeric_prefix(line, iaratio)
            changed["iaratio"] = iaratio
    para_path.write_text("".join(lines), encoding="utf-8")
    return changed


def clean_case(case_dir: Path) -> None:
    for pattern in CLEAN_FILE_PATTERNS:
        for path in case_dir.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()
    for name in CLEAN_DIR_NAMES:
        path = case_dir / name
        if path.exists():
            shutil.rmtree(path)


def prepare_case(src: Path, dst: Path, *, force: bool, para_overrides: dict[str, Any]) -> dict[str, Any]:
    if dst.exists():
        if not force:
            raise FileExistsError(f"{dst} exists; pass --force or use a new --out-root")
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    clean_case(dst)
    changed = rewrite_para_in(dst / "para.in", **para_overrides)
    (dst / "para_overrides.json").write_text(
        json.dumps(changed, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return changed


def run_profile(args: argparse.Namespace, run_dir: Path, profile: str) -> dict[str, Any]:
    log_path = run_dir / f"run_{profile}.log"
    metadata_path = run_dir / f"metadata_{profile}.json"
    runner_log_path = run_dir / f"runner_{profile}.log"
    cmd = [
        sys.executable,
        str(RUN_CONFIGURED),
        str(run_dir),
        "--profile",
        profile,
        "--threads",
        str(args.threads),
        "--omp-threads",
        str(args.omp_threads),
        "--ti-workers",
        str(args.ti_workers),
        "--pdsurf-tile-columns",
        str(args.tile_columns),
        "--ti-tile-columns",
        str(args.tile_columns),
        "--binary",
        str(args.binary),
        "--python",
        str(args.python),
        "--log",
        str(log_path),
        "--metadata",
        str(metadata_path),
    ]
    if args.voro_timer:
        cmd.extend(["--voro-timer", "1"])
    completed = subprocess.run(cmd, text=True, capture_output=True, check=False)
    runner_log_path.write_text(
        completed.stdout + completed.stderr,
        encoding="utf-8",
    )

    parsed = parse_log(log_path)
    metadata: dict[str, Any] = {}
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return {
        "profile": profile,
        "run_dir": str(run_dir),
        "command": cmd,
        "runner_returncode": completed.returncode,
        "runner_log": str(runner_log_path),
        "log": str(log_path),
        "metadata": metadata,
        "parsed_log": parsed,
    }


def load_existing_profile(run_dir: Path, profile: str) -> dict[str, Any]:
    log_path = run_dir / f"run_{profile}.log"
    metadata_path = run_dir / f"metadata_{profile}.json"
    metadata: dict[str, Any] = {}
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return {
        "profile": profile,
        "run_dir": str(run_dir),
        "command": metadata.get("command"),
        "runner_returncode": metadata.get("returncode"),
        "runner_log": str(run_dir / f"runner_{profile}.log"),
        "log": str(log_path),
        "metadata": metadata,
        "parsed_log": parse_log(log_path),
    }


def load_numeric_flat(path: Path) -> tuple[str, list[float], str]:
    if not path.exists():
        return "missing", [], ""
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        vals = []
        for token in line.split():
            try:
                vals.append(float(token))
            except ValueError:
                pass
        if vals:
            rows.append(vals)
    if not rows:
        return "empty", [], "0x0"
    width = min(len(row) for row in rows)
    shape = f"{len(rows)}x{width}"
    return "ok", [value for row in rows for value in row[:width]], shape


def compare_file_pure(base_dir: Path, test_dir: Path, name: str) -> dict[str, Any]:
    base_status, base_values, base_shape = load_numeric_flat(base_dir / name)
    test_status, test_values, test_shape = load_numeric_flat(test_dir / name)
    row: dict[str, Any] = {
        "file": name,
        "base_shape": base_shape,
        "test_shape": test_shape,
    }
    if base_status != "ok" or test_status != "ok":
        row["status"] = base_status if base_status != "ok" else test_status
        return row
    if base_shape != test_shape or len(base_values) != len(test_values):
        row["status"] = "shape_mismatch"
        return row
    if not base_values:
        row.update({"status": "empty", "max_abs": 0.0, "mean_abs": 0.0, "rms_abs": 0.0, "rel_l2": 0.0})
        return row
    diffs = [test - base for base, test in zip(base_values, test_values)]
    abs_diffs = [abs(value) for value in diffs]
    norm_base = math.sqrt(sum(value * value for value in base_values))
    norm_diff = math.sqrt(sum(value * value for value in diffs))
    row.update(
        {
            "status": "ok",
            "n": len(diffs),
            "max_abs": max(abs_diffs),
            "mean_abs": sum(abs_diffs) / len(abs_diffs),
            "rms_abs": math.sqrt(sum(value * value for value in diffs) / len(diffs)),
            "rel_l2": norm_diff / max(norm_base, 1.0e-30),
            "base_min": min(base_values),
            "base_max": max(base_values),
            "test_min": min(test_values),
            "test_max": max(test_values),
        }
    )
    return row


def compare_generic(base_dir: Path, test_dir: Path) -> list[dict[str, Any]]:
    return [compare_file_pure(base_dir, test_dir, name) for name in GENERIC_COMPARE_FILES]


def read_para_overrides(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "para_overrides.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(",".join(fields) + "\n")
        for row in rows:
            values = []
            for field in fields:
                text = str(row.get(field, ""))
                values.append('"' + text.replace('"', '""') + '"')
            handle.write(",".join(values) + "\n")


def write_summary(
    out_root: Path,
    *,
    case: Path,
    base_profile: str,
    test_profile: str,
    base_run: dict[str, Any],
    test_run: dict[str, Any],
    generic_rows: list[dict[str, Any]],
    anisotropy_report: dict[str, Any],
    para_overrides: dict[str, Any],
) -> None:
    summary = {
        "case": str(case),
        "out_root": str(out_root),
        "base_profile": base_profile,
        "test_profile": test_profile,
        "para_overrides": para_overrides,
        "base": base_run,
        "test": test_run,
        "generic_file_diffs": generic_rows,
        "anisotropy": json_safe(anisotropy_report),
    }
    (out_root / "summary_real_backend_result_ab.json").write_text(
        json.dumps(json_safe(summary), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_csv(out_root / "generic_file_diff_summary.csv", generic_rows)

    base_meta = base_run.get("metadata", {})
    test_meta = test_run.get("metadata", {})
    base_log = base_run.get("parsed_log", {})
    test_log = test_run.get("parsed_log", {})
    ani = json_safe(anisotropy_report)
    metrics = ani.get("metrics", {}) if ani.get("status") == "ok" else {}

    lines = [
        "# Real Inversion Backend Result A/B",
        "",
        f"- case: `{case}`",
        f"- base profile: `{base_profile}`",
        f"- test profile: `{test_profile}`",
        f"- output root: `{out_root}`",
        f"- para overrides: `{json.dumps(para_overrides, ensure_ascii=False)}`",
        "",
        "## Runtime And Residual",
        "",
        "| profile | runner rc | program finished | wall_s | max_rss_gb | nar | after abs mean | after RMS | all_time_s |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for run in [base_run, test_run]:
        meta = run.get("metadata", {})
        log = run.get("parsed_log", {})
        lines.append(
            f"| `{run['profile']}` | {run.get('runner_returncode')} | "
            f"{log.get('finished')} | {fmt(meta.get('wall_s'))} | "
            f"{fmt(meta.get('max_rss_gb'))} | {fmt(log.get('nar'))} | "
            f"{fmt(log.get('after_abs_mean_s'))} | {fmt(log.get('after_rms_s'))} | "
            f"{fmt(log.get('all_time_s'))} |"
        )

    speedup = None
    base_wall = base_meta.get("wall_s")
    test_wall = test_meta.get("wall_s")
    if isinstance(base_wall, (int, float)) and isinstance(test_wall, (int, float)) and test_wall > 0:
        speedup = base_wall / test_wall
    lines.extend(
        [
            "",
            f"- wall speedup `{test_profile}` vs `{base_profile}`: `{fmt(speedup)}x`",
            "",
            "## Anisotropy-Aware Metrics",
            "",
            f"- status: `{ani.get('status')}`",
            f"- Gc/Gs vector relative L2: `{fmt(ani.get('vector_rel_l2'))}`",
            "",
            "| metric | mean_abs | rms | max_abs | p95_abs | corr | unit |",
            "|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for name in [
        "vs",
        "azim_amp_file_fraction",
        "amp_from_gcgs_pct",
        "fast_axis_angle_from_gcgs_deg",
        "fast_axis_angle_from_gcgs_deg_amp_ge_threshold",
        "gc_pct",
        "gs_pct",
        "gcgs_vector_diff_pct",
    ]:
        row = metrics.get(name, {})
        lines.append(
            f"| `{name}` | {fmt(row.get('mean_abs'))} | {fmt(row.get('rms'))} | "
            f"{fmt(row.get('max_abs'))} | {fmt(row.get('p95_abs'))} | "
            f"{fmt(row.get('corr'))} | {row.get('unit', '')} |"
        )

    lines.extend(
        [
            "",
            "## Generic File Differences",
            "",
            "| file | status | max_abs | rms_abs | rel_l2 |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in generic_rows:
        lines.append(
            f"| `{row.get('file')}` | {row.get('status')} | {fmt(row.get('max_abs'))} | "
            f"{fmt(row.get('rms_abs'))} | {fmt(row.get('rel_l2'))} |"
        )

    lines.extend(
        [
            "",
            "## Run Directories",
            "",
            f"- `{base_profile}`: `{base_run['run_dir']}`",
            f"- `{test_profile}`: `{test_run['run_dir']}`",
        ]
    )
    (out_root / "summary_real_backend_result_ab.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--base-profile", default="strict_ani_voro")
    parser.add_argument("--test-profile", default="fast_tight_root")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--threads", default="4")
    parser.add_argument("--omp-threads", default="4")
    parser.add_argument("--ti-workers", default="4")
    parser.add_argument("--tile-columns", default="64")
    parser.add_argument("--maxiter", type=int, default=None)
    parser.add_argument("--nrealizations", type=int, default=None)
    parser.add_argument("--sparse-fraction", type=float, default=None)
    parser.add_argument("--ncell", type=int, default=None)
    parser.add_argument("--acell", type=int, default=None)
    parser.add_argument("--iaratio", type=float, default=None)
    parser.add_argument("--amp-threshold", type=float, default=0.001)
    parser.add_argument("--voro-timer", action="store_true")
    parser.add_argument("--compare-only", action="store_true", help="reuse existing profile directories under --out-root")
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    case = args.case.expanduser().resolve()
    out_root = args.out_root.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    para_overrides = {
        "maxiter": args.maxiter,
        "nrealizations": args.nrealizations,
        "sparse_fraction": args.sparse_fraction,
        "ncell": args.ncell,
        "acell": args.acell,
        "iaratio": args.iaratio,
    }
    para_overrides = {key: value for key, value in para_overrides.items() if value is not None}

    base_dir = out_root / args.base_profile
    test_dir = out_root / args.test_profile
    if args.compare_only:
        changed_base = read_para_overrides(base_dir)
        changed_test = read_para_overrides(test_dir)
        if changed_base != changed_test:
            raise RuntimeError(f"para overrides differ between existing runs: {changed_base} vs {changed_test}")
        base_run = load_existing_profile(base_dir, args.base_profile)
        test_run = load_existing_profile(test_dir, args.test_profile)
    else:
        changed_base = prepare_case(case, base_dir, force=args.force, para_overrides=para_overrides)
        changed_test = prepare_case(case, test_dir, force=args.force, para_overrides=para_overrides)
        if changed_base != changed_test:
            raise RuntimeError(f"para overrides differ between runs: {changed_base} vs {changed_test}")

        print(f"running {args.base_profile} in {base_dir}", flush=True)
        base_run = run_profile(args, base_dir, args.base_profile)
        print(
            f"{args.base_profile}: rc={base_run['runner_returncode']} "
            f"wall={fmt(base_run.get('metadata', {}).get('wall_s'))}s",
            flush=True,
        )

        print(f"running {args.test_profile} in {test_dir}", flush=True)
        test_run = run_profile(args, test_dir, args.test_profile)
        print(
            f"{args.test_profile}: rc={test_run['runner_returncode']} "
            f"wall={fmt(test_run.get('metadata', {}).get('wall_s'))}s",
            flush=True,
        )

    generic_rows = compare_generic(base_dir, test_dir)
    anisotropy_report = compare_anisotropy_files(base_dir, test_dir, amp_threshold=args.amp_threshold)
    write_summary(
        out_root,
        case=case,
        base_profile=args.base_profile,
        test_profile=args.test_profile,
        base_run=base_run,
        test_run=test_run,
        generic_rows=generic_rows,
        anisotropy_report=anisotropy_report,
        para_overrides=changed_base,
    )
    print(f"summary={out_root / 'summary_real_backend_result_ab.md'}")
    return max(int(base_run["runner_returncode"]), int(test_run["runner_returncode"]))


if __name__ == "__main__":
    raise SystemExit(main())
