#!/usr/bin/env python3
"""Preflight checks for DAzi/Vorotomo case size and memory risk.

The estimates are intentionally conservative upper bounds.  They are meant to
catch obviously dangerous parameter combinations before a long server run, not
to reproduce the exact resident set size of every compiler/runtime.
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union


REAL4 = 4
REAL8 = 8
INT4 = 4
UNSET_VALUES = {"", "default", "unset", "none"}


def gib(value: Union[float, int]) -> float:
    return float(value) / (1024.0**3)


def is_fortran_comment_line(line: str) -> bool:
    lowered = line.lower()
    if re.fullmatch(r"c+", lowered):
        return True
    if lowered.startswith("cccc"):
        return True
    return lowered.startswith("c ") or lowered.startswith("c\t") or lowered.startswith("c:")


def active_lines(path: Path) -> List[str]:
    lines: List[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if is_fortran_comment_line(line):
            continue
        lines.append(line)
    return lines


def values_before_comment(line: str) -> List[str]:
    return re.split(r"\bc\s*:", line, maxsplit=1, flags=re.IGNORECASE)[0].split()


def collect_tokens(lines: List[str], start: int, count: int) -> Tuple[List[str], int]:
    values: List[str] = []
    index = start
    while len(values) < count and index < len(lines):
        values.extend(values_before_comment(lines[index]))
        index += 1
    if len(values) < count:
        raise ValueError(f"expected {count} tokens from para.in line {start + 1}, got {len(values)}")
    return values[:count], index


def parse_para_runtime_controls(lines: List[str]) -> Dict[str, str]:
    values: Dict[str, str] = {}
    label_map = {
        "profile": "profile",
        "gp_spfra": "gp_spfra",
        "ray_backend": "ray_backend",
        "lsmr_backend": "lsmr_backend",
        "threads": "threads",
        "pdsurf_tile_columns": "pdsurf_tile_columns",
        "ti_tile_columns": "ti_tile_columns",
        "voro_timer": "voro_timer",
    }
    for line in lines:
        if "c:" not in line.lower():
            continue
        value_part, comment = re.split(r"\bc\s*:", line, maxsplit=1, flags=re.IGNORECASE)
        tokens = value_part.split()
        if not tokens:
            continue
        key = comment.split(",", 1)[0].split("(", 1)[0].strip().lower().replace("-", "_")
        for label, mapped in label_map.items():
            if key == label or key.startswith(label):
                values[mapped] = tokens[0]
                break
    return values


def parse_para(path: Path) -> Dict[str, Any]:
    lines = active_lines(path)
    if len(lines) < 18:
        raise ValueError(f"{path} has too few active lines")
    data_file = values_before_comment(lines[0])[0]
    nx, ny, nz = [int(value) for value in values_before_comment(lines[1])[:3]]
    max_source_receiver = int(values_before_comment(lines[6])[0])
    spfra = float(values_before_comment(lines[7])[0])
    maxiter = int(values_before_comment(lines[8])[0])
    iso_mode = values_before_comment(lines[9])[0].strip().upper().startswith("T")
    kmax = int(values_before_comment(lines[13])[0])
    periods, base = collect_tokens(lines, 14, kmax)
    if base + 2 >= len(lines):
        raise ValueError(f"{path} ends before output/Vorotomo controls")
    output_intermediate = int(values_before_comment(lines[base])[0])
    output_raypaths = int(values_before_comment(lines[base + 1])[0])
    voro = int(values_before_comment(lines[base + 2])[0])
    runtime_start = base + 3
    ncell = acell = 0
    nrealizations = 1
    nzrand = 0
    if voro:
        if base + 6 >= len(lines):
            raise ValueError(f"{path} has voronoi mode but incomplete vorotomo controls")
        ncell = int(values_before_comment(lines[base + 3])[0])
        acell = int(values_before_comment(lines[base + 4])[0])
        nrealizations = int(values_before_comment(lines[base + 5])[0])
        nzrand = int(values_before_comment(lines[base + 6])[0])
        runtime_start = base + 7
    runtime = parse_para_runtime_controls(lines[runtime_start:])
    return {
        "data_file": data_file,
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "max_source_receiver": max_source_receiver,
        "spfra": spfra,
        "maxiter": maxiter,
        "iso_mode": iso_mode,
        "kmax": kmax,
        "period_count": len(periods),
        "period_min": float(periods[0]) if periods else None,
        "period_max": float(periods[-1]) if periods else None,
        "output_intermediate": output_intermediate,
        "output_raypaths": output_raypaths,
        "voro": bool(voro),
        "ncell": ncell,
        "acell": acell,
        "nrealizations": nrealizations,
        "nzrand": nzrand,
        "runtime_controls": runtime,
    }


def count_measurements(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            stripped = raw.strip()
            if stripped and not stripped.startswith("#"):
                count += 1
    return count


def parse_float_or_none(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    stripped = str(value).strip().lower()
    if stripped in UNSET_VALUES:
        return None
    return float(stripped)


def estimate_case(
    case_dir: Path,
    infile: str = "para.in",
    *,
    ray_backend: Optional[str] = None,
    threads: Optional[Union[int, str]] = None,
    gp_spfra: Optional[Union[float, str]] = None,
    warn_gb: Optional[float] = 64.0,
    hard_gb: Optional[float] = None,
) -> Dict[str, Any]:
    case_dir = case_dir.expanduser().resolve()
    para_path = case_dir / infile
    para = parse_para(para_path)
    runtime = para["runtime_controls"]
    data_path = case_dir / para["data_file"]
    if not data_path.exists():
        raise FileNotFoundError(data_path)
    measurements = count_measurements(data_path)

    ray = (ray_backend or runtime.get("ray_backend") or "serial").strip().lower()
    if ray not in {"serial", "fmm_parallel"}:
        ray = "serial"
    if threads is None:
        threads_value = runtime.get("threads") or "1"
    else:
        threads_value = str(threads)
    nthreads = max(1, int(float(threads_value)))
    gp_value = parse_float_or_none(str(gp_spfra) if gp_spfra is not None else runtime.get("gp_spfra"))
    if gp_value is None:
        gp_value = 0.8

    nx = int(para["nx"])
    ny = int(para["ny"])
    nz = int(para["nz"])
    kmax = int(para["kmax"])
    max_source_receiver = int(para["max_source_receiver"])
    iso_mode = bool(para["iso_mode"])
    maxvp = (nx - 2) * (ny - 2) * (nz - 1)
    maxm = maxvp if iso_mode else 3 * maxvp
    model_multiplier = 1 if iso_mode else 3
    maxnar = int(float(para["spfra"]) * measurements * nx * ny * nz * model_multiplier)
    dense_gdiag = ray != "fmm_parallel"

    dense_diag_bytes = 3 * REAL4
    if dense_gdiag and iso_mode:
        dense_diag_bytes = measurements * maxvp * REAL4 + 2 * REAL4
    elif dense_gdiag:
        dense_diag_bytes = 3 * measurements * maxvp * REAL4

    sparse_triplet_bytes = maxnar * (REAL4 + 2 * INT4 + INT4)
    nslots = max_source_receiver * max_source_receiver * kmax
    receiver_source_work_bytes = nslots * 3 * REAL4
    data_vector_bytes = measurements * 12 * REAL4
    model_vector_bytes = (
        (4 * maxm + maxm * max(1, int(para["nrealizations"])) + 3 * maxm) * REAL4
        + 3 * maxvp * REAL4
        + maxvp * REAL4
    )
    solver_work_bytes = (measurements + maxm * 3) * REAL4
    lsen_bytes = REAL4 if iso_mode else nx * ny * kmax * (nz - 1) * REAL4
    anis_model_bytes = 2 * REAL4 if iso_mode else 2 * maxvp * REAL4
    trcv_bytes = (nx - 2) * (ny - 2) * kmax * REAL8
    core_other_bytes = (
        receiver_source_work_bytes
        + data_vector_bytes
        + model_vector_bytes
        + solver_work_bytes
        + lsen_bytes
        + anis_model_bytes
        + trcv_bytes
        + measurements * INT4
    )

    voro_projection_per_realization_bytes = 0
    active_realization_threads = 0
    ncells = 0
    nunk = 0
    gp_maxnar = 0
    if para["voro"]:
        ncells = (int(para["ncell"]) + int(para["acell"])) * int(para["nzrand"])
        nunk = ncells if iso_mode else 3 * ncells
        gp_maxnar = int(gp_value * measurements * max(1, nunk))
        grid_work = maxvp * (7 * REAL4) + measurements * REAL4
        cell_work = max(1, ncells) * 10 * REAL4 + max(1, nunk) * 2 * REAL4
        gp_sparse = gp_maxnar * (2 * INT4 + INT4 + REAL4)
        voro_projection_per_realization_bytes = grid_work + cell_work + gp_sparse + maxm * REAL4
        active_realization_threads = min(nthreads, max(1, int(para["nrealizations"])))

    base_upper_bytes = sparse_triplet_bytes + dense_diag_bytes + core_other_bytes
    voro_parallel_upper_bytes = voro_projection_per_realization_bytes * active_realization_threads
    peak_upper_bytes = base_upper_bytes + voro_parallel_upper_bytes

    warnings: List[str] = []
    notes: List[str] = []
    if iso_mode and ray == "fmm_parallel":
        notes.append("ISO+fmm_parallel uses dummy dense diagnostic G arrays; ANI dense diagnostic arrays are not expected.")
    if dense_gdiag and dense_diag_bytes > 32 * 1024**3:
        warnings.append("serial dense diagnostic G arrays are very large; use fmm_parallel for large data.")
    if para["voro"] and active_realization_threads > 1 and voro_projection_per_realization_bytes > 2 * 1024**3:
        warnings.append("Vorotomo projection memory is multiplied by active realization threads; reduce threads for high-cell tests.")
    if warn_gb is not None and gib(peak_upper_bytes) > float(warn_gb):
        warnings.append(f"estimated peak upper bound exceeds warn threshold {float(warn_gb):.1f} GB.")
    hard_exceeded = False
    if hard_gb is not None and gib(peak_upper_bytes) > float(hard_gb):
        hard_exceeded = True
        warnings.append(f"estimated peak upper bound exceeds hard threshold {float(hard_gb):.1f} GB.")

    return {
        "case_dir": str(case_dir),
        "para": para,
        "runtime": {
            "ray_backend": ray,
            "threads": nthreads,
            "gp_spfra": gp_value,
            "dense_gdiag": dense_gdiag,
        },
        "dimensions": {
            "measurements": measurements,
            "maxvp": maxvp,
            "maxm": maxm,
            "maxnar_main_sparse_upper": maxnar,
            "ncells_total": ncells,
            "voro_unknowns": nunk,
            "gp_maxnar_projection_upper": gp_maxnar,
            "active_realization_threads": active_realization_threads,
        },
        "memory_bytes": {
            "main_sparse_triplet_upper": sparse_triplet_bytes,
            "dense_diagnostic_g": dense_diag_bytes,
            "core_other_upper": core_other_bytes,
            "base_upper": base_upper_bytes,
            "voro_projection_per_realization_upper": voro_projection_per_realization_bytes,
            "voro_projection_parallel_upper": voro_parallel_upper_bytes,
            "peak_upper": peak_upper_bytes,
        },
        "memory_gb": {
            "main_sparse_triplet_upper": gib(sparse_triplet_bytes),
            "dense_diagnostic_g": gib(dense_diag_bytes),
            "core_other_upper": gib(core_other_bytes),
            "base_upper": gib(base_upper_bytes),
            "voro_projection_per_realization_upper": gib(voro_projection_per_realization_bytes),
            "voro_projection_parallel_upper": gib(voro_parallel_upper_bytes),
            "peak_upper": gib(peak_upper_bytes),
        },
        "warnings": warnings,
        "notes": notes,
        "hard_exceeded": hard_exceeded,
    }


def format_report(report: Dict[str, Any]) -> str:
    para = report["para"]
    dims = report["dimensions"]
    mem = report["memory_gb"]
    runtime = report["runtime"]
    lines = [
        "Preflight memory estimate:",
        f"  iso_mode={para['iso_mode']} voro={para['voro']} ray_backend={runtime['ray_backend']} threads={runtime['threads']}",
        f"  measurements={dims['measurements']} maxvp={dims['maxvp']} maxm={dims['maxm']} main_maxnar={dims['maxnar_main_sparse_upper']}",
        f"  main_sparse_upper={mem['main_sparse_triplet_upper']:.2f} GB dense_diag={mem['dense_diagnostic_g']:.2f} GB base_upper={mem['base_upper']:.2f} GB",
    ]
    if para["voro"]:
        lines.append(
            "  "
            f"voro_cells={dims['ncells_total']} gp_spfra={runtime['gp_spfra']:.3g} "
            f"active_realization_threads={dims['active_realization_threads']}"
        )
        lines.append(
            "  "
            f"voro_projection_per_realization={mem['voro_projection_per_realization_upper']:.2f} GB "
            f"voro_projection_parallel={mem['voro_projection_parallel_upper']:.2f} GB"
        )
    lines.append(f"  estimated_peak_upper={mem['peak_upper']:.2f} GB")
    for note in report["notes"]:
        lines.append(f"  NOTE: {note}")
    for warning in report["warnings"]:
        lines.append(f"  WARN: {warning}")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Estimate DAzi/Vorotomo case memory before running")
    parser.add_argument("case_dir", nargs="?", type=Path, default=None)
    parser.add_argument("--case-dir", dest="case_dir_opt", type=Path, default=None)
    parser.add_argument("--infile", default="para.in")
    parser.add_argument("--ray-backend", default=None)
    parser.add_argument("--threads", default=None)
    parser.add_argument("--gp-spfra", default=None)
    parser.add_argument("--warn-gb", type=float, default=64.0)
    parser.add_argument("--hard-gb", type=float, default=None)
    parser.add_argument("--json", type=Path, default=None, help="write JSON report to this file")
    parser.add_argument("--json-only", action="store_true")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    case_dir = args.case_dir_opt or args.case_dir or Path(".")
    try:
        report = estimate_case(
            case_dir,
            args.infile,
            ray_backend=args.ray_backend,
            threads=args.threads,
            gp_spfra=args.gp_spfra,
            warn_gb=args.warn_gb,
            hard_gb=args.hard_gb,
        )
    except Exception as exc:
        print(f"dazi_preflight: {exc}", file=sys.stderr)
        return 2
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if args.json_only:
        print(json.dumps(report, indent=2, ensure_ascii=False))
    else:
        print(format_report(report))
    return 3 if report["hard_exceeded"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
