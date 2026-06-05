#!/usr/bin/env python3
"""Compare DAzi anisotropy model files with angle-aware vector metrics."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterable


MODEL_NAME = "Gc_Gs_model.inv"


def resolve_model_path(path: Path) -> Path:
    return path / MODEL_NAME if path.is_dir() else path


def parse_numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for lineno, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                rows.append([float(token) for token in stripped.split()])
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: non-numeric row: {stripped}") from exc
    return rows


def column(rows: list[list[float]], one_based_col: int, nrows: int) -> list[float]:
    idx = one_based_col - 1
    return [rows[row][idx] for row in range(nrows)]


def angle_diff_180(test: float, base: float) -> float:
    return (test - base + 90.0) % 180.0 - 90.0


def angle_from_gcgs(gc_pct: Iterable[float], gs_pct: Iterable[float]) -> list[float]:
    return [0.5 * math.degrees(math.atan2(gs, gc)) % 180.0 for gc, gs in zip(gc_pct, gs_pct)]


def amp_percent_from_gcgs(gc_pct: Iterable[float], gs_pct: Iterable[float]) -> list[float]:
    return [0.5 * math.hypot(gc, gs) for gc, gs in zip(gc_pct, gs_pct)]


def correlation(left: list[float], right: list[float]) -> float | None:
    n = min(len(left), len(right))
    if n == 0:
        return None
    lx = left[:n]
    rx = right[:n]
    ml = sum(lx) / n
    mr = sum(rx) / n
    vl = sum((value - ml) ** 2 for value in lx)
    vr = sum((value - mr) ** 2 for value in rx)
    if vl <= 0.0 or vr <= 0.0:
        return None
    cov = sum((x - ml) * (y - mr) for x, y in zip(lx, rx))
    return cov / math.sqrt(vl * vr)


def percentile_abs(values: list[float], percentile: float) -> float | None:
    if not values:
        return None
    vals = sorted(abs(value) for value in values)
    if len(vals) == 1:
        return vals[0]
    position = (len(vals) - 1) * percentile / 100.0
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return vals[lo]
    weight = position - lo
    return vals[lo] * (1.0 - weight) + vals[hi] * weight


def signed_stats(base: list[float], test: list[float], *, angle: bool = False) -> dict[str, Any]:
    n = min(len(base), len(test))
    if n == 0:
        return empty_stats()
    diffs = [
        angle_diff_180(test[i], base[i]) if angle else test[i] - base[i]
        for i in range(n)
    ]
    abs_diffs = [abs(value) for value in diffs]
    return {
        "n": n,
        "bias": sum(diffs) / n,
        "mean_abs": sum(abs_diffs) / n,
        "rms": math.sqrt(sum(value * value for value in diffs) / n),
        "max_abs": max(abs_diffs),
        "p95_abs": percentile_abs(diffs, 95.0),
        "corr": correlation(base[:n], test[:n]),
    }


def magnitude_stats(values: list[float]) -> dict[str, Any]:
    n = len(values)
    if n == 0:
        return empty_stats()
    abs_values = [abs(value) for value in values]
    return {
        "n": n,
        "bias": None,
        "mean_abs": sum(abs_values) / n,
        "rms": math.sqrt(sum(value * value for value in values) / n),
        "max_abs": max(abs_values),
        "p95_abs": percentile_abs(values, 95.0),
        "corr": None,
    }


def empty_stats() -> dict[str, Any]:
    return {
        "n": 0,
        "bias": None,
        "mean_abs": None,
        "rms": None,
        "max_abs": None,
        "p95_abs": None,
        "corr": None,
    }


def mask_values(values: list[float], mask: list[bool]) -> list[float]:
    return [value for value, keep in zip(values, mask) if keep]


def compare_anisotropy_files(base_path: Path, test_path: Path, amp_threshold: float = 0.001) -> dict[str, Any]:
    base_file = resolve_model_path(base_path)
    test_file = resolve_model_path(test_path)
    if not base_file.exists() or not test_file.exists():
        return {
            "status": "missing",
            "base": str(base_file),
            "test": str(test_file),
            "base_exists": base_file.exists(),
            "test_exists": test_file.exists(),
        }

    base_rows = parse_numeric_rows(base_file)
    test_rows = parse_numeric_rows(test_file)
    nrows = min(len(base_rows), len(test_rows))
    if nrows == 0:
        return {"status": "empty", "base": str(base_file), "test": str(test_file)}

    base_width = min(len(row) for row in base_rows[:nrows])
    test_width = min(len(row) for row in test_rows[:nrows])
    ncols = min(base_width, test_width)
    if ncols < 8:
        return {
            "status": "insufficient_columns",
            "base": str(base_file),
            "test": str(test_file),
            "base_rows": len(base_rows),
            "test_rows": len(test_rows),
            "base_cols": base_width,
            "test_cols": test_width,
        }

    base_vs = column(base_rows, 4, nrows)
    test_vs = column(test_rows, 4, nrows)
    base_angle = column(base_rows, 5, nrows)
    test_angle = column(test_rows, 5, nrows)
    base_amp = column(base_rows, 6, nrows)
    test_amp = column(test_rows, 6, nrows)
    base_gc = column(base_rows, 7, nrows)
    test_gc = column(test_rows, 7, nrows)
    base_gs = column(base_rows, 8, nrows)
    test_gs = column(test_rows, 8, nrows)

    base_angle_gcgs = angle_from_gcgs(base_gc, base_gs)
    test_angle_gcgs = angle_from_gcgs(test_gc, test_gs)
    base_amp_gcgs = amp_percent_from_gcgs(base_gc, base_gs)
    test_amp_gcgs = amp_percent_from_gcgs(test_gc, test_gs)

    amp_mask = [max(abs(base_amp[i]), abs(test_amp[i])) >= amp_threshold for i in range(nrows)]
    amp_pct_threshold = amp_threshold * 100.0
    amp_gcgs_mask = [
        max(abs(base_amp_gcgs[i]), abs(test_amp_gcgs[i])) >= amp_pct_threshold
        for i in range(nrows)
    ]
    vector_diff = [
        math.hypot(test_gc[i] - base_gc[i], test_gs[i] - base_gs[i])
        for i in range(nrows)
    ]
    vector_ref_norm_sq = sum(base_gc[i] * base_gc[i] + base_gs[i] * base_gs[i] for i in range(nrows))
    vector_diff_norm_sq = sum(value * value for value in vector_diff)
    if vector_ref_norm_sq > 0.0:
        vector_rel_l2 = math.sqrt(vector_diff_norm_sq / vector_ref_norm_sq)
    elif vector_diff_norm_sq == 0.0:
        vector_rel_l2 = 0.0
    else:
        vector_rel_l2 = None

    metrics: dict[str, dict[str, Any]] = {
        "vs": {**signed_stats(base_vs, test_vs), "unit": "km/s"},
        "azim_angle_file_deg": {**signed_stats(base_angle, test_angle, angle=True), "unit": "deg_mod_180"},
        "azim_angle_file_deg_amp_ge_threshold": {
            **signed_stats(mask_values(base_angle, amp_mask), mask_values(test_angle, amp_mask), angle=True),
            "unit": f"deg_mod_180_amp_ge_{amp_threshold:g}",
        },
        "azim_amp_file_fraction": {**signed_stats(base_amp, test_amp), "unit": "fraction"},
        "gc_pct": {**signed_stats(base_gc, test_gc), "unit": "percent"},
        "gs_pct": {**signed_stats(base_gs, test_gs), "unit": "percent"},
        "fast_axis_angle_from_gcgs_deg": {
            **signed_stats(base_angle_gcgs, test_angle_gcgs, angle=True),
            "unit": "deg_mod_180",
        },
        "fast_axis_angle_from_gcgs_deg_amp_ge_threshold": {
            **signed_stats(
                mask_values(base_angle_gcgs, amp_gcgs_mask),
                mask_values(test_angle_gcgs, amp_gcgs_mask),
                angle=True,
            ),
            "unit": f"deg_mod_180_amp_pct_ge_{amp_pct_threshold:g}",
        },
        "amp_from_gcgs_pct": {**signed_stats(base_amp_gcgs, test_amp_gcgs), "unit": "percent"},
        "gcgs_vector_diff_pct": {**magnitude_stats(vector_diff), "unit": "percent_vector_norm"},
    }

    if ncols >= 14:
        names = {
            9: "std_vs",
            10: "std_gc_pct",
            11: "std_gs_pct",
            12: "stdall_vs",
            13: "stdall_gc_pct",
            14: "stdall_gs_pct",
        }
        for one_based_col, name in names.items():
            metrics[name] = {
                **signed_stats(column(base_rows, one_based_col, nrows), column(test_rows, one_based_col, nrows)),
                "unit": "reported",
            }

    return {
        "status": "ok",
        "base": str(base_file),
        "test": str(test_file),
        "base_rows": len(base_rows),
        "test_rows": len(test_rows),
        "compared_rows": nrows,
        "base_cols": base_width,
        "test_cols": test_width,
        "compared_cols": ncols,
        "amp_threshold": amp_threshold,
        "vector_rel_l2": vector_rel_l2,
        "metrics": metrics,
    }


def json_safe(value: Any) -> Any:
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)


def format_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Anisotropy-Aware Model Comparison",
        "",
        f"- status: `{report.get('status')}`",
        f"- base: `{report.get('base')}`",
        f"- test: `{report.get('test')}`",
    ]
    if report.get("status") != "ok":
        return "\n".join(lines)

    lines.extend(
        [
            f"- compared rows: `{report.get('compared_rows')}`",
            f"- compared columns: `{report.get('compared_cols')}`",
            f"- angle amplitude threshold: `{report.get('amp_threshold')}`",
            f"- Gc/Gs vector relative L2: `{fmt(report.get('vector_rel_l2'))}`",
            "",
            "| metric | n | bias | mean_abs | rms | max_abs | p95_abs | corr | unit |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for name, row in report["metrics"].items():
        lines.append(
            f"| `{name}` | {row.get('n', '')} | {fmt(row.get('bias'))} | "
            f"{fmt(row.get('mean_abs'))} | {fmt(row.get('rms'))} | "
            f"{fmt(row.get('max_abs'))} | {fmt(row.get('p95_abs'))} | "
            f"{fmt(row.get('corr'))} | {row.get('unit', '')} |"
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path, help="Baseline Gc_Gs_model.inv or output directory")
    parser.add_argument("test", type=Path, help="Test Gc_Gs_model.inv or output directory")
    parser.add_argument("--amp-threshold", type=float, default=0.001)
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--out-md", type=Path, default=None)
    args = parser.parse_args()

    report = compare_anisotropy_files(args.base, args.test, args.amp_threshold)
    safe_report = json_safe(report)
    text = format_markdown(safe_report)
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(safe_report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if args.out_md is not None:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
