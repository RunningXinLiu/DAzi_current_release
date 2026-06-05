#!/usr/bin/env python3
"""Detailed post-process diagnostics for a real DAzi backend A/B."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Iterable


MODEL_FILE = "Gc_Gs_model.inv"
PERIOD_FILE = "period_Azm_tomo.inv"
TRAVELTIME_FILE = "Traveltime_use_05th.dat"

AMP_THRESHOLDS_PCT = [0.1, 0.5, 1.0]
TOP_FRACTIONS = [0.20, 0.10]


def parse_numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            values: list[float] = []
            for token in line.split():
                try:
                    values.append(float(token))
                except ValueError:
                    pass
            if values:
                rows.append(values)
    return rows


def percentile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    vals = sorted(values)
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * p / 100.0
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    w = pos - lo
    return vals[lo] * (1.0 - w) + vals[hi] * w


def correlation(left: list[float], right: list[float]) -> float | None:
    n = min(len(left), len(right))
    if n == 0:
        return None
    a = left[:n]
    b = right[:n]
    ma = sum(a) / n
    mb = sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((x - mb) ** 2 for x in b)
    if va <= 0.0 or vb <= 0.0:
        return None
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return cov / math.sqrt(va * vb)


def signed_stats(base: list[float], test: list[float]) -> dict[str, Any]:
    n = min(len(base), len(test))
    if n == 0:
        return empty_stats()
    diffs = [test[i] - base[i] for i in range(n)]
    abs_diffs = [abs(value) for value in diffs]
    return {
        "n": n,
        "bias": sum(diffs) / n,
        "mean_abs": sum(abs_diffs) / n,
        "rms": math.sqrt(sum(value * value for value in diffs) / n),
        "p50_abs": percentile(abs_diffs, 50.0),
        "p90_abs": percentile(abs_diffs, 90.0),
        "p95_abs": percentile(abs_diffs, 95.0),
        "p99_abs": percentile(abs_diffs, 99.0),
        "max_abs": max(abs_diffs),
        "corr": correlation(base[:n], test[:n]),
    }


def scalar_stats(values: list[float]) -> dict[str, Any]:
    if not values:
        return empty_stats()
    abs_values = [abs(value) for value in values]
    n = len(values)
    return {
        "n": n,
        "bias": None,
        "mean_abs": sum(abs_values) / n,
        "rms": math.sqrt(sum(value * value for value in values) / n),
        "p50_abs": percentile(abs_values, 50.0),
        "p90_abs": percentile(abs_values, 90.0),
        "p95_abs": percentile(abs_values, 95.0),
        "p99_abs": percentile(abs_values, 99.0),
        "max_abs": max(abs_values),
        "corr": None,
    }


def empty_stats() -> dict[str, Any]:
    return {
        "n": 0,
        "bias": None,
        "mean_abs": None,
        "rms": None,
        "p50_abs": None,
        "p90_abs": None,
        "p95_abs": None,
        "p99_abs": None,
        "max_abs": None,
        "corr": None,
    }


def angle_diff_180(test: float, base: float) -> float:
    return (test - base + 90.0) % 180.0 - 90.0


def axis_alignment_from_angle_diffs(angle_diffs_deg: list[float]) -> float | None:
    if not angle_diffs_deg:
        return None
    vals = [math.cos(math.radians(2.0 * value)) for value in angle_diffs_deg]
    return sum(vals) / len(vals)


def angle_from_gcgs(gc_pct: float, gs_pct: float) -> float:
    return 0.5 * math.degrees(math.atan2(gs_pct, gc_pct)) % 180.0


def amp_from_gcgs_pct(gc_pct: float, gs_pct: float) -> float:
    return 0.5 * math.hypot(gc_pct, gs_pct)


def model_record(row: list[float], idx: int) -> dict[str, float | int]:
    if len(row) < 8:
        raise ValueError("Gc_Gs_model.inv needs at least 8 numeric columns")
    rec: dict[str, float | int] = {
        "idx": idx,
        "lon": row[0],
        "lat": row[1],
        "depth": row[2],
        "vs": row[3],
        "angle_file": row[4],
        "amp_file_pct": row[5] * 100.0,
        "gc_pct": row[6],
        "gs_pct": row[7],
        "amp_pct": amp_from_gcgs_pct(row[6], row[7]),
        "angle_gcgs": angle_from_gcgs(row[6], row[7]),
    }
    names = [
        "std_vs",
        "std_gc_pct",
        "std_gs_pct",
        "stdall_vs",
        "stdall_gc_pct",
        "stdall_gs_pct",
    ]
    for offset, name in enumerate(names, start=8):
        rec[name] = row[offset] if len(row) > offset else None  # type: ignore[assignment]
    return rec


def load_model(path: Path) -> list[dict[str, Any]]:
    return [model_record(row, idx) for idx, row in enumerate(parse_numeric_rows(path))]


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.8g}"
    return str(value)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def by_depth(records: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    groups["ALL"] = list(records)
    for rec in records:
        groups[f"{float(rec['depth']):.4g}"].append(rec)
    return dict(groups)


def values(records: list[dict[str, Any]], key: str) -> list[float]:
    return [float(rec[key]) for rec in records if rec.get(key) is not None]


def matching_records(
    base: list[dict[str, Any]],
    test: list[dict[str, Any]],
    keep: Callable[[dict[str, Any], dict[str, Any]], bool] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    out_base: list[dict[str, Any]] = []
    out_test: list[dict[str, Any]] = []
    for b, t in zip(base, test):
        if keep is None or keep(b, t):
            out_base.append(b)
            out_test.append(t)
    return out_base, out_test


def depthwise_rows(base: list[dict[str, Any]], test: list[dict[str, Any]]) -> list[dict[str, Any]]:
    base_groups = by_depth(base)
    test_groups = by_depth(test)
    rows: list[dict[str, Any]] = []
    metric_specs = [
        ("vs_kms", "vs", "vs", "km/s"),
        ("amp_pct", "amp_pct", "amp_pct", "percent"),
        ("gc_pct", "gc_pct", "gc_pct", "percent"),
        ("gs_pct", "gs_pct", "gs_pct", "percent"),
    ]
    for depth in sorted(base_groups, key=lambda item: -1.0 if item == "ALL" else float(item)):
        bgroup = base_groups[depth]
        tgroup = test_groups.get(depth, [])
        for metric, bkey, tkey, unit in metric_specs:
            stats = signed_stats(values(bgroup, bkey), values(tgroup, tkey))
            rows.append({"depth_km": depth, "metric": metric, "unit": unit, **stats})

        angle_diffs = [
            angle_diff_180(float(t["angle_gcgs"]), float(b["angle_gcgs"]))
            for b, t in zip(bgroup, tgroup)
        ]
        rows.append(
            {
                "depth_km": depth,
                "metric": "fast_axis_deg",
                "unit": "deg_mod_180",
                **scalar_stats(angle_diffs),
                "axis_alignment_mean": axis_alignment_from_angle_diffs(angle_diffs),
            }
        )
        vector_diffs = [
            math.hypot(float(t["gc_pct"]) - float(b["gc_pct"]), float(t["gs_pct"]) - float(b["gs_pct"]))
            for b, t in zip(bgroup, tgroup)
        ]
        rows.append(
            {
                "depth_km": depth,
                "metric": "gcgs_vector_pct",
                "unit": "percent_vector_norm",
                **scalar_stats(vector_diffs),
            }
        )
    return rows


def amplitude_mask_rows(base: list[dict[str, Any]], test: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    base_groups = by_depth(base)
    test_groups = by_depth(test)
    for depth in sorted(base_groups, key=lambda item: -1.0 if item == "ALL" else float(item)):
        bgroup = base_groups[depth]
        tgroup = test_groups.get(depth, [])
        mask_specs: list[tuple[str, Callable[[dict[str, Any], dict[str, Any]], bool]]] = []
        for threshold in AMP_THRESHOLDS_PCT:
            mask_specs.append(
                (
                    f"amp_ge_{threshold:g}pct",
                    lambda b, t, threshold=threshold: max(float(b["amp_pct"]), float(t["amp_pct"])) >= threshold,
                )
            )
        max_amp = [max(float(b["amp_pct"]), float(t["amp_pct"])) for b, t in zip(bgroup, tgroup)]
        for frac in TOP_FRACTIONS:
            cutoff = percentile(max_amp, (1.0 - frac) * 100.0)
            if cutoff is None:
                continue
            mask_specs.append(
                (
                    f"top_{int(frac * 100)}pct_amp",
                    lambda b, t, cutoff=cutoff: max(float(b["amp_pct"]), float(t["amp_pct"])) >= cutoff,
                )
            )

        for mask_name, predicate in mask_specs:
            mb, mt = matching_records(bgroup, tgroup, predicate)
            angle_diffs = [
                angle_diff_180(float(t["angle_gcgs"]), float(b["angle_gcgs"]))
                for b, t in zip(mb, mt)
            ]
            rows.append(
                {
                    "depth_km": depth,
                    "mask": mask_name,
                    "n": len(angle_diffs),
                    **scalar_stats(angle_diffs),
                    "axis_alignment_mean": axis_alignment_from_angle_diffs(angle_diffs),
                    "base_amp_mean_pct": sum(values(mb, "amp_pct")) / len(mb) if mb else None,
                    "test_amp_mean_pct": sum(values(mt, "amp_pct")) / len(mt) if mt else None,
                }
            )
    return rows


def z_stats(values_in: list[float]) -> dict[str, Any]:
    row = scalar_stats(values_in)
    n = len(values_in)
    row["frac_le_1"] = sum(1 for value in values_in if abs(value) <= 1.0) / n if n else None
    row["frac_le_2"] = sum(1 for value in values_in if abs(value) <= 2.0) / n if n else None
    return row


def uncertainty_rows(base: list[dict[str, Any]], test: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    base_groups = by_depth(base)
    test_groups = by_depth(test)
    specs = [
        ("vs_selected_std", "vs", ("std_vs",), "km/s"),
        ("vs_all_std", "vs", ("stdall_vs",), "km/s"),
        ("gcgs_vector_selected_std", "gcgs_vector", ("std_gc_pct", "std_gs_pct"), "percent"),
        ("gcgs_vector_all_std", "gcgs_vector", ("stdall_gc_pct", "stdall_gs_pct"), "percent"),
    ]
    eps = 1.0e-12
    for depth in sorted(base_groups, key=lambda item: -1.0 if item == "ALL" else float(item)):
        bgroup = base_groups[depth]
        tgroup = test_groups.get(depth, [])
        for metric, kind, std_keys, unit in specs:
            zvals: list[float] = []
            for b, t in zip(bgroup, tgroup):
                if kind == "vs":
                    diff = abs(float(t["vs"]) - float(b["vs"]))
                    sigma = math.sqrt(float(b[std_keys[0]]) ** 2 + float(t[std_keys[0]]) ** 2)
                else:
                    diff = math.hypot(float(t["gc_pct"]) - float(b["gc_pct"]), float(t["gs_pct"]) - float(b["gs_pct"]))
                    sigma = math.sqrt(
                        float(b[std_keys[0]]) ** 2
                        + float(b[std_keys[1]]) ** 2
                        + float(t[std_keys[0]]) ** 2
                        + float(t[std_keys[1]]) ** 2
                    )
                if sigma > eps:
                    zvals.append(diff / sigma)
            rows.append(
                {
                    "depth_km": depth,
                    "metric": metric,
                    "unit": unit,
                    **z_stats(zvals),
                }
            )
    return rows


def top_set(records: list[dict[str, Any]], key: str, frac: float, high: bool) -> set[int]:
    vals = values(records, key)
    if not vals:
        return set()
    cutoff = percentile(vals, (1.0 - frac) * 100.0 if high else frac * 100.0)
    if cutoff is None:
        return set()
    if high:
        return {int(rec["idx"]) for rec in records if float(rec[key]) >= cutoff}
    return {int(rec["idx"]) for rec in records if float(rec[key]) <= cutoff}


def overlap_rows(base: list[dict[str, Any]], test: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    base_groups = by_depth(base)
    test_groups = by_depth(test)
    specs = [
        ("amp_top", "amp_pct", True),
        ("vs_high", "vs", True),
        ("vs_low", "vs", False),
    ]
    for depth in sorted(base_groups, key=lambda item: -1.0 if item == "ALL" else float(item)):
        bgroup = base_groups[depth]
        tgroup = test_groups.get(depth, [])
        for frac in TOP_FRACTIONS:
            for name, key, high in specs:
                bset = top_set(bgroup, key, frac, high)
                tset = top_set(tgroup, key, frac, high)
                inter = bset & tset
                union = bset | tset
                rows.append(
                    {
                        "depth_km": depth,
                        "kind": f"{name}_{int(frac * 100)}pct",
                        "base_count": len(bset),
                        "test_count": len(tset),
                        "intersection": len(inter),
                        "jaccard": len(inter) / len(union) if union else None,
                        "base_recovered_frac": len(inter) / len(bset) if bset else None,
                        "test_recovered_frac": len(inter) / len(tset) if tset else None,
                    }
                )
    return rows


def residual_rows(base_dir: Path, test_dir: Path) -> list[dict[str, Any]]:
    base_rows = parse_numeric_rows(base_dir / TRAVELTIME_FILE)
    test_rows = parse_numeric_rows(test_dir / TRAVELTIME_FILE)
    n = min(len(base_rows), len(test_rows))
    specs = [
        ("forward_time_s", 2),
        ("residual_s", 3),
        ("weighted_residual", 5),
    ]
    rows: list[dict[str, Any]] = []
    for name, col in specs:
        base_vals = [base_rows[i][col] for i in range(n) if len(base_rows[i]) > col and len(test_rows[i]) > col]
        test_vals = [test_rows[i][col] for i in range(n) if len(base_rows[i]) > col and len(test_rows[i]) > col]
        rows.append({"metric": name, **signed_stats(base_vals, test_vals)})
    for label, rows_in in (("base", base_rows[:n]), ("test", test_rows[:n])):
        residuals = [row[3] for row in rows_in if len(row) > 3]
        weighted = [row[5] for row in rows_in if len(row) > 5]
        rows.append(
            {
                "metric": f"{label}_residual_abs_distribution",
                **scalar_stats(residuals),
                "mean_signed": sum(residuals) / len(residuals) if residuals else None,
            }
        )
        rows.append(
            {
                "metric": f"{label}_weighted_residual_abs_distribution",
                **scalar_stats(weighted),
                "mean_signed": sum(weighted) / len(weighted) if weighted else None,
            }
        )
    return rows


def period_rows(base_dir: Path, test_dir: Path) -> list[dict[str, Any]]:
    base_rows = parse_numeric_rows(base_dir / PERIOD_FILE)
    test_rows = parse_numeric_rows(test_dir / PERIOD_FILE)
    grouped: dict[str, tuple[list[list[float]], list[list[float]]]] = {}
    periods = sorted({row[2] for row in base_rows if len(row) >= 9})
    for period in periods:
        b = [row for row in base_rows if len(row) >= 9 and abs(row[2] - period) < 1.0e-6]
        t = [row for row in test_rows if len(row) >= 9 and abs(row[2] - period) < 1.0e-6]
        grouped[f"{period:.4g}"] = (b, t)
    rows: list[dict[str, Any]] = []
    for period, (bgroup, tgroup) in [("ALL", (base_rows, test_rows)), *grouped.items()]:
        n = min(len(bgroup), len(tgroup))
        if n == 0:
            continue
        phase_base = [bgroup[i][3] for i in range(n) if len(bgroup[i]) >= 9 and len(tgroup[i]) >= 9]
        phase_test = [tgroup[i][3] for i in range(n) if len(bgroup[i]) >= 9 and len(tgroup[i]) >= 9]
        amp_base = [bgroup[i][5] * 100.0 for i in range(n) if len(bgroup[i]) >= 9 and len(tgroup[i]) >= 9]
        amp_test = [tgroup[i][5] * 100.0 for i in range(n) if len(bgroup[i]) >= 9 and len(tgroup[i]) >= 9]
        angle_diffs = [
            angle_diff_180(tgroup[i][4], bgroup[i][4])
            for i in range(n)
            if len(bgroup[i]) >= 9 and len(tgroup[i]) >= 9
        ]
        rows.append({"period_s": period, "metric": "phase_velocity", "unit": "km/s", **signed_stats(phase_base, phase_test)})
        rows.append({"period_s": period, "metric": "azim_amp", "unit": "percent", **signed_stats(amp_base, amp_test)})
        rows.append(
            {
                "period_s": period,
                "metric": "fast_axis",
                "unit": "deg_mod_180",
                **scalar_stats(angle_diffs),
                "axis_alignment_mean": axis_alignment_from_angle_diffs(angle_diffs),
            }
        )
    return rows


def best_row(rows: list[dict[str, Any]], *, metric: str, field: str, skip_all: bool = True) -> dict[str, Any] | None:
    candidates = [
        row
        for row in rows
        if row.get("metric") == metric and row.get(field) is not None and (not skip_all or row.get("depth_km") != "ALL")
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda row: float(row[field]))


def format_summary(
    base_dir: Path,
    test_dir: Path,
    depth_rows: list[dict[str, Any]],
    mask_rows: list[dict[str, Any]],
    uncertainty: list[dict[str, Any]],
    overlap: list[dict[str, Any]],
    residual: list[dict[str, Any]],
    period: list[dict[str, Any]],
) -> str:
    all_depth = [row for row in depth_rows if row.get("depth_km") == "ALL"]
    all_masks = [row for row in mask_rows if row.get("depth_km") == "ALL"]
    all_unc = [row for row in uncertainty if row.get("depth_km") == "ALL"]
    all_overlap = [row for row in overlap if row.get("depth_km") == "ALL"]
    all_period = [row for row in period if row.get("period_s") == "ALL"]
    lines = [
        "# Detailed Real A/B Diagnostics",
        "",
        f"- base: `{base_dir}`",
        f"- test: `{test_dir}`",
        "",
        "## Overall Model Differences",
        "",
        "| metric | mean_abs | rms | p95_abs | p99_abs | max_abs | corr | unit |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in all_depth:
        lines.append(
            f"| `{row['metric']}` | {fmt(row.get('mean_abs'))} | {fmt(row.get('rms'))} | "
            f"{fmt(row.get('p95_abs'))} | {fmt(row.get('p99_abs'))} | {fmt(row.get('max_abs'))} | "
            f"{fmt(row.get('corr'))} | {row.get('unit', '')} |"
        )
    lines.extend(["", "## Amplitude-Masked Fast-Axis Differences", ""])
    lines.append("| mask | n | mean_abs | p95_abs | p99_abs | max_abs | axis_alignment |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for row in all_masks:
        lines.append(
            f"| `{row['mask']}` | {fmt(row.get('n'))} | {fmt(row.get('mean_abs'))} | "
            f"{fmt(row.get('p95_abs'))} | {fmt(row.get('p99_abs'))} | {fmt(row.get('max_abs'))} | "
            f"{fmt(row.get('axis_alignment_mean'))} |"
        )
    lines.extend(["", "## Uncertainty-Normalized Differences", ""])
    lines.append("| metric | mean_abs_z | p95_abs_z | p99_abs_z | max_abs_z | frac<=1 | frac<=2 |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for row in all_unc:
        lines.append(
            f"| `{row['metric']}` | {fmt(row.get('mean_abs'))} | {fmt(row.get('p95_abs'))} | "
            f"{fmt(row.get('p99_abs'))} | {fmt(row.get('max_abs'))} | "
            f"{fmt(row.get('frac_le_1'))} | {fmt(row.get('frac_le_2'))} |"
        )
    lines.extend(["", "## Top-Anomaly Overlap", ""])
    lines.append("| kind | base_count | test_count | jaccard | base_recovered | test_recovered |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for row in all_overlap:
        lines.append(
            f"| `{row['kind']}` | {fmt(row.get('base_count'))} | {fmt(row.get('test_count'))} | "
            f"{fmt(row.get('jaccard'))} | {fmt(row.get('base_recovered_frac'))} | "
            f"{fmt(row.get('test_recovered_frac'))} |"
        )
    lines.extend(["", "## Residual And Period Outputs", ""])
    lines.append("| metric | mean_abs | rms | p95_abs | max_abs | corr |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for row in residual:
        lines.append(
            f"| `{row['metric']}` | {fmt(row.get('mean_abs'))} | {fmt(row.get('rms'))} | "
            f"{fmt(row.get('p95_abs'))} | {fmt(row.get('max_abs'))} | {fmt(row.get('corr'))} |"
        )
    for row in all_period:
        lines.append(
            f"| `period_{row['metric']}` | {fmt(row.get('mean_abs'))} | {fmt(row.get('rms'))} | "
            f"{fmt(row.get('p95_abs'))} | {fmt(row.get('max_abs'))} | {fmt(row.get('corr'))} |"
        )
    lines.extend(["", "## Depths With Largest Differences", ""])
    lines.append("| diagnostic | depth_km | value |")
    lines.append("|---|---:|---:|")
    for label, metric, field in [
        ("Vs p95 abs", "vs_kms", "p95_abs"),
        ("Amplitude p95 abs", "amp_pct", "p95_abs"),
        ("Fast-axis p95 abs", "fast_axis_deg", "p95_abs"),
        ("Gc/Gs vector p95 abs", "gcgs_vector_pct", "p95_abs"),
    ]:
        row = best_row(depth_rows, metric=metric, field=field)
        if row:
            lines.append(f"| `{label}` | {row.get('depth_km')} | {fmt(row.get(field))} |")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--test", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    base_dir = args.base.expanduser().resolve()
    test_dir = args.test.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    base_model = load_model(base_dir / MODEL_FILE)
    test_model = load_model(test_dir / MODEL_FILE)
    if len(base_model) != len(test_model):
        raise ValueError(f"model row count differs: {len(base_model)} vs {len(test_model)}")
    for b, t in zip(base_model, test_model):
        if (
            abs(float(b["lon"]) - float(t["lon"])) > 1.0e-5
            or abs(float(b["lat"]) - float(t["lat"])) > 1.0e-5
            or abs(float(b["depth"]) - float(t["depth"])) > 1.0e-5
        ):
            raise ValueError("model grids differ; interpolation would be required")

    depth_rows_out = depthwise_rows(base_model, test_model)
    mask_rows_out = amplitude_mask_rows(base_model, test_model)
    uncertainty_out = uncertainty_rows(base_model, test_model)
    overlap_out = overlap_rows(base_model, test_model)
    residual_out = residual_rows(base_dir, test_dir)
    period_out = period_rows(base_dir, test_dir)

    write_csv(outdir / "depthwise_model_diff.csv", depth_rows_out)
    write_csv(outdir / "amplitude_masked_angle_diff.csv", mask_rows_out)
    write_csv(outdir / "uncertainty_normalized_diff.csv", uncertainty_out)
    write_csv(outdir / "top_anomaly_overlap.csv", overlap_out)
    write_csv(outdir / "residual_diff_stats.csv", residual_out)
    write_csv(outdir / "periodwise_output_diff.csv", period_out)
    summary_text = format_summary(
        base_dir,
        test_dir,
        depth_rows_out,
        mask_rows_out,
        uncertainty_out,
        overlap_out,
        residual_out,
        period_out,
    )
    (outdir / "DETAILED_AB_DIAGNOSTICS.md").write_text(summary_text, encoding="utf-8")
    (outdir / "detailed_ab_diagnostics.json").write_text(
        json.dumps(
            {
                "base": str(base_dir),
                "test": str(test_dir),
                "outputs": [
                    "depthwise_model_diff.csv",
                    "amplitude_masked_angle_diff.csv",
                    "uncertainty_normalized_diff.csv",
                    "top_anomaly_overlap.csv",
                    "residual_diff_stats.csv",
                    "periodwise_output_diff.csv",
                    "DETAILED_AB_DIAGNOSTICS.md",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(outdir / "DETAILED_AB_DIAGNOSTICS.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
