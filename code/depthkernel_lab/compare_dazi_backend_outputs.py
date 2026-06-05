#!/usr/bin/env python3
"""Compare DAzi output files from two depth-kernel backends."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


FILES = [
    "DSurfTomo.inv",
    "Gc_Gs_model.inv",
    "period_Azm_tomo.inv",
    "period_phaseVMOD.dat",
    "phaseV_FWD.dat",
    "Traveltime_statis_01th.dat",
    "Traveltime_use_01th.dat",
    "cost_terms.dat",
    "costVSiter.dat",
    "misfit_sorted.dat",
    "misfit_voro_raw_sanitized.dat",
]


def parse_log(path: Path) -> dict[str, float | int | str]:
    text = path.read_text(errors="replace")
    out: dict[str, float | int | str] = {"log": str(path)}

    patterns = {
        "nar": r"Number of non-zero element \(nar\) of G is\s+([0-9]+)",
        "before_abs_mean_s": r"Before Inversion: abs mean, std, RMS of Res:\s+([0-9.Ee+-]+)",
        "before_std_s": r"Before Inversion: abs mean, std, RMS of Res:\s+[0-9.Ee+-]+\s+s\s+([0-9.Ee+-]+)",
        "before_rms_s": r"Before Inversion: abs mean, std, RMS of Res:\s+[0-9.Ee+-]+\s+s\s+[0-9.Ee+-]+\s+s\s+([0-9.Ee+-]+)",
        "after_abs_mean_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+([0-9.Ee+-]+)",
        "after_std_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+[0-9.Ee+-]+\s+s\s+([0-9.Ee+-]+)",
        "after_rms_s": r"After Inversion: abs mean, std, RMS of Res\s+:\s+[0-9.Ee+-]+\s+s\s+[0-9.Ee+-]+\s+s\s+([0-9.Ee+-]+)",
        "all_time_s": r"All time cost=\s+([0-9.Ee+-]+)s",
        "pdsurf_batch_s": r"batch_seconds=([0-9.Ee+-]+)",
        "pdsurf_wall_s": r"batch_seconds=[0-9.Ee+-]+\s+wall_seconds=([0-9.Ee+-]+)",
        "ti_batch_wall_s": r"workers=[0-9]+\s+tasks=[0-9]+\s+wall_seconds=([0-9.Ee+-]+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if not match:
            out[key] = ""
            continue
        value = match.group(1)
        out[key] = int(value) if key == "nar" else float(value)

    out["finished"] = "Program finishes successfully" in text
    return out


def load_numeric(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for line in path.read_text(errors="replace").splitlines():
        vals = []
        for token in line.split():
            try:
                vals.append(float(token))
            except ValueError:
                pass
        if vals:
            rows.append(vals)
    if not rows:
        return np.array([], dtype=float)
    width = min(len(row) for row in rows)
    return np.array([row[:width] for row in rows], dtype=float)


def compare_file(base_dir: Path, test_dir: Path, name: str) -> dict[str, str | float | int]:
    base_path = base_dir / name
    test_path = test_dir / name
    row: dict[str, str | float | int] = {"file": name}
    if not base_path.exists() or not test_path.exists():
        row.update({"status": "missing"})
        return row

    a = load_numeric(base_path)
    b = load_numeric(test_path)
    row["base_shape"] = "x".join(map(str, a.shape))
    row["test_shape"] = "x".join(map(str, b.shape))
    if a.shape != b.shape:
        row.update({"status": "shape_mismatch"})
        return row
    if a.size == 0:
        row.update({"status": "empty", "max_abs": 0.0, "rms_abs": 0.0, "rel_l2": 0.0})
        return row

    diff = b - a
    finite = np.isfinite(a) & np.isfinite(b)
    if not np.all(finite):
        a = a[finite]
        b = b[finite]
        diff = b - a
    if diff.size == 0:
        row.update({"status": "nonfinite_only"})
        return row

    norm_a = float(np.linalg.norm(a.ravel()))
    row.update(
        {
            "status": "ok",
            "n": int(diff.size),
            "max_abs": float(np.max(np.abs(diff))),
            "mean_abs": float(np.mean(np.abs(diff))),
            "rms_abs": float(math.sqrt(float(np.mean(diff * diff)))),
            "rel_l2": float(np.linalg.norm(diff.ravel()) / max(norm_a, 1.0e-30)),
            "base_min": float(np.min(a)),
            "base_max": float(np.max(a)),
            "test_min": float(np.min(b)),
            "test_max": float(np.max(b)),
        }
    )
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--test", required=True, type=Path)
    parser.add_argument("--base-log", default="run_default.log")
    parser.add_argument("--test-log", default="run_pdsurf_numba.log")
    parser.add_argument("--base-label", default="default")
    parser.add_argument("--test-label", default="test")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    base_log = parse_log(args.base / args.base_log)
    test_log = parse_log(args.test / args.test_log)

    log_fields = sorted(set(base_log) | set(test_log))
    with (args.outdir / "log_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["backend", *log_fields])
        writer.writeheader()
        writer.writerow({"backend": args.base_label, **base_log})
        writer.writerow({"backend": args.test_label, **test_log})

    rows = [compare_file(args.base, args.test, name) for name in FILES]
    fields = ["file", "status", "base_shape", "test_shape", "n", "max_abs", "mean_abs", "rms_abs", "rel_l2", "base_min", "base_max", "test_min", "test_max"]
    with (args.outdir / "file_diff_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print("log_summary")
    for label, data in [(args.base_label, base_log), (args.test_label, test_log)]:
        print(
            f"{label}: finished={data.get('finished')} nar={data.get('nar')} "
            f"before_abs={data.get('before_abs_mean_s')} after_abs={data.get('after_abs_mean_s')} "
            f"all_time={data.get('all_time_s')} pdsurf_wall={data.get('pdsurf_wall_s')}"
        )

    print("file_diff_summary")
    for row in rows:
        if row.get("status") != "ok":
            print(f"{row['file']}: {row.get('status')}")
        else:
            print(
                f"{row['file']}: max_abs={row['max_abs']:.6g} "
                f"rms_abs={row['rms_abs']:.6g} rel_l2={row['rel_l2']:.6g}"
            )


if __name__ == "__main__":
    main()
