#!/usr/bin/env python3
"""Diagnose sparse G triplet pattern and value differences."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any

import numpy as np


def load_triplets(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rw = np.fromfile(path / "rw_triplets.bin", dtype=np.float32)
    iw = np.fromfile(path / "iw_triplets.bin", dtype=np.int32)
    col = np.fromfile(path / "col_triplets.bin", dtype=np.int32)
    if not (rw.size == iw.size == col.size):
        raise ValueError(f"Triplet size mismatch under {path}: rw={rw.size} iw={iw.size} col={col.size}")
    return iw, col, rw.astype(np.float64)


def make_key(iw: np.ndarray, col: np.ndarray, max_col: int) -> np.ndarray:
    if np.any(col < 1) or np.any(col > max_col):
        raise ValueError(f"column outside expected 1..{max_col}")
    return iw.astype(np.int64) * np.int64(max_col + 1) + col.astype(np.int64)


def block_depth(col: np.ndarray, nparpi: int, cells_per_depth: int) -> tuple[np.ndarray, np.ndarray]:
    local0 = col.astype(np.int64) - 1
    block = local0 // nparpi
    within = local0 % nparpi
    depth0 = within // cells_per_depth
    return block.astype(np.int16), depth0.astype(np.int16)


def block_name(block: int) -> str:
    return ("Vs", "Gc", "Gs")[block] if 0 <= block <= 2 else f"block_{block}"


def ftol_bucket(values: np.ndarray, ftol: float) -> np.ndarray:
    av = np.abs(values)
    buckets = np.full(values.size, ">10x", dtype=object)
    bucket_defs = [
        (1.0, "<=1x"),
        (1.1, "<=1.1x"),
        (1.5, "<=1.5x"),
        (2.0, "<=2x"),
        (5.0, "<=5x"),
        (10.0, "<=10x"),
    ]
    previous = np.zeros(values.size, dtype=bool)
    for multiplier, label in bucket_defs:
        mask = (av <= multiplier * ftol) & ~previous
        buckets[mask] = label
        previous |= mask
    return buckets


def stats(values: np.ndarray, ref: np.ndarray | None = None) -> dict[str, Any]:
    if values.size == 0:
        out: dict[str, Any] = {
            "n": 0,
            "min_abs": 0.0,
            "max_abs": 0.0,
            "mean_abs": 0.0,
            "rms": 0.0,
        }
        if ref is not None:
            out["rel_l2"] = 0.0
        return out
    av = np.abs(values)
    out = {
        "n": int(values.size),
        "min_abs": float(np.min(av)),
        "max_abs": float(np.max(av)),
        "mean_abs": float(np.mean(av)),
        "rms": float(math.sqrt(float(np.mean(values * values)))),
    }
    if ref is not None:
        out["rel_l2"] = float(np.linalg.norm(values) / max(np.linalg.norm(ref), 1.0e-30))
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarize_by_group(
    label: str,
    values: np.ndarray,
    iw: np.ndarray,
    col: np.ndarray,
    *,
    nparpi: int,
    cells_per_depth: int,
    ftol: float,
    ref: np.ndarray | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    block, depth0 = block_depth(col, nparpi, cells_per_depth)
    buckets = ftol_bucket(values, ftol)
    rows.append({
        "set": label,
        "group": "all",
        "block": "all",
        "bucket": "all",
        "row_min": int(np.min(iw)) if iw.size else "",
        "row_max": int(np.max(iw)) if iw.size else "",
        "near_1x": int(np.count_nonzero(np.abs(values) <= ftol)),
        "near_1p5x": int(np.count_nonzero(np.abs(values) <= 1.5 * ftol)),
        **stats(values, ref),
    })
    for blk in sorted(np.unique(block).tolist()):
        mask = block == blk
        rows.append({
            "set": label,
            "group": "block",
            "block": block_name(int(blk)),
            "bucket": "all",
            "row_min": int(np.min(iw[mask])),
            "row_max": int(np.max(iw[mask])),
            "near_1x": int(np.count_nonzero(np.abs(values[mask]) <= ftol)),
            "near_1p5x": int(np.count_nonzero(np.abs(values[mask]) <= 1.5 * ftol)),
            **stats(values[mask], ref[mask] if ref is not None else None),
        })
    for bucket in ["<=1x", "<=1.1x", "<=1.5x", "<=2x", "<=5x", "<=10x", ">10x"]:
        mask = buckets == bucket
        if not np.any(mask):
            continue
        rows.append({
            "set": label,
            "group": "ftol_bucket",
            "block": "all",
            "bucket": bucket,
            "row_min": int(np.min(iw[mask])),
            "row_max": int(np.max(iw[mask])),
            **stats(values[mask], ref[mask] if ref is not None else None),
        })
    for blk in sorted(np.unique(block).tolist()):
        for depth in sorted(np.unique(depth0[block == blk]).tolist()):
            mask = (block == blk) & (depth0 == depth)
            rows.append({
                "set": label,
                "group": "block_depth",
                "block": block_name(int(blk)),
                "depth_index": int(depth + 1),
                "bucket": "all",
                "row_min": int(np.min(iw[mask])),
                "row_max": int(np.max(iw[mask])),
                "near_1x": int(np.count_nonzero(np.abs(values[mask]) <= ftol)),
                "near_1p5x": int(np.count_nonzero(np.abs(values[mask]) <= 1.5 * ftol)),
                **stats(values[mask], ref[mask] if ref is not None else None),
            })
    return rows


def top_entries(
    label: str,
    values: np.ndarray,
    iw: np.ndarray,
    col: np.ndarray,
    *,
    nparpi: int,
    cells_per_depth: int,
    count: int,
) -> list[dict[str, Any]]:
    if values.size == 0:
        return []
    block, depth0 = block_depth(col, nparpi, cells_per_depth)
    order = np.argsort(np.abs(values))[::-1][:count]
    rows: list[dict[str, Any]] = []
    for idx in order:
        rows.append({
            "set": label,
            "row": int(iw[idx]),
            "col": int(col[idx]),
            "block": block_name(int(block[idx])),
            "depth_index": int(depth0[idx] + 1),
            "value": float(values[idx]),
            "abs_value": float(abs(values[idx])),
        })
    return rows


def format_float(value: Any) -> str:
    if value == "":
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_markdown(path: Path, rows: list[dict[str, Any]], top_rows: list[dict[str, Any]], headline: dict[str, Any]) -> None:
    def pick(set_name: str, group: str, block: str = "all", bucket: str = "all") -> dict[str, Any]:
        for row in rows:
            if row.get("set") == set_name and row.get("group") == group and row.get("block") == block and row.get("bucket") == bucket:
                return row
        return {}

    lines = [
        "# G Triplet Pattern Diagnostics",
        "",
        "## Headline",
        "",
        "| metric | value |",
        "|---|---:|",
    ]
    for key, value in headline.items():
        lines.append(f"| `{key}` | {format_float(value)} |")

    lines.extend([
        "",
        "## Main Sets",
        "",
        "| set | n | rel_l2 | min_abs | mean_abs | max_abs | near_1x | near_1.5x | row_min | row_max |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for set_name in ["common_diff", "base_only", "test_only"]:
        row = pick(set_name, "all")
        lines.append(
            f"| `{set_name}` | {row.get('n', '')} | {format_float(row.get('rel_l2', ''))} | "
            f"{format_float(row.get('min_abs', ''))} | {format_float(row.get('mean_abs', ''))} | "
            f"{format_float(row.get('max_abs', ''))} | {row.get('near_1x', '')} | "
            f"{row.get('near_1p5x', '')} | {row.get('row_min', '')} | {row.get('row_max', '')} |"
        )

    lines.extend([
        "",
        "## Block Summary",
        "",
        "| set | block | n | rel_l2 | mean_abs | max_abs | near_1x | near_1.5x |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for set_name in ["common_diff", "base_only", "test_only"]:
        for block in ["Vs", "Gc", "Gs"]:
            row = pick(set_name, "block", block=block)
            if not row:
                continue
            lines.append(
                f"| `{set_name}` | `{block}` | {row.get('n', '')} | {format_float(row.get('rel_l2', ''))} | "
                f"{format_float(row.get('mean_abs', ''))} | {format_float(row.get('max_abs', ''))} | "
                f"{row.get('near_1x', '')} | {row.get('near_1p5x', '')} |"
            )

    lines.extend([
        "",
        "## Pattern-Only ftol Buckets",
        "",
        "| set | bucket | n | min_abs | mean_abs | max_abs |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for set_name in ["base_only", "test_only"]:
        for bucket in ["<=1x", "<=1.1x", "<=1.5x", "<=2x", "<=5x", "<=10x", ">10x"]:
            row = pick(set_name, "ftol_bucket", bucket=bucket)
            if not row:
                continue
            lines.append(
                f"| `{set_name}` | `{bucket}` | {row.get('n', '')} | {format_float(row.get('min_abs', ''))} | "
                f"{format_float(row.get('mean_abs', ''))} | {format_float(row.get('max_abs', ''))} |"
            )

    lines.extend([
        "",
        "## Top Pattern/Common Entries",
        "",
        "| set | row | col | block | depth_index | value | abs_value |",
        "|---|---:|---:|---|---:|---:|---:|",
    ])
    for row in top_rows:
        lines.append(
            f"| `{row.get('set','')}` | {row.get('row','')} | {row.get('col','')} | `{row.get('block','')}` | "
            f"{row.get('depth_index','')} | {format_float(row.get('value',''))} | {format_float(row.get('abs_value',''))} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--test", required=True, type=Path)
    parser.add_argument("--nparpi", required=True, type=int)
    parser.add_argument("--cells-per-depth", required=True, type=int)
    parser.add_argument("--ftol", type=float, default=1.0e-4)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    max_col = args.nparpi * 3
    biw, bcol, brw = load_triplets(args.base)
    tiw, tcol, trw = load_triplets(args.test)
    bkey = make_key(biw, bcol, max_col)
    tkey = make_key(tiw, tcol, max_col)

    common, bidx, tidx = np.intersect1d(bkey, tkey, return_indices=True)
    base_only_mask = ~np.isin(bkey, tkey)
    test_only_mask = ~np.isin(tkey, bkey)
    b_only = np.nonzero(base_only_mask)[0]
    t_only = np.nonzero(test_only_mask)[0]

    common_diff = trw[tidx] - brw[bidx]
    common_base = brw[bidx]
    headline = {
        "base_n": int(brw.size),
        "test_n": int(trw.size),
        "common_n": int(common.size),
        "base_only_n": int(b_only.size),
        "test_only_n": int(t_only.size),
        "base_only_unique_rows": int(np.unique(biw[b_only]).size),
        "test_only_unique_rows": int(np.unique(tiw[t_only]).size),
        "base_only_unique_cols": int(np.unique(bcol[b_only]).size),
        "test_only_unique_cols": int(np.unique(tcol[t_only]).size),
        "common_rel_l2": float(np.linalg.norm(common_diff) / max(np.linalg.norm(common_base), 1.0e-30)),
        "base_only_near_1p5x": int(np.count_nonzero(np.abs(brw[b_only]) <= 1.5 * args.ftol)),
        "test_only_near_1p5x": int(np.count_nonzero(np.abs(trw[t_only]) <= 1.5 * args.ftol)),
    }

    rows: list[dict[str, Any]] = []
    rows.extend(summarize_by_group(
        "common_diff",
        common_diff,
        biw[bidx],
        bcol[bidx],
        nparpi=args.nparpi,
        cells_per_depth=args.cells_per_depth,
        ftol=args.ftol,
        ref=common_base,
    ))
    rows.extend(summarize_by_group(
        "base_only",
        brw[b_only],
        biw[b_only],
        bcol[b_only],
        nparpi=args.nparpi,
        cells_per_depth=args.cells_per_depth,
        ftol=args.ftol,
    ))
    rows.extend(summarize_by_group(
        "test_only",
        trw[t_only],
        tiw[t_only],
        tcol[t_only],
        nparpi=args.nparpi,
        cells_per_depth=args.cells_per_depth,
        ftol=args.ftol,
    ))

    top_rows = []
    top_rows.extend(top_entries(
        "common_diff",
        common_diff,
        biw[bidx],
        bcol[bidx],
        nparpi=args.nparpi,
        cells_per_depth=args.cells_per_depth,
        count=args.top_n,
    ))
    top_rows.extend(top_entries(
        "base_only",
        brw[b_only],
        biw[b_only],
        bcol[b_only],
        nparpi=args.nparpi,
        cells_per_depth=args.cells_per_depth,
        count=args.top_n,
    ))
    top_rows.extend(top_entries(
        "test_only",
        trw[t_only],
        tiw[t_only],
        tcol[t_only],
        nparpi=args.nparpi,
        cells_per_depth=args.cells_per_depth,
        count=args.top_n,
    ))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "g_triplet_pattern_summary.csv", rows)
    write_csv(args.out_dir / "g_triplet_top_entries.csv", top_rows)
    write_markdown(args.out_dir / "G_TRIPLET_PATTERN_DIAGNOSTICS.md", rows, top_rows, headline)

    print(
        f"base_n={headline['base_n']} test_n={headline['test_n']} "
        f"common={headline['common_n']} base_only={headline['base_only_n']} "
        f"test_only={headline['test_only_n']} common_rel_l2={headline['common_rel_l2']:.6g}"
    )
    print(f"wrote {args.out_dir / 'G_TRIPLET_PATTERN_DIAGNOSTICS.md'}")


if __name__ == "__main__":
    main()
