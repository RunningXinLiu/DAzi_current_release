#!/usr/bin/env python3
"""Summarize depth-kernel sweep CSV files and optionally make heatmaps."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


METRICS = (
    "max_abs_velocity_diff",
    "kernel_rel_l2_all",
    "kernel_rel_l2_vs",
    "kernel_rel_l2_vp",
    "kernel_rel_l2_rho",
    "kernel_sign_mismatch_all",
)


def read_rows(path: Path) -> list[dict[str, object]]:
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        rows: list[dict[str, object]] = []
        for row in reader:
            parsed: dict[str, object] = dict(row)
            for key, value in row.items():
                if key in ("label", "column", "compare_file"):
                    continue
                try:
                    parsed[key] = float(value)
                except (TypeError, ValueError):
                    parsed[key] = value
            rows.append(parsed)
    return rows


def write_best_tables(rows: list[dict[str, object]], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    by_label: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_label[str(row.get("label", row.get("column")))].append(row)

    with (outdir / "best_by_column.csv").open("w", newline="") as f:
        fieldnames = [
            "metric",
            "label",
            "best_dln",
            "best_minthk",
            "best_value",
            "worst_dln",
            "worst_minthk",
            "worst_value",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for metric in METRICS:
            for label, label_rows in sorted(by_label.items()):
                best = min(label_rows, key=lambda row: float(row[metric]))
                worst = max(label_rows, key=lambda row: float(row[metric]))
                writer.writerow(
                    {
                        "metric": metric,
                        "label": label,
                        "best_dln": best["dln"],
                        "best_minthk": best["minthk"],
                        "best_value": best[metric],
                        "worst_dln": worst["dln"],
                        "worst_minthk": worst["minthk"],
                        "worst_value": worst[metric],
                    }
                )

    with (outdir / "best_overall.txt").open("w") as f:
        for metric in METRICS:
            best = min(rows, key=lambda row: float(row[metric]))
            worst = max(rows, key=lambda row: float(row[metric]))
            f.write(f"{metric}\n")
            f.write(
                "  best:  "
                f"label={best.get('label', best.get('column'))} dln={best['dln']} "
                f"minthk={best['minthk']} value={best[metric]}\n"
            )
            f.write(
                "  worst: "
                f"label={worst.get('label', worst.get('column'))} dln={worst['dln']} "
                f"minthk={worst['minthk']} value={worst[metric]}\n\n"
            )


def plot_heatmaps(rows: list[dict[str, object]], outdir: Path, metric: str) -> bool:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    by_label: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_label[str(row.get("label", row.get("column")))].append(row)

    for label, label_rows in sorted(by_label.items()):
        dlns = sorted({float(row["dln"]) for row in label_rows}, reverse=True)
        minthks = sorted({float(row["minthk"]) for row in label_rows})
        grid = [[float("nan") for _ in minthks] for _ in dlns]
        index_dln = {value: i for i, value in enumerate(dlns)}
        index_minthk = {value: i for i, value in enumerate(minthks)}
        for row in label_rows:
            i = index_dln[float(row["dln"])]
            j = index_minthk[float(row["minthk"])]
            grid[i][j] = float(row[metric])

        fig, ax = plt.subplots(figsize=(6.0, 4.2), constrained_layout=True)
        image = ax.imshow(grid, aspect="auto", cmap="viridis")
        ax.set_title(f"{label}: {metric}")
        ax.set_xlabel("minthk / nsublayer")
        ax.set_ylabel("finite-difference dln")
        ax.set_xticks(range(len(minthks)), [f"{v:g}" for v in minthks])
        ax.set_yticks(range(len(dlns)), [f"{v:g}" for v in dlns])
        fig.colorbar(image, ax=ax)
        safe_label = label.replace("/", "_").replace(" ", "_")
        fig.savefig(outdir / f"heatmap_{safe_label}_{metric}.png", dpi=180)
        plt.close(fig)
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", default="real_column_sweep_modref/real_columns_sweep_summary.csv")
    parser.add_argument("--outdir", default="real_column_sweep_modref/summary_plots")
    parser.add_argument("--metric", default="kernel_rel_l2_all", choices=METRICS)
    args = parser.parse_args()

    summary = Path(args.summary)
    outdir = Path(args.outdir)
    rows = read_rows(summary)
    if not rows:
        raise SystemExit(f"No rows found in {summary}")

    write_best_tables(rows, outdir)
    plotted = plot_heatmaps(rows, outdir, args.metric)
    print(f"wrote summary files to {outdir}")
    if plotted:
        print(f"wrote heatmaps for metric {args.metric}")
    else:
        print("matplotlib is not available; skipped heatmap PNG output")


if __name__ == "__main__":
    main()
