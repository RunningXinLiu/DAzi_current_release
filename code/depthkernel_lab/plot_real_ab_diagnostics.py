#!/usr/bin/env python3
"""Plot interpretation-level diagnostics for two DAzi inversion results."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


MODEL_FILE = "Gc_Gs_model.inv"
PERIOD_FILE = "period_Azm_tomo.inv"
TRAVELTIME_FILE = "Traveltime_use_05th.dat"
DEFAULT_DEPTHS = [1.0, 3.0, 4.0, 5.0, 8.0, 12.0]
AMP_MASK_PCT = 0.5


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 8.5
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["legend.frameon"] = False
plt.rcParams["figure.dpi"] = 120


def read_numeric_table(path: Path) -> np.ndarray:
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
    if not rows:
        raise ValueError(f"no numeric rows found in {path}")
    width = max(len(row) for row in rows)
    out = np.full((len(rows), width), np.nan, dtype=float)
    for i, row in enumerate(rows):
        out[i, : len(row)] = row
    return out


def load_model(run_dir: Path) -> dict[str, np.ndarray]:
    arr = read_numeric_table(run_dir / MODEL_FILE)
    if arr.shape[1] < 8:
        raise ValueError(f"{run_dir / MODEL_FILE} needs at least 8 columns")
    gc = arr[:, 6]
    gs = arr[:, 7]
    return {
        "raw": arr,
        "lon": arr[:, 0],
        "lat": arr[:, 1],
        "depth": arr[:, 2],
        "vs": arr[:, 3],
        "angle": angle_from_gcgs(gc, gs),
        "amp_pct": 0.5 * np.hypot(gc, gs),
        "gc_pct": gc,
        "gs_pct": gs,
    }


def angle_from_gcgs(gc_pct: np.ndarray, gs_pct: np.ndarray) -> np.ndarray:
    return np.mod(0.5 * np.degrees(np.arctan2(gs_pct, gc_pct)), 180.0)


def angle_diff_180(test: np.ndarray, base: np.ndarray) -> np.ndarray:
    return np.mod(test - base + 90.0, 180.0) - 90.0


def check_same_grid(base: dict[str, np.ndarray], test: dict[str, np.ndarray]) -> None:
    for key in ("lon", "lat", "depth"):
        if base[key].shape != test[key].shape or not np.allclose(base[key], test[key], atol=1.0e-5):
            raise ValueError(f"base/test model grids differ in {key}; interpolation would be required")


def safe_percentile(values: np.ndarray, p: float, default: float = 0.0) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return default
    return float(np.percentile(finite, p))


def signed_abs_stats(values: np.ndarray) -> dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {"n": 0, "mean_abs": math.nan, "p95_abs": math.nan, "p99_abs": math.nan, "max_abs": math.nan}
    abs_values = np.abs(finite)
    return {
        "n": int(finite.size),
        "mean_abs": float(np.mean(abs_values)),
        "p95_abs": float(np.percentile(abs_values, 95.0)),
        "p99_abs": float(np.percentile(abs_values, 99.0)),
        "max_abs": float(np.max(abs_values)),
    }


def nearest_depths(available: np.ndarray, requested: list[float]) -> list[float]:
    unique = np.unique(available)
    depths: list[float] = []
    for depth in requested:
        nearest = float(unique[np.argmin(np.abs(unique - depth))])
        if nearest not in depths:
            depths.append(nearest)
    return depths


def depth_mask(model: dict[str, np.ndarray], depth: float) -> np.ndarray:
    return np.isclose(model["depth"], depth, atol=1.0e-5)


def set_map_axes(ax: plt.Axes, lon: np.ndarray, lat: np.ndarray) -> None:
    ax.set_xlim(float(np.min(lon)), float(np.max(lon)))
    ax.set_ylim(float(np.min(lat)), float(np.max(lat)))
    mean_lat = float(np.mean(lat))
    ax.set_aspect(1.0 / math.cos(math.radians(mean_lat)))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.tick_params(labelsize=7, length=2.5)


def scatter_panel(
    ax: plt.Axes,
    lon: np.ndarray,
    lat: np.ndarray,
    values: np.ndarray,
    *,
    title: str,
    cmap: str,
    vmin: float | None = None,
    vmax: float | None = None,
    label: str,
    mask: np.ndarray | None = None,
) -> Any:
    if mask is None:
        mask = np.isfinite(values)
    else:
        mask = mask & np.isfinite(values)
        ax.scatter(lon[~mask], lat[~mask], s=9, marker="s", c="#dedede", linewidths=0.0, alpha=0.65)
    sc = ax.scatter(
        lon[mask],
        lat[mask],
        c=values[mask],
        s=13,
        marker="s",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.0,
    )
    ax.set_title(title, pad=4)
    set_map_axes(ax, lon, lat)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.86, pad=0.02)
    cbar.set_label(label)
    cbar.ax.tick_params(labelsize=7, length=2.5)
    return sc


def draw_sticks(
    ax: plt.Axes,
    lon: np.ndarray,
    lat: np.ndarray,
    angle_deg: np.ndarray,
    amp_pct: np.ndarray,
    *,
    stride: int,
    max_amp_pct: float,
) -> None:
    if lon.size == 0:
        return
    order = np.lexsort((lat, lon))
    order = order[:: max(1, stride)]
    lon_s = lon[order]
    lat_s = lat[order]
    angle_s = np.deg2rad(angle_deg[order])
    amp_s = amp_pct[order]
    max_amp_pct = max(max_amp_pct, 1.0e-6)
    scale = np.clip(amp_s / max_amp_pct, 0.15, 1.0)
    base_len = min(float(np.max(lon) - np.min(lon)), float(np.max(lat) - np.min(lat))) * 0.060
    half_len = 0.5 * base_len * scale
    dx = half_len * np.cos(angle_s)
    dy = half_len * np.sin(angle_s)
    for x0, y0, xh, yh in zip(lon_s, lat_s, dx, dy):
        ax.plot([x0 - xh, x0 + xh], [y0 - yh, y0 + yh], color="#202020", lw=0.42, solid_capstyle="round")


def plot_depth_slice(
    base: dict[str, np.ndarray],
    test: dict[str, np.ndarray],
    depth: float,
    outdir: Path,
    *,
    stride: int,
    base_label: str,
    test_label: str,
    comparison_title: str,
) -> dict[str, float]:
    mask = depth_mask(base, depth)
    lon = base["lon"][mask]
    lat = base["lat"][mask]
    base_vs = base["vs"][mask]
    test_vs = test["vs"][mask]
    base_amp = base["amp_pct"][mask]
    test_amp = test["amp_pct"][mask]
    base_angle = base["angle"][mask]
    test_angle = test["angle"][mask]

    mean_vs = float(np.mean(base_vs))
    base_vs_anom = (base_vs - mean_vs) / mean_vs * 100.0
    test_vs_anom = (test_vs - mean_vs) / mean_vs * 100.0
    diff_vs_ms = (test_vs - base_vs) * 1000.0
    diff_amp = test_amp - base_amp
    diff_angle = angle_diff_180(test_angle, base_angle)
    high_amp = np.maximum(base_amp, test_amp) >= AMP_MASK_PCT

    vs_vlim = max(0.5, safe_percentile(np.abs(np.r_[base_vs_anom, test_vs_anom]), 98.0))
    diff_vs_vlim = max(1.0, safe_percentile(np.abs(diff_vs_ms), 98.0))
    amp_vmax = max(0.5, safe_percentile(np.r_[base_amp, test_amp], 98.0))
    angle_vlim = max(1.0, min(12.0, safe_percentile(np.abs(diff_angle[high_amp]), 98.0, default=1.0)))
    max_amp = max(float(np.max(base_amp)), float(np.max(test_amp)), 1.0e-6)

    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.6), constrained_layout=True)
    scatter_panel(
        axes[0, 0],
        lon,
        lat,
        base_vs_anom,
        title=f"{base_label} Vs anomaly, {depth:g} km",
        cmap="RdBu_r",
        vmin=-vs_vlim,
        vmax=vs_vlim,
        label="dVs from base depth mean (%)",
    )
    draw_sticks(axes[0, 0], lon, lat, base_angle, base_amp, stride=stride, max_amp_pct=max_amp)
    scatter_panel(
        axes[0, 1],
        lon,
        lat,
        test_vs_anom,
        title=f"{test_label} Vs anomaly, {depth:g} km",
        cmap="RdBu_r",
        vmin=-vs_vlim,
        vmax=vs_vlim,
        label="dVs from base depth mean (%)",
    )
    draw_sticks(axes[0, 1], lon, lat, test_angle, test_amp, stride=stride, max_amp_pct=max_amp)
    scatter_panel(
        axes[0, 2],
        lon,
        lat,
        diff_vs_ms,
        title=f"{test_label} - {base_label} Vs",
        cmap="RdBu_r",
        vmin=-diff_vs_vlim,
        vmax=diff_vs_vlim,
        label="m/s",
    )
    scatter_panel(
        axes[1, 0],
        lon,
        lat,
        base_amp,
        title=f"{base_label} anisotropy amplitude",
        cmap="viridis",
        vmin=0.0,
        vmax=amp_vmax,
        label="%",
    )
    draw_sticks(axes[1, 0], lon, lat, base_angle, base_amp, stride=stride, max_amp_pct=max_amp)
    scatter_panel(
        axes[1, 1],
        lon,
        lat,
        test_amp,
        title=f"{test_label} anisotropy amplitude",
        cmap="viridis",
        vmin=0.0,
        vmax=amp_vmax,
        label="%",
    )
    draw_sticks(axes[1, 1], lon, lat, test_angle, test_amp, stride=stride, max_amp_pct=max_amp)
    scatter_panel(
        axes[1, 2],
        lon,
        lat,
        diff_angle,
        title=f"fast-axis diff, amp >= {AMP_MASK_PCT:g}%",
        cmap="RdBu_r",
        vmin=-angle_vlim,
        vmax=angle_vlim,
        label="deg, signed modulo 180",
        mask=high_amp,
    )
    fig.suptitle(f"{comparison_title} at {depth:g} km", fontsize=11)
    stem = f"depth_slice_{depth:g}km".replace(".", "p")
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    return {
        "depth_km": depth,
        "vs_p95_ms": signed_abs_stats(diff_vs_ms)["p95_abs"],
        "amp_p95_pct": signed_abs_stats(diff_amp)["p95_abs"],
        "angle_p95_deg_amp_ge_0p5": signed_abs_stats(diff_angle[high_amp])["p95_abs"],
        "angle_n_amp_ge_0p5": int(np.count_nonzero(high_amp)),
    }


def depthwise_summary(
    base: dict[str, np.ndarray],
    test: dict[str, np.ndarray],
    outdir: Path,
    *,
    base_label: str,
    test_label: str,
    comparison_title: str,
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for depth in np.unique(base["depth"]):
        mask = depth_mask(base, float(depth))
        diff_vs_ms = (test["vs"][mask] - base["vs"][mask]) * 1000.0
        diff_amp = test["amp_pct"][mask] - base["amp_pct"][mask]
        diff_gcgs = np.hypot(test["gc_pct"][mask] - base["gc_pct"][mask], test["gs_pct"][mask] - base["gs_pct"][mask])
        diff_angle = angle_diff_180(test["angle"][mask], base["angle"][mask])
        max_amp = np.maximum(base["amp_pct"][mask], test["amp_pct"][mask])
        high_0p5 = max_amp >= 0.5
        high_1p0 = max_amp >= 1.0
        top20_cut = safe_percentile(max_amp, 80.0)
        top20 = max_amp >= top20_cut
        rows.append(
            {
                "depth_km": float(depth),
                "vs_p95_ms": signed_abs_stats(diff_vs_ms)["p95_abs"],
                "amp_p95_pct": signed_abs_stats(diff_amp)["p95_abs"],
                "gcgs_p95_pct": signed_abs_stats(diff_gcgs)["p95_abs"],
                "angle_p95_all_deg": signed_abs_stats(diff_angle)["p95_abs"],
                "angle_p95_amp_ge_0p5_deg": signed_abs_stats(diff_angle[high_0p5])["p95_abs"],
                "angle_p95_amp_ge_1p0_deg": signed_abs_stats(diff_angle[high_1p0])["p95_abs"],
                "angle_p95_top20_amp_deg": signed_abs_stats(diff_angle[top20])["p95_abs"],
            }
        )

    depths = np.array([row["depth_km"] for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.0), sharex=True, constrained_layout=True)
    axes[0, 0].plot(depths, [row["vs_p95_ms"] for row in rows], marker="o", lw=1.4)
    axes[0, 0].set_ylabel("p95 |dVs| (m/s)")
    axes[0, 0].set_title("Vs difference")
    axes[0, 1].plot(depths, [row["amp_p95_pct"] for row in rows], marker="o", lw=1.4, color="#2c7fb8")
    axes[0, 1].set_ylabel("p95 |dAmp| (%)")
    axes[0, 1].set_title("Anisotropy amplitude difference")
    axes[1, 0].plot(depths, [row["gcgs_p95_pct"] for row in rows], marker="o", lw=1.4, color="#7b3294")
    axes[1, 0].set_ylabel("p95 |d(Gc,Gs)| (%)")
    axes[1, 0].set_title("Gc/Gs vector difference")
    axes[1, 1].plot(depths, [row["angle_p95_all_deg"] for row in rows], marker="o", lw=1.1, label="all")
    axes[1, 1].plot(depths, [row["angle_p95_amp_ge_0p5_deg"] for row in rows], marker="o", lw=1.1, label="amp >= 0.5%")
    axes[1, 1].plot(depths, [row["angle_p95_amp_ge_1p0_deg"] for row in rows], marker="o", lw=1.1, label="amp >= 1.0%")
    axes[1, 1].plot(depths, [row["angle_p95_top20_amp_deg"] for row in rows], marker="o", lw=1.1, label="top 20% amp")
    axes[1, 1].set_ylabel("p95 fast-axis diff (deg)")
    axes[1, 1].set_title("Fast-axis difference by mask")
    axes[1, 1].legend(fontsize=7)
    for ax in axes.ravel():
        ax.grid(True, alpha=0.25, lw=0.5)
        ax.set_xlabel("Depth (km)")
    fig.suptitle(f"Depthwise summary: {test_label} vs {base_label}", fontsize=11)
    fig.savefig(outdir / "depthwise_stability_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "depthwise_stability_summary.pdf", bbox_inches="tight")
    plt.close(fig)
    return rows


def group_period_rows(arr: np.ndarray) -> dict[float, np.ndarray]:
    periods = np.unique(arr[:, 2])
    return {float(period): arr[np.isclose(arr[:, 2], period)] for period in periods}


def plot_period_diagnostics(
    base_dir: Path,
    test_dir: Path,
    outdir: Path,
    *,
    base_label: str,
    test_label: str,
) -> list[dict[str, float]]:
    base = read_numeric_table(base_dir / PERIOD_FILE)
    test = read_numeric_table(test_dir / PERIOD_FILE)
    base_groups = group_period_rows(base)
    test_groups = group_period_rows(test)
    rows: list[dict[str, float]] = []
    for period in sorted(base_groups):
        b = base_groups[period]
        t = test_groups.get(period)
        if t is None:
            continue
        n = min(len(b), len(t))
        b = b[:n]
        t = t[:n]
        diff_phase_ms = (t[:, 3] - b[:, 3]) * 1000.0
        diff_amp_pct = (t[:, 5] - b[:, 5]) * 100.0
        diff_angle = angle_diff_180(t[:, 4], b[:, 4])
        rows.append(
            {
                "period_s": period,
                "phase_mean_ms": signed_abs_stats(diff_phase_ms)["mean_abs"],
                "phase_p95_ms": signed_abs_stats(diff_phase_ms)["p95_abs"],
                "amp_mean_pct": signed_abs_stats(diff_amp_pct)["mean_abs"],
                "amp_p95_pct": signed_abs_stats(diff_amp_pct)["p95_abs"],
                "angle_mean_deg": signed_abs_stats(diff_angle)["mean_abs"],
                "angle_p95_deg": signed_abs_stats(diff_angle)["p95_abs"],
            }
        )
    periods = np.array([row["period_s"] for row in rows])
    fig, axes = plt.subplots(3, 1, figsize=(8.8, 7.2), sharex=True, constrained_layout=True)
    axes[0].plot(periods, [row["phase_mean_ms"] for row in rows], marker="o", lw=1.2, label="mean")
    axes[0].plot(periods, [row["phase_p95_ms"] for row in rows], marker="o", lw=1.2, label="p95")
    axes[0].set_ylabel("|d phase V| (m/s)")
    axes[0].set_title("Period-level phase velocity difference")
    axes[0].legend()
    axes[1].plot(periods, [row["amp_mean_pct"] for row in rows], marker="o", lw=1.2, label="mean", color="#2c7fb8")
    axes[1].plot(periods, [row["amp_p95_pct"] for row in rows], marker="o", lw=1.2, label="p95", color="#41b6c4")
    axes[1].set_ylabel("|d amp| (%)")
    axes[1].set_title("Period-level anisotropy amplitude difference")
    axes[1].legend()
    axes[2].plot(periods, [row["angle_mean_deg"] for row in rows], marker="o", lw=1.2, label="mean", color="#7b3294")
    axes[2].plot(periods, [row["angle_p95_deg"] for row in rows], marker="o", lw=1.2, label="p95", color="#c2a5cf")
    axes[2].set_ylabel("fast-axis diff (deg)")
    axes[2].set_xlabel("Period (s)")
    axes[2].set_title("Period-level fast-axis difference")
    axes[2].legend()
    for ax in axes:
        ax.grid(True, alpha=0.25, lw=0.5)
    fig.suptitle(f"Period diagnostics: {test_label} vs {base_label}", fontsize=11)
    fig.savefig(outdir / "period_diagnostics.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "period_diagnostics.pdf", bbox_inches="tight")
    plt.close(fig)
    return rows


def plot_residual_diagnostics(
    base_dir: Path,
    test_dir: Path,
    outdir: Path,
    *,
    base_label: str,
    test_label: str,
) -> dict[str, float]:
    base = read_numeric_table(base_dir / TRAVELTIME_FILE)
    test = read_numeric_table(test_dir / TRAVELTIME_FILE)
    n = min(len(base), len(test))
    base = base[:n]
    test = test[:n]
    distance = base[:, 0]
    base_res = base[:, 3]
    test_res = test[:, 3]
    diff_res_ms = (test_res - base_res) * 1000.0
    base_wres = base[:, 5]
    test_wres = test[:, 5]
    diff_wres = test_wres - base_wres

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.5), constrained_layout=True)
    bins = np.linspace(
        safe_percentile(np.r_[base_res, test_res], 0.5),
        safe_percentile(np.r_[base_res, test_res], 99.5),
        80,
    )
    axes[0].hist(base_res, bins=bins, histtype="step", lw=1.3, label=base_label)
    axes[0].hist(test_res, bins=bins, histtype="step", lw=1.3, label=test_label)
    axes[0].set_xlabel("Residual (s)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Final residual distribution")
    axes[0].legend()
    axes[1].hist(diff_res_ms, bins=80, color="#525252", alpha=0.85)
    axes[1].set_xlabel(f"{test_label} - {base_label} residual (ms)")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Residual difference")
    axes[2].scatter(distance, diff_res_ms, s=4, alpha=0.22, linewidths=0.0)
    axes[2].axhline(0.0, color="#202020", lw=0.8)
    axes[2].set_xlabel("Distance (km)")
    axes[2].set_ylabel(f"{test_label} - {base_label} residual (ms)")
    axes[2].set_title("Residual difference vs distance")
    for ax in axes:
        ax.grid(True, alpha=0.25, lw=0.5)
    fig.suptitle("Final traveltime residual diagnostics", fontsize=11)
    fig.savefig(outdir / "residual_diagnostics.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "residual_diagnostics.pdf", bbox_inches="tight")
    plt.close(fig)

    out = signed_abs_stats(diff_res_ms)
    out["weighted_residual_p95"] = signed_abs_stats(diff_wres)["p95_abs"]
    return out


def write_markdown(
    outdir: Path,
    base_dir: Path,
    test_dir: Path,
    depth_rows: list[dict[str, float]],
    selected_depth_rows: list[dict[str, float]],
    period_rows: list[dict[str, float]],
    residual_stats: dict[str, float],
    *,
    base_label: str,
    test_label: str,
    comparison_title: str,
    comparison_note: str,
) -> None:
    all_depth = {
        "vs_p95_ms_max": max(row["vs_p95_ms"] for row in depth_rows),
        "amp_p95_pct_max": max(row["amp_p95_pct"] for row in depth_rows),
        "gcgs_p95_pct_max": max(row["gcgs_p95_pct"] for row in depth_rows),
        "angle_p95_top20_max": max(row["angle_p95_top20_amp_deg"] for row in depth_rows),
        "angle_p95_amp_ge_1_max": max(row["angle_p95_amp_ge_1p0_deg"] for row in depth_rows),
    }
    phase_p95_max = max(row["phase_p95_ms"] for row in period_rows)
    amp_period_p95_max = max(row["amp_p95_pct"] for row in period_rows)
    angle_period_p95_max = max(row["angle_p95_deg"] for row in period_rows)
    lines = [
        f"# {comparison_title}",
        "",
        "Comparison:",
        "",
        f"- base ({base_label}): `{base_dir}`",
        f"- test ({test_label}): `{test_dir}`",
        "",
        "This package is an interpretation-level result check. It compares complete",
        "3-D inversion outputs by depth slice instead of reducing the result to one",
        "normalized scalar.",
        "",
        comparison_note,
        "",
        "## Key Figures",
        "",
        "| file | content |",
        "|---|---|",
        "| `depthwise_stability_summary.png` | depth-wise p95 model differences |",
        "| `period_diagnostics.png` | period-level phase velocity, amplitude, and fast-axis differences |",
        "| `residual_diagnostics.png` | final residual distribution and residual differences |",
    ]
    for row in selected_depth_rows:
        depth = row["depth_km"]
        stem = f"depth_slice_{depth:g}km".replace(".", "p")
        lines.append(f"| `{stem}.png` | {base_label}/{test_label}/difference maps at {depth:g} km |")
    lines.extend(
        [
            "",
            "## Summary Metrics",
            "",
            "| diagnostic | value |",
            "|---|---:|",
            f"| max depth-wise p95 `|dVs|` | {all_depth['vs_p95_ms_max']:.3f} m/s |",
            f"| max depth-wise p95 `|dAmp|` | {all_depth['amp_p95_pct_max']:.4f}% |",
            f"| max depth-wise p95 `|d(Gc,Gs)|` | {all_depth['gcgs_p95_pct_max']:.4f}% |",
            f"| max depth-wise p95 fast-axis diff, amp >= 1% | {all_depth['angle_p95_amp_ge_1_max']:.3f} deg |",
            f"| max depth-wise p95 fast-axis diff, top 20% amp | {all_depth['angle_p95_top20_max']:.3f} deg |",
            f"| max period-wise p95 phase velocity diff | {phase_p95_max:.3f} m/s |",
            f"| max period-wise p95 azimuthal amplitude diff | {amp_period_p95_max:.4f}% |",
            f"| max period-wise p95 fast-axis diff | {angle_period_p95_max:.3f} deg |",
            f"| final residual diff p95 | {residual_stats['p95_abs']:.3f} ms |",
            f"| final residual diff max | {residual_stats['max_abs']:.3f} ms |",
            "",
            "## Interpretation",
            "",
            "Use the depth-slice maps first: they show whether the coherent Vs anomaly,",
            "anisotropy amplitude, and fast-axis patterns remain interpretation-level",
            "consistent through depth. The summary curves are supporting diagnostics,",
            "not the primary acceptance criterion.",
            "",
            "For this production comparison, the generated semi-auto run reaches nearly",
            "the same final RMS residual as the original hand-written run, while the",
            "model differences also include changes in the generated initial MOD/MOD_Ref",
            "and runtime backend profile. Therefore these figures should be interpreted",
            "as an original-vs-generated workflow comparison, not a pure backend-only",
            "numerical equivalence test.",
            "",
        ]
    )
    (outdir / "RESULT_COMPARISON.md").write_text("\n".join(lines), encoding="utf-8")
    (outdir / "result_comparison_summary.json").write_text(
        json.dumps(
            {
                "base": str(base_dir),
                "test": str(test_dir),
                "summary": all_depth,
                "period": {
                    "phase_p95_ms_max": phase_p95_max,
                    "amp_p95_pct_max": amp_period_p95_max,
                    "angle_p95_deg_max": angle_period_p95_max,
                },
                "residual": residual_stats,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True, help="base/reference run directory")
    parser.add_argument("--test", type=Path, required=True, help="test/comparison run directory")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--depths", type=float, nargs="+", default=DEFAULT_DEPTHS)
    parser.add_argument("--stride", type=int, default=4)
    parser.add_argument("--base-label", default="base")
    parser.add_argument("--test-label", default="test")
    parser.add_argument("--comparison-title", default="Real DAzi Result Comparison")
    parser.add_argument(
        "--comparison-note",
        default=(
            "If base and test use different initial models or parameter-generation "
            "rules, the differences shown here include those workflow changes."
        ),
    )
    args = parser.parse_args()

    base_dir = args.base.expanduser().resolve()
    test_dir = args.test.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    base = load_model(base_dir)
    test = load_model(test_dir)
    check_same_grid(base, test)

    selected_depths = nearest_depths(base["depth"], args.depths)
    selected_rows = [
        plot_depth_slice(
            base,
            test,
            depth,
            outdir,
            stride=max(1, args.stride),
            base_label=args.base_label,
            test_label=args.test_label,
            comparison_title=args.comparison_title,
        )
        for depth in selected_depths
    ]
    depth_rows = depthwise_summary(
        base,
        test,
        outdir,
        base_label=args.base_label,
        test_label=args.test_label,
        comparison_title=args.comparison_title,
    )
    period_rows = plot_period_diagnostics(
        base_dir,
        test_dir,
        outdir,
        base_label=args.base_label,
        test_label=args.test_label,
    )
    residual_stats = plot_residual_diagnostics(
        base_dir,
        test_dir,
        outdir,
        base_label=args.base_label,
        test_label=args.test_label,
    )
    write_markdown(
        outdir,
        base_dir,
        test_dir,
        depth_rows,
        selected_rows,
        period_rows,
        residual_stats,
        base_label=args.base_label,
        test_label=args.test_label,
        comparison_title=args.comparison_title,
        comparison_note=args.comparison_note,
    )
    print(outdir / "RESULT_COMPARISON.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
