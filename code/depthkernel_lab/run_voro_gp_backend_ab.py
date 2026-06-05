#!/usr/bin/env python3
"""A/B test legacy vs direct Vorotomo GP assembly."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from run_lsmr_backend_ab import (
    DEFAULT_BINARY,
    DEFAULT_PYTHON,
    compare_outputs,
    parse_log,
    prepare_case,
    run,
)


def write_summary(out_root: Path, summary: dict[str, Any]) -> None:
    (out_root / "summary_voro_gp_backend_ab.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    legacy = summary["legacy"]
    direct = summary["direct"]
    lines = [
        "# Vorotomo GP Backend A/B",
        "",
        f"- case: `{summary['case']}`",
        f"- legacy wall: `{legacy.get('wall_s'):.3f} s`",
        f"- direct wall: `{direct.get('wall_s'):.3f} s`",
        f"- legacy solver/projection: `{legacy.get('solver_s')} s`",
        f"- direct solver/projection: `{direct.get('solver_s')} s`",
        "",
        "## Output Differences",
        "",
        "| file | status | max_abs | rms | mean_abs |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, row in summary["diffs"].items():
        lines.append(
            f"| `{name}` | {row.get('status')} | "
            f"{row.get('max_abs', '')} | {row.get('rms', '')} | {row.get('mean_abs', '')} |"
        )
    if legacy.get("voro_timers") or direct.get("voro_timers"):
        lines.extend(["", "## Vorotomo Timer Rows", ""])
        lines.append("| run | kind | iter | real | cells | map | gp | lsmr | back |")
        lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for label, run_data in (("legacy", legacy), ("direct", direct)):
            for row in run_data.get("voro_timers", []):
                lines.append(
                    f"| {label} | {row['kind']} | {row['iter']} | {row['real']} | "
                    f"{row['cells_s']:.4f} | {row['map_s']:.4f} | {row['gp_s']:.4f} | "
                    f"{row['lsmr_s']:.4f} | {row['back_s']:.4f} |"
                )
    (out_root / "summary_voro_gp_backend_ab.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="A/B test Vorotomo GP assembly backend")
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--case-label", default=None)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--omp-threads", default="10")
    parser.add_argument("--ti-workers", default="4")
    parser.add_argument("--tile-columns", default="64")
    parser.add_argument("--voro-timer", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.out_root.mkdir(parents=True, exist_ok=True)
    legacy_dir = args.out_root / "legacy_gp"
    direct_dir = args.out_root / "direct_gp"
    prepare_case(args.case, legacy_dir, args.force)
    prepare_case(args.case, direct_dir, args.force)

    common_env = {
        "OMP_NUM_THREADS": str(args.omp_threads),
        "DAZI_PYTHON": str(args.python),
        "DAZI_DEPTHKERNEL_BACKEND": "pdsurf_numba_fortran_pv",
        "DAZI_DEPTHKERNEL_TI_BACKEND": "process_batch",
        "DAZI_DEPTHKERNEL_TI_WORKERS": str(args.ti_workers),
        "DAZI_PDSURF_TILE_COLUMNS": str(args.tile_columns),
        "DAZI_DEPTHKERNEL_TI_TILE_COLUMNS": str(args.tile_columns),
    }
    if args.voro_timer:
        common_env["DAZI_VORO_TIMER"] = "1"

    legacy_wall = run(
        [str(args.binary), "para.in"],
        cwd=legacy_dir,
        env={**common_env, "DAZI_VORO_GP_BACKEND": "legacy"},
        stdout=legacy_dir / "run_gp_legacy.log",
    )
    direct_wall = run(
        [str(args.binary), "para.in"],
        cwd=direct_dir,
        env={**common_env, "DAZI_VORO_GP_BACKEND": "direct"},
        stdout=direct_dir / "run_gp_direct.log",
    )

    legacy = parse_log(legacy_dir / "run_gp_legacy.log")
    direct = parse_log(direct_dir / "run_gp_direct.log")
    legacy["wall_s"] = legacy_wall
    direct["wall_s"] = direct_wall
    summary = {
        "case": args.case_label or args.case.name,
        "legacy": legacy,
        "direct": direct,
        "diffs": compare_outputs(legacy_dir, direct_dir),
    }
    write_summary(args.out_root, summary)

    print(f"case={summary['case']}")
    print(f"legacy_wall={legacy_wall:.3f}s direct_wall={direct_wall:.3f}s")
    print(f"legacy_solver={legacy.get('solver_s')} direct_solver={direct.get('solver_s')}")
    print(f"summary={args.out_root / 'summary_voro_gp_backend_ab.md'}")


if __name__ == "__main__":
    main()
