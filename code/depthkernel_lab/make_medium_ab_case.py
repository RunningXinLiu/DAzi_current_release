#!/usr/bin/env python3
"""Build a bounded multi-period DAzi A/B test case from the Anning real case."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


PERIODS = [
    "1",
    "1.20000000000000",
    "1.40000000000000",
    "1.60000000000000",
    "1.80000000000000",
    "2",
    "2.20000000000000",
    "2.40000000000000",
    "2.60000000000000",
    "2.80000000000000",
]


def parse_header(line: str) -> tuple[int, int]:
    parts = line.split()
    if len(parts) < 5 or parts[0] != "#":
        raise ValueError(f"Not a CDTB header: {line.rstrip()}")
    period_index = int(parts[3])
    wave_type = int(parts[4])
    return period_index, wave_type


def filter_cdtb(src: Path, dst: Path, max_period: int, sources_per_period: int) -> dict[int, int]:
    counts = {period: 0 for period in range(1, max_period + 1)}
    kept: list[str] = []
    keep_block = False

    for line in src.read_text().splitlines(keepends=True):
        if line.startswith("#"):
            period, _wave_type = parse_header(line)
            keep_block = period <= max_period and counts[period] < sources_per_period
            if keep_block:
                counts[period] += 1
                kept.append(line)
        elif keep_block:
            kept.append(line)

    dst.write_text("".join(kept))
    return counts


def rewrite_para(src: Path, dst: Path, data_name: str, max_period: int) -> None:
    lines = src.read_text().splitlines()
    out: list[str] = []
    i = 0
    period_block_done = False
    while i < len(lines):
        line = lines[i]
        if i == 3:
            out.append(f"{data_name}\t     	     c: traveltime data file")
        elif "c: maximum of interation" in line:
            out.append("1                                    c: maximum of interation")
        elif "c: kmaxRc" in line and "followed by periods" in line:
            out.append(f"{max_period}                                    c: kmaxRc (followed by periods)")
            out.append(" ".join(PERIODS[:max_period]))
            i += 1
            period_block_done = True
        elif "c: normal voro cell number" in line:
            out.append("10                                   c: normal voro cell number per layer, usually around (nx*ny)/4")
        elif "c: adaptive voro cell number" in line:
            out.append("5                                    c: adaptive voro cell number per layer (adaptive to the ray path density), usually around (nx*ny)/10")
        elif "c: number of realizations" in line:
            out.append("1                                    c: number of realizations (30~100), Total cell number is (nzrand * (ncell + acell))")
        else:
            out.append(line)
        i += 1

    if not period_block_done:
        raise ValueError(f"Could not find period block in {src}")
    dst.write_text("\n".join(out) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-case", required=True, type=Path)
    parser.add_argument("--out-case", required=True, type=Path)
    parser.add_argument("--max-period", type=int, default=5)
    parser.add_argument("--sources-per-period", type=int, default=30)
    args = parser.parse_args()

    if args.out_case.exists():
        raise FileExistsError(args.out_case)
    shutil.copytree(args.source_case, args.out_case)

    data_name = f"CDTB_medium_p{args.max_period}_s{args.sources_per_period}.dat"
    counts = filter_cdtb(
        args.source_case / "CDTB_clustered_6.dat",
        args.out_case / data_name,
        args.max_period,
        args.sources_per_period,
    )
    rewrite_para(args.source_case / "para.in", args.out_case / "para.in", data_name, args.max_period)
    print(f"wrote {args.out_case}")
    print("period_source_counts=" + ",".join(f"{k}:{v}" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    main()
