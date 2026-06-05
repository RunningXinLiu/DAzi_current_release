#!/usr/bin/env python3
"""Convert DAzi inversion outputs to MOD-format true files for forward modeling.

Inputs:
- DSurfTomo.inv  -> MODVs.true (full grid)
- Gc_Gs_model.inv -> MODGc.true, MODGs.true (inner grid, converted to fraction)

This tool is intended for creating a synthetic forward truth case from inversion
results without re-running inversion.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from math import isnan


Field = tuple[float, float, float, float]

def _coord_key(v: float) -> str:
    # Robust key for float coordinates with small tolerance around decimal text outputs.
    return f"{v:.10f}"


def _load_unique_sorted(values: list[float]) -> list[float]:
    return sorted({round(float(v), 10) for v in values})


def _index_maps(coords: list[float]) -> dict[str, int]:
    return {_coord_key(v): i for i, v in enumerate(coords)}


def _parse_iso_inv(path: Path) -> tuple[list[float], list[float], list[float], list[list[list[float]]]]:
    rows: list[tuple[float, float, float, float]] = []
    with path.open(encoding="utf-8", errors="replace") as f:
        for lineno, raw in enumerate(f, 1):
            s = raw.split()
            if not s:
                continue
            if len(s) < 4:
                raise ValueError(f"{path}:{lineno} has fewer than 4 columns")
            try:
                lon, lat, dep, vs = (float(s[0]), float(s[1]), float(s[2]), float(s[3]))
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno} parse error: {raw!r}") from exc
            rows.append((lon, lat, dep, vs))

    if not rows:
        raise ValueError(f"{path} is empty")

    lons = _load_unique_sorted([r[0] for r in rows])
    lats = _load_unique_sorted([r[1] for r in rows])
    deps = _load_unique_sorted([r[2] for r in rows])

    lon_map = _index_maps(lons)
    lat_map = _index_maps(lats)
    dep_map = _index_maps(deps)

    nlon, nlat, ndep = len(lons), len(lats), len(deps)
    cube = [[[float("nan") for _ in range(ndep)] for _ in range(nlon)] for _ in range(nlat)]

    for lon, lat, dep, vs in rows:
        i = lon_map[_coord_key(lon)]
        j = lat_map[_coord_key(lat)]
        k = dep_map[_coord_key(dep)]
        if not isnan(cube[j][i][k]):
            raise ValueError(f"Duplicate grid point: lon={lon}, lat={lat}, dep={dep} in {path}")
        cube[j][i][k] = vs

    for j in range(nlat):
        for i in range(nlon):
            for k in range(ndep):
                if isnan(cube[j][i][k]):
                    raise ValueError(
                        f"Missing value in {path}: lon={lons[i]}, lat={lats[j]}, dep={deps[k]}"
                    )

    return lons, lats, deps, cube


def _parse_ani_inv(
    path: Path,
) -> tuple[list[float], list[float], list[float], list[list[list[float]]], list[list[list[float]]]]:
    rows: list[tuple[float, float, float, float, float]] = []
    with path.open(encoding="utf-8", errors="replace") as f:
        for lineno, raw in enumerate(f, 1):
            s = raw.split()
            if not s:
                continue
            if len(s) < 8:
                raise ValueError(f"{path}:{lineno} has fewer than 8 columns")
            try:
                lon = float(s[0])
                lat = float(s[1])
                dep = float(s[2])
                gc = float(s[6]) / 100.0
                gs = float(s[7]) / 100.0
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno} parse error: {raw!r}") from exc
            rows.append((lon, lat, dep, gc, gs))

    if not rows:
        raise ValueError(f"{path} is empty")

    lons = _load_unique_sorted([r[0] for r in rows])
    lats = _load_unique_sorted([r[1] for r in rows])
    deps = _load_unique_sorted([r[2] for r in rows])

    lon_map = _index_maps(lons)
    lat_map = _index_maps(lats)
    dep_map = _index_maps(deps)

    nlon, nlat, ndep = len(lons), len(lats), len(deps)
    gc_cube = [[[float("nan") for _ in range(ndep)] for _ in range(nlon)] for _ in range(nlat)]
    gs_cube = [[[float("nan") for _ in range(ndep)] for _ in range(nlon)] for _ in range(nlat)]

    for lon, lat, dep, gc, gs in rows:
        i = lon_map[_coord_key(lon)]
        j = lat_map[_coord_key(lat)]
        k = dep_map[_coord_key(dep)]
        if not isnan(gc_cube[j][i][k]) or not isnan(gs_cube[j][i][k]):
            raise ValueError(
                f"Duplicate anisotropic grid point: lon={lon}, lat={lat}, dep={dep} in {path}"
            )
        gc_cube[j][i][k] = gc
        gs_cube[j][i][k] = gs

    for j in range(nlat):
        for i in range(nlon):
            for k in range(ndep):
                if isnan(gc_cube[j][i][k]) or isnan(gs_cube[j][i][k]):
                    raise ValueError(
                        f"Missing anisotropic value in {path}: lon={lons[i]}, lat={lats[j]}, dep={deps[k]}"
                    )

    return lons, lats, deps, gc_cube, gs_cube


def _write_mod_true(path: Path, lons: list[float], lats: list[float], deps: list[float], cube: list[list[list[float]]]) -> None:
    nlat = len(lats)
    nlon = len(lons)
    ndep = len(deps)
    with path.open("w", encoding="utf-8") as f:
        f.write(" ".join(f"{d:g}" for d in deps) + "\n")
        # MOD file stores per depth, all lons are rows; each row is (north->south) lat.
        for k in range(ndep):
            for j in range(nlon):
                row = [cube[nlat - 1 - i][j][k] for i in range(nlat)]
                f.write("".join(f"{v:8.4f}" for v in row) + "\n")


def _write_inner_g(path: Path, lons: list[float], lats: list[float], deps: list[float], cube: list[list[list[float]]]) -> None:
    # cube shape = (nlat, nlon, ndep), with nlat,nlon from inner anisotropic grid.
    nlat = len(lats)
    nlon = len(lons)
    ndep = len(deps)
    with path.open("w", encoding="utf-8") as f:
        for k in range(ndep):
            for j in range(nlon):
                row = [cube[nlat - 1 - i][j][k] for i in range(nlat)]
                f.write("".join(f"{v:8.4f}" for v in row) + "\n")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build MODVs.true/MODGc.true/MODGs.true from inversion outputs")
    parser.add_argument("case_dir", type=Path, help="case directory containing DSurfTomo.inv and Gc_Gs_model.inv")
    parser.add_argument("--iso", default="DSurfTomo.inv", help="path to isotropic inversion output")
    parser.add_argument("--ani", default="Gc_Gs_model.inv", help="path to anisotropic inversion output")
    parser.add_argument("--vs-true", default="MODVs.true", help="output MODVs.true name")
    parser.add_argument("--gc-true", default="MODGc.true", help="output MODGc.true name")
    parser.add_argument("--gs-true", default="MODGs.true", help="output MODGs.true name")
    parser.add_argument("--force", action="store_true", help="overwrite outputs if they already exist")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    case_dir: Path = args.case_dir
    if not case_dir.exists():
        raise FileNotFoundError(f"case_dir not found: {case_dir}")

    iso_path = case_dir / args.iso
    ani_path = case_dir / args.ani
    if not iso_path.exists():
        raise FileNotFoundError(f"missing inversion iso file: {iso_path}")
    if not ani_path.exists():
        raise FileNotFoundError(f"missing inversion ani file: {ani_path}")

    out_vs = case_dir / args.vs_true
    out_gc = case_dir / args.gc_true
    out_gs = case_dir / args.gs_true

    for p in [out_vs, out_gc, out_gs]:
        if p.exists() and not args.force:
            raise FileExistsError(f"{p} exists; re-run with --force")

    lons_iso, lats_iso, deps_iso, vs_cube = _parse_iso_inv(iso_path)
    lons_ani, lats_ani, deps_ani, gc_cube, gs_cube = _parse_ani_inv(ani_path)

    if len(deps_iso) < 2:
        raise ValueError(f"unexpected depth size in {iso_path}: {len(deps_iso)}")

    if len(lons_iso) <= 1 or len(lats_iso) <= 1:
        raise ValueError(f"invalid iso grid size in {iso_path}")

    if len(deps_ani) != len(deps_iso) - 1:
        raise ValueError(
            f"anisotropic depth count mismatch: {ani_path} has {len(deps_ani)} layers, "
            f"but iso has {len(deps_iso)} nodes"
        )
    if len(lons_ani) != len(lons_iso) - 2:
        raise ValueError(
            f"anisotropic longitude count mismatch: {len(lons_ani)} vs inner of iso {len(lons_iso)}"
        )
    if len(lats_ani) != len(lats_iso) - 2:
        raise ValueError(
            f"anisotropic latitude count mismatch: {len(lats_ani)} vs inner of iso {len(lats_iso)}"
        )

    _write_mod_true(out_vs, lons_iso, lats_iso, deps_iso, vs_cube)
    _write_inner_g(out_gc, lons_ani, lats_ani, deps_ani, gc_cube)
    _write_inner_g(out_gs, lons_ani, lats_ani, deps_ani, gs_cube)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
