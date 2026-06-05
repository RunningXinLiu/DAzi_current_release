#!/usr/bin/env python3
"""Run DAzi with one user-facing thread config and explicit backend metadata."""

from __future__ import annotations

import argparse
import json
import os
import re
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from dazi_preflight import estimate_case, format_report

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = REPO_ROOT / "bin/DAzimSurfTomo"
DEFAULT_PYTHON = Path("/opt/miniconda3/envs/seisloc/bin/python")
DEFAULT_PDSURF_SENK_SCRIPT = Path(__file__).resolve().parent / "pdsurf_tile_senk_solver.py"
DEFAULT_TI_BATCH_SCRIPT = Path(__file__).resolve().parent / "depthkernel_ti_process_batch.py"
DEFAULT_TI_WORKER = Path(__file__).resolve().parent / "depthkernel_ti_worker"
DEFAULT_LSMR_SCRIPT = Path(__file__).resolve().parent / "scipy_lsmr_triplet_solver.py"

PROFILES = {
    "legacy_baseline": {
        "depthkernel_backend": "default",
        "depthkernel_ti_backend": "default",
        "voro_gp_backend": "legacy",
        "ray_backend": "serial",
        "lsmr_backend": "unset",
        "description": "original Fortran depthkernel/TI and legacy Vorotomo GP",
    },
    "strict_ani_voro": {
        "depthkernel_backend": "default",
        "depthkernel_ti_backend": "process_batch",
        "voro_gp_backend": "direct",
        "ray_backend": "serial",
        "lsmr_backend": "unset",
        "description": "baseline-equivalent ani_voro profile: exact TI batch and direct GP only",
    },
    "fast_exploratory": {
        "depthkernel_backend": "pdsurf_numba_fortran_pv",
        "depthkernel_ti_backend": "process_batch",
        "voro_gp_backend": "direct",
        "ray_backend": "fmm_parallel",
        "lsmr_backend": "scipy",
        "description": "fast profile: pDSurf sensitivities, TI batch, and direct GP",
    },
    "fast_tight_root": {
        "depthkernel_backend": "pdsurf_numba_tight_root_fortran_pv",
        "depthkernel_ti_backend": "process_batch",
        "voro_gp_backend": "direct",
        "ray_backend": "fmm_parallel",
        "lsmr_backend": "scipy",
        "description": "fast profile with tight-root pDSurf sensitivities, TI batch, and direct GP",
    },
}
DEFAULT_PROFILE = "fast_tight_root"

ENV_KEYS = [
    "DAZI_THREADS",
    "DAZI_OMP_THREADS",
    "OMP_NUM_THREADS",
    "OMP_PROC_BIND",
    "OMP_PLACES",
    "DAZI_PYTHON",
    "DAZI_RAY_BACKEND",
    "DAZI_FMM_THREADS",
    "DAZI_DEPTHKERNEL_BACKEND",
    "DAZI_PDSURF_ENGINE",
    "DAZI_PDSURF_SENK_SCRIPT",
    "DAZI_PDSURF_TILE_COLUMNS",
    "DAZI_DEPTHKERNEL_TI_BACKEND",
    "DAZI_DEPTHKERNEL_TI_BATCH_SCRIPT",
    "DAZI_DEPTHKERNEL_TI_WORKER",
    "DAZI_DEPTHKERNEL_TI_WORKERS",
    "DAZI_DEPTHKERNEL_TI_TILE_COLUMNS",
    "DAZI_VORO_GP_BACKEND",
    "DAZI_VORO_GP_SPFRA",
    "DAZI_VORO_TIMER",
    "DAZI_LSMR_BACKEND",
    "DAZI_MAIN_LSMR_BACKEND",
    "DAZI_VORO_LSMR_BACKEND",
    "DAZI_LSMR_SCRIPT",
    "DAZI_LSMR_KEEP_BINARIES",
    "DAZI_AZI_VS_MODE",
]

ALIASES = {
    "dazi_threads": "threads",
    "threadnum": "threads",
    "thread_num": "threads",
    "omp_num_threads": "omp_threads",
    "openmp_threads": "omp_threads",
    "depthkernel_threads": "omp_threads",
    "depthkernel_ti_workers": "ti_workers",
    "ti_threads": "ti_workers",
    "depthkernel_ti_batch_script": "ti_batch_script",
    "ti_batch": "ti_batch_script",
    "depthkernel_ti_worker": "ti_worker",
    "iso_tile_columns": "pdsurf_tile_columns",
    "tile_columns": "pdsurf_tile_columns",
    "depthkernel_ti_tile_columns": "ti_tile_columns",
    "fmm_threads": "threads",
    "gp_backend": "voro_gp_backend",
    "gp_spfra": "voro_gp_spfra",
    "voro_gp_sparse_fraction": "voro_gp_spfra",
}

UNSET_VALUES = {"", "default", "unset", "none"}
TRUE_VALUES = {"1", "true", "yes", "on", "y"}
FALSE_VALUES = {"0", "false", "no", "off", "n"}


def normalize_key(key: str) -> str:
    key = key.strip().lower().replace("-", "_")
    return ALIASES.get(key, key)


def parse_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for lineno, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{lineno}: expected key=value")
        key, value = line.split("=", 1)
        value = value.split("#", 1)[0].strip().strip('"').strip("'")
        values[normalize_key(key)] = value
    return values


def strip_inline_comment(line: str) -> str:
    return re.split(r"\bc\s*:", line, maxsplit=1, flags=re.IGNORECASE)[0].strip()


def parse_para_runtime_controls(path: Path) -> dict[str, str]:
    """Read the user-facing optional runtime block appended to para.in.

    DAzimSurfTomo itself ignores this tail block after the required ISO/Voro
    parameters.  The wrapper reads it so users can keep case parameters in one
    visible file while still setting OpenMP/backend environment variables
    before launching the Fortran executable.
    """
    values: dict[str, str] = {}
    if not path.exists():
        return values
    label_map = {
        "profile": "profile",
        "gp_spfra": "voro_gp_spfra",
        "ray_backend": "ray_backend",
        "lsmr_backend": "lsmr_backend",
        "threads": "threads",
        "fmm_threads": "threads",
        "pdsurf_tile_columns": "pdsurf_tile_columns",
        "ti_tile_columns": "ti_tile_columns",
        "voro_timer": "voro_timer",
    }
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "c:" not in raw_line.lower():
            continue
        value_part, comment = re.split(r"\bc\s*:", raw_line, maxsplit=1, flags=re.IGNORECASE)
        tokens = value_part.split()
        if not tokens:
            continue
        comment_key = comment.split(",", 1)[0].split("(", 1)[0].strip().lower().replace("-", "_")
        for label, key in label_map.items():
            if comment_key == label or comment_key.startswith(label):
                values[key] = tokens[0]
                break
    return values


def config_get(config: dict[str, str], key: str, default: str | None = None) -> str | None:
    return config.get(normalize_key(key), default)


def cli_or_config(value: str | None, config: dict[str, str], key: str, default: str) -> str:
    if value is not None:
        return value
    return config_get(config, key, default) or default


def parse_bool(value: str | None) -> bool:
    if value is None:
        return False
    lowered = value.strip().lower()
    if lowered in TRUE_VALUES:
        return True
    if lowered in FALSE_VALUES or lowered in UNSET_VALUES:
        return False
    raise ValueError(f"invalid boolean value: {value}")


def set_optional_env(env: dict[str, str], name: str, value: str | None) -> None:
    if value is None or value.strip().lower() in UNSET_VALUES:
        env.pop(name, None)
    else:
        env[name] = value.strip()


def parse_extra_env(items: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"--extra-env expects KEY=VALUE, got {item}")
        key, value = item.split("=", 1)
        if not key:
            raise ValueError(f"--extra-env has empty key in {item}")
        out[key] = value
    return out


def rss_to_bytes(value: int) -> int:
    if sys.platform == "darwin":
        return value
    return value * 1024


def parse_positive_int(value: str, name: str) -> int:
    try:
        parsed = int(str(value).strip())
    except ValueError as exc:
        raise ValueError(f"{name} must be a positive integer, got {value!r}") from exc
    if parsed <= 0:
        raise ValueError(f"{name} must be a positive integer, got {value!r}")
    return parsed


def resolve_optional_thread_count(value: str | None, fallback: int, name: str) -> int:
    if value is None or value.strip().lower() in UNSET_VALUES | {"auto"}:
        return fallback
    return parse_positive_int(value, name)


def build_runtime(args: argparse.Namespace, config: dict[str, str]) -> tuple[list[str], dict[str, str], dict[str, str]]:
    profile_name = cli_or_config(args.profile, config, "profile", DEFAULT_PROFILE)
    if profile_name not in PROFILES:
        raise ValueError(f"unknown profile {profile_name!r}; choose one of {', '.join(sorted(PROFILES))}")
    profile = PROFILES[profile_name]

    threads_raw = args.threads or os.environ.get("DAZI_THREADS") or config_get(config, "threads", "4") or "4"
    threads = parse_positive_int(threads_raw, "threads")
    omp_threads = resolve_optional_thread_count(
        args.omp_threads or os.environ.get("DAZI_OMP_THREADS") or config_get(config, "omp_threads"),
        threads,
        "omp_threads",
    )
    ti_workers = resolve_optional_thread_count(
        args.ti_workers or os.environ.get("DAZI_DEPTHKERNEL_TI_WORKERS") or config_get(config, "ti_workers"),
        threads,
        "ti_workers",
    )
    pdsurf_tile_columns = cli_or_config(args.pdsurf_tile_columns, config, "pdsurf_tile_columns", "64")
    ti_tile_columns = cli_or_config(args.ti_tile_columns, config, "ti_tile_columns", pdsurf_tile_columns)

    binary = Path(cli_or_config(args.binary, config, "binary", str(DEFAULT_BINARY))).expanduser()
    infile = cli_or_config(args.infile, config, "infile", "para.in")
    python = Path(cli_or_config(args.python, config, "python", str(DEFAULT_PYTHON))).expanduser()
    pdsurf_senk_script = Path(
        cli_or_config(args.pdsurf_senk_script, config, "pdsurf_senk_script", str(DEFAULT_PDSURF_SENK_SCRIPT))
    ).expanduser()
    ti_batch_script = Path(
        cli_or_config(args.ti_batch_script, config, "ti_batch_script", str(DEFAULT_TI_BATCH_SCRIPT))
    ).expanduser()
    ti_worker = Path(cli_or_config(args.ti_worker, config, "ti_worker", str(DEFAULT_TI_WORKER))).expanduser()
    lsmr_script = Path(cli_or_config(args.lsmr_script, config, "lsmr_script", str(DEFAULT_LSMR_SCRIPT))).expanduser()
    ray_backend = (
        args.ray_backend
        or os.environ.get("DAZI_RAY_BACKEND")
        or config_get(config, "ray_backend", profile["ray_backend"])
        or profile["ray_backend"]
    )
    depthkernel_backend = cli_or_config(
        args.depthkernel_backend,
        config,
        "depthkernel_backend",
        profile["depthkernel_backend"],
    )
    depthkernel_ti_backend = cli_or_config(
        args.depthkernel_ti_backend,
        config,
        "depthkernel_ti_backend",
        profile["depthkernel_ti_backend"],
    )
    voro_gp_backend = cli_or_config(args.voro_gp_backend, config, "voro_gp_backend", profile["voro_gp_backend"])
    voro_gp_spfra = cli_or_config(args.voro_gp_spfra, config, "voro_gp_spfra", "unset")
    lsmr_backend = (
        args.lsmr_backend
        or os.environ.get("DAZI_LSMR_BACKEND")
        or config_get(config, "lsmr_backend", profile["lsmr_backend"])
        or profile["lsmr_backend"]
    )
    main_lsmr_backend = (
        args.main_lsmr_backend
        or os.environ.get("DAZI_MAIN_LSMR_BACKEND")
        or config_get(config, "main_lsmr_backend", "unset")
        or "unset"
    )
    voro_lsmr_backend = (
        args.voro_lsmr_backend
        or os.environ.get("DAZI_VORO_LSMR_BACKEND")
        or config_get(config, "voro_lsmr_backend", "unset")
        or "unset"
    )
    ray_backend = ray_backend.strip().lower()
    lsmr_backend = lsmr_backend.strip().lower()
    main_lsmr_backend = main_lsmr_backend.strip().lower()
    voro_lsmr_backend = voro_lsmr_backend.strip().lower()
    if ray_backend not in {"serial", "fmm_parallel"}:
        raise ValueError("ray_backend must be serial or fmm_parallel")
    for key, value in {
        "lsmr_backend": lsmr_backend,
        "main_lsmr_backend": main_lsmr_backend,
        "voro_lsmr_backend": voro_lsmr_backend,
    }.items():
        if value not in UNSET_VALUES and value not in {"fortran", "scipy"}:
            raise ValueError(f"{key} must be fortran, scipy, or unset")
    pdsurf_engine = cli_or_config(args.pdsurf_engine, config, "pdsurf_engine", "auto")
    if pdsurf_engine.strip().lower() == "auto":
        if "tight_root" in depthkernel_backend:
            pdsurf_engine = "tight_root"
        else:
            pdsurf_engine = "default"
    omp_proc_bind = cli_or_config(args.omp_proc_bind, config, "omp_proc_bind", "unset")
    omp_places = cli_or_config(args.omp_places, config, "omp_places", "unset")

    voro_timer_raw = args.voro_timer
    if voro_timer_raw is None:
        voro_timer_raw = config_get(config, "voro_timer", "0")
    voro_timer = parse_bool(voro_timer_raw)

    env = os.environ.copy()
    env["DAZI_THREADS"] = str(threads)
    env["DAZI_OMP_THREADS"] = str(omp_threads)
    env["OMP_NUM_THREADS"] = str(omp_threads)
    env["DAZI_PYTHON"] = str(python)
    env["DAZI_FMM_THREADS"] = str(omp_threads)
    env["DAZI_PDSURF_SENK_SCRIPT"] = str(pdsurf_senk_script)
    env["DAZI_PDSURF_ENGINE"] = str(pdsurf_engine)
    env["DAZI_PDSURF_TILE_COLUMNS"] = str(pdsurf_tile_columns)
    env["DAZI_DEPTHKERNEL_TI_BATCH_SCRIPT"] = str(ti_batch_script)
    env["DAZI_DEPTHKERNEL_TI_WORKER"] = str(ti_worker)
    env["DAZI_DEPTHKERNEL_TI_WORKERS"] = str(ti_workers)
    env["DAZI_DEPTHKERNEL_TI_TILE_COLUMNS"] = str(ti_tile_columns)
    env["DAZI_LSMR_SCRIPT"] = str(lsmr_script)
    set_optional_env(env, "OMP_PROC_BIND", omp_proc_bind)
    set_optional_env(env, "OMP_PLACES", omp_places)
    set_optional_env(env, "DAZI_RAY_BACKEND", ray_backend)
    set_optional_env(env, "DAZI_DEPTHKERNEL_BACKEND", depthkernel_backend)
    set_optional_env(env, "DAZI_DEPTHKERNEL_TI_BACKEND", depthkernel_ti_backend)
    set_optional_env(env, "DAZI_VORO_GP_BACKEND", voro_gp_backend)
    set_optional_env(env, "DAZI_VORO_GP_SPFRA", voro_gp_spfra)
    set_optional_env(env, "DAZI_LSMR_BACKEND", lsmr_backend)
    set_optional_env(env, "DAZI_MAIN_LSMR_BACKEND", main_lsmr_backend)
    set_optional_env(env, "DAZI_VORO_LSMR_BACKEND", voro_lsmr_backend)
    if voro_timer:
        env["DAZI_VORO_TIMER"] = "1"
    else:
        env.pop("DAZI_VORO_TIMER", None)
    env.update(parse_extra_env(args.extra_env))

    command = [str(binary), infile]
    selected_env = {key: env[key] for key in ENV_KEYS if key in env}
    runtime = {
        "binary": str(binary),
        "infile": infile,
        "profile": profile_name,
        "profile_description": profile["description"],
        "threads": str(threads),
        "thread_policy": "unified_threads_budget; omp_threads and ti_workers default to threads",
        "ray_tracing_backend": str(ray_backend),
        "ray_tracing_threads": str(omp_threads) if ray_backend == "fmm_parallel" else "serial",
        "depthkernel_backend": str(depthkernel_backend),
        "depthkernel_ti_backend": str(depthkernel_ti_backend),
        "depthkernel_openmp_threads": str(omp_threads),
        "voro_realization_threads": str(omp_threads),
        "omp_threads": str(omp_threads),
        "fmm_threads": str(omp_threads),
        "ti_workers": str(ti_workers),
        "pdsurf_engine": str(pdsurf_engine),
        "pdsurf_senk_script": str(pdsurf_senk_script),
        "ti_batch_script": str(ti_batch_script),
        "ti_worker": str(ti_worker),
        "lsmr_script": str(lsmr_script),
        "voro_gp_backend": str(voro_gp_backend),
        "lsmr_backend": str(lsmr_backend),
        "main_lsmr_backend": str(main_lsmr_backend),
        "voro_lsmr_backend": str(voro_lsmr_backend),
        "pdsurf_tile_columns": str(pdsurf_tile_columns),
        "ti_tile_columns": str(ti_tile_columns),
        "voro_gp_spfra": str(voro_gp_spfra),
    }
    return command, env, {"runtime": runtime, "selected_env": selected_env}


def write_metadata(path: Path | None, metadata: dict[str, object]) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def float_config_value(config: dict[str, str], key: str, env_key: str, default: float | None) -> float | None:
    raw = os.environ.get(env_key) or config_get(config, key)
    if raw is None:
        return default
    lowered = raw.strip().lower()
    if lowered in UNSET_VALUES:
        return None
    return float(raw)


def runtime_header(metadata: dict[str, object]) -> str:
    runtime = metadata.get("runtime", {})
    preflight = metadata.get("preflight")
    lines = [
        "DAzi configured run",
        f"case_dir={metadata.get('case_dir')}",
        f"profile={runtime.get('profile')}",
        f"ray_backend={runtime.get('ray_tracing_backend')} fmm_threads={runtime.get('fmm_threads')}",
        f"depthkernel_backend={runtime.get('depthkernel_backend')} depthkernel_ti_backend={runtime.get('depthkernel_ti_backend')}",
        f"lsmr_backend={runtime.get('lsmr_backend')} main_lsmr_backend={runtime.get('main_lsmr_backend')} voro_lsmr_backend={runtime.get('voro_lsmr_backend')}",
        f"threads={runtime.get('threads')} omp_threads={runtime.get('omp_threads')} ti_workers={runtime.get('ti_workers')}",
        f"voro_gp_backend={runtime.get('voro_gp_backend')} voro_gp_spfra={runtime.get('voro_gp_spfra')}",
    ]
    if isinstance(preflight, dict):
        lines.append(format_report(preflight))
    return "\n".join(lines) + "\n" + "=" * 72 + "\n"


def run_case(args: argparse.Namespace) -> int:
    case_dir = args.case_dir.expanduser().resolve()
    if not case_dir.is_dir():
        raise NotADirectoryError(case_dir)

    if args.config is not None:
        config_path = args.config.expanduser().resolve()
        config = parse_config(config_path)
    else:
        config_path = None
        config = {}
    infile_for_para_config = args.infile or config_get(config, "infile", "para.in") or "para.in"
    para_runtime_config = parse_para_runtime_controls(case_dir / infile_for_para_config)
    config = {**config, **para_runtime_config}
    command, env, runtime_meta = build_runtime(args, config)

    metadata: dict[str, object] = {
        "case_dir": str(case_dir),
        "config": str(config_path) if config_path is not None else None,
        "command": command,
        **runtime_meta,
        "started_at": datetime.now(timezone.utc).isoformat(),
    }
    preflight_report: dict[str, object] | None = None
    if not args.skip_preflight:
        runtime = runtime_meta["runtime"]
        warn_gb = float_config_value(config, "preflight_warn_gb", "DAZI_PREFLIGHT_WARN_GB", 64.0)
        hard_gb = args.preflight_hard_gb
        if hard_gb is None:
            hard_gb = float_config_value(config, "preflight_hard_gb", "DAZI_PREFLIGHT_HARD_GB", None)
        preflight_report = estimate_case(
            case_dir,
            runtime["infile"],
            ray_backend=runtime["ray_tracing_backend"],
            threads=runtime["omp_threads"],
            gp_spfra=runtime["voro_gp_spfra"],
            warn_gb=warn_gb,
            hard_gb=hard_gb,
        )
        metadata["preflight"] = preflight_report

    if args.dry_run:
        print(json.dumps(metadata, indent=2, ensure_ascii=False))
        write_metadata(args.metadata, {**metadata, "dry_run": True})
        return 3 if preflight_report and preflight_report.get("hard_exceeded") else 0

    if preflight_report is not None:
        print(format_report(preflight_report))
        if preflight_report.get("hard_exceeded"):
            metadata.update(
                {
                    "finished_at": datetime.now(timezone.utc).isoformat(),
                    "returncode": 3,
                    "process_returncode": None,
                    "program_finished": False,
                    "fortran_stop_detected": None,
                    "wall_s": 0.0,
                    "max_rss_bytes": 0,
                    "max_rss_gb": 0.0,
                    "log": str(args.log) if args.log is not None else None,
                }
            )
            write_metadata(args.metadata, metadata)
            print("preflight hard memory threshold exceeded; run was not started", file=sys.stderr)
            return 3

    if args.log is not None:
        args.log.parent.mkdir(parents=True, exist_ok=True)
    stdout_target = args.log.open("w") if args.log is not None else None
    if stdout_target is not None:
        stdout_target.write(runtime_header(metadata))
        stdout_target.flush()

    t0 = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=case_dir,
            env=env,
            stdout=stdout_target if stdout_target is not None else None,
            stderr=subprocess.STDOUT if stdout_target is not None else None,
            check=False,
        )
    finally:
        if stdout_target is not None:
            stdout_target.close()
    wall_s = time.perf_counter() - t0
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    max_rss_bytes = rss_to_bytes(int(usage.ru_maxrss))
    process_returncode = completed.returncode
    effective_returncode = process_returncode
    program_finished: bool | None = None
    fortran_stop_detected: bool | None = None
    if args.log is not None and args.log.exists():
        log_text = args.log.read_text(encoding="utf-8", errors="replace")
        program_finished = "Program finishes successfully" in log_text
        fortran_stop_detected = re.search(r"(^|\n)\s*STOP\b", log_text) is not None
        if not args.allow_incomplete_log and process_returncode == 0:
            if fortran_stop_detected or not program_finished:
                effective_returncode = 1

    metadata.update(
        {
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "returncode": effective_returncode,
            "process_returncode": process_returncode,
            "program_finished": program_finished,
            "fortran_stop_detected": fortran_stop_detected,
            "wall_s": wall_s,
            "max_rss_bytes": max_rss_bytes,
            "max_rss_gb": max_rss_bytes / (1024**3),
            "log": str(args.log) if args.log is not None else None,
        }
    )
    write_metadata(args.metadata, metadata)
    print(
        f"returncode={effective_returncode} process_returncode={process_returncode} wall={wall_s:.3f}s "
        f"max_rss={metadata['max_rss_gb']:.3f}GB"
    )
    if args.metadata is not None:
        print(f"metadata={args.metadata}")
    if args.log is not None:
        print(f"log={args.log}")
    return effective_returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run DAzi with unified CPU/backend settings")
    parser.add_argument("case_dir", type=Path)
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="optional legacy key=value config; by default runtime controls are read from CASE_DIR/para.in",
    )
    parser.add_argument("--binary", default=None)
    parser.add_argument("--python", default=None)
    parser.add_argument("--infile", default=None)
    parser.add_argument(
        "--profile",
        default=None,
        choices=sorted(PROFILES),
        help=f"backend profile; default is {DEFAULT_PROFILE}",
    )
    parser.add_argument(
        "--threads",
        default=None,
        help="single CPU budget; used for FMM/OpenMP, Vorotomo realizations, and TI workers unless overridden",
    )
    parser.add_argument("--omp-threads", default=None, help="advanced override for OpenMP/FMM/depthkernel threads")
    parser.add_argument("--ti-workers", default=None, help="advanced override for TI process-batch worker count")
    parser.add_argument("--pdsurf-tile-columns", default=None)
    parser.add_argument("--pdsurf-engine", default=None, help="pDSurf solver engine: auto/default/tight_root/etc.")
    parser.add_argument("--pdsurf-senk-script", default=None, help="path to pdsurf_tile_senk_solver.py")
    parser.add_argument("--ti-batch-script", default=None, help="path to depthkernel_ti_process_batch.py")
    parser.add_argument("--ti-worker", default=None, help="path to compiled depthkernel_ti_worker")
    parser.add_argument("--lsmr-script", default=None, help="path to scipy_lsmr_triplet_solver.py")
    parser.add_argument("--ti-tile-columns", default=None)
    parser.add_argument("--ray-backend", default=None, choices=("serial", "fmm_parallel"))
    parser.add_argument("--depthkernel-backend", default=None)
    parser.add_argument("--depthkernel-ti-backend", default=None)
    parser.add_argument("--voro-gp-backend", default=None)
    parser.add_argument(
        "--voro-gp-spfra",
        default=None,
        help="Voronoi projection sparse allocation fraction; separate from para.in spfra",
    )
    parser.add_argument("--lsmr-backend", default=None)
    parser.add_argument("--main-lsmr-backend", default=None)
    parser.add_argument("--voro-lsmr-backend", default=None)
    parser.add_argument("--omp-proc-bind", default=None)
    parser.add_argument("--omp-places", default=None)
    parser.add_argument("--voro-timer", default=None, nargs="?", const="1")
    parser.add_argument("--extra-env", action="append", default=[], help="additional KEY=VALUE environment override")
    parser.add_argument("--log", type=Path, default=None)
    parser.add_argument("--metadata", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-preflight", action="store_true", help="skip case-size and memory preflight")
    parser.add_argument("--preflight-hard-gb", type=float, default=None, help="abort before running if estimated peak upper bound exceeds this many GB")
    parser.add_argument(
        "--allow-incomplete-log",
        action="store_true",
        help="do not convert a missing success marker or Fortran STOP in the log into a failed return code",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        raise SystemExit(run_case(args))
    except Exception as exc:
        print(f"run_dazi_configured: {exc}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
